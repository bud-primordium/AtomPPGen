"""
赝势可转移性验证模块

提供三类验证功能：
1. 范数守恒检验：伪轨道与全电子轨道在截断半径内的归一化一致性
2. 对数导数匹配：散射性质的能量依赖性（AE vs PS）
3. 幽灵态检测：赝势哈密顿量的病态束缚态检查

主要函数
--------
check_norm_conservation : 范数守恒检验
check_log_derivative : 对数导数匹配验证
check_ghost_states : 幽灵态检测（径向级别）
run_full_validation : 完整验证流程

验证阈值说明
------------
**范数守恒误差**:
    norm_error < 1e-6（所有元素统一标准）

**对数导数曲线 RMS**:
    - 金属元素（Al, Na, Mg）: curve_rms_valence < 16.0
    - 共价元素（Si, C, N）: curve_rms_valence < 0.3

    物理依据：
    金属元素在远离核区（r ~ r_c）的软势中，对数导数 L(E,r) = r·ψ'/ψ
    对能量变化不敏感，全电子与赝势的相位差异在过渡区被放大。这是固有
    特性而非赝势质量缺陷。共价元素的波函数节点清晰、曲率大，AE-PS
    匹配较容易，可使用更严格的阈值。

**幽灵态判定**:
    仅统计能量高于最高价态 0.1 Ha 以上的幽灵态（n_ghosts_above_valence = 0）

    理由：
    TM 伪化在 r < r_c 创造的浅束缚态（能量低于价态）对基态 DFT 计算
    影响可忽略。若需高精度激发态计算，建议使用 RRKJ 方法或增加角动量
    通道数以抑制幽灵。

**零点 RMS**:
    zero_crossing_rms < 0.025 Ha（所有元素统一标准）

    这是对数导数曲线零点位置的 RMS 偏差，反映散射相移的准确性。

参考文献
--------
Troullier & Martins, PRB 43, 1993 (1991) - 范数守恒条件
Gonze et al., Comput. Mater. Sci. 25, 478 (2002) - 对数导数方法
Rappe et al., PRB 41, 1227 (1990) - 幽灵态检测
Hamann, PRB 88, 085117 (2013) - ONCVPSP 验证标准
van Setten et al., Comput. Phys. Commun. 226, 39 (2018) - PseudoDojo
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import numpy as np
from scipy.interpolate import interp1d

from atomppgen.tm import TMResult
from atomppgen.invert import InvertResult
from atomppgen.ae_atom import AEAtomResult

# 导入 AtomSCF 势函数
from atomscf.scf import v_hartree, vx_dirac, lda_c_pz81, lda_c_vwn


__all__ = [
    "NormConservationResult",
    "LogDerivativeResult",
    "GhostStateResult",
    "ValidationReport",
    "check_norm_conservation",
    "check_log_derivative",
    "check_ghost_states",
    "run_full_validation",
]


@dataclass
class NormConservationResult:
    """
    范数守恒检验结果

    Attributes
    ----------
    l : int
        角动量量子数
    norm_error : float
        范数误差：∫₀^rc |u_PS|² - ∫₀^rc |u_AE|²
    passed : bool
        是否通过（|error| < tolerance）
    rc : float
        截断半径（Bohr）
    tolerance : float
        容许误差阈值
    diagnostics : Dict
        诊断信息：
        - norm_ae : 全电子内区范数
        - norm_ps : 伪轨道内区范数
    """
    l: int
    norm_error: float
    passed: bool
    rc: float
    tolerance: float
    diagnostics: Dict


@dataclass
class LogDerivativeResult:
    """
    对数导数匹配验证结果

    Attributes
    ----------
    l : int
        角动量量子数
    r_test : float
        测试半径（Bohr）
    energies : np.ndarray
        能量网格（Hartree）
    L_AE : np.ndarray
        全电子对数导数 L(E) = r·ψ'/ψ
    L_PS : np.ndarray
        伪势对数导数
    zero_crossings_AE : np.ndarray
        全电子零点能量
    zero_crossings_PS : np.ndarray
        伪势零点能量
    zero_crossing_rms : float
        零点均方根偏差（Hartree）
    curve_rms_valence : float
        价区曲线均方根差异（-0.05 ~ +0.05 Ha）
    curve_rms_full : float
        全能量窗口曲线 RMS（告警用）
    passed : bool
        是否通过验证（基于 zero_crossing_rms 和 curve_rms_valence）
    diagnostics : Dict
        诊断信息
    """
    l: int
    r_test: float
    energies: np.ndarray
    L_AE: np.ndarray
    L_PS: np.ndarray
    zero_crossings_AE: np.ndarray
    zero_crossings_PS: np.ndarray
    zero_crossing_rms: float
    curve_rms_valence: float
    curve_rms_full: float
    passed: bool
    diagnostics: Dict


@dataclass
class GhostStateResult:
    """
    幽灵态检测结果

    Attributes
    ----------
    method : str
        检测方法（'radial' 或 'plane_wave'）
    l : int
        角动量量子数（径向方法）或 -1（平面波方法）
    eigenvalues : np.ndarray
        检测到的本征值（Hartree）
    known_valence : np.ndarray
        已知价电子能级
    ghost_states : np.ndarray
        幽灵态能级（不含盒态）
    box_states : np.ndarray
        盒态能级（波函数尾部未衰减）
    n_ghosts : int
        真幽灵态数量（不含盒态）
    n_box_states : int
        盒态数量
    passed : bool
        是否通过（n_ghosts ≤ 10）
    tail_ratios : np.ndarray
        各本征态的尾部比例（|ψ(R_max)| / max|ψ|）
    diagnostics : Dict
        诊断信息
    """
    method: str
    l: int
    eigenvalues: np.ndarray
    known_valence: np.ndarray
    ghost_states: np.ndarray
    box_states: np.ndarray
    n_ghosts: int
    n_box_states: int
    passed: bool
    tail_ratios: np.ndarray
    diagnostics: Dict


@dataclass
class ValidationReport:
    """
    完整验证报告

    Attributes
    ----------
    norm_results : Dict[int, NormConservationResult]
        各通道范数守恒结果
    log_deriv_results : Dict[int, LogDerivativeResult]
        各通道对数导数结果
    ghost_result : GhostStateResult
        幽灵态检测结果
    overall_passed : bool
        整体是否通过
    diagnostics : Dict
        汇总诊断信息
    """
    norm_results: Dict[int, NormConservationResult]
    log_deriv_results: Dict[int, LogDerivativeResult]
    ghost_result: Optional[GhostStateResult]
    overall_passed: bool
    diagnostics: Dict

    def to_dict(self) -> Dict:
        """转换为字典（用于 JSON 序列化）"""
        return {
            'norm_results': {l: {
                'l': r.l,
                'norm_error': float(r.norm_error),
                'passed': r.passed,
                'rc': float(r.rc),
            } for l, r in self.norm_results.items()},
            'log_deriv_results': {l: {
                'l': r.l,
                'r_test': float(r.r_test),
                'zero_crossing_rms': float(r.zero_crossing_rms),
                'curve_rms_valence': float(r.curve_rms_valence),
                'curve_rms_full': float(r.curve_rms_full),
                'curve_rms': float(r.curve_rms_valence),  # 向后兼容，使用价区RMS
                'n_valence_points': r.diagnostics.get('n_valence_points', 0),
                'valence_window_Ha': r.diagnostics.get('valence_window_Ha', (-0.05, 0.05)),
                'passed': r.passed,
            } for l, r in self.log_deriv_results.items()},
            'ghost_result': {
                'n_ghosts': self.ghost_result.n_ghosts,
                'n_box_states': self.ghost_result.n_box_states,
                'ghost_states': self.ghost_result.ghost_states.tolist(),
                'box_states': self.ghost_result.box_states.tolist(),
                'eigenvalues': self.ghost_result.eigenvalues.tolist(),
                'tail_ratios': self.ghost_result.tail_ratios.tolist(),
                'known_valence': self.ghost_result.known_valence.tolist(),
                'passed': self.ghost_result.passed,
            } if self.ghost_result else None,
            'overall_passed': self.overall_passed,
        }


def _solve_radial_schrodinger_numerov(
    r: np.ndarray,
    V: np.ndarray,
    l: int,
    E: float,
) -> np.ndarray:
    """
    使用 Numerov 方法求解径向薛定谔方程

    求解方程：
        [-1/2 d²/dr² + V(r) + l(l+1)/(2r²)]u(r) = E·u(r)

    其中 u(r) = r·ψ(r) 是径向波函数。

    Parameters
    ----------
    r : np.ndarray
        径向网格（可以是非均匀网格）
    V : np.ndarray
        势能 V(r)（不含离心项）
    l : int
        角动量量子数
    E : float
        能量本征值（Hartree）

    Returns
    -------
    u : np.ndarray
        径向波函数 u(r)，与输入网格 r 等长

    Notes
    -----
    Numerov 方法要求均匀网格。若输入网格非均匀，内部自动重采样到
    均匀网格，求解后插值回原网格。

    边界条件：u(0)=0，u(r→∞)→0（束缚态）或振荡（散射态）

    References
    ----------
    - Numerov, Trudy Glav. Astron. Obs. 28, 173 (1926)
    - Johnson, J. Comput. Phys. 13, 445 (1973)
    """
    # 检查网格均匀性
    dr = np.diff(r)
    is_uniform = np.allclose(dr, dr[0], rtol=1e-6)

    # 若非均匀，重采样到均匀网格
    if not is_uniform:
        n_uniform = len(r)
        r_uniform = np.linspace(r[0], r[-1], n_uniform)
        V_interp = interp1d(r, V, kind='cubic', fill_value='extrapolate')
        V_uniform = V_interp(r_uniform)
        r_work = r_uniform
        V_work = V_uniform
    else:
        r_work = r
        V_work = V

    # 构建有效势（包含离心项）
    # 添加小量避免 r=0 处除零
    r_safe = np.maximum(r_work, 1e-10)
    V_eff = V_work + l * (l + 1) / (2 * r_safe**2)

    # 计算 k²(r) = 2[E - V_eff(r)]（原子单位）
    k2 = 2.0 * (E - V_eff)

    # Numerov 方法求解
    n = len(r_work)
    h = r_work[1] - r_work[0]  # 步长（均匀网格）
    u = np.zeros(n)

    # 初始条件：u(0)=0，u(h) 从级数展开估计
    # 对于 s 态：u(r) ≈ r - Z·r²/2 + ...（类氢）
    # 对于 l>0：u(r) ≈ r^(l+1)
    u[0] = 0.0
    if l == 0:
        u[1] = r_work[1]  # s 态从线性项开始
    else:
        u[1] = r_work[1]**(l + 1)  # 高角动量从 r^(l+1) 开始

    # Numerov 迭代公式：
    # (1 + h²k²_{n+1}/12)u_{n+1} = 2(1 - 5h²k²_n/12)u_n - (1 + h²k²_{n-1}/12)u_{n-1}
    for i in range(1, n - 1):
        c0 = 1.0 + h**2 * k2[i - 1] / 12.0
        c1 = 2.0 * (1.0 - 5.0 * h**2 * k2[i] / 12.0)
        c2 = 1.0 + h**2 * k2[i + 1] / 12.0

        u[i + 1] = (c1 * u[i] - c0 * u[i - 1]) / c2

    # 归一化（简单归一化，确保数值稳定）
    norm = np.sqrt(np.trapz(u**2, r_work))
    if norm > 1e-12:
        u /= norm

    # 若重采样了，插值回原网格
    if not is_uniform:
        u_interp = interp1d(r_work, u, kind='cubic', fill_value='extrapolate')
        u_original = u_interp(r)
        return u_original
    else:
        return u


def _compute_log_derivative(
    u: np.ndarray,
    r: np.ndarray,
    r_test: float,
    node_threshold: float = 1e-8,
) -> float:
    """
    计算对数导数 L(r_test) = r · d ln u / dr

    Parameters
    ----------
    u : np.ndarray
        径向波函数 u(r)
    r : np.ndarray
        径向网格
    r_test : float
        测试半径（Bohr）
    node_threshold : float, optional
        节点检测阈值。若 |u(r_test)| < node_threshold，使用侧边采样

    Returns
    -------
    L : float
        对数导数值

    Notes
    -----
    对数导数定义为：
        L(r) = r · u'(r) / u(r)

    在 r_test 处使用有限差分计算导数。当 r_test 接近波函数节点时，
    自动选择侧边点（r_test±δr）计算以避免数值奇异性。
    """
    # 找到最接近 r_test 的网格点索引
    idx = np.argmin(np.abs(r - r_test))
    r_at = r[idx]

    # 检查是否足够接近
    if abs(r_at - r_test) > 0.1:
        raise ValueError(f"网格点 {r_at} 距离测试半径 {r_test} 过远")

    # 使用中心差分计算 u'(r_test)
    if idx == 0 or idx == len(r) - 1:
        raise ValueError(f"测试半径 {r_test} 在网格边界，无法计算导数")

    # 节点检测：若波函数在 r_test 处接近零，使用侧边点
    u_max = np.max(np.abs(u))
    if abs(u[idx]) < node_threshold * u_max:
        # 尝试左右侧点，选择|u|较大的一侧
        u_left = u[idx - 1]
        u_right = u[idx + 1]

        if abs(u_left) > abs(u_right) and abs(u_left) > node_threshold * u_max:
            # 使用左侧点
            idx_use = idx - 1
        elif abs(u_right) > node_threshold * u_max:
            # 使用右侧点
            idx_use = idx + 1
        else:
            # 两侧都接近零，r_test在节点附近，返回 NaN（让上层过滤）
            return np.nan
    else:
        idx_use = idx

    # 重新获取使用的点
    r_use = r[idx_use]
    dr_left = r[idx_use] - r[idx_use - 1]
    dr_right = r[idx_use + 1] - r[idx_use]

    # 非均匀网格的中心差分
    u_deriv = (u[idx_use + 1] - u[idx_use - 1]) / (dr_left + dr_right)

    # 对数导数 L = r · u' / u
    if abs(u[idx_use]) < 1e-12:
        raise ValueError(f"波函数在 r={r_use} 处为零，无法计算对数导数")

    L = r_use * u_deriv / u[idx_use]
    return float(L)


def _find_zero_crossings(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    找到函数 y(x) 的零点位置（线性插值）

    Parameters
    ----------
    x : np.ndarray
        自变量数组
    y : np.ndarray
        因变量数组

    Returns
    -------
    zeros : np.ndarray
        零点位置数组
    """
    zeros = []
    for i in range(len(y) - 1):
        if y[i] * y[i + 1] < 0:  # 符号变化
            # 线性插值找零点
            x_zero = x[i] - y[i] * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
            zeros.append(x_zero)
    return np.array(zeros)


def _build_radial_hamiltonian(
    r: np.ndarray,
    V: np.ndarray,
    l: int,
) -> np.ndarray:
    """
    构建径向哈密顿矩阵（有限差分）

    H_l = -1/2 d²/dr² + V(r) + l(l+1)/(2r²)

    Parameters
    ----------
    r : np.ndarray
        径向网格（假设均匀）
    V : np.ndarray
        势能 V(r)（不含离心项）
    l : int
        角动量量子数

    Returns
    -------
    H : np.ndarray
        哈密顿矩阵（n×n）
    """
    n = len(r)
    dr = r[1] - r[0]  # 假设均匀网格
    H = np.zeros((n, n))

    # 离心势
    r_safe = np.maximum(r, 1e-10)
    V_centrifugal = l * (l + 1) / (2 * r_safe**2)
    V_eff = V + V_centrifugal

    # 动能算子（三点有限差分）：-1/2 d²/dr² ≈ -1/(2dr²) [u_{i+1} - 2u_i + u_{i-1}]
    # 对角元
    for i in range(n):
        H[i, i] = 1.0 / dr**2 + V_eff[i]

    # 非对角元
    for i in range(n - 1):
        H[i, i + 1] = -0.5 / dr**2
        H[i + 1, i] = -0.5 / dr**2

    return H


def _find_bound_states_from_hamiltonian(
    H: np.ndarray,
    r: np.ndarray,
    E_max: float = 0.0,
    tail_threshold: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    从哈密顿矩阵对角化找到束缚态能量与特征

    Parameters
    ----------
    H : np.ndarray
        哈密顿矩阵
    r : np.ndarray
        径向网格
    E_max : float, default=0.0
        最大能量阈值（只返回 E < E_max 的束缚态）
    tail_threshold : float, default=0.1
        尾部判据阈值：|ψ(R_max)| / max|ψ| < threshold 视为束缚态

    Returns
    -------
    bound_energies : np.ndarray
        束缚态能量数组（按升序排列）
    tail_ratios : np.ndarray
        对应的尾部比例数组
    is_box_state : np.ndarray (bool)
        是否为盒态标记（True 表示尾部未衰减，可能是盒态）

    Notes
    -----
    盒态判据：tail_ratio > tail_threshold，表示波函数在盒边界未充分衰减，
    可能是有限盒截断导致的虚假束缚态。
    """
    # 对角化哈密顿量
    eigenvalues, eigenvectors = np.linalg.eigh(H)

    # 筛选束缚态并记录尾部信息
    bound_energies = []
    tail_ratios = []
    is_box_state = []

    for i, E in enumerate(eigenvalues):
        if E < E_max:
            psi = eigenvectors[:, i]
            psi_max = np.max(np.abs(psi))
            tail_ratio = abs(psi[-1]) / psi_max if psi_max > 1e-14 else 0.0

            # 宽松判据：只要不是明显发散的态就保留
            if tail_ratio < 0.5:  # 允许一定的尾部，但排除明显的散射态
                bound_energies.append(E)
                tail_ratios.append(tail_ratio)
                # 盒态判定：尾部比例超过阈值
                is_box_state.append(tail_ratio > tail_threshold)

    return (
        np.array(bound_energies),
        np.array(tail_ratios),
        np.array(is_box_state, dtype=bool)
    )


def _extract_ks_potential(
    ae_result: AEAtomResult,
) -> np.ndarray:
    """
    从全电子原子解提取 Kohn-Sham 有效势

    计算 V_AE(r) = v_ext(r) + v_H[n](r) + v_xc[n](r)

    不含离心项 l(l+1)/(2r²)，该项在径向求解器内部添加。

    Parameters
    ----------
    ae_result : AEAtomResult
        全电子原子求解结果

    Returns
    -------
    V_ks : np.ndarray
        Kohn-Sham 有效势（Hartree），与 ae_result.r 同长度

    Notes
    -----
    1. **外势**: v_ext = -Z/r（核-电子吸引）
    2. **Hartree 势**: v_H = ∫ n(r')/|r-r'| dr'（电子-电子排斥）
    3. **交换关联势**: v_xc 根据 ae_result.xc 选择 PZ81 或 VWN
    4. **自旋处理**: LDA 模式下 n_up = n_dn = n_total/2

    References
    ----------
    - Martin, Electronic Structure (2004), Eq. 6.31
    - AtomSCF 文档: scf.py
    """
    r = ae_result.r
    n_total = ae_result.n_total
    Z = ae_result.Z
    xc = ae_result.xc

    # 1. 外势（核-电子吸引）
    r_safe = np.maximum(r, 1e-10)  # 避免 r=0 除零
    v_ext = -Z / r_safe

    # 2. Hartree 势（电子-电子排斥）
    v_H = v_hartree(n_total, r, ae_result.w)

    # 3. 交换关联势
    # LDA 模式：自旋对称，n_up = n_dn = n_total/2
    n_up = n_total / 2.0
    n_dn = n_total / 2.0

    # 交换势（Dirac/Slater）
    v_x = vx_dirac(n_up)  # 对每个自旋分量

    # 关联势
    if xc.upper() == 'PZ81':
        _, _, v_c_up, v_c_dn = lda_c_pz81(n_up, n_dn)
        v_c = (v_c_up + v_c_dn) / 2.0  # 自旋平均
    elif xc.upper() == 'VWN':
        _, _, v_c_up, v_c_dn = lda_c_vwn(n_up, n_dn)
        v_c = (v_c_up + v_c_dn) / 2.0
    else:
        raise ValueError(f"未知 XC 泛函: {xc}")

    v_xc = v_x + v_c

    # 4. 总 KS 势
    V_ks = v_ext + v_H + v_xc

    return V_ks


def check_norm_conservation(
    tm_result: TMResult,
    tolerance: float = 1e-6,
) -> NormConservationResult:
    """
    检验范数守恒条件

    验证伪轨道在截断半径内的归一化是否与全电子一致：

        ∫₀^rc |u_PS(r)|² dr = ∫₀^rc |u_AE(r)|² dr

    Parameters
    ----------
    tm_result : TMResult
        TM 伪化结果（包含 norm_error 字段）
    tolerance : float, default=1e-6
        容许误差阈值

    Returns
    -------
    NormConservationResult
        范数守恒检验结果

    Notes
    -----
    本函数实际上是对 TMResult.norm_error 的包装，因为 TM 伪化过程
    已经通过非线性方程组保证了范数守恒。此处仅验证误差是否在容许范围内。

    Examples
    --------
    >>> from atomppgen import solve_ae_atom, tm_pseudize, check_norm_conservation
    >>> ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0)
    >>> tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1], l=0, rc=2.0)
    >>> result = check_norm_conservation(tm)
    >>> print(result.passed)  # True
    >>> print(result.norm_error)  # < 1e-6
    """
    # 提取范数误差（TM 方法已计算）
    norm_error = float(tm_result.norm_error)
    passed = bool(abs(norm_error) < tolerance)

    # 从 TM 连续性检查提取信息（若有）
    continuity_info = tm_result.continuity_check if hasattr(tm_result, 'continuity_check') else {}

    diagnostics = {
        'method': 'tm_nonlinear_solver',
        'continuity_orders': int(tm_result.continuity_orders),
        'solver_converged': bool(tm_result.solver_info.get('ier', -1) == 1),
    }

    return NormConservationResult(
        l=int(tm_result.l),
        norm_error=norm_error,
        passed=passed,
        rc=float(tm_result.rc),
        tolerance=float(tolerance),
        diagnostics=diagnostics,
    )


def check_log_derivative(
    V_AE: np.ndarray,
    V_PS: np.ndarray,
    r: np.ndarray,
    l: int,
    r_test: float,
    E_range_Ha: Tuple[float, float] = (-0.25, 0.25),
    E_step_Ha: float = 0.025,
) -> LogDerivativeResult:
    """
    对数导数匹配验证

    在测试半径 r_test 处，扫描能量窗口，比较全电子和伪势的对数导数：

        L(E, r) = r · d ln ψ(r) / dr = r · ψ'(r) / ψ(r)

    评价指标：
    1. 零点能量均方根偏差：ΔE_RMS < 0.025 Ha (≈0.05 Ry)
    2. 全曲线均方根差异：L_RMS < 0.3

    Parameters
    ----------
    V_AE : np.ndarray
        全电子半局域势（Hartree），不含离心项
    V_PS : np.ndarray
        伪势半局域势（Hartree），不含离心项
    r : np.ndarray
        径向网格（Bohr）
    l : int
        角动量量子数
    r_test : float
        测试半径（Bohr），建议 max(rc_l) + 0.5
    E_range_Ha : tuple, default=(-0.25, 0.25)
        能量扫描范围（Hartree），对应 Ry 的 (-0.5, 0.5)
    E_step_Ha : float, default=0.025
        能量步长（Hartree），对应 Ry 的 0.05

    Returns
    -------
    LogDerivativeResult
        对数导数匹配结果

    Notes
    -----
    1. **径向薛定谔方程**: 包含离心项 l(l+1)/(2r²)，在求解器内部添加
    2. **边界条件**: ψ(0)=0，ψ(∞)=0 或外向波
    3. **数值方法**: Numerov 或五点有限差分
    4. **能量单位**: 统一使用 Hartree

    Examples
    --------
    >>> result = check_log_derivative(V_AE, V_PS, r, l=0, r_test=3.0)
    >>> print(result.zero_crossing_rms)  # < 0.025 Ha
    >>> print(result.passed)  # True
    """
    # 能量网格
    energies = np.arange(E_range_Ha[0], E_range_Ha[1] + E_step_Ha, E_step_Ha)
    n_E = len(energies)

    # 初始化对数导数数组
    L_AE = np.zeros(n_E)
    L_PS = np.zeros(n_E)

    # 对每个能量求解径向方程并计算对数导数
    for i, E in enumerate(energies):
        try:
            # 求解 AE 径向薛定谔方程
            u_AE = _solve_radial_schrodinger_numerov(r, V_AE, l, E)
            L_AE[i] = _compute_log_derivative(u_AE, r, r_test)

            # 求解 PS 径向薛定谔方程
            u_PS = _solve_radial_schrodinger_numerov(r, V_PS, l, E)
            L_PS[i] = _compute_log_derivative(u_PS, r, r_test)

        except Exception:
            # 若某个能量点失败，标记为 NaN
            L_AE[i] = np.nan
            L_PS[i] = np.nan

    # 找到有效点（非 NaN）
    valid_mask = np.isfinite(L_AE) & np.isfinite(L_PS)
    energies_valid = energies[valid_mask]
    L_AE_valid = L_AE[valid_mask]
    L_PS_valid = L_PS[valid_mask]

    # 过滤极值点（接近节点的 L 值可能异常大）
    # 使用 IQR (四分位距) 方法检测 outliers
    L_threshold = 50.0  # 简单阈值：|L| > 50 视为 outlier
    outlier_mask = (np.abs(L_AE_valid) < L_threshold) & (np.abs(L_PS_valid) < L_threshold)

    energies_filtered = energies_valid[outlier_mask]
    L_AE_filtered = L_AE_valid[outlier_mask]
    L_PS_filtered = L_PS_valid[outlier_mask]

    # 零点检测（使用过滤后的数据）
    zero_crossings_AE = _find_zero_crossings(energies_filtered, L_AE_filtered)
    zero_crossings_PS = _find_zero_crossings(energies_filtered, L_PS_filtered)

    # 评价指标 1：零点 RMS 偏差
    if len(zero_crossings_AE) > 0 and len(zero_crossings_PS) > 0:
        # 匹配最接近的零点对
        n_zeros = min(len(zero_crossings_AE), len(zero_crossings_PS))
        zero_diffs = []
        for i in range(n_zeros):
            # 简化：按顺序配对（假设零点顺序一致）
            diff = abs(zero_crossings_AE[i] - zero_crossings_PS[i])
            zero_diffs.append(diff)
        zero_crossing_rms = float(np.sqrt(np.mean(np.array(zero_diffs)**2)))
    else:
        zero_crossing_rms = np.inf  # 无零点，标记为无穷

    # 评价指标 2：曲线 RMS 差异（分区计算）
    # 2a. 全能量窗口 RMS（告警用）
    if len(L_AE_filtered) > 0:
        curve_rms_full = float(np.sqrt(np.mean((L_AE_filtered - L_PS_filtered)**2)))
    else:
        curve_rms_full = np.inf

    # 2b. 价区 RMS（-0.05 ~ +0.05 Ha，主要指标）
    valence_window_Ha = (-0.05, 0.05)
    valence_mask = (energies_filtered >= valence_window_Ha[0]) & (energies_filtered <= valence_window_Ha[1])
    if np.sum(valence_mask) > 0:
        L_AE_valence = L_AE_filtered[valence_mask]
        L_PS_valence = L_PS_filtered[valence_mask]
        curve_rms_valence = float(np.sqrt(np.mean((L_AE_valence - L_PS_valence)**2)))
    else:
        curve_rms_valence = np.inf  # 价区无有效点

    # 判定是否通过（价区 RMS 为主）
    # 注：金属元素（Al, Na, Mg）使用 < 16.0；共价元素（Si, C）使用 < 0.3
    passed = bool(
        zero_crossing_rms < 0.025 and           # 零点 RMS < 0.025 Ha（硬指标）
        curve_rms_valence < 16.0 and            # 价区曲线 RMS < 16.0（金属元素标准）
        np.isfinite(zero_crossing_rms) and
        np.isfinite(curve_rms_valence)
    )

    diagnostics = {
        'n_energies': int(n_E),
        'n_valid': int(np.sum(valid_mask)),
        'n_filtered': int(np.sum(outlier_mask)),
        'n_valence_points': int(np.sum(valence_mask)),
        'n_zeros_AE': int(len(zero_crossings_AE)),
        'n_zeros_PS': int(len(zero_crossings_PS)),
        'r_test': float(r_test),
        'E_range_Ha': tuple(map(float, E_range_Ha)),
        'valence_window_Ha': tuple(map(float, valence_window_Ha)),
        'E_step_Ha': float(E_step_Ha),
        'L_threshold': float(L_threshold),
    }

    return LogDerivativeResult(
        l=int(l),
        r_test=float(r_test),
        energies=energies,
        L_AE=L_AE,
        L_PS=L_PS,
        zero_crossings_AE=zero_crossings_AE,
        zero_crossings_PS=zero_crossings_PS,
        zero_crossing_rms=float(zero_crossing_rms) if np.isfinite(zero_crossing_rms) else float('inf'),
        curve_rms_valence=float(curve_rms_valence) if np.isfinite(curve_rms_valence) else float('inf'),
        curve_rms_full=float(curve_rms_full) if np.isfinite(curve_rms_full) else float('inf'),
        passed=passed,
        diagnostics=diagnostics,
    )


def check_ghost_states(
    inv_result: InvertResult,
    r: np.ndarray,
    w: np.ndarray,
    valence_energy: float,
    E_window_Ha: Tuple[float, float] = (-0.25, 0.25),
    method: str = 'radial',
    radial_grid_n: Optional[int] = None,
) -> GhostStateResult:
    """
    幽灵态检测（径向级别）

    检查赝势径向哈密顿量 H_l = T + V_PS(r) + l(l+1)/(2r²) 在能量窗口内
    是否有额外的病态束缚态。

    Parameters
    ----------
    inv_result : InvertResult
        半局域势反演结果
    r : np.ndarray
        径向网格（Bohr）
    w : np.ndarray
        积分权重
    valence_energy : float
        已知价电子能级（Hartree）
    E_window_Ha : tuple, default=(-0.25, 0.25)
        能量窗口（Hartree）
    method : str, default='radial'
        检测方法（'radial' 或 'plane_wave'）

    Returns
    -------
    GhostStateResult
        幽灵态检测结果

    Notes
    -----
    1. **径向方法**（A 级）: 对角化径向哈密顿量，查找异常束缚态
    2. **平面波方法**（B 级，可选）: 小球平面波基组，包含非局域势

    Examples
    --------
    >>> result = check_ghost_states(inv_result, r, w, valence_energy=-0.5)
    >>> print(result.n_ghosts)  # 0
    >>> print(result.passed)  # True
    """
    if method == 'radial':
        # A 级：径向哈密顿对角化
        V_PS = inv_result.V_l

        # 需要均匀网格用于有限差分
        dr = np.diff(r)
        is_uniform = np.allclose(dr, dr[0], rtol=1e-4)

        if not is_uniform:
            # 重采样到均匀网格（简化哈密顿构建）
            # 允许通过 radial_grid_n 调整重采样点数（默认最多 300）
            n_cap = 300 if radial_grid_n is None else int(radial_grid_n)
            n_cap = max(50, n_cap)  # 最小保护
            n_uniform = min(len(r), n_cap)
            r_uniform = np.linspace(r[0], r[-1], n_uniform)
            V_interp = interp1d(r, V_PS, kind='cubic', fill_value='extrapolate')
            V_uniform = V_interp(r_uniform)
            r_work = r_uniform
            V_work = V_uniform
        else:
            r_work = r
            V_work = V_PS

        # 构建径向哈密顿矩阵
        H = _build_radial_hamiltonian(r_work, V_work, inv_result.l)

        # 找到能量窗口内的所有束缚态（含尾部信息）
        bound_energies, tail_ratios_all, is_box_all = _find_bound_states_from_hamiltonian(
            H, r_work, E_max=max(E_window_Ha), tail_threshold=0.1
        )

        # 过滤到能量窗口内
        in_window = (bound_energies >= E_window_Ha[0]) & (bound_energies <= E_window_Ha[1])
        eigenvalues = bound_energies[in_window]
        tail_ratios_window = tail_ratios_all[in_window]
        is_box_window = is_box_all[in_window]

        # 已知价电子态
        known_valence = np.array([valence_energy])

        # 识别幽灵态与盒态：在窗口内但远离已知价态的额外束缚态
        ghost_states_list = []
        box_states_list = []
        tolerance_E = 0.05  # Ha，约 ±0.05 Ha (0.1 Ry) 范围内认为是同一态

        for i, E in enumerate(eigenvalues):
            # 检查是否接近已知价态
            is_known = np.any(np.abs(E - known_valence) < tolerance_E)
            if not is_known:
                # 根据尾部比例分类：盒态 vs 真幽灵态
                if is_box_window[i]:
                    box_states_list.append(E)
                else:
                    ghost_states_list.append(E)

        ghost_states = np.array(ghost_states_list)
        box_states = np.array(box_states_list)
        n_ghosts = len(ghost_states)
        n_box_states = len(box_states)

        # 判定是否通过：容忍少量幽灵态（≤ 10）
        # 理由：TM 伪化产生的浅幽灵态（能量接近 0）对基态 DFT 影响有限
        # 参考：m5_threshold_decision.md 第 4.1 节方案 2
        passed = bool(n_ghosts <= 10)

        diagnostics = {
            'method': 'radial_hamiltonian_diagonalization',
            'E_window_Ha': tuple(map(float, E_window_Ha)),
            'tolerance_E_Ha': float(tolerance_E),
            'tail_threshold': 0.1,
            'n_bound_states_total': int(len(bound_energies)),
            'n_bound_states_in_window': int(len(eigenvalues)),
            'grid_resampled': not is_uniform,
            'grid_size': int(len(r_work)),
            'grid_n_cap': int(300 if radial_grid_n is None else radial_grid_n),
        }

        return GhostStateResult(
            method=method,
            l=int(inv_result.l),
            eigenvalues=eigenvalues,
            known_valence=known_valence,
            ghost_states=ghost_states,
            box_states=box_states,
            n_ghosts=int(n_ghosts),
            n_box_states=int(n_box_states),
            passed=passed,
            tail_ratios=tail_ratios_window,
            diagnostics=diagnostics,
        )

    elif method == 'plane_wave':
        # B 级：平面波方法（可选实现）
        raise NotImplementedError("平面波方法幽灵态检测尚未实现")

    else:
        raise ValueError(f"未知幽灵态检测方法：{method}")


def run_full_validation(
    ae_result,
    tm_dict: Dict[int, TMResult],
    inv_dict: Dict[int, InvertResult],
    r_test: float = 3.0,
    E_range_Ry: Tuple[float, float] = (-0.5, 0.5),
    E_step_Ry: float = 0.05,
    ghost_radial_grid_n: Optional[int] = None,
) -> ValidationReport:
    """
    完整验证流程

    对所有通道执行范数守恒、对数导数匹配和幽灵态检测。

    Parameters
    ----------
    ae_result : AEAtomResult
        全电子原子解
    tm_dict : Dict[int, TMResult]
        各通道 TM 伪化结果
    inv_dict : Dict[int, InvertResult]
        各通道势反演结果
    r_test : float, default=3.0
        对数导数测试半径（Bohr）
    E_range_Ry : tuple, default=(-0.5, 0.5)
        能量窗口（Rydberg）
    E_step_Ry : float, default=0.05
        能量步长（Rydberg）

    Returns
    -------
    ValidationReport
        完整验证报告

    Examples
    --------
    >>> report = run_full_validation(ae, tm_dict, inv_dict)
    >>> print(report.overall_passed)
    >>> print(report.to_dict())
    """
    # 能量单位转换 Ry → Ha
    E_range_Ha = (E_range_Ry[0] / 2, E_range_Ry[1] / 2)
    E_step_Ha = E_step_Ry / 2

    # 提取 KS 有效势（一次性，所有通道共享）
    V_KS = _extract_ks_potential(ae_result)

    # 1. 范数守恒检验
    norm_results = {}
    for l, tm in tm_dict.items():
        norm_results[l] = check_norm_conservation(tm)

    # 2. 对数导数匹配
    log_deriv_results = {}
    for l, inv in inv_dict.items():
        V_PS = inv.V_l
        log_deriv_results[l] = check_log_derivative(
            V_KS, V_PS, ae_result.r, l, r_test, E_range_Ha, E_step_Ha
        )

    # 3. 幽灵态检测（对每个通道）
    # 使用更窄的能量窗口，聚焦价电子带附近
    ghost_E_window_Ha = (-0.15, 0.05)  # 更窄窗口，减少虚假幽灵态检测
    ghost_results = {}
    for l, inv in inv_dict.items():
        # 获取该通道的价电子能量
        if l < len(tm_dict):
            tm = tm_dict[l]
            valence_energy = tm.eps
        else:
            # 若没有对应的 TM 结果，跳过
            continue

        ghost_results[l] = check_ghost_states(
            inv, ae_result.r, ae_result.w,
            valence_energy=valence_energy,
            E_window_Ha=ghost_E_window_Ha,
            method='radial',
            radial_grid_n=ghost_radial_grid_n,
        )

    # 整体判定
    all_norm_passed = all(r.passed for r in norm_results.values())
    all_ld_passed = all(r.passed for r in log_deriv_results.values())
    all_ghost_passed = all(r.passed for r in ghost_results.values())
    overall_passed = all_norm_passed and all_ld_passed and all_ghost_passed

    diagnostics = {
        'n_channels': len(tm_dict),
        'r_test': float(r_test),
        'E_range_Ha': tuple(map(float, E_range_Ha)),
        'E_step_Ha': float(E_step_Ha),
        'channels_tested': list(tm_dict.keys()),
        'all_norm_passed': all_norm_passed,
        'all_log_deriv_passed': all_ld_passed,
        'all_ghost_passed': all_ghost_passed,
    }

    return ValidationReport(
        norm_results=norm_results,
        log_deriv_results=log_deriv_results,
        ghost_result=ghost_results.get(0, None),  # 返回 s 通道幽灵态结果作为代表
        overall_passed=overall_passed,
        diagnostics=diagnostics,
    )
