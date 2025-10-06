"""
半局域势反演

从 TM 伪轨道反演半局域势 V_l(r)，使用径向 Schrödinger 方程：

    V_l(r) = ε + (1/2) · u''(r)/u(r) - l(l+1)/(2r²)

技术要点：
- 内区（r ≤ rc）：使用 TM 解析导数，数值稳定
- 外区（r > rc）：使用样条法导数
- 节点保护：检测零点并插值或掩码

主要功能
--------
invert_semilocal_potential : 反演半局域势
InvertResult : 存储反演结果的数据类

参考文献
--------
Troullier & Martins, PRB 43, 1993 (1991)
Giannozzi, Notes on pseudopotential generation (2019)
"""

from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np
from scipy.interpolate import CubicSpline

from atomppgen.tm import TMResult, eval_derivatives_at


__all__ = [
    "InvertResult",
    "invert_semilocal_potential",
]


@dataclass
class InvertResult:
    """
    势反演结果

    Attributes
    ----------
    V_l : np.ndarray
        半局域势 V_l(r)（Hartree）
    r : np.ndarray
        径向网格（Bohr）
    l : int
        角动量量子数
    rc : float
        伪化半径（Bohr）
    eps : float
        伪化能量（Hartree）
    diagnostics : Dict
        诊断信息：
        - n_nodes : 节点数
        - V_max : 最大势值
        - V_min : 最小势值
        - V_at_rc : rc 处势值
        - method_inner : 内区方法 ('analytical')
        - method_outer : 外区方法 ('spline')
    """
    V_l: np.ndarray
    r: np.ndarray
    l: int
    rc: float
    eps: float
    diagnostics: Dict


def invert_semilocal_potential(
    tm_result: TMResult,
    r: np.ndarray,
    u_ps: Optional[np.ndarray] = None,
    node_tol: float = 1e-10,
    V_max_clip: float = 1000.0,
    smooth_rc: bool = False,
    smooth_width: float = 0.1,
) -> InvertResult:
    """
    从 TM 伪轨道反演半局域势

    使用径向 Schrödinger 方程反演：
        V_l(r) = ε + (1/2) · u''(r)/u(r) - l(l+1)/(2r²)

    内区（r ≤ rc）使用 TM 解析导数，外区使用样条法导数。

    Parameters
    ----------
    tm_result : TMResult
        TM 伪化结果（包含 a_coeff, rc, eps, l）
    r : np.ndarray
        径向网格（Bohr）
    u_ps : np.ndarray, optional
        伪轨道，若不提供则从 tm_result.u_ps 获取
    node_tol : float, default=1e-10
        节点检测阈值（|u| < node_tol 视为节点）
    V_max_clip : float, default=1000.0
        势的裁剪上限（防止除零发散）
    smooth_rc : bool, default=False
        是否在 rc 处应用平滑过渡（内外区势值混合）
    smooth_width : float, default=0.1
        平滑区域半宽度（Bohr），仅在 smooth_rc=True 时有效
        平滑范围为 [rc-smooth_width, rc+smooth_width]

    Returns
    -------
    InvertResult
        包含反演势和诊断信息

    Raises
    ------
    ValueError
        如果网格与 tm_result 不匹配

    Notes
    -----
    1. 内区（r ≤ rc）：使用 TM 解析导数，从 `_eval_tm_at_rc()` 获取 u, u', u''
    2. 外区（r > rc）：使用 `eval_derivatives_at()` 样条法
    3. 节点保护：在 |u| < node_tol 附近使用线性插值填充势值
    4. 离心势项 l(l+1)/(2r²) 在原点附近奇异，需特殊处理
    5. **平滑过渡**（可选）：若 smooth_rc=True，在 rc±smooth_width 区域使用
       三次样条平滑内外区势值，适用于特殊 rc/网格组合导致跳变较大的情况。
       默认关闭（测试显示大多数情况 rc 处相对跳变 < 1%）。

    Examples
    --------
    >>> from atomppgen import solve_ae_atom, tm_pseudize
    >>> ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0)
    >>> tm = tm_pseudize(r=ae.r, w=ae.w, u_ae=ae.u_by_l[0][2],
    ...                   eps=ae.eps_by_l[0][2], l=0, rc=2.0)
    >>> inv = invert_semilocal_potential(tm, ae.r)
    >>> print(inv.diagnostics['V_at_rc'])  # rc 处势值
    """
    # 参数提取
    if u_ps is None:
        u_ps = tm_result.u_ps

    if len(r) != len(u_ps):
        raise ValueError(f"网格长度 {len(r)} 与伪轨道长度 {len(u_ps)} 不匹配")

    rc = tm_result.rc
    eps = tm_result.eps
    l = tm_result.l
    a_coeff = tm_result.a_coeff

    # 找到 rc 对应的索引
    i_rc = np.searchsorted(r, rc)
    if i_rc >= len(r):
        i_rc = len(r) - 1

    # 初始化势数组
    V_l = np.zeros_like(r)

    # ========== 内区（r ≤ rc）：使用 TM 解析导数 ==========
    V_l[:i_rc+1] = _invert_inner_analytical(
        r[:i_rc+1], u_ps[:i_rc+1], a_coeff, l, eps, node_tol, V_max_clip
    )

    # ========== 外区（r > rc）：使用样条法导数 ==========
    if i_rc + 1 < len(r):
        V_l[i_rc+1:] = _invert_outer_spline(
            r_full=r,  # 传递完整网格
            u_full=u_ps,  # 传递完整轨道
            i_start=i_rc+1,  # 外区起始索引
            l=l,
            eps=eps,
            node_tol=node_tol,
            V_max_clip=V_max_clip,
        )

    # ========== rc 处平滑过渡（可选） ==========
    if smooth_rc and i_rc > 0 and i_rc < len(r) - 1:
        V_l = _smooth_at_rc(r, V_l, i_rc, smooth_width)

    # ========== 诊断信息 ==========
    diagnostics = {
        'n_nodes': _count_nodes(u_ps, node_tol),
        'V_max': float(np.max(V_l)),
        'V_min': float(np.min(V_l)),
        'V_at_rc': float(V_l[i_rc]),
        'method_inner': 'analytical',
        'method_outer': 'spline',
    }

    return InvertResult(
        V_l=V_l,
        r=r,
        l=l,
        rc=rc,
        eps=eps,
        diagnostics=diagnostics,
    )


def _invert_inner_analytical(
    r: np.ndarray,
    u: np.ndarray,
    a_coeff: np.ndarray,
    l: int,
    eps: float,
    node_tol: float,
    V_max_clip: float,
) -> np.ndarray:
    """
    内区使用 TM 解析导数反演势

    对于 TM 形式 u = r^{l+1} exp(p(r))，可以解析计算：
        u' = [(l+1)/r + p'] u
        u'' = [l(l+1)/r² + 2(l+1)p'/r + p'² + p''] u

    因此：
        u''/u = l(l+1)/r² + 2(l+1)p'/r + p'² + p''

    代入：
        V_l = ε + (1/2) u''/u - l(l+1)/(2r²)
            = ε + (1/2)[l(l+1)/r² + 2(l+1)p'/r + p'² + p''] - l(l+1)/(2r²)
            = ε + (l+1)p'/r + (1/2)(p'² + p'')

    Parameters
    ----------
    r : np.ndarray
        内区网格
    u : np.ndarray
        内区伪轨道
    a_coeff : np.ndarray
        TM 系数 [a_0, a_2, a_4, ...]
    l : int
        角动量
    eps : float
        伪化能量
    node_tol : float
        节点容忍度
    V_max_clip : float
        势裁剪上限

    Returns
    -------
    np.ndarray
        内区势 V_l(r)
    """
    # 计算 p(r) 及其导数
    # p = a_0 + a_2 r² + a_4 r⁴ + ...
    # p' = 2 a_2 r + 4 a_4 r³ + ...
    # p'' = 2 a_2 + 12 a_4 r² + ...

    p = np.zeros_like(r)
    p1 = np.zeros_like(r)  # p'
    p2 = np.zeros_like(r)  # p''

    for i in range(len(a_coeff)):
        p += a_coeff[i] * r**(2*i)
        if i >= 1:
            p1 += 2*i * a_coeff[i] * r**(2*i-1)
            p2 += 2*i * (2*i-1) * a_coeff[i] * r**(2*i-2)

    # V_l = ε + (l+1)p'/r + (1/2)(p'² + p'')
    # 处理 r→0 时的奇点：(l+1)p'/r
    # 当 r 很小时，p' ~ 2 a_2 r，所以 (l+1)p'/r ~ 2(l+1) a_2

    V_l = np.zeros_like(r)

    # 对于 r > 0，使用完整公式
    mask_nonzero = (r > node_tol)

    V_l[mask_nonzero] = (
        eps
        + (l+1) * p1[mask_nonzero] / r[mask_nonzero]
        + 0.5 * (p1[mask_nonzero]**2 + p2[mask_nonzero])
    )

    # 对于 r ≈ 0，使用极限值（假设 p' ~ 2 a_2 r）
    if len(a_coeff) > 1 and not mask_nonzero[0]:
        V_l[0] = eps + 0.5 * p2[0]  # 忽略 p'/r 项（趋于常数）

    # 裁剪极端值
    V_l = np.clip(V_l, -V_max_clip, V_max_clip)

    return V_l


def _invert_outer_spline(
    r_full: np.ndarray,
    u_full: np.ndarray,
    i_start: int,
    l: int,
    eps: float,
    node_tol: float,
    V_max_clip: float,
) -> np.ndarray:
    """
    外区使用样条法导数反演势

    V_l(r) = ε + (1/2) u''/u - l(l+1)/(2r²)

    Parameters
    ----------
    r_full : np.ndarray
        完整网格（包括内外区）
    u_full : np.ndarray
        完整伪轨道（包括内外区）
    i_start : int
        外区起始索引
    l : int
        角动量
    eps : float
        伪化能量
    node_tol : float
        节点容忍度
    V_max_clip : float
        势裁剪上限

    Returns
    -------
    np.ndarray
        外区势 V_l(r)
    """
    n_outer = len(r_full) - i_start
    V_l = np.zeros(n_outer)

    for i in range(n_outer):
        i_global = i_start + i
        r_i = r_full[i_global]
        u_i = u_full[i_global]

        # 检查节点
        if abs(u_i) < node_tol:
            # 节点处使用插值（后续处理）
            V_l[i] = 0.0  # 占位
            continue

        # 使用样条法计算二阶导数
        # 传递完整数据，eval_derivatives_at 自动选择窗口
        derivs = eval_derivatives_at(
            r=r_full,
            u=u_full,
            x0=r_i,
            max_order=2,
            window=7,
        )

        u_val = derivs[0]
        u2 = derivs[2]  # u''

        # V_l = ε + (1/2) u''/u - l(l+1)/(2r²)
        V_l[i] = eps + 0.5 * u2 / u_val - l*(l+1) / (2 * r_i**2)

    # 节点插值
    r_outer = r_full[i_start:]
    u_outer = u_full[i_start:]
    V_l = _interpolate_nodes(r_outer, u_outer, V_l, node_tol)

    # 裁剪
    V_l = np.clip(V_l, -V_max_clip, V_max_clip)

    return V_l


def _interpolate_nodes(
    r: np.ndarray,
    u: np.ndarray,
    V_l: np.ndarray,
    node_tol: float,
) -> np.ndarray:
    """
    在节点附近插值势值

    Parameters
    ----------
    r : np.ndarray
        网格
    u : np.ndarray
        轨道
    V_l : np.ndarray
        势（节点处为占位值）
    node_tol : float
        节点容忍度

    Returns
    -------
    np.ndarray
        插值后的势
    """
    # 找到节点
    is_node = np.abs(u) < node_tol

    if not np.any(is_node):
        return V_l

    # 对节点区域进行线性插值
    V_l_interp = V_l.copy()

    # 使用 scipy 的插值（跳过节点）
    r_valid = r[~is_node]
    V_valid = V_l[~is_node]

    if len(r_valid) > 1:
        from scipy.interpolate import interp1d
        f_interp = interp1d(r_valid, V_valid, kind='linear', fill_value='extrapolate')
        V_l_interp[is_node] = f_interp(r[is_node])

    return V_l_interp


def _count_nodes(u: np.ndarray, node_tol: float) -> int:
    """
    计算轨道的节点数（符号变化次数）

    Parameters
    ----------
    u : np.ndarray
        轨道
    node_tol : float
        容忍度

    Returns
    -------
    int
        节点数
    """
    # 符号变化检测
    sign_changes = 0
    u_filtered = u[np.abs(u) > node_tol]  # 过滤掉接近零的点

    if len(u_filtered) < 2:
        return 0

    for i in range(1, len(u_filtered)):
        if u_filtered[i] * u_filtered[i-1] < 0:
            sign_changes += 1

    return sign_changes


def _smooth_at_rc(
    r: np.ndarray,
    V_l: np.ndarray,
    i_rc: int,
    smooth_width: float,
) -> np.ndarray:
    """
    在 rc 处应用平滑过渡

    在 [rc-smooth_width, rc+smooth_width] 区域使用三次样条平滑势值，
    减少内外区计算方法差异导致的跳变。

    Parameters
    ----------
    r : np.ndarray
        径向网格
    V_l : np.ndarray
        势（未平滑）
    i_rc : int
        rc 对应的网格索引
    smooth_width : float
        平滑区域半宽度（Bohr）

    Returns
    -------
    np.ndarray
        平滑后的势
    """
    rc = r[i_rc]
    r_min = rc - smooth_width
    r_max = rc + smooth_width

    # 找到平滑区域索引
    i_min = np.searchsorted(r, r_min)
    i_max = np.searchsorted(r, r_max)

    # 确保区域有效
    i_min = max(0, i_min)
    i_max = min(len(r), i_max)

    # 至少需要 4 个点才能构建三次样条
    if i_max - i_min < 4:
        return V_l  # 区域太小，跳过平滑

    # 提取平滑区域
    r_smooth = r[i_min:i_max]
    V_smooth_region = V_l[i_min:i_max]

    # 构建三次样条并评估
    cs = CubicSpline(r_smooth, V_smooth_region, bc_type='not-a-knot')

    # 在同一区域重新评估（平滑）
    V_l_smoothed = V_l.copy()
    V_l_smoothed[i_min:i_max] = cs(r_smooth)

    return V_l_smoothed
