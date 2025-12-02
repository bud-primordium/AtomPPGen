"""
Troullier-Martins (TM) 伪化器

实现 TM 方法生成赝轨道，满足范数守恒和连续性约束。

主要功能
--------
tm_pseudize : 对给定 AE 轨道进行 TM 伪化
TMResult : 存储伪化结果的数据类

参考文献
--------
Troullier & Martins, PRB 43, 1993 (1991)

算法文档
--------
详见 docs/source/algorithm/tm_method.rst
"""

from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np
from scipy.optimize import fsolve, least_squares
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline


__all__ = [
    "TMResult",
    "tm_pseudize",
]


@dataclass
class TMResult:
    """
    TM 伪化结果

    Attributes
    ----------
    u_ps : np.ndarray
        伪轨道 u(r)（径向波函数，已拼接内外区）
    a_coeff : np.ndarray
        TM 多项式系数 [a_0, a_2, a_4, ...]
        内区形式：u = r^{l+1} exp(Σ a_{2i} r^{2i})
    rc : float
        伪化半径（Bohr）
    eps : float
        伪化能量（Hartree）
    l : int
        角动量量子数
    norm_error : float
        范数守恒相对误差 |Q_PS - Q_AE| / Q_AE
    continuity_orders : int
        连续性阶数（匹配到几阶导数，2 或 4）
    continuity_check : Dict
        rc 处的连续性检查结果 {'u': {...}, 'du': {...}, 'd2u': {...}}
    solver_info : Dict
        求解器收敛信息 {'ier': int, 'nfev': int, 'mesg': str, 'fallback': bool}
    """
    u_ps: np.ndarray
    a_coeff: np.ndarray
    rc: float
    eps: float
    l: int
    norm_error: float
    continuity_orders: int
    continuity_check: Dict
    solver_info: Dict


def tm_pseudize(
    r: np.ndarray,
    w: np.ndarray,
    u_ae: np.ndarray,
    eps: float,
    l: int,
    rc: float,
    continuity_orders: int = 2,
) -> TMResult:
    """
    使用 Troullier-Martins 方法构造伪轨道

    在 r ≤ rc 内区，伪轨道形式为：
        u(r) = r^{l+1} · exp(Σ_{i=0}^{N} a_{2i} r^{2i})

    约束条件：
    1. 在 r = rc 处，u 及其导数（到 continuity_orders 阶）与 AE 轨道连续
    2. 内区范数守恒：∫_0^{rc} |u|^2 dr = ∫_0^{rc} |u_ae|^2 dr

    Parameters
    ----------
    r : np.ndarray
        径向网格（Bohr）
    w : np.ndarray
        积分权重
    u_ae : np.ndarray
        全电子径向波函数 u(r)（已归一化）
    eps : float
        伪化能量（Hartree）
    l : int
        角动量量子数
    rc : float
        伪化半径（Bohr）
    continuity_orders : int, optional
        导数连续性阶数，默认 2（匹配 u, u', u'' + 范数）
        可选 4（匹配到 u'''' + 范数）

    Returns
    -------
    TMResult
        包含伪轨道、系数、范数误差等信息

    Raises
    ------
    ValueError
        如果 rc 不在网格范围内
        如果 continuity_orders 不是 2 或 4

    Notes
    -----
    TM 方法关键点：
    1. 系数数量由连续性条件决定：
       - continuity_orders = 2: 需要 4 个系数 [a_0, a_2, a_4, a_6]
       - continuity_orders = 4: 需要 6 个系数 [a_0, ..., a_10]
    2. 数值稳定性：使用对数形式避免指数溢出
    3. 范数守恒通过非线性方程组保证（scipy.optimize.fsolve）
    4. **网格推荐**：强烈建议使用 `exp_transformed` 或 `log` 网格类型。
       `linear` 网格的导数精度不足，可能导致 TM 非线性方程组不收敛。

    References
    ----------
    Troullier & Martins, PRB 43, 1993 (1991), 方程 (14)-(18)
    """
    # 验证输入
    if rc <= r[0] or rc > r[-1]:
        raise ValueError(f"rc={rc:.3f} 超出网格范围 [{r[0]:.3f}, {r[-1]:.3f}]")
    if continuity_orders not in (2, 4):
        raise ValueError(f"continuity_orders 必须是 2 或 4，当前为 {continuity_orders}")

    # 找到 rc 对应的网格索引
    i_rc = np.searchsorted(r, rc)
    if i_rc >= len(r):
        i_rc = len(r) - 1

    # 确定系数数量
    # continuity_orders=2: 3个导数约束(u, u', u'') + 1个范数 = 4个方程 → 4个系数
    # continuity_orders=4: 5个导数约束(u, ..., u'''') + 1个范数 = 6个方程 → 6个系数
    n_coeffs = continuity_orders + 2

    # 计算 AE 轨道在 rc 处的导数（使用样条法，适用于非均匀网格）
    derivs_ae = eval_derivatives_at(r, u_ae, rc, max_order=continuity_orders)

    # 获取 rc 处的符号，用于处理存在节点的价层轨道
    sign_prefactor = np.sign(derivs_ae[0])
    if sign_prefactor == 0:
        sign_prefactor = np.sign(u_ae[i_rc])
    if sign_prefactor == 0:
        sign_prefactor = 1.0
    sign_prefactor = float(sign_prefactor)

    # 计算 AE 轨道的内区范数
    norm_ae = _compute_norm(r, w, u_ae, i_rc)

    # 求解 TM 系数
    a_coeff, solver_info = _solve_tm_coefficients(
        rc=rc,
        l=l,
        derivs_ae=derivs_ae,
        norm_ae=norm_ae,
        n_coeffs=n_coeffs,
        continuity_orders=continuity_orders,
        r_inner=r[:i_rc+1],
        w_inner=w[:i_rc+1],
        sign_prefactor=sign_prefactor,
    )

    # 拼接内外区生成完整伪轨道
    u_ps = _splice_orbital(r, u_ae, a_coeff, l, i_rc, sign_prefactor)

    # 计算范数守恒误差
    norm_ps = _compute_norm(r, w, u_ps, i_rc)
    norm_error = abs(norm_ps - norm_ae) / norm_ae

    # 检查连续性（使用解析导数）
    continuity_check = _check_continuity(
        rc=rc,
        l=l,
        a_coeff=a_coeff,
        derivs_ae=derivs_ae,
        order=continuity_orders,
        sign_prefactor=sign_prefactor,
    )

    return TMResult(
        u_ps=u_ps,
        a_coeff=a_coeff,
        rc=rc,
        eps=eps,
        l=l,
        norm_error=norm_error,
        continuity_orders=continuity_orders,
        continuity_check=continuity_check,
        solver_info=solver_info,
    )


def eval_derivatives_at(
    r: np.ndarray,
    u: np.ndarray,
    x0: float,
    max_order: int = 4,
    window: int = 7,
    bc: str = "not-a-knot",
) -> np.ndarray:
    """
    在指定点 x0 处计算导数（使用局部三次样条）

    此方法适用于非均匀网格，通过在 x0 邻域构建三次样条插值，
    在 x0 处直接评估各阶导数。

    Parameters
    ----------
    r : np.ndarray
        径向网格（可以是非均匀的）
    u : np.ndarray
        轨道函数
    x0 : float
        评估点（例如 rc）
    max_order : int, default=4
        最高导数阶数（2 或 4）
    window : int, default=7
        窗口大小（以 x0 为中心取点数），推荐 7（±3 点）
    bc : str, default="not-a-knot"
        边界条件类型，可选 "not-a-knot", "natural", "clamped"

    Returns
    -------
    np.ndarray
        [u, u', u'', u''', u''''] (根据 max_order 返回相应长度)

    Notes
    -----
    - 对于非均匀网格（如 exp_transformed, log），此方法比有限差分更精确
    - 窗口不足时自动退化到 5 点或 3 点
    - 靠近边界时允许不对称窗口

    References
    ----------
    - scipy.interpolate.CubicSpline
    - Press et al., *Numerical Recipes*, Section 3.3: Cubic Spline Interpolation
    """
    # 找到 x0 的索引
    i_center = np.searchsorted(r, x0)
    if i_center >= len(r):
        i_center = len(r) - 1

    # 确定窗口范围（以 i_center 为中心）
    half_win = window // 2
    i_left = max(0, i_center - half_win)
    i_right = min(len(r), i_center + half_win + 1)

    # 如果窗口太小，尝试扩展到至少 5 点
    if i_right - i_left < 5:
        if i_left == 0:
            i_right = min(len(r), 5)
        elif i_right == len(r):
            i_left = max(0, len(r) - 5)
        else:
            # 中间区域窗口太小，强制至少 5 点
            i_left = max(0, i_center - 2)
            i_right = min(len(r), i_center + 3)

    # 提取窗口
    r_win = r[i_left:i_right]
    u_win = u[i_left:i_right]

    # 构建三次样条
    try:
        cs = CubicSpline(r_win, u_win, bc_type=bc)
    except Exception:
        # 如果样条构建失败，回退到更简单的边界条件
        cs = CubicSpline(r_win, u_win, bc_type="natural")

    # 在 x0 处评估导数
    derivs = []

    # u(x0)
    derivs.append(float(cs(x0)))

    # u'(x0)
    derivs.append(float(cs.derivative(1)(x0)))

    # u''(x0)
    derivs.append(float(cs.derivative(2)(x0)))

    if max_order >= 4:
        # u'''(x0)
        derivs.append(float(cs.derivative(3)(x0)))

        # u''''(x0) - 三次样条的四阶导数为 0，这里用数值差分近似
        # 或者在窗口内用五次样条
        # 简化处理：用三阶导的差商
        dx = 1e-6
        d3_plus = cs.derivative(3)(x0 + dx)
        d3_minus = cs.derivative(3)(x0 - dx)
        u4 = (d3_plus - d3_minus) / (2 * dx)
        derivs.append(float(u4))

    return np.array(derivs)


def _eval_derivatives(
    r: np.ndarray,
    u: np.ndarray,
    i_rc: int,
    order: int = 2,
) -> np.ndarray:
    """
    计算轨道在 rc 处的函数值和导数（已废弃，使用 eval_derivatives_at 替代）

    此函数假设等距网格，对非均匀网格会产生系统误差。
    保留用于向后兼容，建议使用 eval_derivatives_at()。

    Parameters
    ----------
    r : np.ndarray
        径向网格
    u : np.ndarray
        轨道函数
    i_rc : int
        rc 对应的网格索引
    order : int
        最高导数阶数（2 或 4）

    Returns
    -------
    np.ndarray
        [u, u', u'', u''', u''''] (根据 order 返回相应长度)

    See Also
    --------
    eval_derivatives_at : 推荐使用的样条法导数计算
    """
    # 使用中心差分公式（5 点模板）
    # 确保不越界
    i = max(2, min(i_rc, len(r) - 3))
    h = r[i] - r[i-1]  # 假设网格接近等距

    derivs = []

    # u
    derivs.append(u[i])

    # u' (一阶导数，5 点)
    u1 = (-u[i+2] + 8*u[i+1] - 8*u[i-1] + u[i-2]) / (12 * h)
    derivs.append(u1)

    # u'' (二阶导数，5 点)
    u2 = (-u[i+2] + 16*u[i+1] - 30*u[i] + 16*u[i-1] - u[i-2]) / (12 * h**2)
    derivs.append(u2)

    if order >= 4:
        # u''' (三阶导数)
        u3 = (u[i+2] - 2*u[i+1] + 2*u[i-1] - u[i-2]) / (2 * h**3)
        derivs.append(u3)

        # u'''' (四阶导数)
        u4 = (u[i+2] - 4*u[i+1] + 6*u[i] - 4*u[i-1] + u[i-2]) / h**4
        derivs.append(u4)

    return np.array(derivs)


def _compute_norm(
    r: np.ndarray,
    w: np.ndarray,
    u: np.ndarray,
    i_rc: int,
) -> float:
    """
    计算内区范数 ∫_0^{rc} |u|^2 dr

    Parameters
    ----------
    r : np.ndarray
        径向网格
    w : np.ndarray
        积分权重
    u : np.ndarray
        轨道函数
    i_rc : int
        rc 对应的网格索引

    Returns
    -------
    float
        内区范数
    """
    return np.sum(u[:i_rc+1]**2 * w[:i_rc+1])


def _eval_tm_at_rc(rc: float, l: int, a: np.ndarray, sign_prefactor: float = 1.0) -> np.ndarray:
    """
    计算 TM 轨道及其导数在 rc 处的值

    u(r) = r^{l+1} exp(p(r))，其中 p(r) = Σ a_{2i} r^{2i}

    Parameters
    ----------
    rc : float
        伪化半径
    l : int
        角动量
    a : np.ndarray
        系数 [a_0, a_2, a_4, ...]

    sign_prefactor : float, optional
        AE 轨道在 rc 处的符号，用于恢复正确符号

    Returns
    -------
    np.ndarray
        [u, u', u'', u''', u''''] (根据 a 的长度返回对应导数)
    """
    r = rc
    N = len(a) - 1

    # 计算 p(r) 及其导数
    # p = a_0 + a_2 r^2 + a_4 r^4 + ...
    p = sum(a[i] * r**(2*i) for i in range(len(a)))

    # 防止指数溢出（裁剪到安全范围）
    p = np.clip(p, -700, 700)

    # p' = 2 a_2 r + 4 a_4 r^3 + ...
    p1 = sum(2*i * a[i] * r**(2*i-1) for i in range(1, len(a)))

    # p'' = 2 a_2 + 12 a_4 r^2 + ...
    p2 = sum(2*i * (2*i-1) * a[i] * r**(2*i-2) for i in range(1, len(a)))

    # p'''
    p3 = sum(2*i * (2*i-1) * (2*i-2) * a[i] * r**(2*i-3) for i in range(1, len(a)))

    # p''''
    p4 = sum(2*i * (2*i-1) * (2*i-2) * (2*i-3) * a[i] * r**(2*i-4) for i in range(1, len(a)))

    # u = r^{l+1} exp(p)
    r_pow = r**(l+1)
    exp_p = np.exp(p)
    u = r_pow * exp_p

    # u' = [(l+1) r^l + r^{l+1} p'] exp(p)
    u1 = ((l+1) * r**l + r_pow * p1) * exp_p

    # u'' = [l(l+1) r^{l-1} + 2(l+1) r^l p' + r^{l+1} (p'^2 + p'')] exp(p)
    u2 = (l*(l+1) * r**(l-1) + 2*(l+1) * r**l * p1 + r_pow * (p1**2 + p2)) * exp_p

    # u'''
    u3 = (
        l*(l+1)*(l-1) * r**(l-2)
        + 3*l*(l+1) * r**(l-1) * p1
        + 3*(l+1) * r**l * (p1**2 + p2)
        + r_pow * (p1**3 + 3*p1*p2 + p3)
    ) * exp_p

    # u''''
    u4 = (
        l*(l+1)*(l-1)*(l-2) * r**(l-3)
        + 4*l*(l+1)*(l-1) * r**(l-2) * p1
        + 6*l*(l+1) * r**(l-1) * (p1**2 + p2)
        + 4*(l+1) * r**l * (p1**3 + 3*p1*p2 + p3)
        + r_pow * (p1**4 + 6*p1**2*p2 + 4*p1*p3 + 3*p2**2 + p4)
    ) * exp_p

    derivs = np.array([u, u1, u2, u3, u4])
    return sign_prefactor * derivs


def _compute_tm_norm(
    r: np.ndarray,
    w: np.ndarray,
    l: int,
    a: np.ndarray,
) -> float:
    """
    计算 TM 轨道的内区范数

    ∫_0^{rc} |u|^2 dr，其中 u = r^{l+1} exp(Σ a_{2i} r^{2i})

    Parameters
    ----------
    r : np.ndarray
        内区网格
    w : np.ndarray
        积分权重
    l : int
        角动量
    a : np.ndarray
        系数

    Returns
    -------
    float
        内区范数
    """
    # 计算 p(r)
    p = np.zeros_like(r)
    for i in range(len(a)):
        p += a[i] * r**(2*i)

    # 防止指数溢出（裁剪到安全范围）
    p = np.clip(p, -700, 700)

    u = r**(l+1) * np.exp(p)
    return np.sum(u**2 * w)


def _solve_tm_coefficients(
    rc: float,
    l: int,
    derivs_ae: np.ndarray,
    norm_ae: float,
    n_coeffs: int,
    continuity_orders: int,
    r_inner: np.ndarray,
    w_inner: np.ndarray,
    sign_prefactor: float,
) -> tuple[np.ndarray, Dict]:
    """
    求解 TM 系数 a_{2i}

    通过非线性方程组求解，约束条件：
    1. rc 处函数值匹配
    2. rc 处各阶导数匹配
    3. 内区范数守恒

    Parameters
    ----------
    rc : float
        伪化半径
    l : int
        角动量
    derivs_ae : np.ndarray
        AE 轨道在 rc 处的导数 [u, u', u'', ...]
    norm_ae : float
        AE 内区范数
    n_coeffs : int
        系数数量
    continuity_orders : int
        连续性阶数
    r_inner, w_inner : np.ndarray
        内区网格和权重
    sign_prefactor : float
        rc 处 AE 轨道的符号，用于保证伪轨道在节点后的符号一致

    Returns
    -------
    a_solution : np.ndarray
        系数 [a_0, a_2, a_4, ...]
    solver_info : Dict
        求解器信息 {'ier': int, 'nfev': int, 'mesg': str, 'fallback': bool}
    """
    # 定义残差函数
    def residuals(a):
        # 计算 TM 轨道在 rc 处的值
        derivs_ps = _eval_tm_at_rc(rc, l, a, sign_prefactor)

        # 匹配导数
        res = []
        n_deriv_constraints = continuity_orders + 1  # u, u', ..., u^{(continuity_orders)}
        for i in range(n_deriv_constraints):
            res.append(derivs_ps[i] - derivs_ae[i])

        # 范数守恒
        norm_ps = _compute_tm_norm(r_inner, w_inner, l, a)
        res.append(norm_ps - norm_ae)

        return res

    # 初值猜测
    a_init = np.zeros(n_coeffs)
    u_ae_rc = derivs_ae[0]
    if u_ae_rc > 0 and rc > 0:
        a_init[0] = np.log(u_ae_rc / rc**(l+1))
    else:
        a_init[0] = 0.0

    # 求解非线性方程组（增强的fallback策略）
    fallback_used = False
    solver_info = {}

    # 第一次尝试：标准初值
    try:
        a_solution, info, ier, mesg = fsolve(
            residuals, a_init, full_output=True, xtol=1e-10, maxfev=2000
        )
        solver_info = {'ier': ier, 'nfev': info['nfev'], 'mesg': mesg, 'fallback': False}

        if ier == 1:
            return a_solution, solver_info
    except Exception as e:
        mesg = f"First attempt failed: {e}"

    # Fallback 策略 1: 缩小 a_0
    fallback_used = True
    try:
        a_init_fallback = a_init.copy()
        a_init_fallback[0] *= 0.5
        a_solution, info, ier, mesg = fsolve(
            residuals, a_init_fallback, full_output=True, xtol=1e-10, maxfev=2000
        )
        solver_info = {'ier': ier, 'nfev': info['nfev'], 'mesg': mesg, 'fallback': True}

        if ier == 1:
            return a_solution, solver_info
    except Exception as e:
        mesg = f"Fallback 1 failed: {e}"

    # Fallback 策略 2: 零初值
    try:
        a_init_zero = np.zeros(n_coeffs)
        a_solution, info, ier, mesg = fsolve(
            residuals, a_init_zero, full_output=True, xtol=1e-8, maxfev=3000
        )
        solver_info = {'ier': ier, 'nfev': info['nfev'], 'mesg': mesg, 'fallback': True}

        if ier == 1:
            return a_solution, solver_info
    except Exception as e:
        mesg = f"Fallback 2 failed: {e}"

    # Fallback 策略 3: 随机扰动
    try:
        np.random.seed(42)
        a_init_random = a_init + 0.1 * np.random.randn(n_coeffs)
        a_solution, info, ier, mesg = fsolve(
            residuals, a_init_random, full_output=True, xtol=1e-8, maxfev=3000
        )
        solver_info = {'ier': ier, 'nfev': info['nfev'], 'mesg': mesg, 'fallback': True}

        if ier == 1:
            return a_solution, solver_info
    except Exception as e:
        mesg = f"Fallback 3 failed: {e}"

    # Fallback 4: 使用 least_squares (更稳健)
    try:
        res_ls = least_squares(
            residuals,
            a_init,
            method='trf',  # Trust Region Reflective
            ftol=1e-8,     # 放宽公差以避免过度迭代
            xtol=1e-8,
            gtol=1e-8,
            max_nfev=3000
        )
        # 即使 success=False，如果 cost 极小也视为成功
        if res_ls.success or res_ls.cost < 1e-15:
            return res_ls.x, {'ier': 1, 'nfev': res_ls.nfev, 'mesg': 'least_squares success (or low cost)', 'fallback': True}
        else:
            mesg = f"Fallback 4 failed: cost={res_ls.cost}"
    except Exception as e:
        mesg = f"Fallback 4 failed: {e}"

    # 所有尝试失败
    raise RuntimeError(
        f"TM 系数求解失败（所有fallback策略均失败）: {mesg}\n"
        f"最后求解器状态: ier={solver_info.get('ier', -1)}, nfev={solver_info.get('nfev', 0)}"
    )


def _splice_orbital(
    r: np.ndarray,
    u_ae: np.ndarray,
    a_coeff: np.ndarray,
    l: int,
    i_rc: int,
    sign_prefactor: float = 1.0,
) -> np.ndarray:
    """
    拼接内外区生成完整伪轨道

    - r ≤ rc: u = r^{l+1} exp(Σ a_{2i} r^{2i})
    - r > rc: u = u_ae

    Parameters
    ----------
    r : np.ndarray
        完整径向网格
    u_ae : np.ndarray
        AE 轨道
    a_coeff : np.ndarray
        TM 系数
    l : int
        角动量
    i_rc : int
        rc 对应的网格索引
    sign_prefactor : float
        rc 处 AE 轨道的符号，确保拼接时符号一致

    Returns
    -------
    np.ndarray
        拼接后的伪轨道
    """
    u_ps = np.copy(u_ae)

    # 内区：使用 TM 形式
    r_inner = r[:i_rc+1]
    p = np.zeros_like(r_inner)
    for i in range(len(a_coeff)):
        p += a_coeff[i] * r_inner**(2*i)

    # 防止指数溢出（裁剪到安全范围）
    p = np.clip(p, -700, 700)

    u_ps[:i_rc+1] = sign_prefactor * r_inner**(l+1) * np.exp(p)

    # 外区：保持 AE 轨道不变（已经复制）

    return u_ps


def _check_continuity(
    rc: float,
    l: int,
    a_coeff: np.ndarray,
    derivs_ae: np.ndarray,
    order: int,
    sign_prefactor: float = 1.0,
) -> Dict:
    """
    检查内外区连续性（使用解析导数）

    Parameters
    ----------
    rc : float
        伪化半径
    l : int
        角动量
    a_coeff : np.ndarray
        TM 系数
    derivs_ae : np.ndarray
        AE 轨道在 rc 处的导数 [u, u', u'', ...]
    order : int
        检查到几阶导数
    sign_prefactor : float
        AE 轨道在 rc 处的符号

    Returns
    -------
    dict
        包含各阶导数的相对误差
    """
    # 计算 TM 轨道在 rc 处的解析导数
    derivs_ps = _eval_tm_at_rc(rc, l, a_coeff, sign_prefactor)

    result = {}
    labels = ['u', 'du', 'd2u', 'd3u', 'd4u']
    n_derivs = min(order + 1, len(derivs_ps), len(derivs_ae))

    for i, label in enumerate(labels[:n_derivs]):
        if abs(derivs_ae[i]) > 1e-12:
            rel_err = abs(derivs_ps[i] - derivs_ae[i]) / abs(derivs_ae[i])
        else:
            rel_err = abs(derivs_ps[i] - derivs_ae[i])

        result[label] = {
            'ps': float(derivs_ps[i]),
            'ae': float(derivs_ae[i]),
            'rel_error': float(rel_err),
        }

    return result
