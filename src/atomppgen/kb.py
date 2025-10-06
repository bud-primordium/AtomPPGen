"""
Kleinman-Bylander 可分离非局域势转换

将半局域势 {V_l(r)} 转换为可分离形式：

    V_NL = Σ_l |β_l⟩ D_l ⟨β_l|

其中投影子 β_l(r) 和耦合系数 D_l 通过以下方式构造：

1. 选择局域道 V_loc(r) = V_l*(r)（通常选择未占据的高角动量通道，如 d）
2. 构造未归一化投影子：χ_l(r) = [V_l(r) - V_loc(r)] · φ_l(r)
3. 归一化：β_l(r) = χ_l(r) / √⟨χ_l|χ_l⟩，使得 ⟨β_l|β_l⟩ = 1
4. 耦合系数：D_l = ⟨χ_l|χ_l⟩ / ⟨φ_l|ΔV|φ_l⟩
   其中 ΔV = V_l - V_loc

技术要点
--------
- 局域道选择：优先选择无占据态的高 l 通道（如 Al 的 d 通道）
- 径向规约波函数：φ_l(r) = u_l(r)/r
- 数值稳定：r→0 时 u_l ~ r^{l+1}，φ_l ~ r^l，β_l 在原点有限
- 耦合系数：采用归一化投影子形式，D = W/Z 保证物理等价性

主要功能
--------
kb_transform : 半局域势到 KB 可分离形式的转换
KBResult : 存储 KB 转换结果的数据类

参考文献
--------
Kleinman & Bylander, PRL 48, 1425 (1982)
Giannozzi, Notes on pseudopotential generation (2019)
"""

from dataclasses import dataclass
from typing import Dict
import numpy as np

from atomppgen.invert import InvertResult


__all__ = [
    "KBResult",
    "kb_transform",
]


@dataclass
class KBResult:
    """
    Kleinman-Bylander 转换结果

    Attributes
    ----------
    V_loc : np.ndarray
        局域势 V_loc(r)（Hartree）
    beta_l : Dict[int, np.ndarray]
        投影子 {l: β_l(r)}（Bohr^{-3/2}，已归一化）
    D_l : Dict[int, float]
        耦合系数 {l: D_l}（Hartree）
    loc_channel : int
        局域通道角动量量子数 l*
    r : np.ndarray
        径向网格（Bohr）
    diagnostics : Dict
        诊断信息：
        - n_channels : 通道数
        - projector_norms : 归一化前的投影子模 {l: ⟨β_l|β_l⟩}
        - coupling_strengths : 耦合强度 {l: D_l}
        - loc_potential_max : 局域势最大值
        - loc_potential_min : 局域势最小值
    """
    V_loc: np.ndarray
    beta_l: Dict[int, np.ndarray]
    D_l: Dict[int, float]
    loc_channel: int
    r: np.ndarray
    diagnostics: Dict


def kb_transform(
    invert_results: Dict[int, InvertResult],
    u_by_l: Dict[int, np.ndarray],
    r: np.ndarray,
    w: np.ndarray,
    loc_channel: int = 2,
) -> KBResult:
    """
    将半局域势转换为 Kleinman-Bylander 可分离形式

    对于每个非局域通道 l ≠ l*，构造投影子：

        β_l(r) = [V_l(r) - V_loc(r)] · φ_l(r)

    其中 φ_l(r) = u_l(r)/r 是径向规约波函数。

    归一化条件：⟨β_l|β_l⟩ = 1，耦合系数 D_l = ⟨β_l|[V_l - V_loc]|β_l⟩。
    在归一化投影子的情况下，D_l = 1/⟨β_l(原始)|β_l(原始)⟩。

    Parameters
    ----------
    invert_results : Dict[int, InvertResult]
        半局域势反演结果 {l: InvertResult}
    u_by_l : Dict[int, np.ndarray]
        径向波函数 {l: u_l(r)}
    r : np.ndarray
        径向网格（Bohr）
    w : np.ndarray
        积分权重
    loc_channel : int, default=2
        局域通道 l*，默认为 d 通道（l=2）

    Returns
    -------
    KBResult
        包含局域势、投影子、耦合系数和诊断信息

    Raises
    ------
    ValueError
        如果局域通道不在 invert_results 中
        如果网格长度不匹配

    Notes
    -----
    1. **局域道选择**：优先选择无占据态的高角动量通道（如 Al 的 d）。
       占据通道作为局域道会引入额外的数值问题。

    2. **径向规约波函数**：φ_l(r) = u_l(r)/r，满足径向薛定谔方程。
       在 r→0 时，u_l ~ r^{l+1}，因此 φ_l ~ r^l 有限。

    3. **耦合系数公式**：
       标准 KB 形式：V_NL = |χ⟩ D ⟨χ|，其中 D = 1/⟨φ|ΔV|φ⟩
       本实现使用归一化投影子 β = χ/√⟨χ|χ⟩，因此：
       D_l = ⟨χ|χ⟩ / ⟨φ|ΔV|φ⟩ = W / Z
       其中 W = ∫ |χ|² w dr，Z = ∫ ΔV · φ² w dr

    4. **数值稳定**：原点附近 r→0 时使用 np.maximum(r, r_min) 保护除法。

    Examples
    --------
    >>> from atomppgen import solve_ae_atom, tm_pseudize, invert_semilocal_potential
    >>> ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=2)
    >>> invert_results = {}
    >>> for l in [0, 1, 2]:
    ...     tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[l][-1], ae.eps_by_l[l][-1], l, rc=2.0+0.1*l)
    ...     invert_results[l] = invert_semilocal_potential(tm, ae.r)
    >>> kb = kb_transform(invert_results, {l: ae.u_by_l[l][-1] for l in [0,1,2]},
    ...                   ae.r, ae.w, loc_channel=2)
    >>> print(kb.diagnostics['n_channels'])  # 3 (s, p, d)
    """
    # 参数验证
    if loc_channel not in invert_results:
        raise ValueError(f"局域通道 l={loc_channel} 不在 invert_results 中")

    for l, inv in invert_results.items():
        if len(inv.r) != len(r):
            raise ValueError(f"通道 l={l} 的网格长度 {len(inv.r)} 与输入网格 {len(r)} 不匹配")
        if l not in u_by_l:
            raise ValueError(f"通道 l={l} 缺少径向波函数 u_l")

    if len(w) != len(r):
        raise ValueError(f"积分权重长度 {len(w)} 与网格长度 {len(r)} 不匹配")

    # 提取局域势
    V_loc = invert_results[loc_channel].V_l.copy()

    # 构造投影子和耦合系数
    beta_l = {}
    D_l = {}
    projector_norms = {}

    for l, inv in invert_results.items():
        if l == loc_channel:
            # 局域通道不需要投影子
            continue

        # 计算径向规约波函数 φ_l = u_l/r
        # 保护除法：r→0 时使用最小半径
        r_safe = np.maximum(r, 1e-10)
        phi_l = u_by_l[l] / r_safe

        # 构造投影子 β_l = [V_l - V_loc] · φ_l
        delta_V = inv.V_l - V_loc
        beta_raw = delta_V * phi_l

        # 计算归一化前的模 W = ⟨χ|χ⟩ = ∫ |χ|² w dr
        # 其中 χ = ΔV · φ_l（未归一化投影子）
        W = np.sum(beta_raw**2 * w)

        if W <= 0:
            raise ValueError(f"通道 l={l} 的投影子模平方 ⟨χ|χ⟩={W:.3e} ≤ 0，无法归一化")

        # 归一化投影子 β = χ / √W
        beta_normalized = beta_raw / np.sqrt(W)

        # 耦合系数 D_l = W / Z
        # 其中 Z = ⟨φ|ΔV|φ⟩ = ∫ ΔV · φ² w dr
        Z = np.sum(delta_V * (phi_l**2) * w)

        if abs(Z) < 1e-12:
            raise ValueError(f"通道 l={l} 的分母 ⟨φ|ΔV|φ⟩={Z:.3e} 接近零，无法计算耦合系数")

        D_l_val = W / Z

        beta_l[l] = beta_normalized
        D_l[l] = D_l_val
        projector_norms[l] = W

    # 诊断信息
    diagnostics = {
        'n_channels': len(invert_results),
        'projector_norms': projector_norms,
        'coupling_strengths': D_l.copy(),
        'loc_potential_max': float(np.max(V_loc)),
        'loc_potential_min': float(np.min(V_loc)),
    }

    return KBResult(
        V_loc=V_loc,
        beta_l=beta_l,
        D_l=D_l,
        loc_channel=loc_channel,
        r=r,
        diagnostics=diagnostics,
    )
