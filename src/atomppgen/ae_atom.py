"""
全电子（All-Electron, AE）原子求解器封装

调用 AtomSCF 进行 LDA 自洽场计算，获取各角动量通道的径向波函数和本征能量。
用于赝势生成流程中的参考 AE 解。

**重要**：优先使用变量变换方法（transformed）+ 指数网格，精度提升约 7 倍。

主要接口
--------
solve_ae_atom : 求解全电子原子的 LDA-DFT
AEAtomResult : 存储 AE 解的结果数据类

参考
----
- AtomSCF 文档: ../AtomSCF/docs/
- 变量变换方法: AtomSCF/docs/source/algorithm/numerical_methods.rst
"""

from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np

from atomscf.grid import radial_grid_linear, radial_grid_log, radial_grid_exp_transformed
from atomscf.scf import SCFConfig, run_lsda_pz81, run_lsda_vwn


__all__ = [
    "AEAtomResult",
    "solve_ae_atom",
]


@dataclass
class AEAtomResult:
    """
    全电子原子求解结果

    Attributes
    ----------
    Z : int
        原子序数
    xc : str
        交换关联泛函（"PZ81" 或 "VWN"）
    r : np.ndarray
        径向网格，单位 Bohr
    w : np.ndarray
        径向积分权重
    eps_by_l : Dict[int, np.ndarray]
        各角动量通道的本征能量，l → [ε₁, ε₂, ...]（Hartree）
    u_by_l : Dict[int, np.ndarray]
        各角动量通道的径向波函数，l → [u₁(r), u₂(r), ...]
        已归一化：∫ u²(r) dr = 1
    n_total : np.ndarray
        总电子密度（自旋求和），单位 e/Bohr³
    energies : Dict[str, float]
        能量分解，包含 "E_total", "E_ext", "E_H", "E_x", "E_c"（Hartree）
    converged : bool
        SCF 是否收敛
    scf_iterations : int
        SCF 迭代次数
    grid_params : Optional[Dict]
        网格参数（如 delta, Rp），用于变量变换方法
    """

    Z: int
    xc: str
    r: np.ndarray
    w: np.ndarray
    eps_by_l: Dict[int, np.ndarray]
    u_by_l: Dict[int, np.ndarray]
    n_total: np.ndarray
    energies: Dict[str, float]
    converged: bool
    scf_iterations: int
    grid_params: Optional[Dict] = None


def solve_ae_atom(
    Z: int,
    xc: str = "PZ81",
    lmax: int = 2,
    grid_type: str = "exp_transformed",
    grid_params: Optional[Dict] = None,
    scf_params: Optional[Dict] = None,
    spin_mode: str = "LDA",
) -> AEAtomResult:
    """
    求解全电子原子的 LDA 自洽场解

    **推荐配置**：使用 `grid_type="exp_transformed"` + `eig_solver="transformed"`
    以获得最高精度（相比 FD5 方法精度提升约 7 倍）。

    **赝势生成注意**：必须使用 `spin_mode="LDA"`（自旋无关势），这是标准
    赝势格式（如 UPF）的要求。

    Parameters
    ----------
    Z : int
        原子序数（例如 Al 为 13）
    xc : str, optional
        交换关联泛函，"PZ81" 或 "VWN"，默认 "PZ81"
    lmax : int, optional
        最大角动量量子数，默认 2（包含 s, p, d）
    grid_type : str, optional
        网格类型：
        - "exp_transformed": 指数变换网格 + 变量变换求解器（**推荐**，最高精度）
        - "linear": 线性等距网格 + FD5 求解器
        - "log": 对数网格 + Numerov/FD5 求解器
    grid_params : dict, optional
        网格参数，例如：
        - exp_transformed: {"n": 2000, "rmin": 0.0, "rmax": 150.0, "total_span": 7.0}
        - linear: {"n": 1200, "rmin": 1e-6, "rmax": 50.0}
        - log: {"n": 1000, "rmin": 1e-6, "rmax": 50.0}
    scf_params : dict, optional
        SCF 参数，例如：
        - {"tol": 1e-6, "maxiter": 150, "mix_alpha": 0.3}
    spin_mode : str, optional
        自旋模式：
        - "LDA": 强制自旋对称（n_up = n_dn），**赝势生成必须用此模式**
        - "LSDA": 自旋极化（n_up ≠ n_dn），适合开壳层原子的全电子计算
        默认 "LDA"

    Returns
    -------
    AEAtomResult
        包含各角动量通道的波函数、能量、网格等

    Raises
    ------
    ValueError
        如果 xc 不是 "PZ81" 或 "VWN"
        如果 grid_type 不支持
        如果 spin_mode 不是 "LDA" 或 "LSDA"

    Examples
    --------
    >>> # 推荐用法：LDA 模式（赝势生成）
    >>> result = solve_ae_atom(Z=13, spin_mode="LDA", grid_type="exp_transformed")
    >>> print(f"3s energy: {result.eps_by_l[0][2]:.6f} Ha")
    >>> print(f"3p energy: {result.eps_by_l[1][2]:.6f} Ha")

    Notes
    -----
    **LDA vs LSDA**：
    - LDA 强制 n_up = n_dn，产生自旋无关势（UPF 格式要求）
    - LSDA 允许 n_up ≠ n_dn，更适合描述开壳层原子
    - 对于闭壳层原子（如 C），两者结果相同

    **与 NIST 参考数据的差异**：
    AtomSCF 当前精度约 1.5-2.5%（相对于 NIST LSD 数据）。
    对于赝势生成，价层轨道相对精度（~0.03 Ha）通常已足够。

    References
    ----------
    - 变量变换方法：AtomSCF/docs/source/algorithm/numerical_methods.rst
    - NIST 原子数据：https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations
    """
    # 验证输入
    if xc not in ("PZ81", "VWN"):
        raise ValueError(f"xc must be 'PZ81' or 'VWN', got '{xc}'")
    if grid_type not in ("linear", "log", "exp_transformed"):
        raise ValueError(f"grid_type must be 'linear', 'log', or 'exp_transformed', got '{grid_type}'")
    if spin_mode not in ("LDA", "LSDA"):
        raise ValueError(f"spin_mode must be 'LDA' or 'LSDA', got '{spin_mode}'")

    # 设置默认网格参数
    if grid_params is None:
        grid_params = {}

    # 生成径向网格
    delta, Rp = None, None  # 变量变换参数（仅 exp_transformed 需要）

    if grid_type == "exp_transformed":
        # 推荐：变量变换方法
        n = grid_params.get("n", 1200)
        rmin = grid_params.get("rmin", 0.0)  # 可以从 0 开始
        rmax = grid_params.get("rmax", 120.0)
        total_span = grid_params.get("total_span", 6.5)
        r, w, delta, Rp = radial_grid_exp_transformed(
            n=n, rmin=rmin, rmax=rmax, total_span=total_span
        )
        eig_solver = "transformed"
    elif grid_type == "linear":
        n = grid_params.get("n", 1200)
        rmin = grid_params.get("rmin", 1e-6)
        rmax = grid_params.get("rmax", 50.0)
        r, w = radial_grid_linear(n=n, rmin=rmin, rmax=rmax)
        eig_solver = "fd5"
    else:  # log
        n = grid_params.get("n", 1000)
        rmin = grid_params.get("rmin", 1e-6)
        rmax = grid_params.get("rmax", 50.0)
        r, w = radial_grid_log(n=n, rmin=rmin, rmax=rmax)
        eig_solver = "fd5_aux"

    # 设置默认 SCF 参数
    if scf_params is None:
        scf_params = {}
    tol = scf_params.get("tol", 1e-7)
    maxiter = scf_params.get("maxiter", 200)
    mix_alpha = scf_params.get("mix_alpha", 0.25)
    eigs_per_l = scf_params.get("eigs_per_l", 3)  # 每个 l 求 3 个态

    # 构造 SCF 配置
    cfg_kwargs = dict(
        Z=Z,
        r=r,
        w=w,
        lmax=lmax,
        eigs_per_l=eigs_per_l,
        eig_solver=eig_solver,
        xc=xc,
        mix_alpha=mix_alpha,
        tol=tol,
        maxiter=maxiter,
        spin_mode=spin_mode,  # 新增：传递自旋模式
    )

    # 变量变换方法需要额外参数
    if grid_type == "exp_transformed":
        cfg_kwargs["delta"] = delta
        cfg_kwargs["Rp"] = Rp

    cfg = SCFConfig(**cfg_kwargs)

    # 根据 XC 泛函选择 SCF 函数
    if xc == "PZ81":
        scf_result = run_lsda_pz81(cfg, verbose=False)
    else:  # VWN
        scf_result = run_lsda_vwn(cfg, verbose=False)

    # 提取自旋分辨的结果并合并
    eps_by_l = {}
    u_by_l = {}

    for l in range(lmax + 1):
        # 合并 up/down 自旋通道的能量和波函数
        # 注意：对于闭壳层或自旋平均，up 和 down 应该相同或相近
        eps_up = scf_result.eps_by_l_sigma.get((l, "up"), np.array([]))
        eps_dn = scf_result.eps_by_l_sigma.get((l, "down"), np.array([]))

        u_up = scf_result.u_by_l_sigma.get((l, "up"), np.array([]))
        u_dn = scf_result.u_by_l_sigma.get((l, "down"), np.array([]))

        # 策略：优先使用 up 通道，如果 up 为空则用 down
        if len(eps_up) > 0:
            eps_by_l[l] = eps_up
            u_by_l[l] = u_up
        elif len(eps_dn) > 0:
            eps_by_l[l] = eps_dn
            u_by_l[l] = u_dn
        else:
            # 空通道，创建空数组
            eps_by_l[l] = np.array([])
            u_by_l[l] = np.array([]).reshape(0, len(r))

    # 计算总电子密度（自旋求和）
    n_total = scf_result.n_up + scf_result.n_dn

    # 保存网格参数
    grid_params_out = None
    if grid_type == "exp_transformed":
        grid_params_out = {"delta": delta, "Rp": Rp}

    # 组装结果
    result = AEAtomResult(
        Z=Z,
        xc=xc,
        r=r,
        w=w,
        eps_by_l=eps_by_l,
        u_by_l=u_by_l,
        n_total=n_total,
        energies=scf_result.energies,
        converged=scf_result.converged,
        scf_iterations=scf_result.iterations,
        grid_params=grid_params_out,
    )

    return result
