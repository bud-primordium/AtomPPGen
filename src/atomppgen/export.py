"""
赝势导出模块

提供多种格式的赝势文件导出功能：
- JSON：结构化元数据（参数、验证结果）
- NPZ：数值数据（径向网格、势能、波函数）
- UPF：Quantum ESPRESSO 兼容格式（实验性）

主要函数
--------
export_pseudopotential : 统一导出接口
"""

import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Union
from datetime import datetime

from atomppgen.ae_atom import AEAtomResult
from atomppgen.tm import TMResult
from atomppgen.invert import InvertResult
from atomppgen.validate import ValidationReport


__all__ = [
    "PseudopotentialData",
    "export_pseudopotential",
]


@dataclass
class PseudopotentialData:
    """
    完整赝势数据包

    聚合全电子解、TM伪化、势反演、验证报告等所有数据源，
    用于统一格式导出。

    Attributes
    ----------
    Z : int
        原子序数
    symbol : str
        元素符号（如 'Al'）
    spin_mode : str
        自旋模式（'LDA' / 'GGA'）
    xc_functional : str
        交换关联泛函（'PZ81' / 'VWN'）
    generation_params : Dict
        生成参数（TM参数、网格参数等）
    radial_grid : np.ndarray
        径向网格，单位 Bohr
    radial_weights : np.ndarray
        积分权重
    ae_eigenvalues_by_l : Dict[int, np.ndarray]
        全电子本征值（按角动量 l 索引），单位 Hartree
    ae_wavefunctions_by_l : Dict[int, List[np.ndarray]]
        全电子径向波函数（u = r·R）
    pseudo_wavefunctions_by_l : Dict[int, np.ndarray]
        赝波函数（伪化后的价电子态）
    semilocal_potentials_by_l : Dict[int, np.ndarray]
        半局域势 V_l(r)，单位 Hartree
    validation_report : ValidationReport
        完整验证报告
    generation_date : str
        生成时间戳（ISO 8601格式）
    code_version : str
        代码版本号
    git_commit : Optional[str]
        Git 提交哈希（可选）
    """
    Z: int
    symbol: str
    spin_mode: str
    xc_functional: str
    generation_params: Dict
    radial_grid: np.ndarray
    radial_weights: np.ndarray
    ae_eigenvalues_by_l: Dict[int, np.ndarray]
    ae_wavefunctions_by_l: Dict[int, List[np.ndarray]]
    pseudo_wavefunctions_by_l: Dict[int, np.ndarray]
    semilocal_potentials_by_l: Dict[int, np.ndarray]
    validation_report: ValidationReport
    generation_date: str
    code_version: str = "0.1.0"
    git_commit: Optional[str] = None


def _export_json(
    data: PseudopotentialData,
    output_path: Path,
    include_arrays: bool = False,
) -> None:
    """
    导出JSON格式（结构化元数据）

    Parameters
    ----------
    data : PseudopotentialData
        完整赝势数据
    output_path : Path
        输出文件路径
    include_arrays : bool, default=False
        是否包含数值数组（大文件警告）。
        若 False，数组数据仅记录形状和范围。

    Notes
    -----
    JSON格式优先存储元数据（参数、验证结果），数值数据建议使用NPZ格式。

    输出单位约定：
    - 能量：Hartree
    - 长度：Bohr
    """
    output_path = Path(output_path)

    # 构建JSON数据结构
    json_data = {
        "metadata": {
            "element": {
                "Z": data.Z,
                "symbol": data.symbol,
            },
            "xc": {
                "functional": data.xc_functional,
                "spin_mode": data.spin_mode,
            },
            "generation_date": data.generation_date,
            "code": {
                "name": "AtomPPGen",
                "version": data.code_version,
            },
            "units": "Hartree_atomic",  # 明确标注单位
        },
        "generation_params": data.generation_params,
        "all_electron_reference": {
            "eigenvalues_Ha": {
                f"l={l}": eps.tolist() for l, eps in data.ae_eigenvalues_by_l.items()
            },
        },
        "pseudopotential": {
            "channels": list(data.semilocal_potentials_by_l.keys()),
            "data_file": str(output_path.with_suffix('.npz').name),
            "note": "Numerical arrays stored in NPZ format for efficiency",
        },
        "validation_report": data.validation_report.to_dict(),
    }

    # 可选：包含git commit
    if data.git_commit:
        json_data["metadata"]["code"]["git_commit"] = data.git_commit

    # 可选：包含数组数据（仅用于调试/小数据集）
    if include_arrays:
        json_data["radial_grid"] = {
            "grid_Bohr": data.radial_grid.tolist(),
            "weights": data.radial_weights.tolist(),
        }
        json_data["pseudopotential"]["semilocal_potentials_Ha"] = {
            f"l={l}": V.tolist() for l, V in data.semilocal_potentials_by_l.items()
        }
    else:
        # 仅记录数组形状和范围
        json_data["radial_grid"] = {
            "n_points": int(len(data.radial_grid)),
            "r_min_Bohr": float(data.radial_grid[0]),
            "r_max_Bohr": float(data.radial_grid[-1]),
        }

    # 确保输出目录存在
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # 写入JSON文件
    with output_path.open('w', encoding='utf-8') as f:
        json.dump(json_data, f, indent=2, ensure_ascii=False)


def _export_npz(
    data: PseudopotentialData,
    output_path: Path,
) -> None:
    """
    导出NPZ格式（压缩数值数据）

    Parameters
    ----------
    data : PseudopotentialData
        完整赝势数据
    output_path : Path
        输出文件路径

    Notes
    -----
    NPZ格式存储所有数值数组，使用压缩节省空间。

    数组命名约定：
    - 'radial_grid'：径向网格 (Bohr)
    - 'ae_eigenvalues_l{l}'：全电子本征值 (Ha)
    - 'ae_wavefunction_l{l}_n{n}'：全电子波函数 u(r)
    - 'ps_wavefunction_l{l}'：赝波函数（价电子）
    - 'semilocal_potential_l{l}'：半局域势 V_l(r) (Ha)

    输出单位约定：
    - 能量：Hartree
    - 长度：Bohr
    """
    output_path = Path(output_path)

    # 构建NPZ数据字典
    npz_data = {
        # 网格
        'radial_grid': data.radial_grid,
        'radial_weights': data.radial_weights,

        # 元数据（标量）
        'Z': data.Z,
        'xc_code': _xc_to_code(data.xc_functional),
    }

    # 全电子本征值和波函数
    for l, eigenvalues in data.ae_eigenvalues_by_l.items():
        npz_data[f'ae_eigenvalues_l{l}'] = eigenvalues

        # 仅保存价电子波函数（最后一个n）
        if l in data.ae_wavefunctions_by_l:
            wavefunctions = data.ae_wavefunctions_by_l[l]
            if len(wavefunctions) > 0:
                n_valence = len(eigenvalues)  # 价电子是最高占据态
                npz_data[f'ae_wavefunction_l{l}_n{n_valence}'] = wavefunctions[-1]

    # 赝势数据
    for l, wf in data.pseudo_wavefunctions_by_l.items():
        npz_data[f'ps_wavefunction_l{l}'] = wf

    for l, V in data.semilocal_potentials_by_l.items():
        npz_data[f'semilocal_potential_l{l}'] = V

    # 确保输出目录存在
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # 使用压缩保存
    np.savez_compressed(output_path, **npz_data)


def _xc_to_code(xc_functional: str) -> int:
    """
    交换关联泛函名称转换为数值代码

    Parameters
    ----------
    xc_functional : str
        泛函名称（如 'PZ81', 'VWN'）

    Returns
    -------
    int
        数值代码（1=PZ81, 2=VWN, ...）
    """
    xc_map = {
        'PZ81': 1,
        'VWN': 2,
        'PW91': 3,
        'PBE': 4,
    }
    return xc_map.get(xc_functional.upper(), 0)


def export_pseudopotential(
    ae_result: AEAtomResult,
    tm_dict: Dict[int, TMResult],
    inv_dict: Dict[int, InvertResult],
    validation_report: ValidationReport,
    output_prefix: str,
    formats: List[str] = ['json', 'npz'],
    metadata: Optional[Dict] = None,
) -> List[Path]:
    """
    导出赝势到多种格式

    Parameters
    ----------
    ae_result : AEAtomResult
        全电子原子解
    tm_dict : Dict[int, TMResult]
        各通道TM伪化结果（按角动量 l 索引）
    inv_dict : Dict[int, InvertResult]
        各通道势反演结果（按角动量 l 索引）
    validation_report : ValidationReport
        完整验证报告
    output_prefix : str
        输出文件名前缀（如 'outputs/al_lda'）
    formats : List[str], default=['json', 'npz']
        导出格式列表，支持 'json', 'npz', 'upf'
    metadata : Optional[Dict], default=None
        额外元数据（如 git_commit, 备注）

    Returns
    -------
    List[Path]
        生成的文件路径列表

    Raises
    ------
    ValueError
        若 tm_dict 和 inv_dict 的 l 通道不一致

    Examples
    --------
    >>> files = export_pseudopotential(
    ...     ae, tm_dict, inv_dict, report,
    ...     output_prefix='outputs/al_lda',
    ...     formats=['json', 'npz']
    ... )
    >>> print(files)
    [PosixPath('outputs/al_lda.json'), PosixPath('outputs/al_lda.npz')]

    Notes
    -----
    输出单位约定：
    - 能量：Hartree 原子单位
    - 长度：Bohr

    JSON格式包含元数据和验证报告，NPZ格式包含数值数组。
    推荐同时导出JSON+NPZ以获得完整数据集。
    """
    # 输入验证：检查l通道一致性
    if set(tm_dict.keys()) != set(inv_dict.keys()):
        raise ValueError(
            f"TM和势反演的l通道不一致: tm={set(tm_dict.keys())}, inv={set(inv_dict.keys())}"
        )

    # 构建完整数据包
    pp_data = PseudopotentialData(
        Z=ae_result.Z,
        symbol=_get_element_symbol(ae_result.Z),
        spin_mode="LDA",  # 赝势生成固定使用 LDA（自旋无关势）
        xc_functional=ae_result.xc,  # 直接使用 AEAtomResult 的 xc 字段
        generation_params=_collect_generation_params(ae_result, tm_dict),
        radial_grid=ae_result.r,
        radial_weights=ae_result.w,
        ae_eigenvalues_by_l=ae_result.eps_by_l,
        ae_wavefunctions_by_l=ae_result.u_by_l,
        pseudo_wavefunctions_by_l={l: tm.u_ps for l, tm in tm_dict.items()},
        semilocal_potentials_by_l={l: inv.V_l for l, inv in inv_dict.items()},
        validation_report=validation_report,
        generation_date=datetime.now().isoformat(),
        code_version="0.1.0",
        git_commit=metadata.get('git_commit') if metadata else None,
    )

    # 导出到各种格式
    output_prefix = Path(output_prefix)
    output_files = []

    for fmt in formats:
        if fmt.lower() == 'json':
            json_path = output_prefix.with_suffix('.json')
            _export_json(pp_data, json_path)
            output_files.append(json_path)

        elif fmt.lower() == 'npz':
            npz_path = output_prefix.with_suffix('.npz')
            _export_npz(pp_data, npz_path)
            output_files.append(npz_path)

        elif fmt.lower() == 'upf':
            # UPF格式暂未实现
            print(f"警告：UPF格式暂未实现，跳过")

        else:
            print(f"警告：未知格式 '{fmt}'，跳过")

    return output_files


def _get_element_symbol(Z: int) -> str:
    """根据原子序数获取元素符号"""
    symbols = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar',
    }
    return symbols.get(Z, f'Z{Z}')


def _infer_xc_functional(ae_result: AEAtomResult) -> str:
    """从AE结果推断XC泛函名称（已弃用，直接使用 ae_result.xc）"""
    return ae_result.xc


def _collect_generation_params(
    ae_result: AEAtomResult,
    tm_dict: Dict[int, TMResult],
) -> Dict:
    """收集生成参数（TM参数、网格参数等）"""
    # TM参数（从各通道提取rc）
    tm_params = {
        'rc_by_l': {l: float(tm.rc) for l, tm in tm_dict.items()},
        'continuity_order': 2,  # TM固定为2阶
    }

    # 网格参数（从AE结果推断）
    grid_params = {
        'type': 'exp_transformed',  # 假设使用exp变换网格
        'n_points': len(ae_result.r),
        'r_max_Bohr': float(ae_result.r[-1]),
    }

    return {
        'tm_pseudization': tm_params,
        'radial_grid': grid_params,
    }
