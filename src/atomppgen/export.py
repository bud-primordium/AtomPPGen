"""
赝势导出模块

提供多种格式的赝势文件导出功能：
- JSON：结构化元数据（参数、验证结果）
- NPZ：数值数据（径向网格、势能、波函数）
- UPF：Quantum ESPRESSO 兼容格式（实验性，当前仅提供最小可解析结构）

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
import xml.etree.ElementTree as ET

from atomppgen.ae_atom import AEAtomResult
from atomppgen.tm import TMResult
from atomppgen.invert import InvertResult
from atomppgen.kb import KBResult
from atomppgen.validate import ValidationReport


__all__ = [
    "PseudopotentialData",
    "export_pseudopotential",
    "export_upf",
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
    kb_loc_channel : Optional[int]
        KB 局域通道角动量量子数 l*（可选）
    kb_V_loc : Optional[np.ndarray]
        KB 局域势 V_loc(r)，单位 Hartree（可选）
    kb_beta_l : Optional[Dict[int, np.ndarray]]
        KB 投影子 {l: β_l(r)}（可选，已归一化）
    kb_D_l : Optional[Dict[int, float]]
        KB 耦合系数 {l: D_l}（可选，Hartree）
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
    kb_loc_channel: Optional[int] = None
    kb_V_loc: Optional[np.ndarray] = None
    kb_beta_l: Optional[Dict[int, np.ndarray]] = None
    kb_D_l: Optional[Dict[int, float]] = None
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

    if data.kb_V_loc is not None:
        json_data["pseudopotential"]["kb"] = {
            "loc_channel": int(data.kb_loc_channel) if data.kb_loc_channel is not None else None,
            "nonlocal_channels": sorted(int(l) for l in (data.kb_beta_l or {}).keys()),
            "note": "KB arrays (V_loc, beta_l, D_l) are stored in NPZ",
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

    # KB 数据（可选）
    if data.kb_V_loc is not None:
        npz_data["kb_V_loc"] = data.kb_V_loc
        if data.kb_loc_channel is not None:
            npz_data["kb_loc_channel"] = int(data.kb_loc_channel)
        if data.kb_beta_l:
            for l, beta in data.kb_beta_l.items():
                npz_data[f"kb_beta_l{int(l)}"] = beta
        if data.kb_D_l:
            for l, D in data.kb_D_l.items():
                npz_data[f"kb_D_l{int(l)}"] = float(D)

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


def _infer_z_valence(Z: int) -> Optional[float]:
    """
    尽量推断常见元素的价电子数（用于 UPF header）。

    Notes
    -----
    这是教学用途的保守映射。更严格的做法应由“芯态选择/价电子配置”显式给出。
    """
    mapping = {
        11: 1.0,  # Na: 3s^1
        13: 3.0,  # Al: 3s^2 3p^1
        14: 4.0,  # Si: 3s^2 3p^2
    }
    return mapping.get(int(Z))


def _format_float_array(values: np.ndarray, per_line: int = 6) -> str:
    values = np.asarray(values).ravel()
    lines = []
    for i in range(0, len(values), per_line):
        chunk = values[i:i + per_line]
        lines.append(" ".join(f"{float(x):.12e}" for x in chunk))
    return "\n".join(lines)


def export_upf(
    data: PseudopotentialData,
    output_path: Union[str, Path],
    metadata: Optional[Dict] = None,
) -> Path:
    """
    导出 UPF v2（实验性）

    当前目标是提供“最小可解析、字段可追溯”的 UPF 结构，便于后续逐步对齐 QE 的严格要求。
    若要获得严格可用于 Quantum ESPRESSO 的 UPF，请关注后续版本更新。

    Parameters
    ----------
    data : PseudopotentialData
        完整赝势数据包
    output_path : str | Path
        输出文件路径（建议后缀为 .upf）
    metadata : dict, optional
        额外信息；UPF 需要的关键字段可在此给出：
        - z_valence : float（建议显式提供）

    Returns
    -------
    Path
        生成的 UPF 文件路径
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _export_upf(data, output_path, metadata=metadata)
    return output_path


def _export_upf(
    data: PseudopotentialData,
    output_path: Path,
    metadata: Optional[Dict] = None,
) -> None:
    z_valence = None
    if metadata and metadata.get("z_valence") is not None:
        z_valence = float(metadata["z_valence"])
    else:
        z_valence = _infer_z_valence(int(data.Z))

    if z_valence is None:
        raise ValueError("UPF 导出需要提供 z_valence（可在 metadata['z_valence'] 指定）")

    # UPF 传统上使用 Ry 能量单位（1 Ha = 2 Ry）
    ha_to_ry = 2.0

    r = np.asarray(data.radial_grid, dtype=float)
    rab = np.gradient(r)

    root = ET.Element("UPF", {"version": "2.0.1"})

    info = ET.SubElement(root, "PP_INFO")
    info.text = (
        "Generated by AtomPPGen (experimental UPF writer).\n"
        "This file is intended to be a traceable intermediate format; QE compatibility is not guaranteed yet."
    )

    ET.SubElement(
        root,
        "PP_HEADER",
        {
            "generated": "AtomPPGen",
            "element": str(data.symbol),
            "pseudo_type": "NC",
            "relativistic": "nonrelativistic",
            "functional": str(data.xc_functional),
            "z_valence": f"{z_valence:.1f}",
            "l_max": str(int(max(data.semilocal_potentials_by_l.keys())) if data.semilocal_potentials_by_l else 0),
            "mesh_size": str(int(len(r))),
            "has_kb": "T" if data.kb_V_loc is not None else "F",
            "units": "bohr_ry",
        },
    )

    mesh = ET.SubElement(root, "PP_MESH")
    pp_r = ET.SubElement(mesh, "PP_R")
    pp_r.text = _format_float_array(r)
    pp_rab = ET.SubElement(mesh, "PP_RAB")
    pp_rab.text = _format_float_array(rab)

    if data.kb_V_loc is not None and data.kb_beta_l is not None and data.kb_D_l is not None:
        pp_local = ET.SubElement(root, "PP_LOCAL")
        pp_local.text = _format_float_array(np.asarray(data.kb_V_loc) * ha_to_ry)

        nonlocal_el = ET.SubElement(root, "PP_NONLOCAL")
        beta_items = sorted((int(l), np.asarray(beta)) for l, beta in data.kb_beta_l.items())
        dij_diag = []
        for idx, (l, beta) in enumerate(beta_items, start=1):
            beta_el = ET.SubElement(
                nonlocal_el,
                "PP_BETA",
                {
                    "index": str(idx),
                    "angular_momentum": str(int(l)),
                    "cutoff_radius_index": "0",
                },
            )
            beta_el.text = _format_float_array(beta)

            D = float(data.kb_D_l[int(l)])
            dij_diag.append(D * ha_to_ry)

        dij_el = ET.SubElement(nonlocal_el, "PP_DIJ")
        nproj = len(dij_diag)
        dij_matrix = np.zeros((nproj, nproj), dtype=float)
        for i in range(nproj):
            dij_matrix[i, i] = dij_diag[i]
        dij_el.text = _format_float_array(dij_matrix.ravel(), per_line=6)

    else:
        semilocal = ET.SubElement(root, "PP_SEMILOCAL")
        for l in sorted(int(k) for k in data.semilocal_potentials_by_l.keys()):
            vl = ET.SubElement(semilocal, "PP_VL", {"angular_momentum": str(int(l))})
            vl.text = _format_float_array(np.asarray(data.semilocal_potentials_by_l[l]) * ha_to_ry)

    tree = ET.ElementTree(root)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)


def export_pseudopotential(
    ae_result: AEAtomResult,
    tm_dict: Dict[int, TMResult],
    inv_dict: Dict[int, InvertResult],
    validation_report: ValidationReport,
    output_prefix: str,
    kb_result: Optional[KBResult] = None,
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
        kb_loc_channel=kb_result.loc_channel if kb_result is not None else None,
        kb_V_loc=kb_result.V_loc if kb_result is not None else None,
        kb_beta_l=kb_result.beta_l if kb_result is not None else None,
        kb_D_l=kb_result.D_l if kb_result is not None else None,
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
            upf_path = output_prefix.with_suffix('.upf')
            export_upf(pp_data, upf_path, metadata=metadata)
            output_files.append(upf_path)

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
