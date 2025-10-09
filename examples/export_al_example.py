"""
Al 赝势完整生成与导出示例

演示从全电子求解到赝势导出的完整流程：
1. 全电子原子解（AE）
2. TM 伪化
3. 势反演
4. 可转移性验证
5. 导出JSON + NPZ格式

运行方式
--------
从 AtomPPGen 目录执行：
    PYTHONPATH=../AtomSCF/src:src python examples/export_al_example.py

输出文件
--------
- outputs/al_lda.json：元数据 + 验证报告
- outputs/al_lda.npz：数值数据（网格、势能、波函数）
"""

from pathlib import Path
import numpy as np

from atomppgen import (
    solve_ae_atom,
    tm_pseudize,
    invert_semilocal_potential,
    export_pseudopotential,
)
from atomppgen.validate import run_full_validation


def main():
    print("=" * 60)
    print("Al LDA 赝势生成与导出示例")
    print("=" * 60)

    # 1. 全电子原子解
    print("\n[1/5] 全电子原子解 (Z=13, LDA)...")
    ae = solve_ae_atom(
        Z=13,
        spin_mode='LDA',
        lmax=0,  # 仅s通道示例
        grid_type='exp_transformed',
        grid_params={'n': 600}
    )
    print(f"  ✓ SCF收敛: {ae.converged}, 迭代: {ae.scf_iterations}")
    print(f"  ✓ 3s能级: {ae.eps_by_l[0][-1]:.6f} Ha")

    # 2. TM 伪化
    print("\n[2/5] TM伪化 (rc=2.5 Bohr)...")
    tm_s = tm_pseudize(
        ae.r, ae.w,
        ae.u_by_l[0][-1],  # 3s价电子
        ae.eps_by_l[0][-1],
        l=0,
        rc=2.5  # 最优参数（基于M5实验）
    )
    print(f"  ✓ 范数误差: {tm_s.norm_error:.3e}")

    # 3. 势反演
    print("\n[3/5] 势反演...")
    inv_s = invert_semilocal_potential(tm_s, ae.r)
    print(f"  ✓ rc处势值: {inv_s.diagnostics['V_at_rc']:.6f} Ha")
    print(f"  ✓ 节点数: {inv_s.diagnostics['n_nodes']}")

    # 4. 验证
    print("\n[4/5] 可转移性验证...")
    report = run_full_validation(
        ae,
        tm_dict={0: tm_s},
        inv_dict={0: inv_s},
        r_test=3.0,
        E_range_Ry=(-0.5, 0.5),
        E_step_Ry=0.05,
    )

    # 显示验证结果
    norm_res = report.norm_results[0]
    ld_res = report.log_deriv_results[0]
    ghost_res = report.ghost_result

    print(f"  范数守恒: {norm_res.norm_error:.2e} {'✓' if norm_res.passed else '✗'}")
    print(f"  对数导数:")
    print(f"    - 零点RMS: {ld_res.zero_crossing_rms:.6f} Ha {'✓' if ld_res.zero_crossing_rms < 0.025 else '✗'}")
    print(f"    - 价区RMS: {ld_res.curve_rms_valence:.2f} {'✓' if ld_res.curve_rms_valence < 16.0 else '✗'}")
    print(f"  幽灵态: {ghost_res.n_ghosts} (真) + {ghost_res.n_box_states} (盒) {'✓' if ghost_res.passed else '✗'}")
    print(f"  整体: {'✅ 通过' if report.overall_passed else '❌ 未通过'}")

    # 5. 导出
    print("\n[5/5] 导出赝势...")
    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)

    files = export_pseudopotential(
        ae_result=ae,
        tm_dict={0: tm_s},
        inv_dict={0: inv_s},
        validation_report=report,
        output_prefix=str(output_dir / "al_lda"),
        formats=['json', 'npz'],
    )

    print(f"  ✓ 已导出 {len(files)} 个文件:")
    for f in files:
        size_kb = f.stat().st_size / 1024
        print(f"    - {f.name} ({size_kb:.1f} KB)")

    # 验证NPZ加载
    print("\n[验证] 加载NPZ数据...")
    npz_data = np.load(output_dir / "al_lda.npz")
    print(f"  包含字段: {list(npz_data.keys())}")
    print(f"  径向网格: {len(npz_data['radial_grid'])} 点")
    print(f"  s通道势能: {npz_data['semilocal_potential_l0'].shape}")

    print("\n" + "=" * 60)
    print("示例完成！")
    print(f"  JSON文件: {output_dir / 'al_lda.json'}")
    print(f"  NPZ文件: {output_dir / 'al_lda.npz'}")
    print("=" * 60)


if __name__ == '__main__':
    main()
