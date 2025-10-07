"""
Al s 通道赝势验证示例

演示完整验证流程：范数守恒 → 对数导数匹配 → 幽灵态检测
生成 JSON 报告和对数导数曲线图
"""

import json
from pathlib import Path
import numpy as np

from atomppgen import solve_ae_atom, tm_pseudize, invert_semilocal_potential
from atomppgen.validate import run_full_validation

# 可选依赖：matplotlib
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def main():
    """运行 Al s 通道完整验证"""

    print("=" * 60)
    print("Al 3s 赝势可转移性验证")
    print("=" * 60)

    # 1. 全电子原子解
    print("\n[1/4] 求解 Al 全电子原子...")
    ae = solve_ae_atom(
        Z=13,
        spin_mode='LDA',
        lmax=0,
        grid_type='exp_transformed',
        grid_params={'n': 600}
    )
    print(f"  ✓ SCF 收敛: {ae.converged}, 迭代: {ae.scf_iterations}")
    print(f"  ✓ 3s 能级: {ae.eps_by_l[0][-1]:.6f} Ha")

    # 2. TM 伪化
    print("\n[2/4] TM 伪化 (s 通道, rc=2.0 Bohr)...")
    tm_s = tm_pseudize(
        ae.r, ae.w,
        ae.u_by_l[0][-1],
        ae.eps_by_l[0][-1],
        l=0,
        rc=2.0
    )
    print(f"  ✓ 范数误差: {tm_s.norm_error:.3e}")
    print(f"  ✓ 连续性阶数: {tm_s.continuity_orders}")

    # 3. 势反演
    print("\n[3/4] 反演半局域势...")
    inv_s = invert_semilocal_potential(tm_s, ae.r)
    print(f"  ✓ rc 处势值: {inv_s.diagnostics['V_at_rc']:.6f} Ha")
    print(f"  ✓ 节点数: {inv_s.diagnostics['n_nodes']}")

    # 4. 完整验证
    print("\n[4/4] 运行完整验证流程...")
    print("  (范数守恒 → 对数导数匹配 → 幽灵态检测)")

    tm_dict = {0: tm_s}
    inv_dict = {0: inv_s}

    report = run_full_validation(
        ae, tm_dict, inv_dict,
        r_test=3.0,
        E_range_Ry=(-0.5, 0.5),  # 完整能量窗口
        E_step_Ry=0.05
    )

    # 5. 显示结果
    print("\n" + "=" * 60)
    print("验证结果")
    print("=" * 60)

    # 范数守恒
    norm_result = report.norm_results[0]
    print(f"\n[范数守恒]")
    print(f"  误差: {norm_result.norm_error:.3e}")
    print(f"  通过: {'✓' if norm_result.passed else '✗'}")

    # 对数导数
    ld_result = report.log_deriv_results[0]
    print(f"\n[对数导数匹配 @ r={ld_result.r_test} Bohr]")
    print(f"  零点 RMS: {ld_result.zero_crossing_rms:.6f} Ha ({ld_result.zero_crossing_rms*2:.6f} Ry)")
    print(f"  曲线 RMS: {ld_result.curve_rms:.6f}")
    print(f"  零点数 (AE/PS): {len(ld_result.zero_crossings_AE)}/{len(ld_result.zero_crossings_PS)}")
    print(f"  通过: {'✓' if ld_result.passed else '✗'}")

    # 幽灵态
    if report.ghost_result:
        ghost = report.ghost_result
        print(f"\n[幽灵态检测]")
        print(f"  方法: {ghost.method}")
        print(f"  窗口内束缚态: {ghost.diagnostics['n_bound_states_in_window']}")
        print(f"  幽灵态数量: {ghost.n_ghosts}")
        print(f"  通过: {'✓' if ghost.passed else '✗'}")

    # 整体判定
    print(f"\n{'=' * 60}")
    print(f"整体验证: {'✓ 通过' if report.overall_passed else '✗ 未通过'}")
    print(f"{'=' * 60}")

    # 6. 导出 JSON 报告
    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)

    report_path = output_dir / "validation_report_al_s.json"
    print(f"\n[导出] 生成 JSON 报告: {report_path}")
    with open(report_path, 'w', encoding='utf-8') as f:
        json.dump(report.to_dict(), f, indent=2, ensure_ascii=False)
    print(f"  ✓ 报告已保存")

    # 7. 绘制对数导数曲线
    if HAS_MPL:
        fig_path = output_dir / "log_derivative_al_s.png"
        print(f"\n[可视化] 生成对数导数曲线: {fig_path}")

        plt.figure(figsize=(10, 6))
        plt.plot(ld_result.energies, ld_result.L_AE,
                 'b-', linewidth=2, label='All-Electron')
        plt.plot(ld_result.energies, ld_result.L_PS,
                 'r--', linewidth=2, label='Pseudopotential')

        # 标记零点
        if len(ld_result.zero_crossings_AE) > 0:
            plt.scatter(ld_result.zero_crossings_AE,
                       np.zeros_like(ld_result.zero_crossings_AE),
                       c='blue', s=100, marker='o', zorder=5,
                       label='AE zeros')
        if len(ld_result.zero_crossings_PS) > 0:
            plt.scatter(ld_result.zero_crossings_PS,
                       np.zeros_like(ld_result.zero_crossings_PS),
                       c='red', s=100, marker='x', zorder=5,
                       label='PS zeros')

        plt.axhline(0, color='k', linewidth=0.5, linestyle=':')
        plt.xlabel('Energy (Ha)', fontsize=12)
        plt.ylabel(f'L(E, r={ld_result.r_test} Bohr)', fontsize=12)
        plt.title('Log Derivative Matching - Al 3s', fontsize=14, fontweight='bold')
        plt.legend(loc='best', fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_path, dpi=150)
        plt.close()
        print(f"  ✓ 图像已保存")
    else:
        print(f"\n[可视化] 跳过（matplotlib 未安装）")

    print("\n" + "=" * 60)
    print("验证完成！")
    print(f"  - JSON 报告: {report_path}")
    if HAS_MPL:
        print(f"  - L(E) 曲线: {output_dir / 'log_derivative_al_s.png'}")
    print("=" * 60)


if __name__ == "__main__":
    main()
