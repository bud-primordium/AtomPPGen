"""
Na（简单金属）LDA-PZ81 模守恒赝势：端到端生成 + 验证 + 导出

目的
----
1) 给出一个 Na 示例，作为金属阈值对照组；
2) 验证对数导数阈值的差异化策略：Na 使用金属阈值（curve_rms_valence < 16.0）。

运行
----
python AtomPPGen/examples/validate_na_lda.py

输出
----
AtomPPGen/outputs/na_lda/
  - na_lda.json
  - na_lda.npz
"""

from __future__ import annotations

from pathlib import Path

from atomppgen import (
    solve_ae_atom,
    tm_pseudize,
    invert_semilocal_potential,
    kb_transform,
    export_pseudopotential,
)
from atomppgen.validate import (
    run_full_validation,
    is_covalent,
    COVALENT_CURVE_RMS_THRESHOLD,
    METAL_CURVE_RMS_THRESHOLD,
)


def _find_repo_root() -> Path:
    cwd = Path.cwd().resolve()
    return next(p for p in [cwd, *cwd.parents] if (p / "AtomPPGen" / "pyproject.toml").is_file())


def main() -> None:
    Z = 11
    symbol = "Na"
    xc = "PZ81"
    lmax = 2  # s, p, d（d 作为散射/局域通道候选）
    loc_channel = 2

    # Na 的价电子更“软”，通常可以选更大的 rc 以降低平面波截断需求
    rc_by_l = {0: 2.6, 1: 2.8, 2: 3.0}

    print(f"[1/6] AE 求解：{symbol} (Z={Z}), xc={xc}, lmax={lmax}")
    ae = solve_ae_atom(
        Z=Z,
        spin_mode="LDA",
        xc=xc,
        lmax=lmax,
        grid_type="exp_transformed",
        grid_params={"n": 900, "rmax": 140.0},
        scf_params={"tol": 1e-6, "maxiter": 200},
    )
    if not ae.converged:
        raise RuntimeError(f"AE SCF 未收敛：iterations={ae.scf_iterations}")

    print("[2/6] TM 伪化（按通道）")
    tm_dict = {}
    for l in range(lmax + 1):
        tm_dict[l] = tm_pseudize(
            ae.r,
            ae.w,
            ae.u_by_l[l][-1],
            ae.eps_by_l[l][-1],
            l=l,
            rc=rc_by_l[l],
        )
        print(f"  l={l} rc={rc_by_l[l]:.2f} norm_error={tm_dict[l].norm_error:.3e}")

    print("[3/6] 势反演（半局域势）")
    inv_dict = {l: invert_semilocal_potential(tm_dict[l], ae.r) for l in tm_dict}

    print(f"[4/6] KB 转换（loc_channel={loc_channel}）")
    kb = kb_transform(inv_dict, {l: ae.u_by_l[l][-1] for l in tm_dict}, ae.r, ae.w, loc_channel=loc_channel)

    print("[5/6] 完整验证（范数/对数导数/幽灵态）")
    report = run_full_validation(
        ae,
        tm_dict=tm_dict,
        inv_dict=inv_dict,
        r_test=max(rc_by_l.values()) + 0.5,
        E_range_Ry=(-0.5, 0.5),
        E_step_Ry=0.05,
    )

    ld_threshold = COVALENT_CURVE_RMS_THRESHOLD if is_covalent(symbol) else METAL_CURVE_RMS_THRESHOLD
    print(f"  对数导数阈值（curve_rms_valence）：{ld_threshold}（共价元素={is_covalent(symbol)}）")
    for l, r in report.log_deriv_results.items():
        print(
            f"  l={l} zero_rms={r.zero_crossing_rms:.4f} Ha, "
            f"valence_rms={r.curve_rms_valence:.4f}, passed={r.passed}"
        )
    print(f"  幽灵态：n_ghosts={report.ghost_result.n_ghosts}, n_box={report.ghost_result.n_box_states}, passed={report.ghost_result.passed}")
    print(f"  overall_passed={report.overall_passed}")

    print("[6/6] 导出（JSON+NPZ，包含 KB 数据）")
    repo_root = _find_repo_root()
    out_dir = repo_root / "AtomPPGen" / "outputs" / "na_lda"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = out_dir / "na_lda"

    files = export_pseudopotential(
        ae_result=ae,
        tm_dict=tm_dict,
        inv_dict=inv_dict,
        validation_report=report,
        output_prefix=str(output_prefix),
        kb_result=kb,
        formats=["json", "npz"],
    )
    for f in files:
        print(f"  wrote: {f}")


if __name__ == "__main__":
    main()
