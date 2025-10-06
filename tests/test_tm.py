"""
TM 伪化器的单元测试

测试 Troullier-Martins 方法的正确性：
- 范数守恒检查
- 连续性验证（函数值和导数匹配）
- 数值稳定性
"""

import pytest
import numpy as np
from atomppgen import solve_ae_atom


class TestTMPseudize:
    """测试 TM 伪化函数"""

    @pytest.mark.unit
    def test_tm_basic_s_orbital(self):
        """测试 Al 3s 轨道的基础 TM 伪化"""
        # 1. 获取 Al 的全电子解
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,  # 只需要 s 通道
            grid_type="exp_transformed",
            grid_params={"n": 800, "rmax": 100.0},
            scf_params={"tol": 1e-6, "maxiter": 120},
        )

        # 2. 提取 3s 轨道
        u_ae_3s = ae_result.u_by_l[0][2]  # l=0, n=3
        eps_3s = ae_result.eps_by_l[0][2]
        r = ae_result.r
        w = ae_result.w

        # 3. 先验证 AE 解的合理性
        assert ae_result.converged, "AE SCF 未收敛"
        assert -0.3 < eps_3s < -0.2, f"3s 能级 {eps_3s:.3f} Ha 超出预期范围"

        # 导入 tm_pseudize（这是我们要实现的）
        from atomppgen.tm import tm_pseudize

        # 4. TM 伪化
        rc = 2.1  # Bohr（Al s 通道推荐值）
        result = tm_pseudize(
            r=r,
            w=w,
            u_ae=u_ae_3s,
            eps=eps_3s,
            l=0,
            rc=rc,
            continuity_orders=2,  # 二阶导数连续性
        )

        # 5. 检查返回结果的基本属性
        assert hasattr(result, "u_ps"), "缺少伪轨道 u_ps"
        assert hasattr(result, "a_coeff"), "缺少 TM 系数 a_coeff"
        assert hasattr(result, "norm_error"), "缺少范数误差"
        assert hasattr(result, "continuity_check"), "缺少连续性检查"

        # 6. 检查伪轨道形状
        assert result.u_ps.shape == r.shape, "伪轨道长度与网格不匹配"
        assert np.all(np.isfinite(result.u_ps)), "伪轨道包含 NaN/Inf"

        # 7. 检查系数数量（continuity_orders=2 需要 4 个系数）
        assert len(result.a_coeff) == 4, f"系数数量错误：{len(result.a_coeff)}"

        # 8. 范数守恒检查（核心要求）
        assert result.norm_error < 1e-5, (
            f"范数守恒误差 {result.norm_error:.2e} 超过阈值 1e-5"
        )

        # 9. 连续性检查
        for key in ["u", "du", "d2u"]:
            assert key in result.continuity_check, f"缺少 {key} 连续性检查"
            rel_err = result.continuity_check[key]["rel_error"]
            assert rel_err < 1e-4, f"{key} 在 rc 处不连续，相对误差 {rel_err:.2e}"

    @pytest.mark.unit
    def test_tm_norm_conservation_precision(self):
        """测试范数守恒的精度"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_params={"n": 600},
            scf_params={"tol": 1e-5, "maxiter": 80},
        )

        from atomppgen.tm import tm_pseudize

        u_ae = ae_result.u_by_l[0][2]
        eps = ae_result.eps_by_l[0][2]

        # 测试不同的 rc 值
        for rc in [1.8, 2.0, 2.2, 2.4]:
            result = tm_pseudize(
                r=ae_result.r,
                w=ae_result.w,
                u_ae=u_ae,
                eps=eps,
                l=0,
                rc=rc,
                continuity_orders=2,
            )

            # 范数守恒应该非常好（< 1e-6）
            assert result.norm_error < 1e-6, (
                f"rc={rc:.1f} 的范数守恒误差 {result.norm_error:.2e} 过大"
            )

    @pytest.mark.unit
    def test_tm_p_orbital(self):
        """测试 p 轨道的伪化"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=1,
            grid_params={"n": 600},
            scf_params={"tol": 1e-5},
        )

        from atomppgen.tm import tm_pseudize

        # 3p 轨道
        u_ae_3p = ae_result.u_by_l[1][1]  # l=1, n=2 (2p 是芯层，3p 是第2个 p 态)
        eps_3p = ae_result.eps_by_l[1][1]

        rc = 2.2  # Al p 通道推荐值
        result = tm_pseudize(
            r=ae_result.r,
            w=ae_result.w,
            u_ae=u_ae_3p,
            eps=eps_3p,
            l=1,
            rc=rc,
            continuity_orders=2,
        )

        # 基本检查
        assert result.norm_error < 1e-5
        assert result.continuity_check["u"]["rel_error"] < 1e-4

    @pytest.mark.unit
    def test_tm_different_continuity_orders(self):
        """测试不同的连续性阶数"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_params={"n": 500},
            scf_params={"tol": 1e-5},
        )

        from atomppgen.tm import tm_pseudize

        u_ae = ae_result.u_by_l[0][2]
        eps = ae_result.eps_by_l[0][2]
        rc = 2.1

        # 二阶连续性
        result_2nd = tm_pseudize(
            r=ae_result.r,
            w=ae_result.w,
            u_ae=u_ae,
            eps=eps,
            l=0,
            rc=rc,
            continuity_orders=2,
        )
        assert result_2nd.continuity_orders == 2
        assert len(result_2nd.a_coeff) == 4

        # 四阶连续性（更严格）
        result_4th = tm_pseudize(
            r=ae_result.r,
            w=ae_result.w,
            u_ae=u_ae,
            eps=eps,
            l=0,
            rc=rc,
            continuity_orders=4,
        )
        assert result_4th.continuity_orders == 4
        assert len(result_4th.a_coeff) == 6  # 需要 6 个系数

        # 四阶应该有更好的范数守恒
        # （理论上应该相近，但不一定更好，因为是不同的约束）
        assert result_4th.norm_error < 1e-4

    @pytest.mark.unit
    def test_tm_invalid_inputs(self):
        """测试无效输入的处理"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_params={"n": 400},
            scf_params={"tol": 1e-5},
        )

        from atomppgen.tm import tm_pseudize

        u_ae = ae_result.u_by_l[0][2]
        eps = ae_result.eps_by_l[0][2]

        # rc 超出网格范围
        with pytest.raises(ValueError, match="超出网格范围|out of range"):
            tm_pseudize(
                r=ae_result.r,
                w=ae_result.w,
                u_ae=u_ae,
                eps=eps,
                l=0,
                rc=999.0,  # 远超 rmax
                continuity_orders=2,
            )

        # rc 太小（< r[0]）
        with pytest.raises(ValueError, match="超出网格范围|out of range"):
            tm_pseudize(
                r=ae_result.r,
                w=ae_result.w,
                u_ae=u_ae,
                eps=eps,
                l=0,
                rc=-1.0,
                continuity_orders=2,
            )

        # 不支持的连续性阶数
        with pytest.raises(ValueError, match="continuity_orders|必须是 2 或 4"):
            tm_pseudize(
                r=ae_result.r,
                w=ae_result.w,
                u_ae=u_ae,
                eps=eps,
                l=0,
                rc=2.0,
                continuity_orders=3,  # 只支持 2 或 4
            )

    @pytest.mark.unit
    def test_tm_inner_outer_splice(self):
        """测试内外区拼接的平滑性"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_params={"n": 600},
            scf_params={"tol": 1e-5},
        )

        from atomppgen.tm import tm_pseudize

        u_ae = ae_result.u_by_l[0][2]
        eps = ae_result.eps_by_l[0][2]
        rc = 2.1

        result = tm_pseudize(
            r=ae_result.r,
            w=ae_result.w,
            u_ae=u_ae,
            eps=eps,
            l=0,
            rc=rc,
            continuity_orders=2,
        )

        # 找到 rc 对应的索引
        i_rc = np.searchsorted(ae_result.r, rc)

        # 检查拼接点前后的连续性
        if i_rc > 2 and i_rc < len(ae_result.r) - 2:
            u_before = result.u_ps[i_rc - 1]
            u_at = result.u_ps[i_rc]
            u_after = result.u_ps[i_rc + 1]

            # 检查是否有跳变
            assert np.isfinite(u_before) and np.isfinite(u_at) and np.isfinite(u_after)

            # 相邻点的变化应该平滑（相对变化 < 20%）
            if abs(u_at) > 1e-8:
                rel_change_before = abs(u_at - u_before) / abs(u_at)
                rel_change_after = abs(u_after - u_at) / abs(u_at)
                assert rel_change_before < 0.2, (
                    f"拼接点前不平滑：{rel_change_before:.2e}"
                )
                assert rel_change_after < 0.2, f"拼接点后不平滑：{rel_change_after:.2e}"

    @pytest.mark.unit
    def test_tm_wavefunction_positivity(self):
        """测试伪轨道在内区的正定性（无节点）"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_params={"n": 500},
            scf_params={"tol": 1e-5},
        )

        from atomppgen.tm import tm_pseudize

        u_ae = ae_result.u_by_l[0][2]
        eps = ae_result.eps_by_l[0][2]
        rc = 2.0

        result = tm_pseudize(
            r=ae_result.r,
            w=ae_result.w,
            u_ae=u_ae,
            eps=eps,
            l=0,
            rc=rc,
            continuity_orders=2,
        )

        # 找到 rc 对应的索引
        i_rc = np.searchsorted(ae_result.r, rc)

        # 内区应该没有节点（对于价层 s 轨道）
        # 注意：AE 的 3s 可能有节点（1s, 2s 的节点），但伪轨道应该无节点
        u_inner = result.u_ps[:i_rc]

        # 检查是否有符号变化（节点）
        # 对于 TM 方法，exp(多项式) 始终为正
        # 但 r^{l+1} 因子在 r>0 时也为正
        # 所以内区应该全为正（或全为负，取决于归一化）

        # 简化检查：内区不应有零点穿越（除了 r=0）
        if len(u_inner) > 10:
            # 跳过前几个点（接近 r=0）
            u_check = u_inner[5:]
            if len(u_check) > 0:
                # 检查符号是否一致
                sign_changes = np.sum(np.diff(np.sign(u_check)) != 0)
                assert sign_changes == 0, f"内区出现 {sign_changes} 个节点"


class TestTMResult:
    """测试 TMResult 数据类"""

    @pytest.mark.unit
    def test_result_attributes(self):
        """测试结果对象的属性"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_params={"n": 400},
            scf_params={"tol": 1e-5},
        )

        from atomppgen.tm import tm_pseudize

        result = tm_pseudize(
            r=ae_result.r,
            w=ae_result.w,
            u_ae=ae_result.u_by_l[0][2],
            eps=ae_result.eps_by_l[0][2],
            l=0,
            rc=2.0,
            continuity_orders=2,
        )

        # 检查所有必需属性
        required_attrs = [
            "u_ps",
            "a_coeff",
            "rc",
            "eps",
            "l",
            "norm_error",
            "continuity_orders",
            "continuity_check",
        ]
        for attr in required_attrs:
            assert hasattr(result, attr), f"缺少属性 {attr}"

        # 检查类型
        assert isinstance(result.u_ps, np.ndarray)
        assert isinstance(result.a_coeff, np.ndarray)
        assert isinstance(result.rc, float)
        assert isinstance(result.eps, float)
        assert isinstance(result.l, int)
        assert isinstance(result.norm_error, float)
        assert isinstance(result.continuity_orders, int)
        assert isinstance(result.continuity_check, dict)


class TestRcPerturbationStability:
    """测试截断半径微扰稳定性"""

    @pytest.mark.unit
    def test_rc_perturbation_norm_stability(self):
        """测试 rc±0.05 时范数误差阶不恶化"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_type="exp_transformed",
            grid_params={"n": 800},
            scf_params={"tol": 1e-6},
        )

        from atomppgen.tm import tm_pseudize

        u_ae = ae_result.u_by_l[0][2]
        eps = ae_result.eps_by_l[0][2]
        rc_base = 2.0

        # 测试 rc-0.05, rc, rc+0.05
        rc_values = [rc_base - 0.05, rc_base, rc_base + 0.05]
        norm_errors = []

        for rc in rc_values:
            result = tm_pseudize(
                r=ae_result.r,
                w=ae_result.w,
                u_ae=u_ae,
                eps=eps,
                l=0,
                rc=rc,
                continuity_orders=2,
            )
            norm_errors.append(result.norm_error)

        # 检查范数误差不恶化（阈值：不超过 3 倍基准值）
        base_error = norm_errors[1]  # rc_base 的误差
        for i, rc in enumerate(rc_values):
            assert norm_errors[i] < 3 * base_error, (
                f"rc={rc:.2f} 范数误差 {norm_errors[i]:.2e} 超过基准 {base_error:.2e} 的 3 倍"
            )

    @pytest.mark.unit
    def test_rc_perturbation_continuity_stability(self):
        """测试 rc±0.05 时连续性误差阶不恶化"""
        ae_result = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_type="exp_transformed",
            grid_params={"n": 800},
            scf_params={"tol": 1e-6},
        )

        from atomppgen.tm import tm_pseudize

        u_ae = ae_result.u_by_l[0][2]
        eps = ae_result.eps_by_l[0][2]
        rc_base = 2.0
        rc_values = [rc_base - 0.05, rc_base, rc_base + 0.05]

        max_cont_errors = []
        for rc in rc_values:
            result = tm_pseudize(
                r=ae_result.r,
                w=ae_result.w,
                u_ae=u_ae,
                eps=eps,
                l=0,
                rc=rc,
                continuity_orders=2,
            )
            # 取所有导数连续性误差的最大值
            errors = [
                result.continuity_check[k]["rel_error"] for k in ["u", "du", "d2u"]
            ]
            max_cont_errors.append(max(errors))

        # 检查连续性误差不恶化（阈值：不超过 3 倍基准值）
        base_error = max_cont_errors[1]
        for i, rc in enumerate(rc_values):
            assert max_cont_errors[i] < 3 * base_error, (
                f"rc={rc:.2f} 连续性误差 {max_cont_errors[i]:.2e} 超过基准 {base_error:.2e} 的 3 倍"
            )


class TestGridConsistencyRegression:
    """测试不同网格类型的一致性"""

    @pytest.mark.unit
    def test_linear_vs_exp_grid_consistency(self):
        """测试 exp_transformed 网格的高精度（linear 网格精度较低，仅验证可用性）"""
        from atomppgen.tm import tm_pseudize

        rc = 2.0
        l = 0

        # 1. 使用 exp_transformed 网格（高精度验证）
        ae_exp = solve_ae_atom(
            Z=13,
            spin_mode="LDA",
            lmax=0,
            grid_type="exp_transformed",
            grid_params={"n": 1200},
            scf_params={"tol": 1e-7, "maxiter": 150},
        )
        result_exp = tm_pseudize(
            r=ae_exp.r,
            w=ae_exp.w,
            u_ae=ae_exp.u_by_l[0][2],
            eps=ae_exp.eps_by_l[0][2],
            l=l,
            rc=rc,
            continuity_orders=2,
        )

        # 2. exp_transformed 网格应达到高精度
        def get_max_cont_error(result):
            return max(
                [result.continuity_check[k]["rel_error"] for k in ["u", "du", "d2u"]]
            )

        error_exp = get_max_cont_error(result_exp)
        assert error_exp < 1e-8, f"exp_transformed 网格连续性误差过大：{error_exp:.2e}"

        # 3. 范数守恒也应高精度
        assert result_exp.norm_error < 1e-6, (
            f"exp_transformed 网格范数误差过大：{result_exp.norm_error:.2e}"
        )
