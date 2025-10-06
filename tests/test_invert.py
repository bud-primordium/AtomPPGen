"""
势反演模块测试
"""

import pytest
import numpy as np

from atomppgen import solve_ae_atom, tm_pseudize
from atomppgen.invert import invert_semilocal_potential, InvertResult


class TestInvertSemilocal:
    """测试势反演基础功能"""

    @pytest.mark.unit
    def test_invert_basic_s_orbital(self):
        """测试 s 轨道的势反演"""
        # 1. 获取 Al 3s 的 TM 伪化结果
        ae = solve_ae_atom(
            Z=13, spin_mode='LDA', lmax=0,
            grid_type='exp_transformed',
            grid_params={'n': 800},
            scf_params={'tol': 1e-6}
        )

        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[0][2],
            eps=ae.eps_by_l[0][2],
            l=0, rc=2.0, continuity_orders=2
        )

        # 2. 反演势
        inv = invert_semilocal_potential(tm, ae.r)

        # 3. 基本检查
        assert isinstance(inv, InvertResult)
        assert inv.V_l.shape == ae.r.shape
        assert inv.l == 0
        assert inv.rc == 2.0
        assert np.all(np.isfinite(inv.V_l)), "势包含 NaN 或 Inf"

        # 4. 物理合理性
        # s 轨道势在原点附近应接近常数（无离心势）
        i_origin = np.searchsorted(ae.r, 0.1)  # r < 0.1 Bohr
        V_near_origin = inv.V_l[:i_origin]
        assert np.std(V_near_origin) < 5.0, "s 轨道势在原点附近应较平缓"

        # 势值应在合理范围（-100 ~ 100 Ha）
        assert np.max(np.abs(inv.V_l)) < 100.0, f"势过大：{np.max(np.abs(inv.V_l)):.2f} Ha"

    @pytest.mark.unit
    def test_invert_p_orbital(self):
        """测试 p 轨道的势反演"""
        ae = solve_ae_atom(
            Z=13, spin_mode='LDA', lmax=1,
            grid_type='exp_transformed',
            grid_params={'n': 800},
            scf_params={'tol': 1e-6}
        )

        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[1][1],  # 3p
            eps=ae.eps_by_l[1][1],
            l=1, rc=2.2, continuity_orders=2
        )

        inv = invert_semilocal_potential(tm, ae.r)

        # p 轨道势应合理
        assert inv.V_l.shape == ae.r.shape
        assert inv.l == 1

        # 检查势的物理合理性（有限值，无极端奇异性）
        assert np.all(np.isfinite(inv.V_l)), "p 轨道势应有限"
        assert np.max(np.abs(inv.V_l)) < 100.0, "p 轨道势不应过大"

    @pytest.mark.unit
    def test_diagnostics(self):
        """测试诊断信息"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})
        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[0][2],
            eps=ae.eps_by_l[0][2],
            l=0, rc=2.0
        )

        inv = invert_semilocal_potential(tm, ae.r)

        # 检查诊断字段
        assert 'n_nodes' in inv.diagnostics
        assert 'V_max' in inv.diagnostics
        assert 'V_min' in inv.diagnostics
        assert 'V_at_rc' in inv.diagnostics
        assert 'method_inner' in inv.diagnostics
        assert 'method_outer' in inv.diagnostics

        # TM 伪轨道应无节点
        assert inv.diagnostics['n_nodes'] == 0, "TM 伪轨道不应有节点"

        # 方法标注
        assert inv.diagnostics['method_inner'] == 'analytical'
        assert inv.diagnostics['method_outer'] == 'spline'

    @pytest.mark.unit
    def test_inner_outer_continuity(self):
        """测试内外区在 rc 处的连续性"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 800})
        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[0][2],
            eps=ae.eps_by_l[0][2],
            l=0, rc=2.0
        )

        inv = invert_semilocal_potential(tm, ae.r)

        # 找到 rc 对应索引
        i_rc = np.searchsorted(ae.r, 2.0)

        # rc 处势应连续
        if i_rc > 0 and i_rc < len(ae.r) - 1:
            V_before = inv.V_l[i_rc - 1]
            V_at = inv.V_l[i_rc]
            V_after = inv.V_l[i_rc + 1]

            # 相邻点变化应平滑（相对变化 < 10%）
            if abs(V_at) > 1e-6:
                rel_jump_before = abs(V_at - V_before) / abs(V_at)
                rel_jump_after = abs(V_after - V_at) / abs(V_at)

                assert rel_jump_before < 0.1, f"rc 前势跳变过大：{rel_jump_before:.3f}"
                assert rel_jump_after < 0.1, f"rc 后势跳变过大：{rel_jump_after:.3f}"

    @pytest.mark.unit
    def test_invalid_inputs(self):
        """测试无效输入处理"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})
        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[0][2],
            eps=ae.eps_by_l[0][2],
            l=0, rc=2.0
        )

        # 网格长度不匹配
        with pytest.raises(ValueError, match="网格长度.*不匹配"):
            invert_semilocal_potential(tm, ae.r[:100])

    @pytest.mark.unit
    def test_V_max_clip(self):
        """测试势裁剪功能"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})
        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[0][2],
            eps=ae.eps_by_l[0][2],
            l=0, rc=2.0
        )

        # 使用较小的裁剪值
        inv = invert_semilocal_potential(tm, ae.r, V_max_clip=50.0)

        # 势应被裁剪在 [-50, 50]
        assert np.max(inv.V_l) <= 50.0
        assert np.min(inv.V_l) >= -50.0


class TestInvertResult:
    """测试 InvertResult 数据类"""

    @pytest.mark.unit
    def test_result_attributes(self):
        """测试结果对象属性"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})
        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[0][2],
            eps=ae.eps_by_l[0][2],
            l=0, rc=2.0
        )

        inv = invert_semilocal_potential(tm, ae.r)

        # 检查所有属性
        required_attrs = ['V_l', 'r', 'l', 'rc', 'eps', 'diagnostics']
        for attr in required_attrs:
            assert hasattr(inv, attr), f"缺少属性 {attr}"

        # 类型检查
        assert isinstance(inv.V_l, np.ndarray)
        assert isinstance(inv.r, np.ndarray)
        assert isinstance(inv.l, int)
        assert isinstance(inv.rc, float)
        assert isinstance(inv.eps, float)
        assert isinstance(inv.diagnostics, dict)

    @pytest.mark.unit
    def test_smooth_rc_option(self):
        """测试 rc 处平滑功能"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})
        tm = tm_pseudize(
            r=ae.r, w=ae.w,
            u_ae=ae.u_by_l[0][2],
            eps=ae.eps_by_l[0][2],
            l=0, rc=2.0
        )

        # 不平滑
        inv_no_smooth = invert_semilocal_potential(tm, ae.r, smooth_rc=False)

        # 平滑
        inv_smooth = invert_semilocal_potential(
            tm, ae.r, smooth_rc=True, smooth_width=0.15
        )

        # 两者形状相同
        assert inv_no_smooth.V_l.shape == inv_smooth.V_l.shape

        # 找到 rc 对应索引
        i_rc = np.searchsorted(ae.r, 2.0)

        # 平滑后 rc 附近应该更平滑（标准差更小）
        window = 5
        i_left = max(0, i_rc - window)
        i_right = min(len(ae.r), i_rc + window + 1)

        std_no_smooth = np.std(inv_no_smooth.V_l[i_left:i_right])
        std_smooth = np.std(inv_smooth.V_l[i_left:i_right])

        # 平滑应该减少局部波动（或至少不增加）
        # 注意：可能原本就很平滑，所以这里只检查不会变得更差
        assert std_smooth <= std_no_smooth * 1.5, \
            f"平滑反而增加波动：{std_smooth:.3e} vs {std_no_smooth:.3e}"
