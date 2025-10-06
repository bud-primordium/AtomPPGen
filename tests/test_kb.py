"""
KB 转换模块测试
"""

import pytest
import numpy as np

from atomppgen import solve_ae_atom, tm_pseudize
from atomppgen.invert import invert_semilocal_potential
from atomppgen.kb import kb_transform, KBResult


class TestKBTransform:
    """测试 KB 可分离形式转换"""

    @pytest.mark.unit
    def test_basic_kb_transform_simple(self):
        """测试基础 KB 转换（简化版，只用 s 通道）"""
        # 1. 生成 Al 全电子解（只用 s 通道）
        ae = solve_ae_atom(
            Z=13, spin_mode='LDA', lmax=0,
            grid_type='exp_transformed',
            grid_params={'n': 600},
            scf_params={'tol': 1e-6}
        )

        # 2. 对 s 通道进行两次不同 rc 的 TM 伪化，模拟两个通道
        # 这样可以避免 p 轨道伪化的收敛问题
        invert_results = {}
        u_by_l = {}

        u_ae_s = ae.u_by_l[0][-1]
        eps_s = ae.eps_by_l[0][-1]

        # 模拟 "s 通道"（rc=2.0）
        tm_0 = tm_pseudize(ae.r, ae.w, u_ae_s, eps_s, l=0, rc=2.0, continuity_orders=2)
        inv_0 = invert_semilocal_potential(tm_0, ae.r)
        invert_results[0] = inv_0
        u_by_l[0] = tm_0.u_ps

        # 模拟 "p 通道"（使用同一个 s 轨道但不同 rc=2.2）
        # 注意：这在物理上不正确，但用于测试 KB 转换逻辑足够
        tm_1 = tm_pseudize(ae.r, ae.w, u_ae_s, eps_s, l=0, rc=2.2, continuity_orders=2)
        inv_1 = invert_semilocal_potential(tm_1, ae.r)
        # 手动标记为 l=1 用于测试
        inv_1.l = 1
        invert_results[1] = inv_1
        u_by_l[1] = tm_1.u_ps

        # 3. KB 转换（"p 通道" 作为局域势）
        kb = kb_transform(invert_results, u_by_l, ae.r, ae.w, loc_channel=1)

        # 4. 基本检查
        assert isinstance(kb, KBResult)
        assert kb.loc_channel == 1
        assert kb.V_loc.shape == ae.r.shape
        assert np.all(np.isfinite(kb.V_loc)), "局域势包含 NaN 或 Inf"

        # 5. 投影子检查
        assert len(kb.beta_l) == 1, "应有 1 个投影子"
        assert 0 in kb.beta_l, "缺少通道 0 投影子"
        assert 1 not in kb.beta_l, "局域通道不应有投影子"

        for l, beta in kb.beta_l.items():
            assert beta.shape == ae.r.shape, f"通道 l={l} 投影子形状错误"
            assert np.all(np.isfinite(beta)), f"通道 l={l} 投影子包含 NaN 或 Inf"

        # 6. 耦合系数检查
        assert len(kb.D_l) == 1, "应有 1 个耦合系数"
        for l, D in kb.D_l.items():
            assert D > 0, f"通道 l={l} 耦合系数应为正数，得到 {D:.3e}"
            # 放宽范围，因为这是人工构造的测试
            assert 0.001 < D < 10000, f"通道 l={l} 耦合系数 {D:.3e} 超出合理范围"

        # 7. 诊断信息检查
        assert kb.diagnostics['n_channels'] == 2
        assert 'projector_norms' in kb.diagnostics
        assert 'coupling_strengths' in kb.diagnostics

    @pytest.mark.unit
    def test_projector_normalization(self):
        """测试投影子归一化"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})

        invert_results = {}
        u_by_l = {}

        u_ae_s = ae.u_by_l[0][-1]
        eps_s = ae.eps_by_l[0][-1]

        # 创建两个通道（不同 rc）
        for i, rc in enumerate([2.0, 2.2]):
            tm = tm_pseudize(ae.r, ae.w, u_ae_s, eps_s, l=0, rc=rc)
            inv = invert_semilocal_potential(tm, ae.r)
            inv.l = i  # 手动设置 l
            invert_results[i] = inv
            u_by_l[i] = tm.u_ps

        # 通道 1 作为局域势
        kb = kb_transform(invert_results, u_by_l, ae.r, ae.w, loc_channel=1)

        # 验证通道 0 投影子归一化
        beta_0 = kb.beta_l[0]
        norm = np.sum(beta_0**2 * ae.w)

        assert abs(norm - 1.0) < 1e-6, f"投影子归一化误差过大：{norm:.3e}"

        # 验证耦合系数 D = W/Z
        # 重新计算 W 和 Z 以验证
        V_loc = invert_results[1].V_l
        delta_V = invert_results[0].V_l - V_loc
        r_safe = np.maximum(ae.r, 1e-10)
        phi_0 = u_by_l[0] / r_safe
        chi_0 = delta_V * phi_0

        W = np.sum(chi_0**2 * ae.w)
        Z = np.sum(delta_V * (phi_0**2) * ae.w)

        D_expected = W / Z
        D_actual = kb.D_l[0]

        assert abs(D_actual - D_expected) < 1e-10, \
            f"耦合系数验证失败：D_actual={D_actual:.6e}, D_expected={D_expected:.6e}"

    @pytest.mark.unit
    def test_different_local_channels(self):
        """测试不同局域道选择"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})

        invert_results = {}
        u_by_l = {}

        u_ae_s = ae.u_by_l[0][-1]
        eps_s = ae.eps_by_l[0][-1]

        # 创建两个通道
        for i, rc in enumerate([2.0, 2.2]):
            tm = tm_pseudize(ae.r, ae.w, u_ae_s, eps_s, l=0, rc=rc)
            inv = invert_semilocal_potential(tm, ae.r)
            inv.l = i
            invert_results[i] = inv
            u_by_l[i] = tm.u_ps

        # 尝试通道 0 作为局域势
        kb_0 = kb_transform(invert_results, u_by_l, ae.r, ae.w, loc_channel=0)
        assert kb_0.loc_channel == 0
        assert 0 not in kb_0.beta_l
        assert 1 in kb_0.beta_l

        # 尝试通道 1 作为局域势
        kb_1 = kb_transform(invert_results, u_by_l, ae.r, ae.w, loc_channel=1)
        assert kb_1.loc_channel == 1
        assert 1 not in kb_1.beta_l
        assert 0 in kb_1.beta_l

        # 两种选择都应物理合理
        assert np.all(np.isfinite(kb_0.V_loc))
        assert np.all(np.isfinite(kb_1.V_loc))

    @pytest.mark.unit
    def test_invalid_inputs(self):
        """测试无效输入处理"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})

        u_ae = ae.u_by_l[0][-1]
        eps = ae.eps_by_l[0][-1]

        tm = tm_pseudize(ae.r, ae.w, u_ae, eps, l=0, rc=2.0)
        inv = invert_semilocal_potential(tm, ae.r)

        invert_results = {0: inv}
        u_by_l = {0: tm.u_ps}

        # 局域通道不存在
        with pytest.raises(ValueError, match="局域通道.*不在"):
            kb_transform(invert_results, u_by_l, ae.r, ae.w, loc_channel=2)

        # 网格长度不匹配
        with pytest.raises(ValueError, match="网格长度.*不匹配"):
            kb_transform(invert_results, u_by_l, ae.r[:100], ae.w, loc_channel=0)

        # 缺少波函数
        with pytest.raises(ValueError, match="缺少径向波函数"):
            kb_transform(invert_results, {}, ae.r, ae.w, loc_channel=0)

    @pytest.mark.unit
    def test_diagnostics(self):
        """测试诊断信息"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})

        invert_results = {}
        u_by_l = {}

        u_ae_s = ae.u_by_l[0][-1]
        eps_s = ae.eps_by_l[0][-1]

        for i, rc in enumerate([2.0, 2.2]):
            tm = tm_pseudize(ae.r, ae.w, u_ae_s, eps_s, l=0, rc=rc)
            inv = invert_semilocal_potential(tm, ae.r)
            inv.l = i
            invert_results[i] = inv
            u_by_l[i] = tm.u_ps

        kb = kb_transform(invert_results, u_by_l, ae.r, ae.w, loc_channel=1)

        # 检查诊断字段
        assert 'n_channels' in kb.diagnostics
        assert 'projector_norms' in kb.diagnostics
        assert 'coupling_strengths' in kb.diagnostics
        assert 'loc_potential_max' in kb.diagnostics
        assert 'loc_potential_min' in kb.diagnostics

        # 检查值
        assert kb.diagnostics['n_channels'] == 2
        assert len(kb.diagnostics['projector_norms']) == 1  # 只有一个投影子
        assert 0 in kb.diagnostics['projector_norms']


class TestKBResult:
    """测试 KBResult 数据类"""

    @pytest.mark.unit
    def test_result_attributes(self):
        """测试结果对象属性"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})

        invert_results = {}
        u_by_l = {}

        u_ae_s = ae.u_by_l[0][-1]
        eps_s = ae.eps_by_l[0][-1]

        for i, rc in enumerate([2.0, 2.2]):
            tm = tm_pseudize(ae.r, ae.w, u_ae_s, eps_s, l=0, rc=rc)
            inv = invert_semilocal_potential(tm, ae.r)
            inv.l = i
            invert_results[i] = inv
            u_by_l[i] = tm.u_ps

        kb = kb_transform(invert_results, u_by_l, ae.r, ae.w, loc_channel=1)

        # 检查所有属性
        required_attrs = ['V_loc', 'beta_l', 'D_l', 'loc_channel', 'r', 'diagnostics']
        for attr in required_attrs:
            assert hasattr(kb, attr), f"缺少属性 {attr}"

        # 类型检查
        assert isinstance(kb.V_loc, np.ndarray)
        assert isinstance(kb.beta_l, dict)
        assert isinstance(kb.D_l, dict)
        assert isinstance(kb.loc_channel, int)
        assert isinstance(kb.r, np.ndarray)
        assert isinstance(kb.diagnostics, dict)
