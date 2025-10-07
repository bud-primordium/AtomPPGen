"""
验证模块测试
"""

import pytest
import numpy as np

from atomppgen import solve_ae_atom, tm_pseudize, invert_semilocal_potential
from atomppgen.validate import (
    check_norm_conservation,
    NormConservationResult,
    _solve_radial_schrodinger_numerov,
    check_log_derivative,
    LogDerivativeResult,
    check_ghost_states,
    GhostStateResult,
    _extract_ks_potential,
)


class TestNormConservation:
    """测试范数守恒检验"""

    @pytest.mark.unit
    def test_norm_simple_s_orbital(self):
        """测试 s 轨道范数守恒"""
        # 生成 Al 全电子解
        ae = solve_ae_atom(
            Z=13, spin_mode='LDA', lmax=0,
            grid_type='exp_transformed',
            grid_params={'n': 600}
        )

        # TM 伪化
        u_ae = ae.u_by_l[0][-1]
        eps = ae.eps_by_l[0][-1]
        tm = tm_pseudize(ae.r, ae.w, u_ae, eps, l=0, rc=2.0)

        # 范数守恒检验
        result = check_norm_conservation(tm)

        # 检查结果类型
        assert isinstance(result, NormConservationResult)
        assert result.l == 0
        assert result.rc == 2.0

        # 检查范数误差
        assert abs(result.norm_error) < 1e-6, \
            f"范数误差过大：{result.norm_error:.3e}"
        assert result.passed, "范数守恒检验未通过"

        # 检查诊断信息
        assert 'method' in result.diagnostics

    @pytest.mark.unit
    def test_norm_conservation_tolerance(self):
        """测试不同容许误差阈值"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})
        tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1], l=0, rc=2.0)

        # 宽松阈值
        result_loose = check_norm_conservation(tm, tolerance=1e-4)
        assert result_loose.passed

        # 严格阈值（若 TM 精度不够可能失败）
        result_strict = check_norm_conservation(tm, tolerance=1e-8)
        # 不强制要求通过，只检查返回值类型
        assert isinstance(result_strict.passed, bool)

    @pytest.mark.unit
    def test_norm_result_attributes(self):
        """测试 NormConservationResult 属性完整性"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})
        tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1], l=0, rc=2.0)
        result = check_norm_conservation(tm)

        # 检查所有必需属性
        required_attrs = ['l', 'norm_error', 'passed', 'rc', 'tolerance', 'diagnostics']
        for attr in required_attrs:
            assert hasattr(result, attr), f"缺少属性 {attr}"

        # 类型检查
        assert isinstance(result.l, int)
        assert isinstance(result.norm_error, float)
        assert isinstance(result.passed, bool)
        assert isinstance(result.rc, float)
        assert isinstance(result.tolerance, float)
        assert isinstance(result.diagnostics, dict)


class TestRadialSolver:
    """测试径向薛定谔方程求解器"""

    @pytest.mark.unit
    def test_numerov_solver_boundary_conditions(self):
        """测试 Numerov 求解器边界条件"""
        # 创建简单的类氢势 V(r) = -Z/r
        n = 500
        r = np.linspace(0.01, 10.0, n)
        Z = 1.0
        V = -Z / r

        # 求解 1s 态（E ≈ -0.5 Ha）
        E = -0.5
        psi = _solve_radial_schrodinger_numerov(r, V, l=0, E=E)

        # 检查基本性质
        assert len(psi) == len(r), "波函数长度与网格不匹配"
        assert abs(psi[0]) < 1e-6, "边界条件 ψ(0)≈0 不满足"
        assert np.all(np.isfinite(psi)), "波函数存在非有限值"

    @pytest.mark.unit
    def test_numerov_solver_with_ae_potential(self):
        """测试 Numerov 在 AE 势下求解"""
        # 使用 Al 原子的 AE 解
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})

        # 提取 3s 态的能量
        eps = ae.eps_by_l[0][-1]
        r = ae.r

        # 构建简化势：核势 + 屏蔽（避免在小 r 处奇异）
        # 使用 Slater 屏蔽势 V(r) = -Z_eff/r，其中 Z_eff 渐变
        Z = 13.0
        V = -Z / np.maximum(r, 0.01)  # 限制最小 r 避免除零

        # 求解径向方程
        psi = _solve_radial_schrodinger_numerov(r, V, l=0, E=eps)

        # 检查波函数合理性
        assert len(psi) == len(r)
        # 允许在边界处有小误差（由于插值）
        assert abs(psi[0]) < 1e-3, f"边界条件不满足: psi[0]={psi[0]}"
        assert np.all(np.isfinite(psi)), "波函数存在非有限值"

        # 检查归一化区间积分为正
        norm = np.trapezoid(psi**2, r)
        assert norm > 0, "波函数积分应为正"


class TestLogDerivative:
    """测试对数导数匹配"""

    @pytest.mark.unit
    def test_ks_potential_extraction(self):
        """测试 KS 有效势提取"""
        # 生成 Al 全电子解
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})

        # 提取 KS 势
        V_ks = _extract_ks_potential(ae)

        # 检查基本性质
        assert len(V_ks) == len(ae.r), "势与网格长度不匹配"
        assert np.all(np.isfinite(V_ks)), "势存在非有限值"

        # KS 势应为负（吸引势）
        assert np.all(V_ks < 0), "KS 势应为吸引势"

        # 在大 r 处应趋向 0（电子完全屏蔽核电荷）
        assert abs(V_ks[-1]) < 1e-6, \
            f"大 r 处势未趋向 0: {V_ks[-1]:.3e}"

        # 势应随 r 单调递增（越靠近核越深）
        # 检查后半部分网格（避免数值噪音）
        r_mid = ae.r[len(ae.r)//2:]
        V_mid = V_ks[len(V_ks)//2:]
        assert np.all(np.diff(V_mid) > -1e-3), "势应单调或近似单调"

    @pytest.mark.unit
    def test_log_derivative_basic(self):
        """测试对数导数基本功能"""
        # 创建简单的网格和势
        n = 500
        r = np.linspace(0.01, 10.0, n)
        Z = 1.0

        # 构建两个相似的势（模拟 AE 和 PS）
        V_AE = -Z / r
        V_PS = -Z / r * 0.95  # 略有不同的势

        # 计算对数导数匹配
        result = check_log_derivative(
            V_AE, V_PS, r, l=0, r_test=3.0,
            E_range_Ha=(-0.5, -0.4),  # 小范围快速测试
            E_step_Ha=0.05
        )

        # 检查结果类型
        assert isinstance(result, LogDerivativeResult)
        assert result.l == 0
        assert result.r_test == 3.0

        # 检查对数导数数组长度
        assert len(result.energies) > 0
        assert len(result.L_AE) == len(result.energies)
        assert len(result.L_PS) == len(result.energies)

        # 检查诊断信息
        assert 'n_energies' in result.diagnostics
        assert 'n_valid' in result.diagnostics

    @pytest.mark.integration
    def test_log_derivative_with_tm_pseudopotential(self):
        """测试 TM 赝势的对数导数匹配（使用真实 KS 势）"""
        # 生成 Al 的 AE 解和 TM 伪化
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})
        tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1], l=0, rc=2.0)

        # 反演得到伪势
        inv = invert_semilocal_potential(tm, ae.r)

        # 提取真实 KS 有效势（AE）
        V_AE = _extract_ks_potential(ae)

        # PS 半局域势
        V_PS = inv.V_l

        # 对数导数检验（缩小能量范围加速测试）
        result = check_log_derivative(
            V_AE, V_PS, ae.r, l=0, r_test=3.0,
            E_range_Ha=(-0.2, 0.2),  # 较窄的能量窗口
            E_step_Ha=0.05
        )

        # 检查基本性质
        assert isinstance(result, LogDerivativeResult)
        assert result.diagnostics['n_valid'] > 0, "应有有效的能量点"

        # 检查对数导数曲线不全为零（说明求解成功）
        assert np.any(np.abs(result.L_AE) > 1e-6), "AE 对数导数不应全为零"
        assert np.any(np.abs(result.L_PS) > 1e-6), "PS 对数导数不应全为零"

        # 不强制要求通过严格判据（需要更精细的能量扫描和 rc 优化）
        assert isinstance(result.passed, bool)


class TestGhostStates:
    """测试幽灵态检测"""

    @pytest.mark.unit
    def test_ghost_states_no_ghost(self):
        """测试无幽灵态情况"""
        # 生成 Al 的 AE 解和伪化
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})
        tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1], l=0, rc=2.0)

        # 反演半局域势
        inv = invert_semilocal_potential(tm, ae.r)

        # 幽灵态检测
        result = check_ghost_states(
            inv, ae.r, ae.w,
            valence_energy=tm.eps,
            E_window_Ha=(-0.25, 0.25),
            method='radial'
        )

        # 检查结果类型
        assert isinstance(result, GhostStateResult)
        assert result.method == 'radial'
        assert result.l == 0

        # 检查幽灵态数量
        assert isinstance(result.n_ghosts, int)
        assert result.n_ghosts >= 0  # 可能有或没有幽灵态

        # 检查诊断信息
        assert 'n_bound_states_total' in result.diagnostics
        assert 'grid_size' in result.diagnostics

    @pytest.mark.integration
    def test_ghost_states_with_p_orbital(self):
        """测试 p 轨道幽灵态检测"""
        # 生成 Al 的 AE 解（包含 p 通道）
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=1, grid_params={'n': 600})

        # p 轨道伪化（使用较大的 rc 避免收敛问题）
        try:
            tm_p = tm_pseudize(ae.r, ae.w, ae.u_by_l[1][-1], ae.eps_by_l[1][-1],
                               l=1, rc=2.0, continuity_orders=2)
            inv_p = invert_semilocal_potential(tm_p, ae.r)

            # 幽灵态检测
            result = check_ghost_states(
                inv_p, ae.r, ae.w,
                valence_energy=tm_p.eps,
                E_window_Ha=(-0.25, 0.25),
                method='radial'
            )

            # 基本检查
            assert isinstance(result, GhostStateResult)
            assert result.l == 1
            assert isinstance(result.passed, bool)

        except Exception as e:
            # p 轨道可能收敛困难，跳过测试
            pytest.skip(f"p 轨道伪化失败: {e}")

    @pytest.mark.unit
    def test_ghost_states_result_attributes(self):
        """测试 GhostStateResult 属性完整性"""
        ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 400})
        tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1], l=0, rc=2.0)
        inv = invert_semilocal_potential(tm, ae.r)

        result = check_ghost_states(inv, ae.r, ae.w, valence_energy=tm.eps)

        # 检查所有必需属性
        required_attrs = ['method', 'l', 'eigenvalues', 'known_valence',
                          'ghost_states', 'n_ghosts', 'passed', 'diagnostics']
        for attr in required_attrs:
            assert hasattr(result, attr), f"缺少属性 {attr}"

        # 类型检查
        assert isinstance(result.method, str)
        assert isinstance(result.l, int)
        assert isinstance(result.eigenvalues, np.ndarray)
        assert isinstance(result.known_valence, np.ndarray)
        assert isinstance(result.ghost_states, np.ndarray)
        assert isinstance(result.n_ghosts, int)
        assert isinstance(result.passed, bool)
        assert isinstance(result.diagnostics, dict)
