"""
导出模块的单元测试

测试赝势导出功能的正确性：
- JSON 格式导出（元数据、验证报告）
- NPZ 格式导出（数值数组、命名约定）
- 统一导出接口
- L-channel 一致性检查
"""

import pytest
import numpy as np
import json
from atomppgen import solve_ae_atom, tm_pseudize, invert_semilocal_potential
from atomppgen.export import export_pseudopotential, PseudopotentialData
from atomppgen.validate import run_full_validation


@pytest.fixture
def al_s_data():
    """
    提供 Al s 通道的测试数据（全电子 + TM + 反演）

    使用参考参数：rc=2.5 Bohr
    """
    # 1. 全电子解（使用较小网格加速测试）
    ae = solve_ae_atom(
        Z=13,
        spin_mode="LDA",
        lmax=0,  # 仅 s 通道
        grid_type="exp_transformed",
        grid_params={"n": 400, "rmax": 80.0},
        scf_params={"tol": 1e-6, "maxiter": 100},
    )

    # 2. TM 伪化
    tm_s = tm_pseudize(
        ae.r, ae.w,
        ae.u_by_l[0][-1],  # 3s 价电子
        ae.eps_by_l[0][-1],
        l=0,
        rc=2.5,
    )

    # 3. 势反演
    inv_s = invert_semilocal_potential(tm_s, ae.r)

    # 4. 验证（简化参数）
    report = run_full_validation(
        ae,
        tm_dict={0: tm_s},
        inv_dict={0: inv_s},
        r_test=3.0,
        E_range_Ry=(-0.5, 0.5),
        E_step_Ry=0.1,  # 粗网格加速
    )

    return {
        'ae': ae,
        'tm_dict': {0: tm_s},
        'inv_dict': {0: inv_s},
        'report': report,
    }


class TestPseudopotentialData:
    """测试 PseudopotentialData 数据类"""

    @pytest.mark.unit
    def test_dataclass_creation(self, al_s_data):
        """测试数据类的基本创建"""
        data = al_s_data

        # 手动构建 PseudopotentialData
        pp_data = PseudopotentialData(
            Z=13,
            symbol='Al',
            spin_mode='LDA',
            xc_functional='PZ81',
            generation_params={'rc_by_l': {0: 2.5}},
            radial_grid=data['ae'].r,
            radial_weights=data['ae'].w,
            ae_eigenvalues_by_l=data['ae'].eps_by_l,
            ae_wavefunctions_by_l=data['ae'].u_by_l,
            pseudo_wavefunctions_by_l={0: data['tm_dict'][0].u_ps},
            semilocal_potentials_by_l={0: data['inv_dict'][0].V_l},
            validation_report=data['report'],
            generation_date='2025-10-09T00:00:00',
            code_version='0.1.0',
        )

        # 基本属性检查
        assert pp_data.Z == 13
        assert pp_data.symbol == 'Al'
        assert pp_data.spin_mode == 'LDA'
        assert pp_data.xc_functional == 'PZ81'
        assert len(pp_data.radial_grid) == len(data['ae'].r)
        assert 0 in pp_data.ae_eigenvalues_by_l
        assert 0 in pp_data.pseudo_wavefunctions_by_l
        assert 0 in pp_data.semilocal_potentials_by_l


class TestJSONExport:
    """测试 JSON 导出功能"""

    @pytest.mark.unit
    def test_json_export_basic(self, al_s_data, tmp_path):
        """测试基本 JSON 导出功能"""
        data = al_s_data

        # 导出到临时目录
        output_prefix = tmp_path / "test_al"
        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['json'],
        )

        # 检查文件生成
        assert len(files) == 1
        json_file = files[0]
        assert json_file.exists()
        assert json_file.suffix == '.json'

        # 检查文件可读取
        with json_file.open('r') as f:
            json_data = json.load(f)

        # 检查顶层结构
        assert 'metadata' in json_data
        assert 'generation_params' in json_data
        assert 'all_electron_reference' in json_data
        assert 'pseudopotential' in json_data
        assert 'validation_report' in json_data
        assert 'radial_grid' in json_data

    @pytest.mark.unit
    def test_json_metadata_content(self, al_s_data, tmp_path):
        """测试 JSON 元数据内容的正确性"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['json'],
        )

        with files[0].open('r') as f:
            json_data = json.load(f)

        # 检查元素信息
        assert json_data['metadata']['element']['Z'] == 13
        assert json_data['metadata']['element']['symbol'] == 'Al'

        # 检查 XC 泛函
        assert json_data['metadata']['xc']['functional'] == 'PZ81'
        assert json_data['metadata']['xc']['spin_mode'] == 'LDA'

        # 检查单位声明
        assert json_data['metadata']['units'] == 'Hartree_atomic'

        # 检查代码信息
        assert json_data['metadata']['code']['name'] == 'AtomPPGen'
        assert json_data['metadata']['code']['version'] == '0.1.0'

    @pytest.mark.unit
    def test_json_validation_report(self, al_s_data, tmp_path):
        """测试 JSON 中的验证报告"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['json'],
        )

        with files[0].open('r') as f:
            json_data = json.load(f)

        report = json_data['validation_report']

        # 检查范数守恒结果
        assert 'norm_results' in report
        assert '0' in report['norm_results']
        assert report['norm_results']['0']['l'] == 0
        assert 'norm_error' in report['norm_results']['0']
        assert 'passed' in report['norm_results']['0']

        # 检查对数导数结果
        assert 'log_deriv_results' in report
        assert '0' in report['log_deriv_results']
        assert 'zero_crossing_rms' in report['log_deriv_results']['0']
        assert 'curve_rms_valence' in report['log_deriv_results']['0']

        # 检查幽灵态结果
        assert 'ghost_result' in report
        assert 'n_ghosts' in report['ghost_result']
        assert 'n_box_states' in report['ghost_result']
        assert 'passed' in report['ghost_result']

        # 检查总体通过状态
        assert 'overall_passed' in report
        assert isinstance(report['overall_passed'], bool)


class TestNPZExport:
    """测试 NPZ 导出功能"""

    @pytest.mark.unit
    def test_npz_export_basic(self, al_s_data, tmp_path):
        """测试基本 NPZ 导出功能"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['npz'],
        )

        # 检查文件生成
        assert len(files) == 1
        npz_file = files[0]
        assert npz_file.exists()
        assert npz_file.suffix == '.npz'

        # 检查文件可加载
        npz_data = np.load(npz_file)

        # 检查基本字段存在
        assert 'radial_grid' in npz_data
        assert 'radial_weights' in npz_data
        assert 'Z' in npz_data
        assert 'xc_code' in npz_data

    @pytest.mark.unit
    def test_npz_naming_convention(self, al_s_data, tmp_path):
        """测试 NPZ 命名约定的正确性"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['npz'],
        )

        npz_data = np.load(files[0])
        keys = list(npz_data.keys())

        # 检查 s 通道（l=0）的命名约定
        assert 'ae_eigenvalues_l0' in keys, "缺少全电子本征值"
        assert 'ps_wavefunction_l0' in keys, "缺少赝波函数"
        assert 'semilocal_potential_l0' in keys, "缺少半局域势"

        # 检查波函数命名（应包含主量子数）
        wf_keys = [k for k in keys if k.startswith('ae_wavefunction_')]
        assert len(wf_keys) > 0, "缺少全电子波函数"

    @pytest.mark.unit
    def test_npz_array_shapes(self, al_s_data, tmp_path):
        """测试 NPZ 数组形状的一致性"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['npz'],
        )

        npz_data = np.load(files[0])
        n_grid = len(data['ae'].r)

        # 检查网格长度
        assert npz_data['radial_grid'].shape == (n_grid,)
        assert npz_data['radial_weights'].shape == (n_grid,)

        # 检查势和波函数长度
        assert npz_data['ps_wavefunction_l0'].shape == (n_grid,)
        assert npz_data['semilocal_potential_l0'].shape == (n_grid,)

        # 检查本征值数量（应为标量数组）
        assert npz_data['ae_eigenvalues_l0'].ndim == 1

    @pytest.mark.unit
    def test_npz_metadata_scalars(self, al_s_data, tmp_path):
        """测试 NPZ 中的元数据标量"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['npz'],
        )

        npz_data = np.load(files[0])

        # 检查原子序数
        assert npz_data['Z'] == 13

        # 检查 XC 代码（PZ81 应为 1）
        assert npz_data['xc_code'] == 1


class TestUnifiedExportInterface:
    """测试统一导出接口"""

    @pytest.mark.unit
    def test_multi_format_export(self, al_s_data, tmp_path):
        """测试同时导出多种格式"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['json', 'npz'],
        )

        # 检查生成了两个文件
        assert len(files) == 2

        # 检查文件类型
        suffixes = {f.suffix for f in files}
        assert '.json' in suffixes
        assert '.npz' in suffixes

        # 检查两个文件都存在
        for f in files:
            assert f.exists()

    @pytest.mark.unit
    def test_l_channel_consistency_check(self, al_s_data, tmp_path):
        """测试 L-channel 一致性检查"""
        data = al_s_data
        output_prefix = tmp_path / "test_al"

        # 故意制造不一致：TM 有 l=0，inv 有 l=1
        tm_dict_bad = {0: data['tm_dict'][0]}
        inv_dict_bad = {1: data['inv_dict'][0]}

        # 应该抛出 ValueError
        with pytest.raises(ValueError, match="l通道不一致"):
            export_pseudopotential(
                ae_result=data['ae'],
                tm_dict=tm_dict_bad,
                inv_dict=inv_dict_bad,
                validation_report=data['report'],
                output_prefix=str(output_prefix),
                formats=['json'],
            )

    @pytest.mark.unit
    def test_output_path_handling(self, al_s_data, tmp_path):
        """测试输出路径处理"""
        data = al_s_data

        # 测试不同的路径格式
        output_prefix = tmp_path / "subdir" / "test_al"

        files = export_pseudopotential(
            ae_result=data['ae'],
            tm_dict=data['tm_dict'],
            inv_dict=data['inv_dict'],
            validation_report=data['report'],
            output_prefix=str(output_prefix),
            formats=['json'],
        )

        # 检查文件名正确
        assert files[0].stem == 'test_al'
        assert files[0].parent.name == 'subdir'


class TestHelperFunctions:
    """测试辅助函数"""

    @pytest.mark.unit
    def test_element_symbol_lookup(self):
        """测试元素符号查找"""
        from atomppgen.export import _get_element_symbol

        # 常见元素
        assert _get_element_symbol(1) == 'H'
        assert _get_element_symbol(6) == 'C'
        assert _get_element_symbol(13) == 'Al'
        assert _get_element_symbol(14) == 'Si'

        # 未定义元素（应返回 Z{n} 格式）
        assert _get_element_symbol(99) == 'Z99'

    @pytest.mark.unit
    def test_xc_code_mapping(self):
        """测试 XC 泛函代码映射"""
        from atomppgen.export import _xc_to_code

        # 标准泛函
        assert _xc_to_code('PZ81') == 1
        assert _xc_to_code('VWN') == 2
        assert _xc_to_code('PW91') == 3
        assert _xc_to_code('PBE') == 4

        # 大小写不敏感
        assert _xc_to_code('pz81') == 1
        assert _xc_to_code('vwn') == 2

        # 未知泛函（应返回 0）
        assert _xc_to_code('UNKNOWN') == 0
