export - 赝势导出模块
=======================

本模块提供多种格式的赝势文件导出功能，支持 JSON（元数据）、NPZ（数值数组）和 UPF（Quantum ESPRESSO）格式。

.. automodule:: atomppgen.export
   :members:
   :undoc-members:
   :show-inheritance:

主要函数
--------

.. autofunction:: atomppgen.export.export_pseudopotential

数据类
------

.. autoclass:: atomppgen.export.PseudopotentialData
   :members:
   :undoc-members:

使用示例
--------

完整工作流
~~~~~~~~~~

.. code-block:: python

   from atomppgen import (
       solve_ae_atom,
       tm_pseudize,
       invert_semilocal_potential,
       export_pseudopotential,
   )
   from atomppgen.validate import run_full_validation

   # 1. 全电子原子解
   ae = solve_ae_atom(
       Z=13,  # Al
       spin_mode='LDA',
       lmax=0,
       grid_type='exp_transformed',
       grid_params={'n': 600}
   )

   # 2. TM 伪化
   tm_s = tm_pseudize(
       ae.r, ae.w,
       ae.u_by_l[0][-1],  # 3s 价电子
       ae.eps_by_l[0][-1],
       l=0,
       rc=2.5
   )

   # 3. 势反演
   inv_s = invert_semilocal_potential(tm_s, ae.r)

   # 4. 验证
   report = run_full_validation(
       ae,
       tm_dict={0: tm_s},
       inv_dict={0: inv_s},
       r_test=3.0,
       E_range_Ry=(-0.5, 0.5),
       E_step_Ry=0.05,
   )

   # 5. 导出到 JSON + NPZ 格式
   files = export_pseudopotential(
       ae_result=ae,
       tm_dict={0: tm_s},
       inv_dict={0: inv_s},
       validation_report=report,
       output_prefix='outputs/al_lda',
       formats=['json', 'npz'],
   )

   print(f"已导出 {len(files)} 个文件:")
   for f in files:
       print(f"  - {f.name} ({f.stat().st_size / 1024:.1f} KB)")

单一格式导出
~~~~~~~~~~~~

.. code-block:: python

   # 仅导出 JSON（元数据 + 验证报告）
   json_files = export_pseudopotential(
       ae_result=ae,
       tm_dict=tm_dict,
       inv_dict=inv_dict,
       validation_report=report,
       output_prefix='outputs/al_lda',
       formats=['json'],
   )

   # 仅导出 NPZ（数值数组）
   npz_files = export_pseudopotential(
       ae_result=ae,
       tm_dict=tm_dict,
       inv_dict=inv_dict,
       validation_report=report,
       output_prefix='outputs/al_lda',
       formats=['npz'],
   )

添加元数据
~~~~~~~~~~

.. code-block:: python

   # 添加 Git 提交信息等元数据
   files = export_pseudopotential(
       ae_result=ae,
       tm_dict=tm_dict,
       inv_dict=inv_dict,
       validation_report=report,
       output_prefix='outputs/al_lda',
       formats=['json', 'npz'],
       metadata={
           'git_commit': 'a1b2c3d',
           'notes': 'Optimized rc parameter',
       }
   )

UPF 导出（实验性）
~~~~~~~~~~~~~~~~~

UPF writer 当前提供“最小可解析、字段可追溯”的结构，用于后续逐步对齐 QE 的严格要求。
建议仍以 JSON/NPZ 作为权威数据源，UPF 视为中间格式。

导出 UPF 时需要提供价电子数 ``z_valence``（教学示例可用 Na/Al/Si 的默认推断，但推荐显式指定）：

.. code-block:: python

   files = export_pseudopotential(
       ae_result=ae,
       tm_dict=tm_dict,
       inv_dict=inv_dict,
       validation_report=report,
       output_prefix='outputs/al_lda',
       kb_result=kb,  # 建议提供，以输出 PP_LOCAL/PP_NONLOCAL
       formats=['upf'],
       metadata={'z_valence': 3.0},
   )

加载导出数据
~~~~~~~~~~~~

.. code-block:: python

   import json
   import numpy as np

   # 加载 JSON 元数据
   with open('outputs/al_lda.json', 'r') as f:
       metadata = json.load(f)

   print(f"元素: {metadata['metadata']['element']['symbol']}")
   print(f"泛函: {metadata['metadata']['xc']['functional']}")
   print(f"验证: {metadata['validation_report']['overall_passed']}")

   # 加载 NPZ 数值数据
   npz_data = np.load('outputs/al_lda.npz')

   r = npz_data['radial_grid']
   V_s = npz_data['semilocal_potential_l0']
   psi_s = npz_data['ps_wavefunction_l0']

   print(f"网格点数: {len(r)}")
   print(f"s 通道势能形状: {V_s.shape}")

   # 绘图
   import matplotlib.pyplot as plt

   plt.plot(r, V_s * r, label='$r V_s(r)$')
   plt.xlabel('r (Bohr)')
   plt.ylabel('$r V_s$ (Ha·Bohr)')
   plt.legend()
   plt.savefig('potential_s.png')

多通道导出
~~~~~~~~~~

.. code-block:: python

   # Al 的 s, p, d 通道
   ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=2)

   # TM 伪化各通道
   rc_by_l = {0: 2.1, 1: 2.2, 2: 2.4}
   tm_dict = {}
   inv_dict = {}

   for l in range(3):
       # 价层索引（Al: 3s, 3p, 3d）
       n_valence = {0: 2, 1: 1, 2: 0}
       u_ae = ae.u_by_l[l][n_valence[l]]
       eps = ae.eps_by_l[l][n_valence[l]]

       # TM 伪化
       tm = tm_pseudize(ae.r, ae.w, u_ae, eps, l=l, rc=rc_by_l[l])
       tm_dict[l] = tm

       # 势反演
       inv = invert_semilocal_potential(tm, ae.r)
       inv_dict[l] = inv

   # 验证
   report = run_full_validation(
       ae, tm_dict, inv_dict,
       r_test=3.0,
       E_range_Ry=(-0.5, 0.5),
       E_step_Ry=0.05,
   )

   # 导出
   files = export_pseudopotential(
       ae_result=ae,
       tm_dict=tm_dict,
       inv_dict=inv_dict,
       validation_report=report,
       output_prefix='outputs/al_lda_spd',
       formats=['json', 'npz'],
   )

技术细节
--------

JSON 格式说明
~~~~~~~~~~~~~

JSON 文件包含以下顶层字段：

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 字段
     - 内容
   * - metadata
     - 元素、泛函、生成日期、代码版本
   * - generation_params
     - TM 参数（rc_by_l）、网格参数
   * - all_electron_reference
     - 全电子本征值（按角动量 l 索引）
   * - pseudopotential
     - 通道列表、数据文件引用
   * - validation_report
     - 范数守恒、对数导数、幽灵态检测结果
   * - radial_grid
     - 网格点数、范围（不包含数组数据）

单位约定：

- 能量：Hartree 原子单位
- 长度：Bohr
- 在 ``metadata.units`` 中明确标注为 ``"Hartree_atomic"``

NPZ 格式说明
~~~~~~~~~~~~

NPZ 文件采用压缩格式（``np.savez_compressed``），包含以下数组：

**网格数据**：

- ``radial_grid``：径向网格 :math:`r` (Bohr)
- ``radial_weights``：积分权重 :math:`w`

**标量元数据**：

- ``Z``：原子序数
- ``xc_code``：泛函代码（1=PZ81, 2=VWN, 3=PW91, 4=PBE）

**全电子数据**（按角动量 l 索引）：

- ``ae_eigenvalues_l{l}``：本征值数组 (Ha)
- ``ae_wavefunction_l{l}_n{n}``：径向波函数 :math:`u(r)`

**赝势数据**（按角动量 l 索引）：

- ``ps_wavefunction_l{l}``：伪波函数 :math:`u_{\mathrm{PS}}(r)`
- ``semilocal_potential_l{l}``：半局域势 :math:`V_l(r)` (Ha)

**KB 可分离形式数据**（可选）：

如果在导出时传入 ``kb_result``，NPZ 还会包含 KB 形式所需的关键数组：

- ``kb_loc_channel``：局域通道 :math:`l^*`
- ``kb_V_loc``：局域势 :math:`V_{\mathrm{loc}}(r)` (Ha)
- ``kb_beta_l{l}``：投影子 :math:`\beta_l(r)`（仅对非局域通道存在）
- ``kb_D_l{l}``：耦合系数 :math:`D_l` (Ha)

命名约定示例：

.. code-block:: python

   # s 通道（l=0）
   'ae_eigenvalues_l0'        # [ε_1s, ε_2s, ε_3s]
   'ae_wavefunction_l0_n3'    # u_3s(r)（价电子）
   'ps_wavefunction_l0'       # u_PS,s(r)
   'semilocal_potential_l0'   # V_s(r)

   # p 通道（l=1）
   'ae_eigenvalues_l1'        # [ε_2p, ε_3p]
   'ae_wavefunction_l1_n2'    # u_3p(r)
   'ps_wavefunction_l1'       # u_PS,p(r)
   'semilocal_potential_l1'   # V_p(r)

数据加载验证
~~~~~~~~~~~~

.. code-block:: python

   import numpy as np

   # 加载 NPZ 数据
   npz_data = np.load('outputs/al_lda.npz')

   # 验证数组形状一致性
   n_grid = len(npz_data['radial_grid'])
   assert npz_data['ps_wavefunction_l0'].shape == (n_grid,)
   assert npz_data['semilocal_potential_l0'].shape == (n_grid,)

   # 验证元数据
   assert npz_data['Z'] == 13  # Al
   assert npz_data['xc_code'] == 1  # PZ81

   # 验证本征值数量
   n_valence_states = len(npz_data['ae_eigenvalues_l0'])
   print(f"s 通道价层态数: {n_valence_states}")

   # 检查数组中无 NaN/Inf
   for key in npz_data.keys():
       arr = npz_data[key]
       if isinstance(arr, np.ndarray):
           assert np.all(np.isfinite(arr)), f"{key} 包含 NaN/Inf"

输出路径处理
~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path

   # 自动创建父目录
   output_prefix = Path('outputs/subdir/test_al')
   files = export_pseudopotential(
       ae_result=ae,
       tm_dict=tm_dict,
       inv_dict=inv_dict,
       validation_report=report,
       output_prefix=str(output_prefix),
       formats=['json', 'npz'],
   )
   # 输出文件:
   # - outputs/subdir/test_al.json
   # - outputs/subdir/test_al.npz

   # 使用字符串路径
   files = export_pseudopotential(
       ...,
       output_prefix='outputs/al_lda',
       ...
   )

   # 使用 Path 对象
   files = export_pseudopotential(
       ...,
       output_prefix=Path.cwd() / 'outputs' / 'al_lda',
       ...
   )

错误处理
~~~~~~~~

.. code-block:: python

   # L-channel 一致性检查
   tm_dict = {0: tm_s, 1: tm_p}  # s, p 通道
   inv_dict = {0: inv_s}          # 仅 s 通道

   try:
       files = export_pseudopotential(
           ae_result=ae,
           tm_dict=tm_dict,
           inv_dict=inv_dict,
           validation_report=report,
           output_prefix='outputs/test',
           formats=['json'],
       )
   except ValueError as e:
       print(f"错误: {e}")
       # 输出: "TM和势反演的l通道不一致: tm={0, 1}, inv={0}"

文件大小估算
~~~~~~~~~~~~

典型文件大小（Al, lmax=0, 600 网格点）：

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - 格式
     - 文件大小
     - 说明
   * - JSON
     - ~3 KB
     - 元数据 + 验证报告（不含数组）
   * - NPZ
     - ~24 KB
     - 压缩数值数组
   * - JSON (含数组)
     - ~1 MB
     - include_arrays=True（不推荐）

多通道情况（lmax=2）：

- JSON: ~5 KB
- NPZ: ~60 KB

建议
~~~~

**推荐实践**：

1. **同时导出 JSON + NPZ**：JSON 存储元数据和验证报告，NPZ 存储数值数据
2. **JSON 不包含数组**：默认 ``include_arrays=False``，避免文件过大
3. **验证输出**：导出后立即加载 NPZ，检查数组形状和数值范围
4. **版本控制**：在 ``metadata`` 中添加 ``git_commit`` 字段

**避免的做法**：

- 不要在 JSON 中包含大数组（``include_arrays=True``）
- 不要混用不同版本的 AtomPPGen 导出的文件
- 不要手动编辑 NPZ 文件（使用 Python 重新导出）

参考示例
--------

完整示例脚本：

- ``examples/export_al_example.py``：Al LDA 赝势生成与导出完整流程

相关模块
--------

- :doc:`ae_atom`：全电子原子解
- :doc:`tm`：TM 伪化器
- :doc:`invert`：势反演
- :doc:`validate`：可转移性验证
