AtomPPGen 文档
================

模守恒赝势生成器 - 基于 Troullier-Martins 方法的教学实现

.. toctree::
   :maxdepth: 2
   :caption: 目录

   introduction
   algorithm/index
   api/index
   examples/index

.. toctree::
   :maxdepth: 2
   :caption: 用户指南

   user_guide/parameter_selection
   user_guide/validation_standards

.. toctree::
   :maxdepth: 2
   :caption: 实验结果

   results/al_lda_validation

简介
----

AtomPPGen 是一个教学用赝势生成工具，实现了：

- **全电子原子解** - 调用 AtomSCF 求解参考态
- **TM 伪化** - Troullier-Martins 方法生成赝轨道
- **势反演** - 从伪轨道反演半局域势
- **KB 形式** - Kleinman-Bylander 可分离非局域投影
- **可转移性检验** - 范数守恒、对数导数曲线、幽灵态检测
- **数据导出** - JSON/NPZ 格式

快速开始
--------

安装
~~~~

.. code-block:: bash

   # 克隆仓库
   cd AtomPPGen

   # 使用 uv（推荐）
   uv venv
   source .venv/bin/activate
   uv pip install -e ../AtomSCF  # 安装依赖
   uv pip install -e ".[dev,docs]"

   # 或使用标准 venv
   python -m venv .venv
   source .venv/bin/activate
   pip install -e ../AtomSCF
   pip install -e ".[dev,docs]"

示例：Al 原子全电子解
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from atomppgen import solve_ae_atom

   # 求解 Al 原子（使用变量变换方法）
   result = solve_ae_atom(
       Z=13,
       xc="PZ81",
       lmax=2,  # s, p, d 通道
       grid_type="exp_transformed",
       grid_params={"n": 800, "rmax": 100.0},
       scf_params={"tol": 1e-6, "maxiter": 150}
   )

   print(f"3s 能级: {result.eps_by_l[0][2]:.6f} Ha")
   print(f"3p 能级: {result.eps_by_l[1][2]:.6f} Ha")
   print(f"总能量: {result.energies['E_total']:.6f} Ha")

设计目标
--------

**主要目标**: 生成 Al（Z=13）的 LDA 模守恒赝势

- 方法：Troullier-Martins (TM) + Kleinman-Bylander (KB)
- 交换关联：LDA-PZ81
- 通道：s, p, d（d 作为散射投影）
- 验收标准：

  - 范数守恒误差 < 1e-6
  - 能级偏差 < 5 mRy
  - 对数导数曲线一致性良好
  - 无幽灵态

文档结构
--------

- **算法原理**: TM 伪化、势反演、KB 转换的数学推导
- **API 参考**: 完整的函数和类接口文档
- **示例**: Al 赝势生成流程和可转移性检验

索引与搜索
----------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
