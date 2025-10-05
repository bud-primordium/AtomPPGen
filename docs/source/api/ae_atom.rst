ae_atom - 全电子原子解
=======================

本模块封装了 AtomSCF，用于获取全电子原子的自洽场解。

.. automodule:: atomppgen.ae_atom
   :members:
   :undoc-members:
   :show-inheritance:

主要函数
--------

.. autofunction:: atomppgen.ae_atom.solve_ae_atom

数据类
------

.. autoclass:: atomppgen.ae_atom.AEAtomResult
   :members:
   :undoc-members:

使用示例
--------

基础用法
~~~~~~~~

.. code-block:: python

   from atomppgen import solve_ae_atom

   # 求解 Al 原子
   result = solve_ae_atom(
       Z=13,
       xc="PZ81",
       lmax=2,
       grid_type="exp_transformed",
       grid_params={"n": 800, "rmax": 100.0},
       scf_params={"tol": 1e-6, "maxiter": 150}
   )

   # 查看结果
   print(f"SCF 收敛: {result.converged}")
   print(f"迭代次数: {result.scf_iterations}")
   print(f"3s 能级: {result.eps_by_l[0][2]:.6f} Ha")
   print(f"3p 能级: {result.eps_by_l[1][2]:.6f} Ha")

不同网格类型
~~~~~~~~~~~~

.. code-block:: python

   # 指数变换网格（推荐，最高精度）
   result_exp = solve_ae_atom(
       Z=13,
       grid_type="exp_transformed",
       grid_params={"n": 1200, "rmax": 120.0, "total_span": 6.5}
   )

   # 线性网格
   result_linear = solve_ae_atom(
       Z=13,
       grid_type="linear",
       grid_params={"n": 1200, "rmin": 1e-6, "rmax": 50.0}
   )

   # 对数网格
   result_log = solve_ae_atom(
       Z=13,
       grid_type="log",
       grid_params={"n": 1000, "rmin": 1e-6, "rmax": 50.0}
   )

访问波函数和密度
~~~~~~~~~~~~~~~~

.. code-block:: python

   result = solve_ae_atom(Z=13, xc="PZ81", lmax=2)

   # 获取 3s 轨道
   r = result.r
   u_3s = result.u_by_l[0][2]  # l=0（s），第3个态（n=3）

   # 总电子密度
   n_total = result.n_total

   # 能量分解
   print(f"总能量: {result.energies['E_total']:.6f} Ha")
   print(f"交换能: {result.energies['E_x']:.6f} Ha")
   print(f"关联能: {result.energies['E_c']:.6f} Ha")

技术细节
--------

网格-求解器兼容性
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - 网格类型
     - 推荐求解器
     - 注意事项
   * - ``exp_transformed``
     - ``transformed``
     - 最高精度（~7x 提升），需要 delta/Rp 参数
   * - ``linear``
     - ``fd5``
     - 等距网格，可用对称 FD 矩阵
   * - ``log``
     - ``fd5_aux``
     - 非等距，ln(r) 等差

与 NIST 数据对比
~~~~~~~~~~~~~~~~

当前实现为非相对论 LSDA，与 NIST LSD 参考数据存在系统性差异：

.. list-table::
   :header-rows: 1

   * - 项目
     - AtomSCF
     - NIST LSD
     - 差异
   * - 总能量 (Ha)
     - -237.30
     - -241.32
     - ~4 Ha
   * - 1s (Ha)
     - -54.27
     - -55.15
     - ~0.88 Ha
   * - 3s (Ha)
     - -0.25
     - -0.30
     - ~0.05 Ha

说明：

- 差异来源可能包括数值方法、泛函实现细节等
- 价层轨道（3s, 3p）相对差异较小（~0.03-0.05 Ha）
- 对于赝势生成，价层相对精度通常已足够
