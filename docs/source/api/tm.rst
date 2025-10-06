tm - Troullier-Martins 伪化器
===============================

本模块实现 Troullier-Martins (TM) 方法，将全电子轨道伪化为光滑无节点的伪轨道，满足范数守恒条件。

.. automodule:: atomppgen.tm
   :members:
   :undoc-members:
   :show-inheritance:

主要函数
--------

.. autofunction:: atomppgen.tm.tm_pseudize

数据类
------

.. autoclass:: atomppgen.tm.TMResult
   :members:
   :undoc-members:

使用示例
--------

基础用法
~~~~~~~~

.. code-block:: python

   from atomppgen import solve_ae_atom, tm_pseudize
   import numpy as np

   # 1. 获取全电子解
   ae_result = solve_ae_atom(
       Z=13,  # Al
       spin_mode="LDA",
       lmax=0,
       grid_type="exp_transformed",
       grid_params={"n": 800, "rmax": 100.0},
       scf_params={"tol": 1e-6}
   )

   # 2. 提取 3s 轨道
   r = ae_result.r
   w = ae_result.w
   u_ae_3s = ae_result.u_by_l[0][2]
   eps_3s = ae_result.eps_by_l[0][2]

   # 3. TM 伪化
   result = tm_pseudize(
       r=r,
       w=w,
       u_ae=u_ae_3s,
       eps=eps_3s,
       l=0,
       rc=2.1,  # Bohr
       continuity_orders=2,
   )

   # 4. 查看结果
   print(f"范数守恒误差: {result.norm_error:.2e}")
   print(f"TM 系数: {result.a_coeff}")

   # 5. 绘图对比
   import matplotlib.pyplot as plt
   plt.plot(r, u_ae_3s, label='AE')
   plt.plot(r, result.u_ps, '--', label='PS')
   plt.axvline(result.rc, color='k', ls=':', alpha=0.5, label=f'rc={result.rc}')
   plt.legend()
   plt.show()

不同连续性阶数
~~~~~~~~~~~~~~

.. code-block:: python

   # 二阶连续性（4 个系数）
   result_2nd = tm_pseudize(
       r=r, w=w, u_ae=u_ae, eps=eps, l=0, rc=2.1,
       continuity_orders=2
   )
   print(f"系数数量: {len(result_2nd.a_coeff)}")  # 4
   print(f"范数误差: {result_2nd.norm_error:.2e}")

   # 四阶连续性（6 个系数，更严格）
   result_4th = tm_pseudize(
       r=r, w=w, u_ae=u_ae, eps=eps, l=0, rc=2.1,
       continuity_orders=4
   )
   print(f"系数数量: {len(result_4th.a_coeff)}")  # 6
   print(f"范数误差: {result_4th.norm_error:.2e}")

多通道伪化
~~~~~~~~~~

.. code-block:: python

   # Al 的 s, p, d 通道
   ae_result = solve_ae_atom(Z=13, spin_mode="LDA", lmax=2)

   rc_by_l = {0: 2.1, 1: 2.2, 2: 2.4}  # 推荐截断半径
   tm_results = {}

   for l in range(3):  # s, p, d
       # 价层轨道索引（Al: 3s, 3p, 3d）
       n_valence = {0: 2, 1: 1, 2: 0}  # 对应 u_by_l 的索引

       u_ae = ae_result.u_by_l[l][n_valence[l]]
       eps = ae_result.eps_by_l[l][n_valence[l]]

       result = tm_pseudize(
           r=ae_result.r, w=ae_result.w,
           u_ae=u_ae, eps=eps, l=l,
           rc=rc_by_l[l],
           continuity_orders=2
       )
       tm_results[l] = result

       print(f"l={l}: norm_error={result.norm_error:.2e}")

连续性检查
~~~~~~~~~~

.. code-block:: python

   result = tm_pseudize(r=r, w=w, u_ae=u_ae, eps=eps, l=0, rc=2.0)

   # 检查各阶导数的连续性
   for key in ['u', 'du', 'd2u']:
       check = result.continuity_check[key]
       print(f"{key}:")
       print(f"  PS = {check['ps']:.6e}")
       print(f"  AE = {check['ae']:.6e}")
       print(f"  相对误差 = {check['rel_error']:.2e}")

技术细节
--------

截断半径选择指南
~~~~~~~~~~~~~~~~

推荐值（单位：Bohr）：

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - 元素
     - 通道
     - rc 推荐值
   * - Al (Z=13)
     - s
     - 2.0-2.2
   * -
     - p
     - 2.1-2.3
   * -
     - d
     - 2.3-2.5

选择原则：

1. **最小值**：rc 必须大于最外层节点半径（确保内区无节点）
2. **最大值**：rc 应小于共价半径的 1.5 倍（避免过软）
3. **平衡点**：较小的 rc 产生更硬的赝势（可转移性好），但计算成本高

连续性阶数对比
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 25 25 20

   * - continuity_orders
     - 系数数量
     - 约束条件
     - 典型精度
   * - 2
     - 4
     - u, u', u'', 范数
     - 1e-6
   * - 4
     - 6
     - u, u', ..., u'''', 范数
     - 1e-7

说明：

- **二阶**（推荐）：平衡精度与求解速度，适合大多数情况
- **四阶**：更平滑的伪势，适合高精度计算或重元素

数值稳定性
~~~~~~~~~~

TM 伪化内部采用以下保护措施：

1. **指数溢出保护**：在 exp(p(r)) 中，如果 p(r) > 100，自动调整初值
2. **范数积分精度**：使用 Simpson 法则，精度 O(h^4)
3. **非线性求解器容差**：xtol=1e-10，自动重试不收敛情况

如果遇到收敛问题：

- 尝试调整 rc（±0.1 Bohr）
- 检查 AE 解是否正确（converged=True）
- 降低 continuity_orders（从 4 降到 2）

范数守恒验收标准
~~~~~~~~~~~~~~~~

**教学标准**：norm_error < 1e-5

**生产标准**：norm_error < 1e-6（推荐用于发表）

典型误差水平：

.. list-table::
   :header-rows: 1

   * - rc (Bohr)
     - 范数误差（Al 3s）
   * - 1.8
     - ~1e-7
   * - 2.0
     - ~5e-7
   * - 2.2
     - ~2e-6
   * - 2.4
     - ~8e-6

说明：较大的 rc 会增加范数守恒误差，但通常仍在可接受范围内。

参考文献
--------

- **原始论文**：N. Troullier and J. L. Martins,
  *Efficient pseudopotentials for plane-wave calculations*,
  **Phys. Rev. B** 43, 1993-2006 (1991).
  DOI: `10.1103/PhysRevB.43.1993 <https://doi.org/10.1103/PhysRevB.43.1993>`_

- **QE 实现指南**：P. Giannozzi, *Notes on pseudopotential generation* (2019),
  https://www.quantum-espresso.org/wp-content/uploads/2022/03/pseudo-gen.pdf

- **算法细节**：:doc:`../algorithm/tm_method`
