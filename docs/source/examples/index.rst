使用示例
========

本节提供完整的使用示例和常见使用模式。

.. toctree::
   :maxdepth: 2

   al_pseudopotential
   testing_transferability
   export_formats

快速开始
--------

基础流程
~~~~~~~~

最简单的使用方式：

.. code-block:: python

   from atomppgen import solve_ae_atom

   # 步骤 1: 求解全电子原子
   ae_result = solve_ae_atom(Z=13, xc="PZ81", lmax=2)

   print(f"3s 能级: {ae_result.eps_by_l[0][2]:.6f} Ha")
   print(f"总能量: {ae_result.energies['E_total']:.6f} Ha")

完整工作流（待实现）
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from atomppgen import (
       solve_ae_atom,
       tm_pseudize,
       invert_potential,
       kb_transform,
       export_pseudopotential,
   )

   # 1. 全电子解
   ae = solve_ae_atom(Z=13, xc="PZ81", lmax=2)

   # 2. TM 伪化（待实现）
   # tm_s = tm_pseudize(r=ae.r, w=ae.w, v_ae=ae.u_by_l[0][2], ...)

   # 3. 势反演（待实现）
   # inv_s = invert_potential(r=ae.r, v_ps=tm_s.v_ps, ...)

   # 4. KB 转换（待实现）
   # kb = kb_transform(V_by_l={...}, ...)

   # 5. 导出（待实现）
   # export_pseudopotential("al_lda.json", ...)

常见使用模式
------------

调整网格参数
~~~~~~~~~~~~

.. code-block:: python

   # 高精度计算
   result_high = solve_ae_atom(
       Z=13,
       grid_type="exp_transformed",
       grid_params={"n": 1500, "rmax": 150.0, "total_span": 7.0},
       scf_params={"tol": 1e-8, "maxiter": 300}
   )

   # 快速测试
   result_fast = solve_ae_atom(
       Z=13,
       grid_params={"n": 500, "rmax": 60.0},
       scf_params={"tol": 1e-5, "maxiter": 80}
   )

不同元素
~~~~~~~~

.. code-block:: python

   # 氢原子
   h_result = solve_ae_atom(Z=1, lmax=0)

   # 碳原子
   c_result = solve_ae_atom(Z=6, lmax=2)

   # 硅原子
   si_result = solve_ae_atom(Z=14, lmax=2)

故障排查
--------

SCF 不收敛
~~~~~~~~~~

如果遇到 SCF 不收敛问题：

1. 增加迭代次数：``scf_params={"maxiter": 300}``
2. 调整混合参数：``scf_params={"mix_alpha": 0.2}``（减小）
3. 增加网格点数：``grid_params={"n": 1200}``
4. 降低收敛阈值：``scf_params={"tol": 1e-5}``（放松）

网格选择
~~~~~~~~

- **高精度需求**：使用 ``exp_transformed`` + ``n >= 1000``
- **快速测试**：使用 ``linear`` + ``n = 600``
- **重原子**：增大 ``rmax``，减小 ``rmin``
