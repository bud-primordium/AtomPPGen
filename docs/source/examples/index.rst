使用示例
========

本节提供完整的使用示例和常见使用模式。

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

完整工作流（Al, LDA）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from atomppgen import (
       solve_ae_atom,
       tm_pseudize,
       invert_semilocal_potential,
       kb_transform,
       export_pseudopotential,
   )
   from atomppgen.validate import run_full_validation

   # 1. 全电子解
   ae = solve_ae_atom(Z=13, xc="PZ81", lmax=2, spin_mode="LDA")

   # 2. TM 伪化（示例：取每个通道能量最高的参考态）
   rc_map = {0: 2.1, 1: 2.2, 2: 2.4}
   tm_results = {}
   for l, rc in rc_map.items():
       tm_results[l] = tm_pseudize(
           r=ae.r,
           w=ae.w,
           u_ae=ae.u_by_l[l][-1],
           eps=ae.eps_by_l[l][-1],
           l=l,
           rc=rc,
           continuity_orders=4,
       )

   # 3. 势反演
   inv_results = {l: invert_semilocal_potential(tm_results[l], ae.r) for l in tm_results}

   # 4. KB 转换（选择 d 通道作为局域道）
   kb = kb_transform(
       invert_results=inv_results,
       u_by_l={l: tm_results[l].u_ps for l in tm_results},
       r=ae.r,
       w=ae.w,
       loc_channel=2,
   )

   # 5. 验证
   report = run_full_validation(ae, tm_results, inv_results, r_test=3.0)

   # 6. 导出
   export_pseudopotential(
       ae_result=ae,
       tm_dict=tm_results,
       inv_dict=inv_results,
       validation_report=report,
       output_prefix="outputs/al_lda_tm",
       formats=["json", "npz"],
   )

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
2. 调整混合参数：``scf_params={"mix_alpha": 0.2}`` (减小)
3. 增加网格点数：``grid_params={"n": 1200}``
4. 降低收敛阈值：``scf_params={"tol": 1e-5}`` (放松)

网格选择
~~~~~~~~

- **高精度需求**：使用 ``exp_transformed`` + ``n >= 1000``
- **快速测试**：使用 ``linear`` + ``n = 600``
- **重原子**：增大 ``rmax``，减小 ``rmin``
