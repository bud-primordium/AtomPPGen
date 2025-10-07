validate - 赝势可转移性验证
==============================

本模块提供赝势质量评估工具，包括范数守恒、对数导数匹配和幽灵态检测。

.. automodule:: atomppgen.validate
   :members:
   :undoc-members:
   :show-inheritance:

主要函数
--------

范数守恒检验
~~~~~~~~~~~~

.. autofunction:: atomppgen.validate.check_norm_conservation

对数导数匹配
~~~~~~~~~~~~

.. autofunction:: atomppgen.validate.check_log_derivative

幽灵态检测
~~~~~~~~~~

.. autofunction:: atomppgen.validate.check_ghost_states

完整验证流程
~~~~~~~~~~~~

.. autofunction:: atomppgen.validate.run_full_validation

数据类
------

.. autoclass:: atomppgen.validate.NormConservationResult
   :members:
   :undoc-members:

.. autoclass:: atomppgen.validate.LogDerivativeResult
   :members:
   :undoc-members:

.. autoclass:: atomppgen.validate.GhostStateResult
   :members:
   :undoc-members:

.. autoclass:: atomppgen.validate.ValidationReport
   :members:
   :undoc-members:

使用示例
--------

基础验证流程
~~~~~~~~~~~~

.. code-block:: python

   from atomppgen import solve_ae_atom, tm_pseudize, invert_semilocal_potential
   from atomppgen.validate import run_full_validation

   # 1. 生成 Al 的全电子解和伪化
   ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0, grid_params={'n': 600})
   tm_s = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1],
                      l=0, rc=2.0)
   inv_s = invert_semilocal_potential(tm_s, ae.r)

   # 2. 构建通道字典
   tm_dict = {0: tm_s}
   inv_dict = {0: inv_s}

   # 3. 完整验证
   report = run_full_validation(
       ae, tm_dict, inv_dict,
       r_test=3.0,                    # 测试半径 (Bohr)
       E_range_Ry=(-0.5, 0.5),       # 能量窗口 (Rydberg)
       E_step_Ry=0.05                # 能量步长 (Rydberg)
   )

   # 4. 查看结果
   print(f"整体通过: {report.overall_passed}")
   print(f"范数守恒 (s): {report.norm_results[0].passed}")
   print(f"对数导数 (s): {report.log_deriv_results[0].passed}")
   if report.ghost_result:
       print(f"幽灵态数量: {report.ghost_result.n_ghosts}")

   # 5. 导出 JSON 报告
   import json
   with open('validation_report.json', 'w') as f:
       json.dump(report.to_dict(), f, indent=2)

范数守恒单项检验
~~~~~~~~~~~~~~~~

.. code-block:: python

   from atomppgen import solve_ae_atom, tm_pseudize
   from atomppgen.validate import check_norm_conservation

   ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=0)
   tm = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1],
                    l=0, rc=2.0)

   # 范数守恒检验
   result = check_norm_conservation(tm, tolerance=1e-6)

   print(f"范数误差: {result.norm_error:.3e}")
   print(f"通过检验: {result.passed}")
   print(f"截断半径: {result.rc} Bohr")
   print(f"求解器收敛: {result.diagnostics['solver_converged']}")

对数导数匹配单项检验
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from atomppgen.validate import check_log_derivative, _extract_ks_potential

   # 提取 KS 有效势和伪势
   V_AE = _extract_ks_potential(ae)
   V_PS = inv.V_l

   # 对数导数匹配（缩小参数加速示例）
   result = check_log_derivative(
       V_AE, V_PS, ae.r,
       l=0,                          # s 轨道
       r_test=3.0,                   # 测试半径 (Bohr)
       E_range_Ha=(-0.2, 0.2),       # 能量窗口 (Hartree)
       E_step_Ha=0.05                # 能量步长 (Hartree)
   )

   print(f"零点 RMS: {result.zero_crossing_rms:.6f} Ha")
   print(f"曲线 RMS: {result.curve_rms:.6f}")
   print(f"通过检验: {result.passed}")

   # 绘制对数导数曲线
   import matplotlib.pyplot as plt
   plt.plot(result.energies, result.L_AE, label='AE')
   plt.plot(result.energies, result.L_PS, label='PS', linestyle='--')
   plt.axhline(0, color='k', linewidth=0.5)
   plt.xlabel('Energy (Ha)')
   plt.ylabel('L(E, r_test)')
   plt.legend()
   plt.savefig('log_derivative.png')

幽灵态检测单项检验
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from atomppgen.validate import check_ghost_states

   # 幽灵态检测 (径向方法)
   result = check_ghost_states(
       inv, ae.r, ae.w,
       valence_energy=tm.eps,
       E_window_Ha=(-0.25, 0.25),
       method='radial'
   )

   print(f"检测到幽灵态数量: {result.n_ghosts}")
   print(f"能量窗口内束缚态数: {result.diagnostics['n_bound_states_in_window']}")
   print(f"通过检验: {result.passed}")

   # 查看所有束缚态
   if len(result.eigenvalues) > 0:
       print("窗口内本征值 (Ha):")
       for E in result.eigenvalues:
           print(f"  {E:.6f}")

多通道验证
~~~~~~~~~~

.. code-block:: python

   # 生成 s, p, d 三通道
   ae = solve_ae_atom(Z=13, spin_mode='LDA', lmax=2)

   # s 通道
   tm_s = tm_pseudize(ae.r, ae.w, ae.u_by_l[0][-1], ae.eps_by_l[0][-1],
                      l=0, rc=2.0)
   inv_s = invert_semilocal_potential(tm_s, ae.r)

   # p 通道 (需调参)
   tm_p = tm_pseudize(ae.r, ae.w, ae.u_by_l[1][-1], ae.eps_by_l[1][-1],
                      l=1, rc=1.9, continuity_orders=2)
   inv_p = invert_semilocal_potential(tm_p, ae.r)

   # d 通道
   tm_d = tm_pseudize(ae.r, ae.w, ae.u_by_l[2][-1], ae.eps_by_l[2][-1],
                      l=2, rc=2.4)
   inv_d = invert_semilocal_potential(tm_d, ae.r)

   # 完整验证
   tm_dict = {0: tm_s, 1: tm_p, 2: tm_d}
   inv_dict = {0: inv_s, 1: inv_p, 2: inv_d}

   report = run_full_validation(ae, tm_dict, inv_dict)

   # 逐通道查看结果
   for l in [0, 1, 2]:
       label = ['s', 'p', 'd'][l]
       norm_pass = report.norm_results[l].passed
       ld_pass = report.log_deriv_results[l].passed
       print(f"{label} 通道: 范数={norm_pass}, 对数导数={ld_pass}")

技术细节
--------

能量单位转换
~~~~~~~~~~~~

所有内部计算使用 **Hartree (Ha)** 原子单位。外部 API 支持 Rydberg (Ry) 输入，
自动转换：

.. code-block:: python

   # 用户指定 Rydberg
   E_range_Ry = (-0.5, 0.5)
   E_step_Ry = 0.05

   # 内部转换为 Hartree
   E_range_Ha = (E_range_Ry[0] / 2, E_range_Ry[1] / 2)  # (-0.25, 0.25) Ha
   E_step_Ha = E_step_Ry / 2                            # 0.025 Ha

Numerov 求解器
~~~~~~~~~~~~~~

对于非均匀网格（如 ``exp_transformed``），自动重采样到等距网格：

.. code-block:: python

   # 检查网格均匀性
   dr = np.diff(r)
   is_uniform = np.allclose(dr, dr[0], rtol=1e-6)

   if not is_uniform:
       # 重采样到均匀网格
       n_uniform = len(r)
       r_uniform = np.linspace(r[0], r[-1], n_uniform)
       V_uniform = interp1d(r, V, kind='cubic')(r_uniform)
       # 在均匀网格上求解...
       # 结果插值回原网格

对数导数零点匹配
~~~~~~~~~~~~~~~~

零点通过符号变化检测，线性插值精确化：

.. code-block:: python

   def find_zero_crossings(E, L):
       zeros = []
       for i in range(len(L) - 1):
           if L[i] * L[i+1] < 0:  # 符号变化
               # 线性插值
               E_zero = E[i] - L[i] * (E[i+1] - E[i]) / (L[i+1] - L[i])
               zeros.append(E_zero)
       return zeros

幽灵态判定准则
~~~~~~~~~~~~~~

束缚态筛选条件：

1. 能量 :math:`E < E_{\\text{max}}` (默认 0 Ha)
2. 边界衰减：:math:`|\\psi(r_{\\text{max}})| < 0.1 \\cdot \\max|\\psi|`

幽灵态判定：

- 在能量窗口 :math:`[-0.25, +0.25]` Ha 内
- 距已知价态 :math:`>0.1` Ha 的额外束缚态

性能优化建议
~~~~~~~~~~~~

1. **对数导数扫描**：粗扫 (0.1 Ry) 定位零点，细扫 (0.025 Ry) 精确化
2. **幽灵态检测**：限制网格点数 (≤300) 加速对角化
3. **多通道并行**：各通道独立，可并行计算（当前串行）

参考算法文档
------------

详细数学推导和原理见 :doc:`../algorithm/validation_methods`。
