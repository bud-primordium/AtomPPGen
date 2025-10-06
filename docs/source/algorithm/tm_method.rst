Troullier-Martins 伪化方法
===========================

本节详细介绍 Troullier-Martins (TM) 方法的数学推导和数值实现。

数学基础
--------

伪化目标
~~~~~~~~

给定全电子（AE）轨道 :math:`\phi_{\text{AE}}(r)` 和本征能量 :math:`\varepsilon_l`，
构造伪轨道 :math:`\phi_{\text{PS}}(r)`，满足：

1. **外区一致性** （:math:`r > r_c`）：

   .. math::

      \phi_{\text{PS}}(r) = \phi_{\text{AE}}(r), \quad r > r_c

2. **内区光滑性** （:math:`r \leq r_c`）：

   - 无节点（避免不必要的振荡）
   - 足够光滑（高阶导数连续）

3. **范数守恒**：

   .. math::

      \int_0^{r_c} |\phi_{\text{PS}}(r)|^2 dr = \int_0^{r_c} |\phi_{\text{AE}}(r)|^2 dr

径向波函数形式
~~~~~~~~~~~~~~

使用径向波函数 :math:`u_{nl}(r) = r R_{nl}(r)`，在 :math:`r \leq r_c` 内采用 TM 参数化：

.. math::

   u_l(r) = r^{l+1} \exp\left( \sum_{i=0}^{N} a_{2i} r^{2i} \right)

**说明**：
- 因子 :math:`r^{l+1}` 保证原点处行为正确（:math:`u_l(0) = 0`）
- 指数形式 :math:`\exp(p(r))` 保证正定性
- 偶次幂 :math:`r^{2i}` 保证光滑性（所有导数在原点连续）

TM 方法
-------

连续性约束
~~~~~~~~~~~

在截断半径 :math:`r = r_c` 处，要求伪轨道与 AE 轨道的函数值及前 :math:`2N-1` 阶导数连续：

.. math::

   u_{\text{PS}}^{(k)}(r_c) = u_{\text{AE}}^{(k)}(r_c), \quad k = 0, 1, 2, \ldots, 2N-1

原始 TM 方法使用 :math:`N=3`（6 个约束），我们简化为 :math:`N=2`（4 个约束）：

- 函数值匹配：:math:`u(r_c)`
- 一阶导数匹配：:math:`u'(r_c)`
- 二阶导数匹配：:math:`u''(r_c)`
- 范数守恒：:math:`\int_0^{r_c} |u|^2 dr`

这给出 4 个方程，求解 4 个系数 :math:`a_0, a_2, a_4, a_6`。

导数计算
~~~~~~~~

设 :math:`p(r) = \sum_{i=0}^N a_{2i} r^{2i}`，则：

.. math::

   u(r) &= r^{l+1} e^{p(r)} \\\\
   u'(r) &= \left[(l+1) r^l + r^{l+1} p'(r)\right] e^{p(r)} \\\\
   u''(r) &= \left[l(l+1) r^{l-1} + 2(l+1) r^l p'(r) + r^{l+1} (p'^2(r) + p''(r))\right] e^{p(r)}

其中：

.. math::

   p'(r) &= \sum_{i=1}^N 2i \, a_{2i} r^{2i-1} \\\\
   p''(r) &= \sum_{i=1}^N 2i(2i-1) a_{2i} r^{2i-2}

范数积分
~~~~~~~~

内区范数：

.. math::

   Q = \int_0^{r_c} u^2(r) dr = \int_0^{r_c} r^{2(l+1)} e^{2p(r)} dr

这是一个超越积分，需要数值求解（使用 Simpson 法则或 Gauss 积分）。

非线性方程组
~~~~~~~~~~~~

定义残差函数 :math:`\mathbf{F}(\mathbf{a})` 如下（:math:`\mathbf{a} = [a_0, a_2, a_4, a_6]^T`）：

.. math::

   F_1 &= u_{\text{PS}}(r_c) - u_{\text{AE}}(r_c) \\\\
   F_2 &= u'_{\text{PS}}(r_c) - u'_{\text{AE}}(r_c) \\\\
   F_3 &= u''_{\text{PS}}(r_c) - u''_{\text{AE}}(r_c) \\\\
   F_4 &= Q_{\text{PS}} - Q_{\text{AE}}

求解 :math:`\mathbf{F}(\mathbf{a}) = \mathbf{0}`，使用 **scipy.optimize.fsolve**（混合 Powell 算法）。

数值实现
--------

初值选择
~~~~~~~~

系数初值对收敛性影响很大。推荐初值：

.. math::

   a_0 &\approx \ln\left(\frac{u_{\text{AE}}(r_c)}{r_c^{l+1}}\right) \\\\
   a_{2i} &= 0, \quad i = 1, 2, 3

内外区拼接
~~~~~~~~~~

完整伪轨道：

.. math::

   u_{\text{PS}}(r) = \begin{cases}
   r^{l+1} \exp\left( \sum_{i=0}^N a_{2i} r^{2i} \right), & r \leq r_c \\\\
   u_{\text{AE}}(r), & r > r_c
   \end{cases}

数值稳定性
~~~~~~~~~~

1. **指数溢出保护**：

   在计算 :math:`e^{p(r)}` 时，如果 :math:`p(r)` 过大（:math:`> 100`），可能溢出。
   使用对数形式处理：

   .. math::

      \ln u(r) = (l+1) \ln r + p(r)

2. **范数积分精度**：

   使用自适应积分（如 :code:`scipy.integrate.quad`）或高阶 Simpson 法则。

3. **求解器容差**：

   :code:`fsolve` 的容差设为 :code:`xtol=1e-8, maxfev=1000`。

数值实现指南
------------

伪代码框架
~~~~~~

.. code-block:: python

   def tm_pseudize(r, w, u_ae, eps, l, rc, N=2):
       # 1. 找到 rc 对应的网格索引
       i_rc = find_index(r, rc)

       # 2. 计算 AE 在 rc 处的导数
       u_rc, du_rc, d2u_rc = compute_derivatives(r, u_ae, i_rc)

       # 3. 计算 AE 内区范数
       Q_ae = simpson_integral(u_ae[:i_rc]**2, r[:i_rc])

       # 4. 定义残差函数
       def residuals(a):
           u_ps, du_ps, d2u_ps = eval_tm_at_rc(rc, l, a)
           Q_ps = compute_tm_norm(r[:i_rc], l, a)
           return [
               u_ps - u_rc,
               du_ps - du_rc,
               d2u_ps - d2u_rc,
               Q_ps - Q_ae
           ]

       # 5. 求解
       a_init = [log(u_rc / rc**(l+1)), 0, 0, 0]
       a_solution = fsolve(residuals, a_init)

       # 6. 拼接
       u_ps = splice(u_ae, a_solution, l, i_rc)

       return u_ps, a_solution

精度检验
--------

范数守恒误差
~~~~~~~~~~~~

定义相对误差：

.. math::

   \delta_{\text{norm}} = \frac{|Q_{\text{PS}} - Q_{\text{AE}}|}{Q_{\text{AE}}}

**验收标准**：:math:`\delta_{\text{norm}} < 10^{-5}`

连续性误差
~~~~~~~~~~

检查在 :math:`r_c` 处的相对误差：

.. math::

   \delta_k = \frac{|u_{\text{PS}}^{(k)}(r_c) - u_{\text{AE}}^{(k)}(r_c)|}{|u_{\text{AE}}^{(k)}(r_c)|}

**验收标准**：:math:`\delta_k < 10^{-4}` for :math:`k = 0, 1, 2`

参考文献
--------

- **原始论文**：N. Troullier and J. L. Martins,
  "Efficient pseudopotentials for plane-wave calculations,"
  *Phys. Rev. B* **43**, 1993-2006 (1991).
  方程 (14)-(18) 定义了 TM 方法的核心公式。

- **QE 实现指南**：P. Giannozzi, *Notes on pseudopotential generation* (2019),
  第 2.2-2.3 节详细说明了 TM 伪化、截断半径选择与实现要点。
  https://www.quantum-espresso.org/wp-content/uploads/2022/03/pseudo-gen.pdf
