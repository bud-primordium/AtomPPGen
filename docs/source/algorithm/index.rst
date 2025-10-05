算法原理
========

本节详细介绍赝势生成的数学原理和数值方法。

.. toctree::
   :maxdepth: 2

   ae_solve
   tm_method
   potential_inversion
   kb_transformation
   transferability

概述
----

赝势生成流程包含以下关键步骤：

1. **全电子原子解** - 获取参考态波函数和能级
2. **TM 伪化** - 在截断半径内构造光滑伪轨道
3. **势反演** - 从伪轨道反演半局域势
4. **KB 转换** - 构造可分离非局域投影
5. **可转移性检验** - 验证赝势质量

理论基础
--------

径向 Schrödinger 方程
~~~~~~~~~~~~~~~~~~~~~

对于球对称势，径向波函数 :math:`u_{nl}(r) = r R_{nl}(r)` 满足：

.. math::

   \left[-\frac{1}{2}\frac{d^2}{dr^2} + V_l(r) + \frac{l(l+1)}{2r^2}\right] u_{nl}(r) = \varepsilon_{nl} u_{nl}(r)

其中 :math:`V_l(r)` 是半局域势（依赖于角动量 :math:`l`）。

范数守恒条件
~~~~~~~~~~~~

在截断半径 :math:`r_c` 内，伪轨道 :math:`\tilde{u}_{nl}` 与全电子轨道 :math:`u_{nl}` 的范数相等：

.. math::

   \int_0^{r_c} |\tilde{u}_{nl}(r)|^2 dr = \int_0^{r_c} |u_{nl}(r)|^2 dr

这保证了电荷分布在 :math:`r \leq r_c` 内的一致性。

可分离非局域形式
~~~~~~~~~~~~~~~~

完整的非局域势可表示为：

.. math::

   V_{\text{NL}}(r, r') = V_{\text{loc}}(r) \delta(r - r') + \sum_l |\beta_l\rangle D_l \langle\beta_l|

其中：

- :math:`V_{\text{loc}}(r)` - 局域势
- :math:`|\beta_l\rangle` - 投影算符
- :math:`D_l` - 耦合系数（单位：Rydberg）
