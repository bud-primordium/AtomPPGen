KB 转换模块 (kb)
==================

概述
----
将半局域势 :math:`\{V_l(r)\}` 转换为 Kleinman-Bylander 可分离形式：

.. math::

   V_{\text{NL}} = \sum_l |\beta_l\rangle D_l \langle\beta_l|

其中 :math:`\beta_l(r)` 是投影子，:math:`D_l` 是耦合系数。

核心功能
--------
.. autofunction:: kb_transform

数据类
------
.. autoclass:: KBResult
   :members:

算法说明
--------

投影子构造
^^^^^^^^^^

对于非局域通道 :math:`l \neq l^*`，构造未归一化投影子：

.. math::

   \chi_l(r) = [V_l(r) - V_{\text{loc}}(r)] \cdot \phi_l(r)

其中：

* :math:`V_{\text{loc}}(r) = V_{l^*}(r)` 是局域势（通常选择未占据的高角动量通道）
* :math:`\phi_l(r) = u_l(r)/r` 是径向规约波函数

归一化条件
^^^^^^^^^^

投影子归一化为 :math:`\langle\beta_l|\beta_l\rangle = 1`：

.. math::

   \beta_l(r) = \frac{\chi_l(r)}{\sqrt{\langle\chi_l|\chi_l\rangle}}

耦合系数
^^^^^^^^

标准 KB 形式的耦合系数为：

.. math::

   D_l^{\text{std}} = \frac{1}{\langle\phi_l|\Delta V|\phi_l\rangle}

其中 :math:`\Delta V = V_l - V_{\text{loc}}`。

本实现采用归一化投影子，因此耦合系数调整为：

.. math::

   D_l = \frac{\langle\chi_l|\chi_l\rangle}{\langle\phi_l|\Delta V|\phi_l\rangle} = \frac{W}{Z}

其中：

* :math:`W = \int |\chi_l|^2 w \, dr`
* :math:`Z = \int \Delta V \cdot \phi_l^2 w \, dr`

这两种形式物理等价，归一化形式在数值上更稳定。

数值稳定性
^^^^^^^^^^

1. **原点保护**：在 :math:`r \to 0` 时，:math:`u_l \sim r^{l+1}`，因此 :math:`\phi_l = u_l/r \sim r^l` 有限。
   使用 ``np.maximum(r, 1e-10)`` 保护除法。

2. **积分权重**：使用提供的积分权重 :math:`w` 计算内积，保证数值精度。

3. **除零检查**：计算 :math:`D_l = W/Z` 时检查 :math:`|Z| < 10^{-12}`，避免分母接近零。

参考文献
--------

* Kleinman & Bylander, *PRL* 48, 1425 (1982)
* Giannozzi, *Notes on pseudopotential generation* (2019)
