半局域势反演
============

本节介绍从 TM 伪轨道反演半局域势 :math:`V_l(r)` 的技术细节。

物理基础
--------

反演的物理含义
~~~~~~~~~~~~~~

势反演并不是“拟合”一个势，而是把径向 Schrödinger 方程 **反过来当作对未知势的代数方程**。
对某个角动量通道 :math:`l`，一旦你已经构造了满足外区匹配条件的伪轨道 :math:`u_l(r)`，
并确定了对应能量 :math:`\varepsilon_l`，那么在每个半径点上，势 :math:`V_l(r)` 必须取某个值，
才能让这条曲线成为该能量下的本征解。

从反演公式

.. math::

   V_l(r) = \varepsilon_l + \frac{1}{2}\frac{u''(r)}{u(r)} - \frac{l(l+1)}{2r^2}

可以直接看到两点直觉：

1. :math:`u''/u` 反映 **局域曲率**。在某个区域如果 :math:`u(r)` “弯得更厉害”（曲率大），
   动能项的局域贡献也更强；为了在同一个 :math:`\varepsilon_l` 下保持方程成立，
   势项必须相应地补偿（这就是反演在数值上对二阶导数很敏感的原因）。
2. :math:`-l(l+1)/(2r^2)` 是 **离心势**。它把高角动量通道在小 :math:`r` 处强烈排斥开来，
   因而即使外区波函数相同，不同 :math:`l` 通道也会得到不同的 :math:`V_l(r)`。

这也解释了“半局域势（semilocal）”这个名字：它在径向坐标上是局域的 :math:`V_l(r)`，
但在角动量空间上分通道。把这种 :math:`l` 依赖的算符放到一般三维基组（尤其是平面波）中时，
就会表现为非局域作用，这正是后续 KB 可分离形式要解决的计算瓶颈。

反演公式
~~~~~~~~

从径向 Schrödinger 方程出发：

.. math::

   \left[-\frac{1}{2}\frac{d^2}{dr^2} + \frac{l(l+1)}{2r^2} + V_l(r)\right] u(r) = \varepsilon u(r)

改写为：

.. math::

   \frac{d^2 u}{dr^2} = \left[\frac{l(l+1)}{r^2} + 2V_l(r) - 2\varepsilon\right] u(r)

解出半局域势：

.. math::

   V_l(r) = \varepsilon + \frac{1}{2}\frac{u''(r)}{u(r)} - \frac{l(l+1)}{2r^2}

数值实现
--------

内区（r ≤ rc）：解析导数
~~~~~~~~~~~~~~~~~~~~~~~~

对于 TM 形式的伪轨道：

.. math::

   u(r) = r^{l+1} \exp\left(\sum_{i=0}^{N} a_{2i} r^{2i}\right)

定义 :math:`p(r) = \sum_{i=0}^{N} a_{2i} r^{2i}`，则：

.. math::

   u'(r) &= \left[\frac{l+1}{r} + p'(r)\right] u(r) \\
   u''(r) &= \left[\frac{l(l+1)}{r^2} + \frac{2(l+1)p'(r)}{r} + p'^2(r) + p''(r)\right] u(r)

因此：

.. math::

   \frac{u''(r)}{u(r)} = \frac{l(l+1)}{r^2} + \frac{2(l+1)p'(r)}{r} + p'^2(r) + p''(r)

代入反演公式：

.. math::

   V_l(r) = \varepsilon + \frac{l+1}{r}p'(r) + \frac{1}{2}\left[p'^2(r) + p''(r)\right]

其中：

.. math::

   p'(r) &= \sum_{i=1}^{N} 2i \, a_{2i} r^{2i-1} \\
   p''(r) &= \sum_{i=1}^{N} 2i(2i-1) a_{2i} r^{2i-2}

**优点**：无需数值微分，精度高，数值稳定。

外区（r > rc）：样条法导数
~~~~~~~~~~~~~~~~~~~~~~~~~~

对于外区（AE 轨道），使用样条插值计算导数：

1. 在每个点 :math:`r_i` 附近构建三次样条（窗口 7 点）
2. 计算 :math:`u(r_i)` 和 :math:`u''(r_i)`
3. 代入公式：

.. math::

   V_l(r_i) = \varepsilon + \frac{1}{2}\frac{u''(r_i)}{u(r_i)} - \frac{l(l+1)}{2r_i^2}

**优点**：适用于任意形式的轨道，无需解析表达式。

节点保护
--------

节点检测
~~~~~~~~

节点定义为 :math:`|u(r)| < \epsilon_{\text{tol}}`（默认 :math:`\epsilon_{\text{tol}} = 10^{-10}`）。

在节点附近，:math:`u''(r)/u(r)` 可能发散，需要特殊处理。

插值策略
~~~~~~~~

1. 标记节点位置：:math:`|u(r_i)| < \epsilon_{\text{tol}}`
2. 在非节点区域正常计算势值
3. 对节点区域使用线性插值填充：

.. math::

   V_l(r_{\text{node}}) = \text{interp1d}(r_{\text{valid}}, V_{\text{valid}}, r_{\text{node}})

势裁剪
------

为防止数值奇异性，对势值进行裁剪：

.. math::

   V_l(r) \in [-V_{\max}, V_{\max}]

默认 :math:`V_{\max} = 1000` Hartree。

数值验证
--------

内外区连续性
~~~~~~~~~~~~

在 :math:`r = r_c` 处，内区（解析）和外区（样条）应给出连续的势值：

.. math::

   \left|\frac{V_l^{\text{inner}}(r_c) - V_l^{\text{outer}}(r_c)}{V_l^{\text{inner}}(r_c)}\right| < 0.1

物理合理性
~~~~~~~~~~

1. **s 轨道** （l=0）：原点附近势应平缓（无离心项）
2. **p/d 轨道**：势应有限，不应有极端奇异性
3. **远场行为**：:math:`V_l(r) \to 0` as :math:`r \to \infty`

参考文献
--------

- **TM 方法**: N. Troullier and J. L. Martins, *PRB* **43**, 1993 (1991)
- **势反演**: P. Giannozzi, *Notes on pseudopotential generation* (2019), Section 3.1
- **样条插值**: Press et al., *Numerical Recipes*, Chapter 3
