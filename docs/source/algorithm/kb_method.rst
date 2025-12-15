Kleinman-Bylander 可分离形式
============================

在 TM 伪化与势反演之后，我们得到各角动量通道的半局域势 :math:`V_l(r)`。本节解释：

- 为什么半局域势在平面波基组中不方便；
- KB 可分离形式如何把它改写成“局域势 + 少量投影子”的结构，从而显著加速；
- 投影子 :math:`\beta_l` 与耦合系数 :math:`D_l` 的物理含义；
- 为什么局域道通常选择未占据的高角动量通道（例如 Al 的 d 通道）。

物理动机
--------

半局域赝势（semilocal pseudopotential）可以写成角动量投影算符的和：

.. math::

   \hat V^{\mathrm{SL}}
   = \sum_{l,m} |Y_{lm}\rangle \, V_l(r) \, \langle Y_{lm}|

这里 :math:`Y_{lm}(\hat{\mathbf r})` 是球谐函数，:math:`V_l(r)` 只依赖半径但**依赖角动量通道**。

**为什么半局域势在平面波基组中是问题？**

平面波基函数是 :math:`\phi_{\mathbf G}(\mathbf r)=e^{i\mathbf G\cdot\mathbf r}`，它没有确定的 :math:`l,m`。
因此如果要在平面波表示中作用 :math:`\hat V^{\mathrm{SL}}`，就需要把平面波在原子附近展开到球谐函数上，
并对每个 :math:`l,m` 分量分别施加不同的 :math:`V_l(r)`。从实现角度看，这会引入昂贵的角动量分解与通道求和，
并让“对每个波函数施加非局域势”的代价接近（甚至退化到）对大矩阵的稠密作用。

KB 的核心思想是：把“按 :math:`l` 选择不同势”的复杂结构，压缩成少数几个**秩-1（可分离）**投影修正。

从半局域到可分离
----------------

第一步选择一个通道 :math:`l^*` 作为**局域势**：

.. math::

   V_{\mathrm{loc}}(r) = V_{l^*}(r), \qquad \Delta V_l(r) = V_l(r) - V_{\mathrm{loc}}(r)

则半局域势可改写为

.. math::

   \hat V^{\mathrm{SL}} = V_{\mathrm{loc}}(r)
   + \sum_{l\ne l^*}\sum_m |Y_{lm}\rangle\,\Delta V_l(r)\,\langle Y_{lm}|

剩下的问题是如何高效处理 :math:`\Delta V_l` 这部分。

KB 的做法是为每个 :math:`l` 选择一个参考径向函数 :math:`\phi_l(r)` (通常取伪化参考态的径向规约波函数)，
并构造

.. math::

   |\chi_l\rangle = \Delta V_l\,|\phi_l\rangle

然后把对应通道的非局域修正写成可分离算符：

.. math::

   \hat V^{\mathrm{KB}}_{\mathrm{NL}}
   = \sum_{l\ne l^*}\sum_m |\beta_{lm}\rangle\, D_l\, \langle\beta_{lm}|

其中 :math:`|\beta_{lm}\rangle` 是投影子（见下一节）。关键点在于，
:math:`|\beta\rangle D\langle\beta|` 是秩-1 结构：

.. math::

   \langle \psi_i |\hat V_{\mathrm{NL}}^{\mathrm{KB}}| \psi_j \rangle
   = \sum_{l\ne l^*}\sum_m
   \bigl(\langle \psi_i | \beta_{lm} \rangle\bigr)\, D_l \, \bigl(\langle \beta_{lm} | \psi_j \rangle\bigr)

**可分离形式为什么能加速计算？**

半局域形式下，平面波并不携带明确的 :math:`l,m`，往往需要做球谐展开并对通道求和；
而可分离形式把非局域势的矩阵元变成“两个标量积的乘积”。在平面波计算中，
对每个波函数只需计算若干个重叠积分 :math:`\langle\beta|\psi\rangle`，再乘以标量 :math:`D_l` 回写，
典型代价可从“近似 :math:`O(N_G^2)` 的通道作用”降低到“近似 :math:`O(N_G)` 的投影累加”，
从而显著加速。

投影子构造
----------

在 AtomPPGen 的径向实现里，我们使用径向规约波函数

.. math::

   \phi_l(r) = \frac{u_l(r)}{r}

并定义未归一化投影子（代码中也称 :math:`\chi_l`）：

.. math::

   \chi_l(r) = \Delta V_l(r)\,\phi_l(r)

令

.. math::

   W_l = \langle\chi_l|\chi_l\rangle = \int_0^{\infty} |\chi_l(r)|^2\, w(r)\,dr
   \qquad
   Z_l = \langle\phi_l|\Delta V_l|\phi_l\rangle
       = \int_0^{\infty} \Delta V_l(r)\,\phi_l(r)^2\, w(r)\,dr

则归一化投影子与耦合系数取为

.. math::

   \beta_l(r) = \frac{\chi_l(r)}{\sqrt{W_l}},
   \qquad
   D_l = \frac{W_l}{Z_l}

这与 `AtomPPGen/src/atomppgen/kb.py` 的实现一致。

**投影子 :math:`\beta_l` 的物理意义是什么？**

:math:`\Delta V_l = V_l - V_{\mathrm{loc}}` 是“该通道相对局域势的差异”，它主要局域在核附近；
因此 :math:`\beta_l \propto \Delta V_l\,\phi_l` 会自动把波函数中**需要非局域修正的散射分量**
“挑出来”。换句话说，KB 非局域项只对那些“在核附近、并且属于该角动量通道”的成分施加修正，
在核外区域几乎不做任何事情。

局域道选择
----------

**为什么局域道常选未占据通道？**

- 占据的价电子通道（例如 s、p）直接决定基态电子结构与化学键，通常希望用非局域投影更精确地再现其散射性质；
- 未占据的高角动量通道（例如许多元素的 d 或 f）对基态贡献小，更适合作为“背景局域势”，
  把主要误差留给非局域修正去补偿；
- 在实际经验中，选择高 :math:`l` 通道作为 :math:`V_{\mathrm{loc}}` 往往还能降低某些幽灵态风险。

**对 Al 为什么常选 d 通道？**

Al 的基态主要由 3s/3p 决定，而 d 通道通常作为散射投影通道使用。
把 d 选作局域势意味着：s/p 的差异通过非局域投影显式校正，而 d 作为“较不敏感”的背景势，
在效率与稳定性上往往更合适。

数值实现
--------

AtomPPGen 的 `kb_transform()` 做的是径向层面的 KB 构造（为后续平面波程序的非局域项准备投影子数据）。
其核心步骤可以概括为：

.. code-block:: python

   # 输入：各通道半局域势 V_by_l、参考径向波函数 u_by_l、网格 r 与权重 w
   # 选择局域通道 l*（例如 d 通道）
   loc_channel = 2
   V_loc = V_by_l[loc_channel]

   for l, V_l in V_by_l.items():
       if l == loc_channel:
           continue

       r_safe = np.maximum(r, 1e-10)
       phi_l = u_by_l[l] / r_safe
       delta_V = V_l - V_loc
       chi = delta_V * phi_l

       W = np.sum(chi**2 * w)
       Z = np.sum(delta_V * (phi_l**2) * w)

       beta_l = chi / np.sqrt(W)
       D_l = W / Z

实现中还包含：

- :math:`r\to 0` 处的除法保护（用 `np.maximum(r, 1e-10)`）；
- 对 :math:`W_l\le 0` 或 :math:`Z_l\approx 0` 的显式异常提示，避免静默产生不可用投影子。

参考文献
--------

- L. Kleinman and D. M. Bylander, *Phys. Rev. Lett.* **48**, 1425 (1982)
- P. Giannozzi, *Notes on pseudopotential generation* (2019)
