赝势可转移性验证方法
========================

本节介绍赝势质量评估的三类验证方法：范数守恒、对数导数匹配和幽灵态检测。

范数守恒验证
------------

物理意义
~~~~~~~~

范数守恒（Norm Conservation）要求伪轨道在截断半径 :math:`r_c` 内的归一化
与全电子轨道一致，保证赝势保留了真实电荷分布的 **局域信息**：

.. math::

   \int_0^{r_c} |u_{\text{PS}}(r)|^2 dr = \int_0^{r_c} |u_{\text{AE}}(r)|^2 dr

其中 :math:`u(r) = r R(r)` 为径向波函数。

数值实现
~~~~~~~~

TM 伪化方法通过非线性方程组隐式保证范数守恒。验证时提取 TM 结果
中的 ``norm_error`` 字段：

.. math::

   \Delta_{\text{norm}} = \frac{Q_{\text{PS}} - Q_{\text{AE}}}{Q_{\text{AE}}}

验收标准：:math:`|\Delta_{\text{norm}}| < 10^{-6}`。

对数导数匹配
------------

对数导数的物理意义
~~~~~~~~~~~~~~~~~~

对数导数（Logarithmic Derivative）:math:`L(E, r)` 描述散射性质在能量窗口内的
变化，是可转移性的核心指标：

.. math::

   L(E, r) = r \frac{\partial \ln \psi(E, r)}{\partial r} = r \frac{\psi'(E, r)}{\psi(E, r)}

为什么它能反映散射相移？

当 :math:`r` 选在核外（对赝势而言通常取 :math:`r \gtrsim r_c`）时，
径向散射解可以写成带相移 :math:`\delta_l(E)` 的渐近形式（忽略长程库仑修正时）：

.. math::

   u_l(r) \propto \sin\bigl(kr - l\pi/2 + \delta_l(E)\bigr), \quad k = \sqrt{2E}

代入 :math:`L(E,r) = r u'(r)/u(r)` 得到

.. math::

   L_l(E, r) = kr\,\cot\bigl(kr - l\pi/2 + \delta_l(E)\bigr)

因此在固定 :math:`r=r_{\text{test}}` 下，:math:`L_l(E)` 与 :math:`\delta_l(E)` 一一对应
（本质上是把“相位信息”映射成“斜率/幅度比”的能量函数）。这就是为什么比较
AE 与 PS 的 :math:`L(E)` 曲线，等价于比较它们在该能量窗口内的散射相移。

为什么用 :math:`L(E)` 而不是直接算 :math:`\delta_l(E)`？

- :math:`L(E)` 只需要在一个有限半径处取 :math:`u'/u`，不依赖波函数归一化，数值上更稳健。
- 直接求 :math:`\delta_l(E)` 往往需要更“远场”的渐近拟合，并且要处理相位的 :math:`\pi` 分支问题；
  在教学实现里，用 :math:`L(E)` 作为中间量更直接。

“零点”为什么有用？

AtomPPGen 在评价指标中使用的是 :math:`L(E)` 的 **过零点** (:math:`L=0`) 位置。
由上式可知，:math:`L=0 \iff \cot(\cdots)=0`，也就是相位满足

.. math::

   kr_{\text{test}} - l\pi/2 + \delta_l(E) = \left(n + \tfrac{1}{2}\right)\pi

这是一条“量子化条件”：它把连续的相移信息离散化成一串能量刻度。
如果赝势在某个能量附近的散射行为（相移）出现偏差，过零点会整体平移；
而在相移随能量快速变化的区域（常见于准束缚/共振特征附近），过零点位置尤其敏感。
因此，比较过零点能量的偏差是一种紧凑而有效的可转移性指标。

在测试半径 :math:`r_{\text{test}}` 处（通常选取 :math:`r_{\text{test}} = \max(r_c^l) + 0.5` a₀），
扫描能量窗口 :math:`E \in [-0.25, +0.25]` Ha（对应 [-0.5, +0.5] Ry），比较
全电子和赝势的 :math:`L(E)` 曲线。

评价指标
~~~~~~~~

1. **零点均方根偏差** （Zero-Crossing RMS）

   对数导数的零点对应固定相位条件（可视为散射相移的离散刻度）。匹配零点位置确保赝势在不同化学环境下
   散射行为一致：

   .. math::

      \Delta E_{\text{RMS}} = \sqrt{\frac{1}{N} \sum_{i=1}^N (E_i^{\text{AE}} - E_i^{\text{PS}})^2}

   **验收阈值**：:math:`\Delta E_{\text{RMS}} < 0.025` Ha (:math:`\approx 0.05` Ry)

   **实现细节（重要）**：

   在当前教学实现中，只有当扫描窗口内 AE/PS 都至少出现 2 个过零点时，才会计算并启用零点 RMS 判据。
   若零点数量不足（常见于某些通道/某些 :math:`r_{\text{test}}` 选择），零点 RMS 记为 N/A，不参与通过判定，
   此时主要依据价区曲线 RMS 进行评估。

2. **曲线均方根差异** （Curve RMS）

   整体形状匹配度：

   .. math::

      L_{\text{RMS}} = \sqrt{\frac{1}{M} \sum_{j=1}^M [L_{\text{AE}}(E_j) - L_{\text{PS}}(E_j)]^2}

   **验收阈值（元素类型差异化）**：

   - **金属元素** （Al, Na, Mg）：:math:`L_{\text{RMS}}^{\text{valence}} < 16.0`
   - **共价元素** （Si, C, N）：:math:`L_{\text{RMS}}^{\text{valence}} < 3.0`

   其中 :math:`L_{\text{RMS}}^{\text{valence}}` 为价区窗口 :math:`E \in [-0.05, +0.05]` Ha 内的曲线 RMS。

   **物理依据**：

   金属元素在远离核区（:math:`r \sim r_c`）的软库仑势中，价电子对数导数
   :math:`L_{\text{AE}}(E, r) \approx 0` 几乎不随能量变化（曲率小）。赝势在此过渡区的
   相位差异被放大，导致较大的曲线 RMS，这是金属元素固有特性而非赝势质量缺陷。

   共价元素的波函数节点清晰、曲率大，AE-PS 匹配较容易，可使用更严格的阈值。

   **价区验证建议**：

   - **能量步长推荐**：:math:`E_{\text{step}} \leq 0.02` Ry（保证价区至少 10 个采样点）
   - **统计下限**：价区点数 :math:`\geq 3`（:math:`E_{\text{step}} = 0.05` Ry 时的临界值）

   **全窗口 RMS**：

   能量窗口 :math:`E \in [-0.25, +0.25]` Ha 的全曲线 RMS 仅作警告参考（金属元素可达 20-25）。

径向薛定谔方程求解（Numerov 方法）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

对数导数计算需要在每个能量点 :math:`E` 求解径向薛定谔方程：

.. math::

   \left[ -\frac{1}{2} \frac{d^2}{dr^2} + V(r) + \frac{l(l+1)}{2r^2} \right] u(r) = E \cdot u(r)

**Numerov 方法特点**：

- **精度**：:math:`O(h^5)` 截断误差
- **网格要求**：等距步长 :math:`h`
- **迭代公式**：

  .. math::

     (1 + \tfrac{h^2}{12} k_{n+1}^2) u_{n+1} =
     2(1 - \tfrac{5h^2}{12} k_n^2) u_n -
     (1 + \tfrac{h^2}{12} k_{n-1}^2) u_{n-1}

  其中 :math:`k^2(r) = 2[E - V_{\text{eff}}(r)]`

**非均匀网格处理**：

若输入网格非均匀（如 ``exp_transformed`` 类型），先用三次样条插值重采样到
等距网格，求解后再插值回原网格。

**初值选择**：

- s 轨道（:math:`l=0`）：:math:`u(h) = h`
- 高角动量（:math:`l > 0`）：:math:`u(h) = h^{l+1}`

**单位约定**：

所有能量和势能使用 **Hartree 原子单位**。外部 API 支持 Rydberg 输入，
内部自动转换（:math:`E_{\text{Ry}} = 2 E_{\text{Ha}}`）。

KS 有效势提取
~~~~~~~~~~~~~

全电子对数导数使用 **Kohn-Sham 有效势**：

.. math::

   V_{\text{AE}}(r) = v_{\text{ext}}(r) + v_H[n](r) + v_{xc}[n](r)

其中：

- :math:`v_{\text{ext}} = -Z/r`：核-电子吸引
- :math:`v_H = \int \frac{n(r')}{|r - r'|} dr'`：Hartree 势（电子-电子排斥）
- :math:`v_{xc}`：交换关联势（PZ81 或 VWN 泛函）

**注意**：离心项 :math:`l(l+1)/(2r^2)` **不包含** 在 :math:`V_{\text{AE}}` 中，
由径向求解器内部添加。

伪势使用反演得到的半局域势 :math:`V_{\text{PS}}^l(r)`。

幽灵态检测
----------

物理意义
~~~~~~~~

幽灵态（Ghost State）是赝势哈密顿量中 **非物理的深束缚态**，
出现在价电子能级附近，导致赝势在某些环境下不可用。

检测原理（A 级：径向方法）
~~~~~~~~~~~~~~~~~~~~~~~~~~~

对每个角动量通道 :math:`l`，构建径向哈密顿矩阵：

.. math::

   H_l = T + V_{\text{PS}}^l(r) + \frac{l(l+1)}{2r^2}

**有限差分构造**：

在均匀网格上，动能算子 :math:`T = -\frac{1}{2} \frac{d^2}{dr^2}`
离散为三点有限差分：

.. math::

   T_{ij} = \begin{cases}
   \frac{1}{\Delta r^2} + V_{\text{eff}}(r_i), & i = j \\
   -\frac{1}{2 \Delta r^2}, & |i - j| = 1 \\
   0, & \text{otherwise}
   \end{cases}

其中 :math:`V_{\text{eff}}(r) = V_{\text{PS}}^l(r) + \frac{l(l+1)}{2r^2}`。

**对角化与筛选**：

1. 对 :math:`H_l` 进行厄米对角化（``scipy.linalg.eigh``）
2. 筛选束缚态：

   - 能量 :math:`E < E_{\text{max}}` （默认 0 Ha）
   - 波函数边界条件：:math:`\left|\psi(r_{\text{max}})\right| < 0.1 \cdot \max_r \left|\psi(r)\right|` （盒态过滤）

3. 识别幽灵态：

   在能量窗口 :math:`E \in [-0.15, +0.05]` Ha 内，采用能量感知分类判定。

**盒态过滤逻辑**：

有限网格在边界处形成"无穷势阱"，产生人工束缚态（**盒态**，Box States）。
盒态特征为波函数在 :math:`r_{\text{max}}` 处未充分衰减。

**尾部比例检测**：

.. math::

   \tau = \frac{|\psi(r_{\text{max}})|}{\max_r |\psi(r)|}

- :math:`\tau < 0.1`：真束缚态（物理或幽灵）
- :math:`\tau > 0.1`：盒态（网格人工产物，不计入幽灵态数）

**能量感知分类**：

仅基于 :math:`\tau` 无法区分危险幽灵态与安全的 Rydberg 激发态。改进判据如下：

1. **正能散射态** （:math:`E > 0`）：

   连续谱在有限盒子中的离散化产物，数学上非真正束缚态，归为盒态。

2. **Rydberg 激发态** （:math:`0 > E > \varepsilon_{\text{valence}} - \delta`）：

   高主量子数束缚态序列（如 Al 的 4s, 5s, 6s...），能量趋于电离阈值 0 Ha。
   这些态位于价电子能级之上，距离费米能级较远，对基态 DFT 计算无影响，归为安全态。

3. **潜在危险幽灵态** （:math:`E < \varepsilon_{\text{valence}} - \delta`）：

   能量显著低于价态的额外束缚态，可能占据基态，导致错误电子结构。
   需进一步用 :math:`\tau` 判据区分真幽灵态与盒态。

能量容差 :math:`\delta = 0.01` Ha 用于排除数值误差附近的态。

**物理意义**：

Rydberg 激发态的存在证明赝势保留了正确的长程库仑行为（:math:`-1/r`）。
将这些态误判为幽灵态会导致对赝势质量的错误评估。能量感知分类确保只标记
真正危险的态（价态下方的非物理束缚态）。

**验收标准**：真幽灵态数量 :math:`N_{\text{ghost}} \leq 10` （TM 伪化产生的浅幽灵态
能量接近 0，对基态 DFT 影响有限）。

**数值优化**：

- 非均匀网格重采样到 300 点均匀网格（加速对角化）
- 矩阵大小限制防止内存溢出

B 级方法（小球平面波）
~~~~~~~~~~~~~~~~~~~~~~

**可选实现**。在小球半径 :math:`R_{\text{cut}}` 内构建平面波基组，
对角化包含非局域 KB 投影子的完整哈密顿量：

.. math::

   H = T + V_{\text{loc}} + \sum_{lm} |\beta_{lm}\rangle D_l \langle \beta_{lm}|

该方法更严格，但计算成本高。

完整验证流程
------------

函数：``run_full_validation()``

输入
~~~~

- ``ae_result``: 全电子原子解（``AEAtomResult``）
- ``tm_dict``: 各通道 TM 伪化结果（``Dict[int, TMResult]``）
- ``inv_dict``: 各通道势反演结果（``Dict[int, InvertResult]``）
- ``r_test``: 对数导数测试半径（默认 3.0 a₀）
- ``E_range_Ry``: 能量窗口（Rydberg，默认 [-0.5, 0.5]）
- ``E_step_Ry``: 能量步长（Rydberg，默认 0.05）

流程
~~~~

1. 提取 KS 有效势（所有通道共享）
2. **范数守恒检验**：对每个通道调用 ``check_norm_conservation()``
3. **对数导数匹配**：对每个通道调用 ``check_log_derivative()``
4. **幽灵态检测**：对每个通道调用 ``check_ghost_states()``
5. 汇总结果，生成 ``ValidationReport``

输出
~~~~

``ValidationReport`` 包含：

- ``norm_results``: 范数守恒结果字典（按 :math:`l` 索引）
- ``log_deriv_results``: 对数导数结果字典
- ``ghost_result``: 幽灵态结果（代表性通道）
- ``overall_passed``: 整体判定（所有检验均通过）
- ``diagnostics``: 诊断信息（通道数、测试参数、分项通过状态）

**JSON 导出**：

.. code-block:: python

   report = run_full_validation(ae, tm_dict, inv_dict)
   with open('outputs/validation_report.json', 'w') as f:
       json.dump(report.to_dict(), f, indent=2)

参考文献
--------

- Troullier & Martins, *PRB* **43**, 1993 (1991) - 范数守恒条件
- Gonze et al., *Comput. Mater. Sci.* **25**, 478 (2002) - 对数导数方法
- Rappe et al., *PRB* **41**, 1227 (1990) - 幽灵态检测
- Numerov, *Trudy Glav. Astron. Obs.* **28**, 173 (1926) - Numerov 方法
