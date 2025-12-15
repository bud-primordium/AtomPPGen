================
验证标准解读
================

本节面向第一次阅读 ``ValidationReport.summary()`` 的同学，按照范数守恒、对数导数与幽灵态检测
三条主线梳理验收阈值、公式来源以及金属/共价元素的差异化标准。

范数守恒
--------

**定义**：截断半径 :math:`r_c` 内伪轨道与全电子轨道的电荷积分必须一致。

.. math::

   Q_l^{\text{AE}}(r_c) = \int_0^{r_c} |u_l^{\text{AE}}(r)|^2 dr,
   \qquad
   Q_l^{\text{PS}}(r_c) = \int_0^{r_c} |u_l^{\text{PS}}(r)|^2 dr

.. math::

   \Delta_{\text{norm}} = \frac{Q_l^{\text{PS}} - Q_l^{\text{AE}}}{Q_l^{\text{AE}}}

**阈值**：

.. math::

   |\Delta_{\text{norm}}| < 10^{-6}

**实践提示**：

- 参数扫描的 9 组 rc 组合中，所有通道的 ``max_norm_error`` 都低至 :math:`10^{-13}`，
  说明 TM 求解器能把范数误差控制在数值舍入范围，若出现 :math:`10^{-4}` 量级误差
  必然是伪化失败（参数或代码层面）。
- 通过 ``ValidationReport.summary()`` 查看 ``norm_error.s/p/d``，若单通道失败，
  可优先减小对应的 :math:`r_c`。

对数导数
--------

**目的**：保持散射相移一致，从而保证赝势在不同化学环境下的转移性。
对数导数定义为：

.. math::

   L_l(E, r) = r \frac{d}{dr} \ln \psi_l(E, r) = r \frac{\psi_l'(E, r)}{\psi_l(E, r)}

**误差度量**：

1. **曲线均方根差 (Curve RMS)**

   .. math::

      L_{\text{RMS}} = \sqrt{ \frac{1}{M} \sum_{j=1}^M \left[L_l^{\text{AE}}(E_j) - L_l^{\text{PS}}(E_j)\right]^2 }

2. **零点 RMS (Zero-Crossing RMS)**

   .. math::

      \Delta E_{\text{RMS}} = \sqrt{ \frac{1}{N} \sum_{i=1}^N \left(E_i^{\text{AE}} - E_i^{\text{PS}}\right)^2 }

**阈值设置**：

- **金属元素** (以 Al 为代表)：:math:`L_{\text{RMS}}^{\text{valence}} < 16.0`。金属价电子在
  软势中近似自由电子，曲线本身趋近零，因此 RMS 可能达到 8-10 的"平台"，可以接受。
- **共价元素** ：:math:`L_{\text{RMS}}^{\text{valence}} < 3.0`；节点明显、相位变化大，需更严格标准（教学实现的实践阈值）。
- 零点 RMS 主要用于确认赝势是否重现全电子的散射共振。若 ``ValidationReport`` 中显示 ``Infinity``，
  表示在当前能量窗口内零点数量不足（实现上要求 AE/PS 都至少出现 2 个零点才启用该判据），
  可通过缩小 :math:`E_{\text{step}}`、扩大能量范围或调整 :math:`r_{\text{test}}` 重跑。

**读数方法**：

- 通过 ``ValidationReport.log_deriv`` 属性获取对数导数数据，可用 matplotlib 绘制 AE/PS 曲线对比。
- Metal 与 Covalent 的阈值在 ``ValidationReport`` 中通过元素分类自动选择。

幽灵态检测
----------

**定义**：赝势哈密顿量产生的非物理束缚态，能量通常位于价态附近，扰乱平面波求解。

检测流程：

1. 在均匀网格上重建径向哈密顿量：

   .. math::

      H_l = -\tfrac{1}{2} \frac{d^2}{dr^2} + V_l^{\text{loc}}(r) + \frac{l(l+1)}{2r^2}

2. 对角化得到本征对：:math:`(E_i, \psi_i)`。
3. 计算尾部比例 :math:`\tau_i = |\psi_i(r_{\max})| / \max_r |\psi_i(r)|`，以区分盒态与真束缚态。
4. 若 :math:`E_i \in [-0.15, +0.05]` Ha 且 :math:`\tau_i < 0.1`，则判为幽灵态；
   :math:`\tau_i > 0.1` 的则是数值盒态，可忽略。

**能量感知分类**：

仅基于尾部比例 :math:`\tau` 无法区分危险幽灵态与安全的 Rydberg 激发态。改进判据如下：

- :math:`E > 0`：正能散射态（盒态），连续谱在有限盒子中的离散化产物。
- :math:`0 > E > \varepsilon_{\text{valence}} - \delta`：Rydberg 激发态，高主量子数束缚态序列
  (如 4s, 5s, 6s...)，对基态 DFT 无影响。
- :math:`E < \varepsilon_{\text{valence}} - \delta`：潜在危险幽灵态，需进一步用 :math:`\tau` 判据区分。

其中 :math:`\delta = 0.01` Ha 为能量容差。

**阈值**：教学任务采用较宽松的 ``ghost_total <= 10``，以允许 s 通道在过大 :math:`r_c` 下暴露问题；
一旦进入生产环境，应要求 ``ghost_total = 0``。

可转移性验证流程
----------------

综合三个指标的推荐顺序：

1. **范数守恒**：失败直接判定伪化无效；无需继续对数导数与幽灵态验证。
2. **对数导数**：在范数通过后，锁定影响材料性质的通道（例如 Al 的 p 通道），根据元素类型应用
   金属/共价阈值。
3. **幽灵态扫描**：最后检查赝势在平面波基中的稳定性。

若三个步骤全部通过，即可在 QE 中配置以下验证策略：

- 以 :math:`E_{\text{cut}} = 35` Ry 启动 ``vc-relax`` 计算，观察能量收敛趋势。
- 对比 ``ValidationReport.summary()`` 输出，确保阈值可以迁移到实战。
