项目介绍
========

背景
----

赝势（Pseudopotential）是固体物理和量子化学计算中的重要工具，用于简化第一性原理计算中的电子-离子相互作用。

模守恒赝势的核心要求
--------------------

1. **范数守恒**

   在截断半径 :math:`r_c` 内，伪轨道与全电子轨道的范数相等：

   .. math::

      \int_0^{r_c} |\psi_{\text{PS}}(r)|^2 dr = \int_0^{r_c} |\psi_{\text{AE}}(r)|^2 dr

2. **能量一致性**

   伪原子的价轨道能级应与全电子原子一致。

3. **可转移性**

   在不同化学环境下，赝势应能给出一致的结果。

技术路线
--------

本项目采用以下技术路线：

1. **全电子原子解**

   - 使用 AtomSCF 求解 LDA 自洽场方程
   - 采用变量变换方法提高精度
   - 获取各角动量通道的波函数和能级

2. **Troullier-Martins 伪化**

   - 在 :math:`r \leq r_c` 内，构造形式：

     .. math::

        v_l(r) = r^{l+1} \exp\left(\sum_{i=0}^N a_{2i} r^{2i}\right)

   - 通过非线性方程组求解系数 :math:`a_{2i}`
   - 约束条件：

     - 在 :math:`r=r_c` 处函数值及导数连续
     - 内区范数守恒

3. **半局域势反演**

   - 从径向 Schrödinger 方程反演：

     .. math::

        V_l(r) = \varepsilon_l + \frac{1}{2}\frac{v''_l}{v_l} - \frac{l(l+1)}{2r^2}

   - 数值保护：节点附近平滑插值

4. **Kleinman-Bylander 可分离形式**

   - 选择局域通道 :math:`V_{\text{loc}}`（推荐 d）
   - 构造投影子：

     .. math::

        \beta_l(r) \propto [V_l(r) - V_{\text{loc}}(r)] \psi_l(r)

   - 计算耦合系数 :math:`D_l`

5. **可转移性检验**

   - 范数守恒检查
   - 对数导数曲线对比
   - 幽灵态检测

6. **数据导出**

   - JSON/NPZ 格式
   - 包含元数据、网格、势、投影子

参考文献
--------

- **TM 方法**: N. Troullier and J. L. Martins, *PRB* **43**, 1993 (1991)
- **KB 形式**: L. Kleinman and D. M. Bylander, *PRL* **48**, 1425 (1982)
- **QE 文档**: P. Giannozzi, *Notes on pseudopotential generation* (2019)
