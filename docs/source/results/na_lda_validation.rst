Na (LDA) 验证结果
=================

本页给出 Na（简单金属）的教学示例验证结果，用于对比 Si（共价）案例，并展示
“金属/共价”阈值策略的差异。

实验环境
--------

- 生成脚本：`AtomPPGen/examples/validate_na_lda.py`
- XC：LDA-PZ81
- 通道：s, p, d（d 作为 KB 局域通道候选）

参数与结论
----------

本示例使用的截断半径（Bohr）：

- :math:`(r_c^s, r_c^p, r_c^d) = (2.6, 2.8, 3.0)`
- 对数导数测试半径：:math:`r_{\text{test}} = 3.5` Bohr（脚本取 :math:`\max(r_c)+0.5`）

验证摘要（典型输出）
-------------------

对数导数价区 RMS 阈值：金属元素使用 :math:`L_{\mathrm{RMS}}^{\mathrm{valence}} < 16.0`。

.. list-table::
   :header-rows: 1

   * - 通道
     - :math:`L_{\mathrm{RMS}}^{\mathrm{valence}}`
     - 零点 RMS
     - 结论
   * - s (l=0)
     - 2.16
     - N/A（零点不足，不纳入判据）
     - PASS
   * - p (l=1)
     - 0.12
     - N/A（零点不足，不纳入判据）
     - PASS
   * - d (l=2)
     - 0.09
     - N/A（零点不足，不纳入判据）
     - PASS

幽灵态检测：真幽灵态 0 个；存在若干盒态（有限边界离散化产物），不影响通过判定。

输出产物
--------

脚本会导出：

- `AtomPPGen/outputs/na_lda/na_lda.json`：元数据 + 验证报告摘要
- `AtomPPGen/outputs/na_lda/na_lda.npz`：数值数组（含 KB 数据）

