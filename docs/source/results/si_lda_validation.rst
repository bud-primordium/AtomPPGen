Si (LDA) 验证结果
=================

本页给出 Si（共价元素）的教学示例验证结果，用于对比 Al/Na（金属）案例。

实验环境
--------

- 生成脚本：`AtomPPGen/examples/validate_si_lda.py`
- XC：LDA-PZ81
- 通道：s, p, d（d 作为 KB 局域通道候选）

参数与结论
----------

本示例使用的截断半径（Bohr）：

- :math:`(r_c^s, r_c^p, r_c^d) = (1.8, 2.0, 2.2)`
- 对数导数测试半径：:math:`r_{\text{test}} = 6.0` Bohr

验证摘要（典型输出）
-------------------

对数导数价区 RMS 阈值：共价元素使用 :math:`L_{\mathrm{RMS}}^{\mathrm{valence}} < 3.0`（见验证标准章节）。

.. list-table::
   :header-rows: 1

   * - 通道
     - :math:`L_{\mathrm{RMS}}^{\mathrm{valence}}`
     - 零点 RMS
     - 结论
   * - s (l=0)
     - 2.31
     - N/A（零点不足，不纳入判据）
     - PASS
   * - p (l=1)
     - 0.49
     - N/A（零点不足，不纳入判据）
     - PASS
   * - d (l=2)
     - 0.31
     - 0.018 Ha
     - PASS

幽灵态检测：真幽灵态 0 个；存在若干盒态（有限边界离散化产物），不影响通过判定。

输出产物
--------

脚本会导出：

- `AtomPPGen/outputs/si_lda/si_lda.json`：元数据 + 验证报告摘要
- `AtomPPGen/outputs/si_lda/si_lda.npz`：数值数组（含 KB 数据）

