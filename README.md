# AtomPPGen

模守恒（Norm-Conserving）赝势生成器，使用 Troullier-Martins 方法。

**在线文档**: https://bud-primordium.github.io/AtomPPGen/

---

## 项目简介

**AtomPPGen** 是一个教学用赝势生成工具，实现以下功能：

- **TM 伪化**：Troullier-Martins 方法生成赝波函数
- **KB 形式**：Kleinman-Bylander 可分离非局域投影子
- **可转移性检验**：范数守恒、对数导数曲线、幽灵态检测
- **输出格式**：JSON/NPZ（后续支持 UPF）

**设计目标**：生成 Al（Z=13）的 LDA 模守恒赝势。

### 开发状态

**核心功能**：
- ✅ TM 伪化算法
- ✅ KB 转换与局域道选择
- ✅ 可转移性验证工具（范数守恒、对数导数、幽灵态检测）
- ✅ JSON/NPZ 导出
- 🚧 CLI 命令行工具（Python API 可用）

**待完善项目**：
- 验证工具的参数优化与可靠性增强
- 多通道验证测试的完整覆盖
- 统计报告与可视化输出

**使用建议**：当前版本适合通过 Python API 使用（参见下方示例），CLI 工具正在开发中。

---

## 安装

### 方式一：使用 uv（推荐）

```bash
cd AtomPPGen
uv venv
source .venv/bin/activate
uv pip install -e ../AtomSCF  # 安装依赖
uv pip install -e ".[dev]"    # 安装本项目
```

### 方式二：使用标准 venv

```bash
cd AtomPPGen
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
# .venv\Scripts\activate   # Windows

pip install -e ../AtomSCF  # 安装 AtomSCF 依赖
pip install -e ".[dev]"    # 安装 AtomPPGen 开发依赖
```

---

## 快速开始

### 生成 Al 赝势（CLI 开发中）

```bash
# 生成 Al 的 LDA 赝势
ppgen generate \
  --element Al \
  --xc LDA \
  --channels s,p,d \
  --rc 2.1,2.2,2.4 \
  --loc d \
  --out out/al_lda_tm.json

# 可转移性检验
ppgen test \
  --input out/al_lda_tm.json \
  --rderiv 3.0 \
  --e-range -0.5:0.5:0.05 \
  --ghost-check

# 绘图
ppgen plot \
  --input out/al_lda_tm.json \
  --out figs/al_lda/
```

### Python API 示例

```python
from atomppgen.ae_atom import solve_ae_atom

# 获取 Al 原子的全电子解（使用变量变换方法）
result = solve_ae_atom(Z=13, xc="PZ81", grid_type="exp_transformed")

print(f"3s 能级: {result.eps_by_l[0][2]:.6f} Ha")
print(f"3p 能级: {result.eps_by_l[1][2]:.6f} Ha")
print(f"总能量: {result.energies['E_total']:.6f} Ha")
```

---

## 开发指南

### 运行测试

```bash
# 运行所有测试
pytest

# 仅单元测试
pytest -m unit

# 详细输出
pytest -v

# 特定模块测试
pytest tests/test_ae_atom.py -v
```

### 项目结构

```
AtomPPGen/
├── src/atomppgen/          # 核心源码
│   ├── ae_atom.py          # 全电子原子求解器（调用 AtomSCF）
│   ├── tm.py               # TM 伪化器
│   ├── invert.py           # 半局域势反演
│   ├── kb.py               # KB 转换
│   ├── validate.py         # 可转移性检验
│   ├── io.py               # 数据导出
│   └── cli.py              # 命令行工具
├── tests/                  # 单元测试
├── examples/               # 使用示例
└── pyproject.toml          # 项目配置
```

---

## 技术细节

### 全电子原子解（AE）

**推荐配置**：使用变量变换方法（精度提升约 7 倍）

```python
from atomppgen.ae_atom import solve_ae_atom

# 变量变换方法（默认）
result = solve_ae_atom(
    Z=13,
    xc="PZ81",
    lmax=2,
    grid_type="exp_transformed",  # 指数网格 + 变量变换
    grid_params={"n": 1200, "rmax": 120.0, "total_span": 6.5},
    scf_params={"tol": 1e-7, "maxiter": 200, "mix_alpha": 0.25}
)
```

### 与 NIST 数据对比

当前实现为**非相对论 LSDA**，与 NIST LSD 参考数据存在系统性差异：

| 项目 | AtomSCF | NIST LSD | 差异 |
|------|---------|----------|------|
| 总能量 (Ha) | -237.30 | -241.32 | ~4 Ha |
| 1s (Ha) | -54.27 | -55.15 | ~0.88 Ha |
| 3s (Ha) | -0.25 | -0.30 | ~0.05 Ha |

**说明**：
- 差异来源可能包括数值方法、泛函实现细节等
- 价层轨道（3s, 3p）相对差异较小（~0.03-0.05 Ha）
- 对于赝势生成，价层相对精度通常已足够

---

## 依赖

- **AtomSCF** (>= 0.1.0): 全电子原子求解器
- **NumPy** (>= 1.20)
- **SciPy** (>= 1.7)

---

## 开发规范

- **语言**：代码注释与文档使用**中文**
- **Docstring**：NumPy 风格
- **公式**：使用 LaTeX（`:math:` 或 `.. math::`）
- **测试**：pytest（`tests/` 目录，`-m unit` 标记）
- **Git 提交**：Conventional Commits 格式（英文，无 AI 痕迹）

---

## 参考资料

- **TM 方法**: Troullier & Martins, *PRB* 43, 1993 (1991)
- **KB 形式**: Kleinman & Bylander, *PRL* 48, 1425 (1982)
- **QE 文档**: Giannozzi, *Notes on pseudopotential generation* (2019)
- **AtomSCF**: 变量变换方法详见 `../AtomSCF/docs/source/algorithm/numerical_methods.rst`

---

## 许可

MIT License - 教学课程作业，供学习参考。
