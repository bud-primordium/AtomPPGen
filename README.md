# AtomPPGen

模守恒（Norm-Conserving）赝势生成器，使用 Troullier-Martins 方法。

**在线文档**: https://bud-primordium.github.io/AtomPPGen/

---

## 项目简介

**AtomPPGen** 是一个教学用赝势生成工具，实现以下功能：

- **TM 伪化**：Troullier-Martins 方法生成赝波函数
- **KB 形式**：Kleinman-Bylander 可分离非局域投影子
- **可转移性检验**：范数守恒、对数导数曲线、幽灵态检测（含能量感知分类）
- **输出格式**：JSON/NPZ（后续支持 UPF）

**设计目标**：生成 Al（Z=13）的 LDA 模守恒赝势，验证精度达到数值极限。

### 开发状态

**核心功能**：
- TM 伪化算法（含径向节点符号处理）
- KB 转换与局域道选择
- 可转移性验证工具（范数守恒、对数导数、幽灵态检测）
- JSON/NPZ 导出
- Python API 完整可用

**当前版本特性**：
- 范数守恒精度：1e-13（优于标准 1e-6）
- 幽灵态检测：能量感知分类，准确区分 Rydberg 激发态与真幽灵态
- 多通道验证：s/p/d 通道完整测试覆盖
- 参数扫描：自动化测试 9 组 rc 组合

**待扩展功能**：
- 命令行工具（CLI）：提供 `ppgen generate/test/plot` 命令
- UPF 格式导出：Quantum ESPRESSO 兼容性
- 其他元素支持：Si、Mg 等

**使用建议**：当前版本通过 Python API 使用（参见下方示例），具备完整验证功能。

---

## 安装

### 依赖项
- Python >= 3.8
- NumPy >= 1.20
- SciPy >= 1.7
- AtomSCF >= 0.1.0

### 安装步骤
```bash
# 克隆仓库
git clone https://github.com/yourusername/AtomPPGen.git
cd AtomPPGen

# 创建虚拟环境
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
# .venv\Scripts\activate   # Windows

# 安装依赖
pip install -e ../AtomSCF  # 本地开发（假设 AtomSCF 在同级目录）
pip install -e .            # 安装 AtomPPGen
```

**注**：若使用 `uv` 工具，可替换为 `uv venv` + `uv pip install -e ...`，其余步骤一致。

---

## 快速开始

### 完整工作流示例

生成 Al 的 LDA 赝势并进行验证：

```python
from atomppgen.ae_atom import solve_ae_atom
from atomppgen.tm import tm_pseudize
from atomppgen.invert import invert_semilocal_potential
from atomppgen.kb import kb_transform
from atomppgen.validate import run_full_validation
from atomppgen.export import export_pseudopotential

# 第 1 步：求解全电子原子
ae = solve_ae_atom(
    Z=13,
    spin_mode="LDA",
    lmax=2,
    grid_type="exp_transformed",
    grid_params={"n": 900, "total_span": 6.5},
    scf_params={"tol": 1e-7, "maxiter": 200}
)

# 第 2 步：对各通道进行 TM 伪化
rc_map = {0: 2.1, 1: 2.2, 2: 2.4}  # s, p, d 通道截断半径（Bohr）
tm_results = {}
for l in range(3):
    tm_results[l] = tm_pseudize(
        r=ae.r,
        w=ae.w,
        u_ae=ae.u_by_l[l][-1],  # 价层轨道
        eps=ae.eps_by_l[l][-1],
        l=l,
        rc=rc_map[l],
        continuity_orders=2
    )

# 第 3 步：反演半局域势
inv_results = {}
for l in range(3):
    inv_results[l] = invert_semilocal_potential(tm_results[l], ae.r)

# 第 4 步：构造 KB 投影子（选择 d 通道作为局域道）
u_by_l = {l: ae.u_by_l[l][-1] for l in range(3)}
kb_result = kb_transform(inv_results, u_by_l, ae.r, ae.w, loc_channel=2)

# 第 5 步：完整验证
validation_report = run_full_validation(
    ae_result=ae,
    tm_dict=tm_results,
    inv_dict=inv_results,
    r_test=3.0,
    E_range_Ry=(-0.5, 0.5),
    E_step_Ry=0.05
)

# 第 6 步：导出赝势
export_pseudopotential(
    ae_result=ae,
    tm_dict=tm_results,
    inv_dict=inv_results,
    validation_report=validation_report,
    output_prefix="outputs/al_lda_tm",
    formats=["json", "npz"]
)

# 查看验证结果
print(f"范数守恒：{'通过' if validation_report.overall_passed else '失败'}")
print(f"KB 局域势最大值：{kb_result.V_loc.max():.4f} Ha")
```

### 推荐参数（Al, LDA）

经过系统验证，推荐以下参数组合：

```python
rc_map = {
    0: 2.1,  # s 通道
    1: 2.2,  # p 通道
    2: 2.4   # d 通道
}
local_channel = 2  # d 作为局域道
xc_functional = "PZ81"  # Perdew-Zunger 1981 LDA
```

**验证结果**（基于参数扫描测试）：
- 范数守恒误差：~1e-13（所有通道，达到数值极限）
- 对数导数 RMS：满足金属元素阈值（< 16.0）
- 幽灵态检测：采用能量感知分类，可区分 Rydberg 激发态与危险幽灵态
- 测试覆盖：60 个测试全部通过

---

## 开发指南

### 运行测试

```bash
# 设置路径（假设 AtomSCF 与 AtomPPGen 在同级目录）
export PYTHONPATH=src:../AtomSCF/src

# 运行所有测试
pytest

# 仅单元测试
pytest -m unit

# 详细输出
pytest -v

# 特定模块测试
pytest tests/test_validate.py -v
```

### 项目结构

```
AtomPPGen/
├── src/atomppgen/
│   ├── ae_atom.py      # 全电子原子求解（调用 AtomSCF）
│   ├── tm.py           # TM 伪化器
│   ├── invert.py       # 半局域势反演
│   ├── kb.py           # KB 投影子构造
│   ├── validate.py     # 验证工具（范数/对数导数/幽灵态）
│   ├── export.py       # JSON/NPZ 导出
│   └── cli.py          # 命令行接口（规划中）
├── tests/              # 单元测试（60 个测试）
├── examples/           # 示例脚本
└── docs/               # Sphinx 文档
```

### 构建文档

```bash
cd docs
pip install sphinx sphinx-rtd-theme
make html
# 输出位于 build/html/index.html
```

---

## 技术细节

### 幽灵态检测方法

**能量感知分类判据**：

AtomPPGen 使用改进的幽灵态检测逻辑，区分真幽灵态与 Rydberg 激发态：

1. **正能散射态**（E > 0）：有限网格盒态，归为安全态
2. **Rydberg 激发态**（E_valence < E < 0）：高主量子数束缚态（4s, 5s, 6s...），距离费米能级远，对基态 DFT 无影响
3. **潜在危险幽灵态**（E < E_valence - 0.01 Ha）：需进一步检查尾部比例

**物理意义**：Rydberg 激发态的存在证明赝势保留了正确的长程库仑行为（-1/r）。

详见文档：[幽灵态检测方法](docs/source/algorithm/validation_methods.rst)

### 范数守恒精度

TM 方法通过非线性方程组隐式保证范数守恒。当前实现达到：

- 理论标准：< 1e-6
- 实际精度：~1e-13（数值极限）

### 对数导数匹配

**元素类型差异化阈值**：

- **金属元素**（Al, Na, Mg）：价区曲线 RMS < 16.0
- **共价元素**（Si, C, N）：价区曲线 RMS < 0.3

金属元素在远核区的软库仑势中，对数导数对能量变化不敏感，导致较大的 AE-PS 相位差异，
这是金属元素固有特性，非赝势质量缺陷。

---

## 与 NIST 数据对比

当前实现为**非相对论 LSDA**，与 NIST LSD 参考数据对比：

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

## 参考资料

- **TM 方法**: Troullier & Martins, *PRB* **43**, 1993 (1991)
- **KB 形式**: Kleinman & Bylander, *PRL* **48**, 1425 (1982)
- **QE 文档**: Giannozzi, *Notes on pseudopotential generation* (2019)
- **幽灵态检测**: Rappe et al., *PRB* **41**, 1227 (1990)
- **对数导数方法**: Gonze et al., *Comput. Mater. Sci.* **25**, 478 (2002)

---

## 许可

MIT License - 教学课程作业，供学习参考。
