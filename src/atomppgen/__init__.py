"""
AtomPPGen: 模守恒赝势生成器

基于 Troullier-Martins (TM) 方法生成 Kleinman-Bylander (KB) 形式的
模守恒（Norm-Conserving）赝势。

主要模块：
- ae_atom: 全电子原子解（调用 AtomSCF）
- tm: TM 伪化器
- invert: 半局域势反演
- kb: KB 可分离形式转换
- validate: 可转移性检验（范数/对数导数/幽灵态）（开发中）
- io: 数据导入导出（开发中）
"""

__version__ = "0.1.0"

# 导入主要模块
from atomppgen.ae_atom import solve_ae_atom, AEAtomResult
from atomppgen.tm import tm_pseudize, TMResult
from atomppgen.invert import invert_semilocal_potential, InvertResult
from atomppgen.kb import kb_transform, KBResult

__all__ = [
    "__version__",
    "solve_ae_atom",
    "AEAtomResult",
    "tm_pseudize",
    "TMResult",
    "invert_semilocal_potential",
    "InvertResult",
    "kb_transform",
    "KBResult",
]
