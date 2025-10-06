"""
AtomPPGen 命令行接口（CLI）

占位实现，完整功能将在后续开发。

基本用法
--------
ppgen generate --element Al --xc LDA --rc 2.0,2.0,2.5
ppgen test --input al_lda.json
ppgen plot --input al_lda.json --out figs/
"""

import sys
import argparse


def main():
    """
    AtomPPGen CLI 主入口

    当前为占位实现，仅提供基本框架。
    完整功能（generate, test, plot）将在后续里程碑实现。
    """
    parser = argparse.ArgumentParser(
        prog="ppgen",
        description="Norm-Conserving Pseudopotential Generator (Troullier-Martins method)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法：
  ppgen generate --element Al --xc LDA --rc 2.0,2.0,2.5
  ppgen test --input al_lda.json
  ppgen plot --input al_lda.json --out figs/

当前版本为占位实现，完整功能开发中。
        """
    )

    subparsers = parser.add_subparsers(dest="command", help="可用子命令")

    # generate 子命令占位
    gen_parser = subparsers.add_parser(
        "generate",
        help="生成赝势（占位）"
    )
    gen_parser.add_argument("--element", help="元素符号，如 Al")
    gen_parser.add_argument("--xc", help="交换关联泛函，如 LDA")
    gen_parser.add_argument("--rc", help="截断半径（逗号分隔），如 2.0,2.0,2.5")

    # test 子命令占位
    test_parser = subparsers.add_parser(
        "test",
        help="可转移性检验（占位）"
    )
    test_parser.add_argument("--input", help="输入赝势文件")

    # plot 子命令占位
    plot_parser = subparsers.add_parser(
        "plot",
        help="绘图（占位）"
    )
    plot_parser.add_argument("--input", help="输入赝势文件")
    plot_parser.add_argument("--out", help="输出目录")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 0

    # 占位实现：仅打印收到的命令
    print(f"[占位] ppgen {args.command} 命令已接收")
    print(f"[占位] 参数: {vars(args)}")
    print()
    print("注意：当前为占位实现，完整功能将在后续里程碑中实现。")
    print("      预计在 M6（导出与示例）阶段完成 CLI 全部功能。")

    return 0


if __name__ == "__main__":
    sys.exit(main())
