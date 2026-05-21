import argparse
from pathlib import Path
from typing import Optional


def repo_root() -> Path:
    """Return repository root as the directory containing this script."""
    return Path(__file__).resolve().parent


def resolve_target_path(cli_path: Optional[str]) -> Path:
    """Resolve target file path for verification.

    Priority:
    1) --file argument
    2) default repo-relative modules/chip_analysis.R
    """
    if cli_path:
        p = Path(cli_path).expanduser()
        return p if p.is_absolute() else (repo_root() / p)

    return repo_root() / "modules" / "chip_analysis.R"


def _is_r_source_file(path: Path) -> bool:
    return path.suffix.lower() == ".r"


def run_check(target: Path) -> int:
    if not target.exists():
        print("=== 检查芯片分析模块代码 ===\n")
        print(f"❌ 目标文件不存在: {target}")
        print("提示: 可用 --file 指定目标文件，例如:\n"
              "  python verify_code.py --file modules/chip_analysis.R")
        return 1

    if not _is_r_source_file(target):
        print("=== 检查芯片分析模块代码 ===\n")
        print(f"❌ 不支持的文件类型: {target}")
        print("提示: 该检查器仅用于 .R 源文件。")
        return 3

    content = target.read_text(encoding="utf-8")
    lines = content.split("\n")

    print("=== 检查芯片分析模块代码 ===\n")
    print(f"目标文件: {target}\n")

    print("检查关键代码:")

    ui_pattern = 'uiOutput("chip_soft_column_selection_panel")'
    server_pattern = "output$chip_soft_column_selection_panel <- renderUI"
    select_id_pattern = 'selectInput("chip_soft_id_col"'
    select_gene_pattern = 'selectInput("chip_soft_gene_col"'

    ui_found = any(ui_pattern in line for line in lines)
    print("✅ UI部分: " + ("找到" if ui_found else "未找到") + " uiOutput")

    server_found = any(server_pattern in line for line in lines)
    print("✅ Server部分: " + ("找到" if server_found else "未找到") + " renderUI定义")

    select_id = any(select_id_pattern in line for line in lines)
    select_gene = any(select_gene_pattern in line for line in lines)
    print("✅ selectInput: " + ("找到" if (select_id and select_gene) else "未找到") + " 直接生成的selectInput")

    print("\n关键代码位置:")
    for i, line in enumerate(lines, 1):
        if ui_pattern in line:
            print(f"  第{i}行 (UI): {line.strip()}")
        if server_pattern in line:
            print(f"  第{i}行 (Server): {line.strip()}")
        if select_id_pattern in line or select_gene_pattern in line:
            print(f"  第{i}行 (selectInput): {line.strip()}")

    ok = ui_found and server_found and select_id and select_gene
    result = "✅ 全部通过" if ok else "❌ 检查失败"
    print(f"\n=== 结果: {result} ===")

    if ok:
        print("\n请完全重启应用后测试！")

    return 0 if ok else 2


def main() -> int:
    parser = argparse.ArgumentParser(description="检查 chip_analysis.R 中关键 UI/Server 片段是否存在")
    parser.add_argument("--file", help="目标 R 文件路径（默认 modules/chip_analysis.R）", default=None)
    args = parser.parse_args()

    target = resolve_target_path(args.file)
    return run_check(target)


if __name__ == "__main__":
    raise SystemExit(main())
