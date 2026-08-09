#!/usr/bin/env python3

import argparse
import os
import shutil
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import run_pipeline

METHOD_TOOL_LABELS = {
    "treepl": "TreePL executable",
    "mdcat": "MD-Cat executable",
    "wlogdate": "wLogDate executable",
    "lsd2": "LSD2 executable",
}

METHOD_DEFAULT_PATHS = {
    "treepl": "treePL",
    "mdcat": "md_cat.py",
    "wlogdate": "launch_wLogDate.py",
    "lsd2": "lsd2",
}

PYTHON_MODULES = ["dendropy"]


def module_available(module_name):
    try:
        __import__(module_name)
        return True
    except Exception:
        return False


def resolve_tool(path, label):
    return run_pipeline.resolve_executable(path, label, required=False)


def main():
    parser = argparse.ArgumentParser(
        description="Validate external dependencies for the coalescent dating pipeline",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=["treepl", "mdcat", "wlogdate", "lsd2", "all"],
        default=["all"],
        help="Methods to validate (default: all)",
    )
    parser.add_argument(
        "--astral4-bin",
        default="bin/astral4",
        help="ASTRAL executable path/name to check (default: bin/astral4)",
    )
    parser.add_argument(
        "--treepl-bin",
        default=METHOD_DEFAULT_PATHS["treepl"],
        help="TreePL executable path/name to check",
    )
    parser.add_argument(
        "--mdcat-bin",
        default=METHOD_DEFAULT_PATHS["mdcat"],
        help="MD-Cat executable/script path/name to check",
    )
    parser.add_argument(
        "--wlogdate-bin",
        default=METHOD_DEFAULT_PATHS["wlogdate"],
        help="wLogDate executable/script path/name to check",
    )
    parser.add_argument(
        "--lsd2-bin",
        default=METHOD_DEFAULT_PATHS["lsd2"],
        help="LSD2 executable path/name to check",
    )

    args = parser.parse_args()

    requested_methods = ["treepl", "mdcat", "wlogdate", "lsd2"]
    if "all" not in args.methods:
        requested_methods = args.methods

    print("Checking pipeline dependencies...\n")

    failures = []

    astral_path = resolve_tool(args.astral4_bin, "ASTRAL/CASTLES-Pro executable")
    if astral_path:
        print(f"[OK] ASTRAL/CASTLES-Pro executable: {astral_path}")
    else:
        failures.append(
            "ASTRAL/CASTLES-Pro executable (expected one of: "
            f"{args.astral4_bin}, astral4, astral)"
        )
        print("[MISSING] ASTRAL/CASTLES-Pro executable")

    tool_paths = {
        "treepl": args.treepl_bin,
        "mdcat": args.mdcat_bin,
        "wlogdate": args.wlogdate_bin,
        "lsd2": args.lsd2_bin,
    }

    for method in requested_methods:
        label = METHOD_TOOL_LABELS[method]
        tool = resolve_tool(tool_paths[method], label)
        if tool:
            print(f"[OK] {label}: {tool}")
        else:
            failures.append(
                f"{label} (expected one of: {tool_paths[method]}, "
                + ", ".join(run_pipeline.TOOL_ALIASES.get(label, ()))
                + ")"
            )
            print(f"[MISSING] {label}")

    for module_name in PYTHON_MODULES:
        if module_available(module_name):
            print(f"[OK] Python module: {module_name}")
        else:
            failures.append(f"Python module '{module_name}'")
            print(f"[MISSING] Python module: {module_name}")

    if failures:
        print("\nDependency check failed. Missing items:")
        for item in failures:
            print(f"  - {item}")

        bootstrap_script = os.path.join(REPO_ROOT, "scripts", "bootstrap_dependencies.sh")
        print("\nSuggested fix:")
        print("  1) Activate your conda environment")
        print(f"  2) Run: bash {bootstrap_script}")
        print("  3) Add tool path to PATH (if not already): export PATH=\"$PWD/.local/bin:$PATH\"")
        return 1

    print("\nAll requested dependencies are available.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
