#!/usr/bin/env python3

import argparse
import os
import shlex
import shutil
import subprocess
import sys


def run(cmd):
    print(f"\n[RUN] {cmd}\n")
    result = subprocess.run(cmd, shell=True)

    if result.returncode != 0:
        sys.exit(f"[ERROR] Command failed:\n{cmd}")


def check_exists(path, label):
    if not os.path.exists(path):
        sys.exit(f"[ERROR] {label} not found: {path}")


def check_executable(path, label):
    if os.path.exists(path):
        return

    if os.sep not in path and shutil.which(path):
        return

    sys.exit(f"[ERROR] {label} not found: {path}")


def mkdir(path):
    os.makedirs(path, exist_ok=True)


def shell_join(args):
    return " ".join(shlex.quote(str(arg)) for arg in args)


def parse_ci(ci_args):
    if ci_args is None:
        return None

    if len(ci_args) == 1:
        ci_args = ci_args[0].split()

    if len(ci_args) != 3:
        sys.exit(
            "[ERROR] --CI expects three values: samples lower_quantile "
            "upper_quantile"
        )

    samples, lower, upper = ci_args

    try:
        samples_int = int(samples)
        lower_float = float(lower)
        upper_float = float(upper)
    except ValueError:
        sys.exit("[ERROR] --CI values must be numeric")

    if samples_int <= 0:
        sys.exit("[ERROR] --CI samples must be positive")

    if not 0 <= lower_float < upper_float <= 1:
        sys.exit("[ERROR] --CI quantiles must satisfy 0 <= lower < upper <= 1")

    return str(samples_int), str(lower_float), str(upper_float)


def run_astral4(
    astral4_bin,
    gene_trees,
    output_tree,
    species_tree=None,
    outgroup=None,
    gene_length=None
):
    cmd = [
        astral4_bin,
        "-i",
        gene_trees,
        "-o",
        output_tree,
    ]

    if species_tree:
        cmd.extend(["-C", "-c", species_tree])

    if outgroup:
        cmd.extend(["--root", outgroup])

    if gene_length is not None:
        cmd.extend(["--genelength", gene_length])

    run(shell_join(cmd))


def run_treepl(tree, calibrations, output_dir):
    config = os.path.join(output_dir, "treepl.config")
    dated_tree = os.path.join(output_dir, "dated_tree.tre")

    with open(config, "w") as f:
        f.write(f"treefile = {tree}\n")
        f.write("smooth = 100\n")
        f.write("numsites = 500000\n")
        f.write(f"outfile = {dated_tree}\n\n")

        with open(calibrations) as c:
            f.write(c.read())

    run(f"treePL {config}")

    return dated_tree


def run_mdcat(
    tree,
    calibrations,
    output_dir,
    ci=None,
    seq_length=None,
    p=10
):
    dated_tree = os.path.join(output_dir, "dated_tree.tre")

    cmd = [
        "python3",
        "md_cat.py",
        "-i",
        tree,
        "-o",
        dated_tree,
        "-p",
        p,
        "-t",
        calibrations,
        "-b",
    ]

    if seq_length is not None:
        cmd.extend(["-l", seq_length])

    if ci is not None:
        cmd.extend(["--CI", " ".join(ci)])

    run(shell_join(cmd))

    return dated_tree


def run_wlogdate(tree, calibrations, output_dir):
    dated_tree = os.path.join(output_dir, "dated_tree.tre")

    cmd = (
        f"python launch_wLogDate.py "
        f"-i {tree} "
        f"-c {calibrations} "
        f"-o {dated_tree}"
    )

    run(cmd)

    return dated_tree


def run_lsd2(tree, calibrations, output_dir):
    prefix = os.path.join(output_dir, "lsd2")

    cmd = (
        f"lsd2 "
        f"-i {tree} "
        f"-d {calibrations} "
        f"-o {prefix}"
    )

    run(cmd)

    dated_tree = prefix + ".date.nwk"

    if os.path.exists(dated_tree):
        shutil.copy(dated_tree, os.path.join(output_dir, "dated_tree.tre"))

    return os.path.join(output_dir, "dated_tree.tre")


def main():
    parser = argparse.ArgumentParser(
        description="Coalescent-aware dating pipeline"
    )

    parser.add_argument(
        "--gene-trees",
        required=True,
        help="Gene trees in Newick format"
    )

    parser.add_argument(
        "--species-tree",
        default=None,
        help="Optional user-provided species tree"
    )

    parser.add_argument(
        "--calibrations",
        required=True,
        help="Calibration file"
    )

    parser.add_argument(
        "--method",
        required=True,
        choices=["treepl", "mdcat", "wlogdate", "lsd2"],
        help="Dating method"
    )

    parser.add_argument(
        "--outgroup",
        default=None,
        help="Optional outgroup taxon passed to ASTRAL/CASTLES-Pro with --root"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "--astral4-bin",
        default="bin/astral4",
        help="Path to the ASTRAL/CASTLES-Pro executable (default: bin/astral4)"
    )

    parser.add_argument(
        "--gene-length",
        type=int,
        default=None,
        help="Optional gene length passed to ASTRAL/CASTLES-Pro with --genelength"
    )

    parser.add_argument(
        "--CI",
        nargs="+",
        default=None,
        metavar="CI",
        help=(
            "Compute MD-Cat confidence intervals, e.g. "
            '--CI "1000 0.025 0.975"'
        )
    )

    parser.add_argument(
        "--seq-length",
        type=int,
        default=None,
        help="Optional sequence length passed to MD-Cat with -l"
    )

    parser.add_argument(
        "--mdcat-p",
        type=int,
        default=10,
        help="MD-Cat -p value (default: 10)"
    )

    args = parser.parse_args()
    ci = parse_ci(args.CI)

    if ci is not None and args.method != "mdcat":
        sys.exit(
            "[ERROR] Confidence intervals are currently supported only with "
            "--method mdcat"
        )

    if args.seq_length is not None and args.method != "mdcat":
        sys.exit("[ERROR] --seq-length is currently supported only with MD-Cat")

    if args.seq_length is not None and args.seq_length <= 0:
        sys.exit("[ERROR] --seq-length must be positive")

    if args.gene_length is not None and args.gene_length <= 0:
        sys.exit("[ERROR] --gene-length must be positive")

    if args.mdcat_p <= 0:
        sys.exit("[ERROR] --mdcat-p must be positive")

    check_exists(args.gene_trees, "Gene trees")
    check_exists(args.calibrations, "Calibration file")
    check_executable(args.astral4_bin, "ASTRAL/CASTLES-Pro executable")

    if args.species_tree:
        check_exists(args.species_tree, "Species tree")

    mkdir(args.output)

    intermediate_dir = os.path.join(args.output, "intermediate")
    logs_dir = os.path.join(args.output, "logs")

    mkdir(intermediate_dir)
    mkdir(logs_dir)

    su_tree = os.path.join(
        args.output,
        "species_tree_su.tre"
    )

    print("\n=== STEP 1: ASTRAL/CASTLES-Pro SU tree estimation ===\n")

    run_astral4(
        args.astral4_bin,
        args.gene_trees,
        su_tree,
        species_tree=args.species_tree,
        outgroup=args.outgroup,
        gene_length=args.gene_length
    )

    dating_input = su_tree

    print("\n=== STEP 2: Molecular dating ===\n")

    if args.method == "treepl":
        run_treepl(
            dating_input,
            args.calibrations,
            args.output
        )

    elif args.method == "mdcat":
        run_mdcat(
            dating_input,
            args.calibrations,
            args.output,
            ci=ci,
            seq_length=args.seq_length,
            p=args.mdcat_p
        )

    elif args.method == "wlogdate":
        run_wlogdate(
            dating_input,
            args.calibrations,
            args.output
        )

    elif args.method == "lsd2":
        run_lsd2(
            dating_input,
            args.calibrations,
            args.output
        )

    print("\nPipeline completed successfully.\n")
    print(f"SU tree: {su_tree}")
    print(f"Dated tree: {os.path.join(args.output, 'dated_tree.tre')}")

    if ci is not None:
        print(
            "Confidence intervals: "
            f"{ci[0]} MD-Cat samples, quantiles {ci[1]} and {ci[2]}"
        )


if __name__ == "__main__":
    main()
