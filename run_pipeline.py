#!/usr/bin/env python3

import argparse
import os
import re
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


def is_numeric_value(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def is_calendar_date(value):
    return re.match(r"^\d{1,4}(?:-\d{1,2}){0,2}$", value) is not None


def normalize_calibration_value(value, line_no):
    cleaned = value.strip().strip("`")

    if is_numeric_value(cleaned) or is_calendar_date(cleaned):
        return cleaned

    sys.exit(
        f"[ERROR] Calibration date on line {line_no} must be numeric or "
        f"year-month-day format: {value}"
    )


def parse_constraint_value(value, line_no):
    stripped = value.strip()
    wrapper_match = re.match(r"^(l|u|b)\s*\((.*)\)$", stripped, re.IGNORECASE)

    if not wrapper_match:
        return {
            "exact": normalize_calibration_value(stripped, line_no),
            "min": None,
            "max": None,
        }

    wrapper = wrapper_match.group(1).lower()
    contents = wrapper_match.group(2)

    if wrapper in {"l", "u"}:
        date = normalize_calibration_value(contents, line_no)
        return {
            "exact": None,
            "min": date if wrapper == "l" else None,
            "max": date if wrapper == "u" else None,
        }

    parts = [part.strip() for part in contents.split(",")]

    if len(parts) != 2:
        sys.exit(
            f"[ERROR] Bounded calibration on line {line_no} must have two "
            "dates"
        )

    return {
        "exact": None,
        "min": normalize_calibration_value(parts[0], line_no),
        "max": normalize_calibration_value(parts[1], line_no),
    }


def apply_operator_constraint(constraints, operator, line_no):
    if operator is None or operator == "=":
        return constraints

    if constraints["exact"] is None:
        sys.exit(
            f"[ERROR] Calibration line {line_no} cannot combine {operator} "
            "with l(...), u(...), or b(...)"
        )

    exact = constraints["exact"]

    if operator == ">=":
        return {"exact": None, "min": exact, "max": None}

    return {"exact": None, "min": None, "max": exact}


def parse_calibration_constraints(value, operator, line_no):
    constraints = parse_constraint_value(value, line_no)
    return apply_operator_constraint(constraints, operator, line_no)


def strip_inline_comment(line):
    return line.split("#", 1)[0].strip()


def parse_calibration_line(line, line_no):
    stripped = strip_inline_comment(line)

    if not stripped:
        return {"type": "skip"}

    if re.match(r"^(mrca|min|max)\s*=", stripped, re.IGNORECASE):
        return {"type": "treepl", "line": line.rstrip("\n")}

    mrca_match = re.match(
        r"^mrca\s*\(\s*([^()]+?)\s*\)\s*(>=|<=|=)?\s*(.+)$",
        stripped,
        re.IGNORECASE
    )

    if mrca_match:
        taxa = tuple(
            taxon.strip()
            for taxon in mrca_match.group(1).split(",")
            if taxon.strip()
        )

        if len(taxa) < 2:
            sys.exit(
                f"[ERROR] Calibration MRCA on line {line_no} must include "
                "at least two terminal taxa"
            )

        operator = mrca_match.group(2)

        return {
            "type": "mrca",
            "taxa": taxa,
            "constraints": parse_calibration_constraints(
                mrca_match.group(3),
                operator,
                line_no
            )
        }

    label_match = re.match(
        r"^([^\s<>=]+)\s*(?:(>=|<=|=)\s*)?(.+)$",
        stripped
    )

    if label_match:
        operator = label_match.group(2)

        return {
            "type": "label",
            "label": label_match.group(1),
            "constraints": parse_calibration_constraints(
                label_match.group(3),
                operator,
                line_no
            )
        }

    sys.exit(
        f"[ERROR] Unsupported calibration format on line {line_no}: "
        f"{line.rstrip()}"
    )


def parse_calibrations(calibrations):
    entries = []

    with open(calibrations) as f:
        for line_no, line in enumerate(f, start=1):
            entries.append(parse_calibration_line(line, line_no))

    return entries


def load_tree(tree_path):
    try:
        import dendropy
    except ImportError:
        sys.exit("[ERROR] DendroPy is required for MRCA calibration handling")

    return dendropy.Tree.get(
        path=tree_path,
        schema="newick",
        rooting="force-rooted",
        preserve_underscores=True
    )


def write_tree(tree, output):
    with open(output, "w") as f:
        f.write(
            tree.as_string(
                schema="newick",
                suppress_rooting=True,
                suppress_annotations=True,
                unquoted_underscores=True
            )
        )


def get_leaf_label(node):
    if node.taxon is not None:
        return node.taxon.label

    return node.label


def find_node_by_label(tree, label):
    for node in tree.preorder_node_iter():
        if node.label == label:
            return node

        if node.taxon is not None and node.taxon.label == label:
            return node

    return None


def first_leaf_label(node):
    for leaf in node.leaf_iter():
        return get_leaf_label(leaf)

    return None


def representative_mrca_taxa(node, label):
    children = list(node.child_node_iter())

    if len(children) < 2:
        sys.exit(
            f"[ERROR] Calibration label {label} does not identify an internal "
            "node with at least two child clades"
        )

    taxa = []

    for child in children:
        taxon = first_leaf_label(child)

        if taxon is not None:
            taxa.append(taxon)

        if len(taxa) == 2:
            return taxa[0], taxa[1]

    sys.exit(f"[ERROR] Could not identify descendant taxa for {label}")


def existing_labels(tree):
    labels = set()

    for node in tree.preorder_node_iter():
        if node.label:
            labels.add(node.label)

        if node.taxon is not None and node.taxon.label:
            labels.add(node.taxon.label)

    return labels


def next_mrca_label(labels, index):
    while True:
        candidate = f"mrca{index}"

        if candidate not in labels:
            labels.add(candidate)
            return candidate, index + 1

        index += 1


def get_mrca_node(tree, taxa):
    try:
        return tree.mrca(taxon_labels=list(taxa))
    except Exception as exc:
        sys.exit(
            "[ERROR] Could not find MRCA for "
            f"{', '.join(taxa)}: {exc}"
        )


def format_mrca_label(taxa):
    return f"mrca({','.join(taxa)})"


def add_bound(bounds, key, label, constraints):
    if key not in bounds:
        bounds[key] = {"label": label, "min": None, "max": None, "exact": None}

    if constraints["exact"] is not None:
        bounds[key]["exact"] = constraints["exact"]
        bounds[key]["min"] = constraints["exact"]
        bounds[key]["max"] = constraints["exact"]
        return

    if constraints["min"] is not None:
        bounds[key]["min"] = constraints["min"]

    if constraints["max"] is not None:
        bounds[key]["max"] = constraints["max"]


def validate_bound(label, lower, upper):
    if (
        lower is not None
        and upper is not None
        and is_numeric_value(lower)
        and is_numeric_value(upper)
        and float(lower) > float(upper)
    ):
        sys.exit(
            f"[ERROR] Calibration minimum is greater than maximum for {label}: "
            f"{lower} > {upper}"
        )


def validate_numeric_calibration(label, bounds, method):
    values = [
        value
        for value in (bounds["exact"], bounds["min"], bounds["max"])
        if value is not None
    ]

    if method == "mdcat" and bounds["exact"] is not None:
        return

    if method != "lsd2" and any(not is_numeric_value(value) for value in values):
        sys.exit(
            f"[ERROR] Calendar-date calibrations are supported only with "
            f"--method lsd2 or exact-date --method mdcat; {label} has a "
            "non-numeric date"
        )


def labeled_calibration_age(bounds, method):
    label = bounds["label"]
    lower = bounds["min"]
    upper = bounds["max"]
    exact = bounds["exact"]

    validate_bound(label, lower, upper)
    validate_numeric_calibration(label, bounds, method)

    if exact is not None:
        return exact

    if method in {"mdcat", "wlogdate"}:
        sys.exit(
            "[ERROR] Minimum/maximum calibration bounds are not supported "
            f"with --method {method}; use exact MRCA ages instead"
        )

    if lower is not None and upper is not None:
        return f"b({lower},{upper})"

    if lower is not None:
        return f"l({lower})"

    return f"u({upper})"


def prepare_treepl_calibrations(entries, tree_path, output):
    tree = None
    out_lines = []
    calibration_bounds = {}
    calibration_index = 1

    for entry in entries:
        if entry["type"] == "skip":
            continue

        if entry["type"] == "treepl":
            out_lines.append(entry["line"])
            continue

        if entry["type"] == "mrca" and len(entry["taxa"]) == 2:
            taxon1, taxon2 = entry["taxa"]
            key = ("taxa", tuple(sorted(entry["taxa"])))
        else:
            if tree is None:
                tree = load_tree(tree_path)

            if entry["type"] == "mrca":
                label = format_mrca_label(entry["taxa"])
                node = get_mrca_node(tree, entry["taxa"])
            else:
                label = entry["label"]
                node = find_node_by_label(tree, label)

            if node is None:
                sys.exit(
                    f"[ERROR] Calibration label not found in tree: "
                    f"{label}"
                )

            taxon1, taxon2 = representative_mrca_taxa(node, label)
            key = ("node", id(node))

        if key not in calibration_bounds:
            calibration_name = f"pipeline_calib{calibration_index}"
            calibration_index += 1
            calibration_bounds[key] = {
                "label": calibration_name,
                "taxon1": taxon1,
                "taxon2": taxon2,
                "min": None,
                "max": None,
                "exact": None,
            }

        add_bound(
            calibration_bounds,
            key,
            calibration_bounds[key]["label"],
            entry["constraints"]
        )

    for bounds in calibration_bounds.values():
        validate_bound(bounds["label"], bounds["min"], bounds["max"])
        validate_numeric_calibration(bounds["label"], bounds, "treepl")
        out_lines.extend([
            f"mrca = {bounds['label']} {bounds['taxon1']} {bounds['taxon2']}",
        ])

        if bounds["min"] is not None:
            out_lines.append(f"min = {bounds['label']} {bounds['min']}")

        if bounds["max"] is not None:
            out_lines.append(f"max = {bounds['label']} {bounds['max']}")

    with open(output, "w") as f:
        f.write("\n".join(out_lines))

        if out_lines:
            f.write("\n")

    return output


def prepare_treepl_unit_calibration(tree_path, output):
    tree = load_tree(tree_path)
    taxon1, taxon2 = representative_mrca_taxa(tree.seed_node, "root")

    with open(output, "w") as f:
        f.write(f"mrca = pipeline_root {taxon1} {taxon2}\n")
        f.write("min = pipeline_root 1\n")
        f.write("max = pipeline_root 1\n")

    return output


def prepare_labeled_calibrations(
    entries,
    tree_path,
    output_tree,
    output_calibrations,
    method
):
    tree = None
    labels = None
    label_index = 1
    mrca_labels = {}
    calibration_bounds = {}
    has_mrca = any(entry["type"] == "mrca" for entry in entries)

    if has_mrca:
        tree = load_tree(tree_path)
        labels = existing_labels(tree)

    for entry in entries:
        if entry["type"] == "skip":
            continue

        if entry["type"] == "treepl":
            sys.exit(
                "[ERROR] TreePL-style calibration lines are supported only "
                "with --method treepl"
            )

        if entry["type"] == "label":
            if tree is not None and find_node_by_label(tree, entry["label"]) is None:
                sys.exit(
                    f"[ERROR] Calibration label not found in tree: "
                    f"{entry['label']}"
                )

            key = ("label", entry["label"])
            add_bound(
                calibration_bounds,
                key,
                entry["label"],
                entry["constraints"]
            )
            continue

        node = get_mrca_node(tree, entry["taxa"])
        node_key = id(node)

        if node_key not in mrca_labels:
            label, label_index = next_mrca_label(labels, label_index)
            node.label = label
            mrca_labels[node_key] = label

        key = ("node", node_key)
        add_bound(
            calibration_bounds,
            key,
            mrca_labels[node_key],
            entry["constraints"]
        )

    out_lines = [
        f"{bounds['label']} {labeled_calibration_age(bounds, method)}"
        for bounds in calibration_bounds.values()
    ]

    if method == "lsd2":
        out_lines.insert(0, str(len(out_lines)))

    with open(output_calibrations, "w") as f:
        f.write("\n".join(out_lines))

        if out_lines:
            f.write("\n")

    if tree is None:
        return tree_path, output_calibrations

    write_tree(tree, output_tree)
    return output_tree, output_calibrations


def prepare_calibrations(tree_path, calibrations, method, intermediate_dir):
    if calibrations is None:
        if method == "treepl":
            output = os.path.join(intermediate_dir, "calibrations.treepl.txt")
            return tree_path, prepare_treepl_unit_calibration(tree_path, output)

        return tree_path, None

    entries = parse_calibrations(calibrations)
    parsed_entries = [entry for entry in entries if entry["type"] != "skip"]

    if not parsed_entries:
        sys.exit("[ERROR] Calibration file is empty")

    has_treepl_entries = any(entry["type"] == "treepl" for entry in parsed_entries)
    has_simple_entries = any(
        entry["type"] in {"label", "mrca"} for entry in parsed_entries
    )

    if has_treepl_entries and method != "treepl":
        sys.exit(
            "[ERROR] TreePL-style calibration lines are supported only with "
            "--method treepl"
        )

    if method == "treepl":
        if has_treepl_entries and not has_simple_entries:
            return tree_path, calibrations

        output = os.path.join(intermediate_dir, "calibrations.treepl.txt")
        return tree_path, prepare_treepl_calibrations(entries, tree_path, output)

    output_tree = os.path.join(intermediate_dir, "species_tree_su.calibrated.tre")
    output_calibrations = os.path.join(
        intermediate_dir,
        f"calibrations.{method}.txt"
    )

    return prepare_labeled_calibrations(
        entries,
        tree_path,
        output_tree,
        output_calibrations,
        method
    )


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


def run_treepl(
    tree,
    calibrations,
    output_dir,
    smooth=100,
    numsites=500000,
    options=None
):
    config = os.path.join(output_dir, "treepl.config")
    dated_tree = os.path.join(output_dir, "dated_tree.tre")

    with open(config, "w") as f:
        f.write(f"treefile = {tree}\n")
        f.write(f"smooth = {smooth}\n")
        f.write(f"numsites = {numsites}\n")
        f.write(f"outfile = {dated_tree}\n\n")

        with open(calibrations) as c:
            f.write(c.read())

        if options:
            f.write("\n")
            f.write("\n".join(options))
            f.write("\n")

    run(shell_join(["treePL", config]))

    return dated_tree


def calibration_file_has_calendar_dates(calibrations):
    with open(calibrations) as f:
        for line in f:
            if re.search(r"\b\d{4}-\d{1,2}(?:-\d{1,2})?\b", line):
                return True

    return False


def treepl_options_from_args(args):
    options = []

    flag_options = {
        "treepl_thorough": "thorough",
        "treepl_prime": "prime",
        "treepl_moredetailcvad": "moredetailcvad",
    }

    value_options = {
        "treepl_opt": "opt",
        "treepl_optad": "optad",
        "treepl_optcvad": "optcvad",
        "treepl_nthreads": "nthreads",
    }

    for attr, option in flag_options.items():
        if getattr(args, attr):
            options.append(option)

    for attr, option in value_options.items():
        value = getattr(args, attr)

        if value is not None:
            options.append(f"{option} = {value}")

    return options


def treepl_args_requested(args):
    if args.treepl_smooth != 100 or args.treepl_numsites != 500000:
        return True

    return bool(treepl_options_from_args(args))


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
    ]

    if calibrations is not None:
        cmd.extend(["-t", calibrations])

        has_calendar_dates = calibration_file_has_calendar_dates(calibrations)

        if has_calendar_dates:
            cmd.append("-d")
        else:
            cmd.append("-b")

    if seq_length is not None:
        cmd.extend(["-l", seq_length])

    if ci is not None:
        cmd.extend(["--CI", " ".join(ci)])

    run(shell_join(cmd))

    return dated_tree


def run_wlogdate(tree, calibrations, output_dir):
    dated_tree = os.path.join(output_dir, "dated_tree.tre")

    cmd = [
        "python",
        "launch_wLogDate.py",
        "-i",
        tree,
        "-o",
        dated_tree,
    ]

    if calibrations is not None:
        cmd.extend(["-t", calibrations, "-b"])

    run(shell_join(cmd))

    return dated_tree


def run_lsd2(
    tree,
    calibrations,
    output_dir,
    seq_length=None,
    min_branch_length=0.001
):
    prefix = os.path.join(output_dir, "lsd2")

    cmd = [
        "lsd2",
        "-i",
        tree,
    ]

    if calibrations is None:
        cmd.extend(["-a", 0, "-z", 1])
    else:
        cmd.extend(["-d", calibrations])

    if seq_length is not None:
        cmd.extend(["-s", seq_length])

    cmd.extend([
        "-u",
        min_branch_length,
        "-o",
        prefix,
    ])

    run(shell_join(cmd))

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
        default=None,
        help="Optional calibration file"
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
        default=".",
        help="Output directory (default: current directory)"
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
        "--mdcat-ci",
        "--CI",
        dest="CI",
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
        dest="seq_length",
        type=int,
        default=None,
        help="Optional sequence length passed to MD-Cat with -l or LSD2 with -s"
    )

    parser.add_argument(
        "--mdcat-seq-length",
        dest="seq_length",
        type=int,
        help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--lsd2-seq-length",
        dest="seq_length",
        type=int,
        help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--lsd2-min-branch-length",
        type=float,
        default=0.001,
        help="LSD2 -u minimum branch length value (default: 0.001)"
    )

    parser.add_argument(
        "--mdcat-p",
        type=int,
        default=10,
        help="MD-Cat -p value (default: 10)"
    )

    parser.add_argument(
        "--treepl-smooth",
        type=float,
        default=100,
        help="TreePL smooth value (default: 100)"
    )

    parser.add_argument(
        "--treepl-numsites",
        type=int,
        default=500000,
        help="TreePL numsites value (default: 500000)"
    )

    parser.add_argument(
        "--treepl-nthreads",
        type=int,
        default=None,
        help="TreePL nthreads value"
    )

    parser.add_argument(
        "--treepl-thorough",
        action="store_true",
        help="Add TreePL thorough option"
    )

    parser.add_argument(
        "--treepl-prime",
        action="store_true",
        help="Add TreePL prime option"
    )

    parser.add_argument(
        "--treepl-moredetailcvad",
        action="store_true",
        help="Add TreePL moredetailcvad option"
    )

    parser.add_argument(
        "--treepl-opt",
        type=int,
        default=None,
        help="TreePL opt value"
    )

    parser.add_argument(
        "--treepl-optad",
        type=int,
        default=None,
        help="TreePL optad value"
    )

    parser.add_argument(
        "--treepl-optcvad",
        type=int,
        default=None,
        help="TreePL optcvad value"
    )

    args = parser.parse_args()
    ci = parse_ci(args.CI)

    if ci is not None and args.method != "mdcat":
        sys.exit(
            "[ERROR] Confidence intervals are currently supported only with "
            "--method mdcat"
        )

    if args.seq_length is not None and args.method not in {"mdcat", "lsd2"}:
        sys.exit("[ERROR] --seq-length is supported only with MD-Cat or LSD2")

    if (
        args.calibrations is None
        and args.method == "lsd2"
        and args.seq_length is None
    ):
        sys.exit(
            "[ERROR] --seq-length is required with --method lsd2 when no "
            "calibration file is provided"
        )

    if args.method != "lsd2" and args.lsd2_min_branch_length != 0.001:
        sys.exit("[ERROR] --lsd2-min-branch-length is supported only with LSD2")

    if args.method != "treepl" and treepl_args_requested(args):
        sys.exit("[ERROR] --treepl-* options are supported only with TreePL")

    if args.seq_length is not None and args.seq_length <= 0:
        sys.exit("[ERROR] --seq-length must be positive")

    if args.lsd2_min_branch_length <= 0:
        sys.exit("[ERROR] --lsd2-min-branch-length must be positive")

    if args.gene_length is not None and args.gene_length <= 0:
        sys.exit("[ERROR] --gene-length must be positive")

    if args.mdcat_p <= 0:
        sys.exit("[ERROR] --mdcat-p must be positive")

    if args.treepl_smooth <= 0:
        sys.exit("[ERROR] --treepl-smooth must be positive")

    if args.treepl_numsites <= 0:
        sys.exit("[ERROR] --treepl-numsites must be positive")

    positive_treepl_ints = {
        "--treepl-nthreads": args.treepl_nthreads,
    }

    for option, value in positive_treepl_ints.items():
        if value is not None and value <= 0:
            sys.exit(f"[ERROR] {option} must be positive")

    check_exists(args.gene_trees, "Gene trees")

    if args.calibrations is not None:
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

    dating_input, dating_calibrations = prepare_calibrations(
        su_tree,
        args.calibrations,
        args.method,
        intermediate_dir
    )

    print("\n=== STEP 2: Molecular dating ===\n")

    if args.method == "treepl":
        run_treepl(
            dating_input,
            dating_calibrations,
            args.output,
            smooth=args.treepl_smooth,
            numsites=args.treepl_numsites,
            options=treepl_options_from_args(args)
        )

    elif args.method == "mdcat":
        run_mdcat(
            dating_input,
            dating_calibrations,
            args.output,
            ci=ci,
            seq_length=args.seq_length,
            p=args.mdcat_p
        )

    elif args.method == "wlogdate":
        run_wlogdate(
            dating_input,
            dating_calibrations,
            args.output
        )

    elif args.method == "lsd2":
        run_lsd2(
            dating_input,
            dating_calibrations,
            args.output,
            seq_length=args.seq_length,
            min_branch_length=args.lsd2_min_branch_length
        )

    print("\nPipeline completed successfully.\n")
    print(f"SU tree: {su_tree}")
    print(f"Dating input tree: {dating_input}")
    print(f"Dating calibrations: {dating_calibrations}")
    print(f"Dated tree: {os.path.join(args.output, 'dated_tree.tre')}")

    if ci is not None:
        print(
            "Confidence intervals: "
            f"{ci[0]} MD-Cat samples, quantiles {ci[1]} and {ci[2]}"
        )


if __name__ == "__main__":
    main()
