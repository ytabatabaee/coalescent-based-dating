import os
import shlex
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import run_pipeline


def capture_commands(monkeypatch):
    commands = []
    monkeypatch.setattr(run_pipeline, "run", commands.append)
    return commands


def write_file(path, contents):
    with open(path, "w") as handle:
        handle.write(contents)


def split_command(command):
    return shlex.split(command)


def test_mdcat_numeric_calibration_uses_paper_options(tmp_path, monkeypatch):
    commands = capture_commands(monkeypatch)
    calibrations = tmp_path / "calibrations.txt"
    write_file(calibrations, "mrca1 1\n")

    run_pipeline.run_mdcat(
        "tree.nwk",
        str(calibrations),
        str(tmp_path),
        ci=("1000", "0.025", "0.975"),
        seq_length=40000,
        p=10,
    )

    args = split_command(commands[0])
    assert args[:6] == [
        "python3",
        "md_cat.py",
        "-i",
        "tree.nwk",
        "-o",
        str(tmp_path / "dated_tree.tre"),
    ]
    assert args[args.index("-p") + 1] == "10"
    assert args[args.index("-t") + 1] == str(calibrations)
    assert "-b" in args
    assert "-d" not in args
    assert args[args.index("-l") + 1] == "40000"
    assert args[args.index("--CI") + 1] == "1000 0.025 0.975"


def test_mdcat_calendar_dates_use_date_mode(tmp_path, monkeypatch):
    commands = capture_commands(monkeypatch)
    calibrations = tmp_path / "calibrations.txt"
    write_file(calibrations, "A 2013-01-01\n")

    run_pipeline.run_mdcat("tree.nwk", str(calibrations), str(tmp_path), p=10)

    args = split_command(commands[0])
    assert args[args.index("-t") + 1] == str(calibrations)
    assert "-d" in args
    assert "-b" not in args


def test_mdcat_without_calibrations_omits_calibration_options(tmp_path, monkeypatch):
    commands = capture_commands(monkeypatch)

    run_pipeline.run_mdcat("tree.nwk", None, str(tmp_path), seq_length=40000, p=10)

    args = split_command(commands[0])
    assert "-t" not in args
    assert "-b" not in args
    assert "-d" not in args
    assert args[args.index("-l") + 1] == "40000"


def test_wlogdate_calibrated_and_uncalibrated_options(tmp_path, monkeypatch):
    commands = capture_commands(monkeypatch)
    calibrations = tmp_path / "calibrations.txt"
    write_file(calibrations, "mrca1 1\n")

    run_pipeline.run_wlogdate("tree.nwk", str(calibrations), str(tmp_path))
    run_pipeline.run_wlogdate("tree.nwk", None, str(tmp_path))

    calibrated = split_command(commands[0])
    assert calibrated[:4] == ["python", "launch_wLogDate.py", "-i", "tree.nwk"]
    assert calibrated[calibrated.index("-t") + 1] == str(calibrations)
    assert "-b" in calibrated

    uncalibrated = split_command(commands[1])
    assert "-t" not in uncalibrated
    assert "-b" not in uncalibrated


def test_lsd2_calibrated_uses_supplement_options(tmp_path, monkeypatch):
    commands = capture_commands(monkeypatch)
    calibrations = tmp_path / "calibrations.txt"
    write_file(calibrations, "1\nmrca1 1\n")

    run_pipeline.run_lsd2(
        "tree.nwk",
        str(calibrations),
        str(tmp_path),
        seq_length=40000,
        min_branch_length=0.001,
    )

    args = split_command(commands[0])
    assert args[:4] == ["lsd2", "-i", "tree.nwk", "-d"]
    assert args[args.index("-d") + 1] == str(calibrations)
    assert args[args.index("-s") + 1] == "40000"
    assert args[args.index("-u") + 1] == "0.001"
    assert args[args.index("-o") + 1] == str(tmp_path / "lsd2")
    assert "-a" not in args
    assert "-z" not in args


def test_lsd2_without_calibrations_uses_unit_ultrametric_options(tmp_path, monkeypatch):
    commands = capture_commands(monkeypatch)

    run_pipeline.run_lsd2(
        "tree.nwk",
        None,
        str(tmp_path),
        seq_length=40000,
        min_branch_length=0.001,
    )

    args = split_command(commands[0])
    assert args[:4] == ["lsd2", "-i", "tree.nwk", "-a"]
    assert args[args.index("-a") + 1] == "0"
    assert args[args.index("-z") + 1] == "1"
    assert args[args.index("-s") + 1] == "40000"
    assert args[args.index("-u") + 1] == "0.001"
    assert "-d" not in args


def test_treepl_options_from_args_are_limited_to_paper_workflow():
    args = SimpleNamespace(
        treepl_thorough=True,
        treepl_prime=True,
        treepl_moredetailcvad=True,
        treepl_opt=5,
        treepl_optad=5,
        treepl_optcvad=1,
        treepl_nthreads=16,
    )

    assert run_pipeline.treepl_options_from_args(args) == [
        "thorough",
        "prime",
        "moredetailcvad",
        "opt = 5",
        "optad = 5",
        "optcvad = 1",
        "nthreads = 16",
    ]


def test_treepl_config_includes_calibrations_and_paper_options(tmp_path, monkeypatch):
    commands = capture_commands(monkeypatch)
    calibrations = tmp_path / "calibrations.treepl.txt"
    write_file(calibrations, "mrca = calib A B\nmin = calib 1\nmax = calib 2\n")

    run_pipeline.run_treepl(
        "tree.nwk",
        str(calibrations),
        str(tmp_path),
        smooth=1000,
        numsites=63430000,
        options=["thorough", "prime", "moredetailcvad", "opt = 5"],
    )

    config = (tmp_path / "treepl.config").read_text()
    assert "treefile = tree.nwk\n" in config
    assert "smooth = 1000\n" in config
    assert "numsites = 63430000\n" in config
    assert f"outfile = {tmp_path / 'dated_tree.tre'}\n" in config
    assert "mrca = calib A B\nmin = calib 1\nmax = calib 2\n" in config
    assert "thorough\nprime\nmoredetailcvad\nopt = 5\n" in config
    assert split_command(commands[0]) == ["treePL", str(tmp_path / "treepl.config")]


def test_treepl_without_user_calibrations_creates_root_unit_calibration(tmp_path):
    tree = tmp_path / "tree.nwk"
    write_file(tree, "((A:1,B:1):1,(C:1,D:1):1);\n")

    dating_tree, calibrations = run_pipeline.prepare_calibrations(
        str(tree),
        None,
        "treepl",
        str(tmp_path),
    )

    assert dating_tree == str(tree)
    generated = (tmp_path / "calibrations.treepl.txt").read_text().splitlines()
    assert generated == [
        "mrca = pipeline_root A C",
        "min = pipeline_root 1",
        "max = pipeline_root 1",
    ]
    assert calibrations == str(tmp_path / "calibrations.treepl.txt")


def test_non_treepl_without_user_calibrations_passes_none(tmp_path):
    tree = tmp_path / "tree.nwk"
    write_file(tree, "((A:1,B:1):1,(C:1,D:1):1);\n")

    dating_tree, calibrations = run_pipeline.prepare_calibrations(
        str(tree),
        None,
        "mdcat",
        str(tmp_path),
    )

    assert dating_tree == str(tree)
    assert calibrations is None
