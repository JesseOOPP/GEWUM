"""Mini workflow regression tests (local, non-SLURM mode).

These tests exercise the file-copy + template-rendering layer of each CLI
command without needing a SLURM cluster or a uMLIP backend: they verify that
`gewum <command> --mode <mode> --dest <dir>` produces the expected scripts in
the destination directory and that SLURM placeholders are rendered.
"""

from types import SimpleNamespace

import pytest

from gewum.commands import ELA, MD, PT, QHA, RD, TC


@pytest.fixture
def work_dir(tmp_path):
    dest = tmp_path / "work"
    dest.mkdir()
    return dest


def _run_copy(module, mode, dest, **extra):
    args = SimpleNamespace(mode=mode, dest=str(dest), **extra)
    module.execute(args)
    return dest


def _assert_copied(dest, names):
    for name in names:
        assert (dest / name).exists(), f"expected copied file: {name}"


@pytest.mark.parametrize(
    "module,mode,expected",
    [
        (RD, "cifgen", ["cifgen.sh"]),
        (RD, "relax", ["relax_umlip.sh"]),
        (RD, "refine", ["refine_umlip.sh"]),
        (RD, "post", ["post_relax.sh"]),
        (RD, "auto", ["cifgen.sh", "relax_umlip.sh", "post_relax.sh", "run_srss.sh"]),
        (PT, "sub", ["replacements.yaml"]),
        (PT, "mutate", ["INPUT"]),
        (PT, "dp", ["doping.yaml"]),
        (PT, "relax", ["relax_umlip_pt.sh"]),
        (ELA, "cal", ["cal_ela.sh", "VPKIT.in1", "VPKIT.in2"]),
        (QHA, "cal", ["cal_qha.sh"]),
        (TC, "fc3", ["tc_cal.sh"]),
        (TC, "post", ["tc_post.sh"]),
        (MD, "nvt", ["run_md_nvt.sh"]),
        (MD, "post", ["post_md.sh"]),
    ],
)
def test_mode_copies_expected_scripts(work_dir, module, mode, expected):
    _run_copy(module, mode, work_dir)
    _assert_copied(work_dir, expected)


def test_rd_script_rendered_with_slurm_config(work_dir):
    import yaml

    # Provide a local slurm_config.yaml so template rendering is applied.
    cfg = {
        "slurm": {"time": "01:00:00", "cpus_per_task": 4, "partition": "tiny"},
        "environment": {"module_purge": False, "modules": []},
    }
    (work_dir / "slurm_config.yaml").write_text(
        yaml.safe_dump(cfg), encoding="utf-8"
    )
    _run_copy(RD, "cifgen", work_dir)
    content = (work_dir / "cifgen.sh").read_text(encoding="utf-8")
    assert "#SBATCH --time=01:00:00" in content
    assert "#SBATCH --cpus-per-task=4" in content
    assert "#SBATCH -p tiny" in content


def test_missing_mode_reports_no_files_copied(work_dir, capsys):
    # A mode with an empty file map copies nothing but must not raise.
    _run_copy(RD, "reorder", work_dir)
    out = capsys.readouterr().out
    assert "reorder" in out
