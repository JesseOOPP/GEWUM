"""Tests for SLURM template rendering utilities (commands.template_utils)."""

import os

from gewum.commands.template_utils import (
    CONFIG_FILENAME,
    copy_with_template,
    find_config_file,
    generate_env_setup,
    generate_slurm_header,
    get_config_info,
    load_config,
    process_shell_script,
)

MINIMAL_CONFIG = {
    "slurm": {
        "time": "12:00:00",
        "cpus_per_task": 8,
        "partition": "test_part",
        "nodes": 2,
    },
    "environment": {
        "module_purge": True,
        "modules": ["gcc/12.2"],
        "conda_path": "/opt/conda",
        "conda_env": "gewum",
    },
    "parallel": {"path": "/usr/bin/parallel"},
}


def _write_config(tmp_path, data=None):
    path = tmp_path / CONFIG_FILENAME
    path.parent.mkdir(parents=True, exist_ok=True)
    import yaml

    path.write_text(yaml.safe_dump(data or MINIMAL_CONFIG), encoding="utf-8")
    return path


# --- find_config_file / load_config -------------------------------------------


def test_find_config_file_prefers_dest_dir(tmp_path, monkeypatch):
    dest_cfg = _write_config(tmp_path / "work")
    root_cfg = _write_config(tmp_path)
    monkeypatch.chdir(tmp_path)
    found = find_config_file(dest_dir=str(tmp_path / "work"))
    assert os.path.abspath(found) == os.path.abspath(str(dest_cfg))
    assert root_cfg.exists()


def test_find_config_file_falls_back_to_installation_dir(tmp_path, monkeypatch):
    # No config in dest dir or cwd: falls back to the GEWUM installation dir,
    # which contains the bundled slurm_config.yaml.
    monkeypatch.chdir(tmp_path)
    found = find_config_file(dest_dir=str(tmp_path))
    assert found is not None
    assert os.path.exists(found)


def test_find_config_file_returns_none_when_missing(monkeypatch):
    monkeypatch.setattr(
        "gewum.commands.template_utils.os.path.exists", lambda p: False
    )
    assert find_config_file(dest_dir=".") is None


def test_load_config_returns_none_for_missing_file(tmp_path):
    assert load_config(config_path=str(tmp_path / "nope.yaml")) is None


def test_load_config_parses_yaml(tmp_path):
    cfg_path = _write_config(tmp_path)
    cfg = load_config(config_path=str(cfg_path))
    assert cfg["slurm"]["cpus_per_task"] == 8
    assert cfg["environment"]["conda_env"] == "gewum"


def test_get_config_info(tmp_path):
    cfg_path = _write_config(tmp_path)
    path, cfg = get_config_info(dest_dir=str(tmp_path))
    assert path == str(cfg_path)
    assert cfg["slurm"]["partition"] == "test_part"


# --- generate_slurm_header ----------------------------------------------------


def test_generate_slurm_header_none_config():
    assert generate_slurm_header(None) is None


def test_generate_slurm_header_defaults():
    header = generate_slurm_header({}, job_name="myjob")
    assert header.startswith("#!/bin/bash")
    assert "#SBATCH --job-name=myjob" in header
    assert "#SBATCH --time=2400:00:00" in header  # default
    assert "#SBATCH --cpus-per-task=64" in header  # default
    assert "#SBATCH -p <partition>" in header  # default


def test_generate_slurm_header_custom():
    header = generate_slurm_header(MINIMAL_CONFIG, job_name="jobA")
    assert "#SBATCH --time=12:00:00" in header
    assert "#SBATCH --cpus-per-task=8" in header
    assert "#SBATCH -p test_part" in header
    assert "#SBATCH -N 2" in header


# --- generate_env_setup -------------------------------------------------------


def test_generate_env_setup_none_config():
    assert generate_env_setup(None) is None


def test_generate_env_setup_full():
    setup = generate_env_setup(MINIMAL_CONFIG)
    assert "module purge" in setup
    assert "module load gcc/12.2" in setup
    assert "source /opt/conda/etc/profile.d/conda.sh" in setup
    assert "conda activate gewum" in setup
    assert "export PATH=/usr/bin/parallel:$PATH" in setup


def test_generate_env_setup_empty():
    assert generate_env_setup({}) == ""


# --- process_shell_script -----------------------------------------------------


def test_process_shell_script_returns_unchanged_without_config():
    content = "echo hello\n{{SLURM_TIME}}\n"
    assert process_shell_script(content, None) == content


def test_process_shell_script_replaces_placeholders():
    script = (
        "#SBATCH --job-name=relax\n"
        "{{SLURM_TIME}} {{SLURM_CPUS}} {{SLURM_PARTITION}} "
        "{{SLURM_NODES}} {{JOB_NAME}}\n"
    )
    out = process_shell_script(script, MINIMAL_CONFIG)
    assert "12:00:00" in out
    assert "8" in out
    assert "test_part" in out
    assert "2" in out
    assert "relax" in out


def test_process_shell_script_injects_header_and_env():
    script = "{{SLURM_HEADER}}\n{{ENV_SETUP}}\nrun_cmd\n"
    out = process_shell_script(script, MINIMAL_CONFIG, job_name="tc")
    assert "#SBATCH --job-name=tc" in out
    assert "#SBATCH --time=12:00:00" in out
    assert "module purge" in out
    assert "conda activate gewum" in out
    assert "run_cmd" in out


# --- copy_with_template -------------------------------------------------------


def test_copy_with_template_renders_sh(tmp_path):
    src = tmp_path / "template.sh"
    dst = tmp_path / "out.sh"
    src.write_text(
        "#!/bin/bash\n#SBATCH --job-name=default\n{{SLURM_TIME}}\n",
        encoding="utf-8",
    )
    ok = copy_with_template(str(src), str(dst), config=MINIMAL_CONFIG)
    assert ok
    content = dst.read_text(encoding="utf-8")
    assert "12:00:00" in content


def test_copy_with_template_plain_copy(tmp_path):
    src = tmp_path / "data.yaml"
    dst = tmp_path / "data_out.yaml"
    src.write_text("a: 1\n", encoding="utf-8")
    ok = copy_with_template(str(src), str(dst), config=MINIMAL_CONFIG)
    assert ok
    assert dst.read_text(encoding="utf-8") == "a: 1\n"


def test_copy_with_template_missing_source(tmp_path):
    ok = copy_with_template(
        str(tmp_path / "nope.sh"), str(tmp_path / "out.sh"), config=MINIMAL_CONFIG
    )
    assert not ok
