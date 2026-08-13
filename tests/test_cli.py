"""CLI smoke tests: run the actual command-line entry point in a subprocess."""

import subprocess
import sys

EXPECTED_COMMANDS = [
    "cifgen",
    "RD",
    "PT",
    "ELA",
    "QHA",
    "TC",
    "MD",
    "DB",
]

VERSION_STRING = "GEWUM Version 1.0.0"


def _run_cli(*args):
    return subprocess.run(
        [sys.executable, "-m", "gewum.main", *args],
        capture_output=True,
        text=True,
        timeout=180,
    )


def test_version_flag():
    result = _run_cli("--version")
    assert result.returncode == 0
    assert VERSION_STRING in result.stdout


def test_help_without_command_lists_all():
    # With no command, main() prints help and exits with code 1.
    result = _run_cli("--help")
    assert result.returncode == 1
    output = result.stdout + result.stderr
    for cmd in EXPECTED_COMMANDS:
        assert cmd in output, f"command {cmd} missing from --help"


def test_each_command_help():
    for cmd in EXPECTED_COMMANDS:
        result = _run_cli(cmd, "-h")
        assert result.returncode == 0, f"{cmd} -h failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert cmd in output, f"{cmd} missing from its own help"
        assert "--mode" in output or cmd == "DB", f"{cmd} -h lacks --mode"
