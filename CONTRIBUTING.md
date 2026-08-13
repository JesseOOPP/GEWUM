# Contributing to GEWUM
# 为 GEWUM 做贡献

Thank you for your interest in contributing to GEWUM (General Exploration Workflow for the Utopia of Materials). This document describes the development workflow and the standards we expect from contributions.

感谢您对 GEWUM 的关注与贡献。本文档说明了开发工作流与贡献标准。

## Table of Contents / 目录

- [Development Setup / 开发环境](#development-setup--开发环境)
- [Code Style / 代码风格](#code-style--代码风格)
- [Testing / 测试](#testing--测试)
- [Submitting a Pull Request / 提交 PR](#submitting-a-pull-request--提交-pr)
- [Reporting Issues / 报告问题](#reporting-issues--报告问题)

## Development Setup / 开发环境

1. Clone the repository:

   ```
   git clone https://github.com/JesseOOPP/GEWUM.git
   cd GEWUM
   ```

2. Create a dedicated environment (Python >= 3.8, <= 3.11) and install in editable mode with development dependencies:

   ```
   conda create -n gewum-dev python=3.10
   conda activate gewum-dev
   pip install -e ".[dev]"
   ```

3. (Optional) Install a uMLIP backend for relaxation-related tests, e.g.:

   ```
   pip install ".[ml]"        # MatterSim
   ```

## Code Style / 代码风格

- Python code must be formatted with [black](https://black.readthedocs.io/) and pass [flake8](https://flake8.pycqa.org/) checks.
- Keep all source code and comments ASCII-compatible; use English for code comments and docstrings.
- Follow the existing modular layout: CLI commands live in `commands/`, workflow scripts in `src/*workflows/`, shared modules in `src/common/`.
- Keep shell scripts (`src/*/**.sh`) bash-compatible (bash, not sh) and POSIX-style for local execution.
- Docstrings follow Google style; every public function should have a one-line summary.

## Testing / 测试

All changes must keep the test suite green:

```
pytest tests/ -v
```

- Unit tests live in `tests/` and use `pytest` + `tmp_path` fixtures; they must not require a SLURM cluster or a uMLIP backend unless explicitly marked.
- Before opening a PR, also run the linters locally:

  ```
  black --check .
  flake8 .
  ```

- The CI workflow (`.github/workflows/ci.yml`) runs the same suite on Python 3.8-3.11.

## Submitting a Pull Request / 提交 PR

1. Create a feature branch from `main`: `git checkout -b feat/your-change`.
2. Make your changes; add or update tests for the changed behaviour.
3. Run `pytest tests/`, `black --check .` and `flake8 .` locally.
4. Push and open a PR with a clear description of the change. Reference the related issue, if any.
5. CI must pass before the PR is merged.

## Reporting Issues / 报告问题

When reporting a bug, include:

- GEWUM version (`gewum --version`) and Python version
- Full command line used and the complete error output
- Relevant configuration (e.g., `slurm_config.yaml` snippets, with secrets removed)
- Minimal steps to reproduce
