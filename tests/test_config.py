"""Tests for the GEWUM configuration module (gewum.config)."""

import os

from gewum import config


def test_source_directory_exists():
    assert os.path.isdir(config.GEWUM_SOURCE_DIR)


def test_all_repositories_exist():
    repos = [
        config.COMMON_REPOSITORY,
        config.RD_REPOSITORY,
        config.PT_REPOSITORY,
        config.CIFGEN_MODULE,
        config.ELA_REPOSITORY,
        config.QHA_REPOSITORY,
        config.MD_REPOSITORY,
    ]
    for repo in repos:
        assert os.path.isdir(repo), f"missing repository: {repo}"


def test_common_file_map_structure():
    assert isinstance(config.COMMON_FILE_MAP, dict)
    assert config.COMMON_FILE_MAP, "COMMON_FILE_MAP must not be empty"
    for mode, files in config.COMMON_FILE_MAP.items():
        assert isinstance(mode, str) and mode
        assert isinstance(files, list)
        for rel in files:
            full = os.path.join(config.COMMON_REPOSITORY, rel)
            assert os.path.exists(full), f"missing common file: {rel}"


def test_ensure_directories_is_idempotent():
    # Must be safe to call repeatedly (module import already calls it once).
    config.ensure_directories()
    config.ensure_directories()
