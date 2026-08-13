"""
GEWUM Setup Script
For backward compatibility with older pip versions
"""

from setuptools import setup, find_packages
import os

# Read the README for long description
readme_path = os.path.join(os.path.dirname(__file__), "GEWUM_Manual.md")
if os.path.exists(readme_path):
    with open(readme_path, encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "GEWUM: General Exploration Workflow for the Utopia of Materials"

# Find all subpackages and add gewum prefix
subpackages = find_packages(where=".")
packages = ["gewum"] + [f"gewum.{pkg}" for pkg in subpackages]

setup(
    name="gewum",
    version="1.0.0",
    author="Jiexi Song",
    author_email="songjx@szlab.ac.cn",
    description="GEWUM: General Exploration Workflow for the Utopia of Materials",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JesseOOPP/GEWUM",
    packages=packages,
    package_dir={"gewum": "."},
    include_package_data=True,
    package_data={
        "gewum": [
            "src/**/*.py",
            "src/**/*.sh",
            "src/**/*.in*",
            "slurm_config.yaml",
            "replacements.yaml",
        ],
    },
    python_requires=">=3.8",
    install_requires=[
        # Core scientific computing
        "numpy>=1.20",
        "pandas>=1.3",
        "matplotlib>=3.4",
        "scipy>=1.7",
        # Atomistic simulation & crystal structure
        "ase>=3.22",
        "pymatgen>=2022.0",
        "pyyaml>=5.4",
        "spglib>=2.0",
        "pyxtal>=0.6",
        # Phonon calculations
        "phonopy>=2.20",
        "phono3py>=2.0",
        # Materials Project API & offline data parsing
        "mp-api>=0.30",
        "ijson>=3.1",
        # Machine learning (structure selection & visualization)
        "scikit-learn>=1.0",
        "dscribe>=1.2",
        "hdbscan>=0.8",
        "umap-learn>=0.5",
        # Utilities
        "tqdm>=4.60",
    ],
    extras_require={
        "ml": ["mattersim>=1.0"],
        "dev": ["pytest", "black", "flake8"],
    },
    entry_points={
        "console_scripts": [
            "gewum=gewum.main:main",
        ],
    },
)
