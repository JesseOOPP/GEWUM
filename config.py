"""GEWUM Configuration Module
Centralized configuration for all workflow paths and file mappings
"""
import os


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GEWUM_SOURCE_DIR = os.path.join(BASE_DIR, "src")

COMMON_REPOSITORY = os.path.join(GEWUM_SOURCE_DIR, "common")
RD_REPOSITORY = os.path.join(GEWUM_SOURCE_DIR, "RDworkflows")
PT_REPOSITORY = os.path.join(GEWUM_SOURCE_DIR, "PTworkflows")
CIFGEN_MODULE = os.path.join(GEWUM_SOURCE_DIR, "cifgen_input")
ELA_REPOSITORY = os.path.join(GEWUM_SOURCE_DIR, "ELAworkflows")  
QHA_REPOSITORY = os.path.join(GEWUM_SOURCE_DIR, "QHAworkflows")
MD_REPOSITORY = os.path.join(GEWUM_SOURCE_DIR, "MDworkflows")

COMMON_FILE_MAP = {
    "relax": [
    ],
    "refine": [
    ],
    "relaxhp": [
    ],
    "post": [
    ],
    "posthp": [
    ],
    "Ehull": [
        "ehull/Ehull_compatibility.py",
        "ehull/Ehull_no_compatibility.py", 
        "ehull/Ehull_post.py",
    ],
    "sym": [
    ],
    "sym2d": [
    ],
    "ph": [
        "phonon/ph_cal.sh",
    ],
    "phpost": [
    ],
    "select": [
        "selection/run_selection.sh",
    ],
    "auto": [
        "selection/run_selection.sh",
    ],
    "viz": [
        "postprocess/run_visualize.sh",
    ],
    "viz2": ["postprocess/run_viz2.sh"]
}


def ensure_directories():
    """Ensure all required directories exist"""
    directories = [
        COMMON_REPOSITORY,
        RD_REPOSITORY, 
        PT_REPOSITORY, 
        CIFGEN_MODULE, 
        ELA_REPOSITORY, 
        QHA_REPOSITORY,
        MD_REPOSITORY,
    ]
    for directory in directories:
        os.makedirs(directory, exist_ok=True)


ensure_directories()
