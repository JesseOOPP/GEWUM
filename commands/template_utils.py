"""
GEWUM Template Utilities
Handles SLURM configuration and shell script template processing
"""
import os
import yaml
import re
from pathlib import Path

CONFIG_FILENAME = "slurm_config.yaml"


def find_config_file(dest_dir="."):
    """
    Find slurm_config.yaml in order of priority:
    1. Current working directory / destination directory
    2. GEWUM installation directory
    
    Args:
        dest_dir: Destination directory for file copy
    
    Returns:
        Path to config file or None
    """
    local_config = os.path.join(dest_dir, CONFIG_FILENAME)
    if os.path.exists(local_config):
        return local_config
    
    cwd_config = os.path.join(os.getcwd(), CONFIG_FILENAME)
    if os.path.exists(cwd_config):
        return cwd_config
    
    gewum_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gewum_config = os.path.join(gewum_dir, CONFIG_FILENAME)
    if os.path.exists(gewum_config):
        return gewum_config
    
    return None


def load_config(config_path=None, dest_dir="."):
    """
    Load SLURM configuration from yaml file
    
    Args:
        config_path: Explicit path to config file (optional)
        dest_dir: Destination directory to search for config
    
    Returns:
        dict: Configuration dictionary or None if not found
    """
    if config_path is None:
        config_path = find_config_file(dest_dir)
    
    if config_path is None or not os.path.exists(config_path):
        return None
    
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Warning: Failed to load config from {config_path}: {e}")
        return None


def generate_slurm_header(config, job_name="gewum_job"):
    """
    Generate SLURM header from configuration
    
    Args:
        config: Configuration dictionary
        job_name: Job name for this script
    
    Returns:
        str: SLURM header block
    """
    if config is None:
        return None
    
    slurm = config.get('slurm', {})
    
    lines = [
        "#!/bin/bash",
        "",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={job_name}.out",
        f"#SBATCH --error={job_name}.err",
        f"#SBATCH --time={slurm.get('time', '2400:00:00')}",
        f"#SBATCH --cpus-per-task={slurm.get('cpus_per_task', 64)}",
        f"#SBATCH -p {slurm.get('partition', '<partition>')}",
        f"#SBATCH -N {slurm.get('nodes', 1)}",
    ]
    
    return "\n".join(lines)


def generate_env_setup(config):
    """
    Generate environment setup commands from configuration
    
    Args:
        config: Configuration dictionary
    
    Returns:
        str: Environment setup block
    """
    if config is None:
        return None
    
    env = config.get('environment', {})
    lines = []
    
    if env.get('module_purge', True):
        lines.append("module purge")
    
    modules = env.get('modules', [])
    for module in modules:
        lines.append(f"module load {module}")
    
    conda_path = env.get('conda_path')
    conda_env = env.get('conda_env')
    if conda_path and conda_env:
        lines.append(f"source {conda_path}/etc/profile.d/conda.sh")
        lines.append(f"conda activate {conda_env}")
    
    extra_env = env.get('extra_env', [])
    for extra in extra_env:
        lines.append(extra)
    
    parallel_cfg = config.get('parallel', {})
    parallel_path = parallel_cfg.get('path')
    if parallel_path:
        lines.append(f"export PATH={parallel_path}:$PATH")
    
    return "\n".join(lines)


def process_shell_script(script_content, config, job_name=None):
    """
    Process shell script content, replacing template placeholders
    
    Recognizes these placeholder patterns:
    - {{SLURM_HEADER}} or lines starting with #SBATCH
    - {{ENV_SETUP}} or module/conda blocks
    - {{SLURM_TIME}}, {{SLURM_CPUS}}, etc. for individual values
    
    Args:
        script_content: Original script content
        config: Configuration dictionary
        job_name: Job name (extracted from script if not provided)
    
    Returns:
        str: Processed script content
    """
    if config is None:
        return script_content
    
    if job_name is None:
        match = re.search(r'#SBATCH\s+--job-name=(\S+)', script_content)
        if match:
            job_name = match.group(1)
        else:
            job_name = "gewum_job"
    
    slurm = config.get('slurm', {})
    env = config.get('environment', {})
    
    replacements = {
        '{{SLURM_TIME}}': slurm.get('time', '2400:00:00'),
        '{{SLURM_CPUS}}': str(slurm.get('cpus_per_task', 64)),
        '{{SLURM_PARTITION}}': slurm.get('partition', '<partition>'),
        '{{SLURM_NODES}}': str(slurm.get('nodes', 1)),
        '{{CONDA_PATH}}': env.get('conda_path', ''),
        '{{CONDA_ENV}}': env.get('conda_env', ''),
        '{{JOB_NAME}}': job_name,
    }
    
    result = script_content
    for placeholder, value in replacements.items():
        result = result.replace(placeholder, value)
    
    if '{{SLURM_HEADER}}' in result:
        header = generate_slurm_header(config, job_name)
        if header:
            result = result.replace('{{SLURM_HEADER}}', header)
    
    if '{{ENV_SETUP}}' in result:
        env_setup = generate_env_setup(config)
        if env_setup:
            result = result.replace('{{ENV_SETUP}}', env_setup)
    
    return result


def copy_with_template(src_path, dst_path, config=None, job_name=None):
    """
    Copy a shell script file with template processing
    
    Args:
        src_path: Source file path
        dst_path: Destination file path
        config: Configuration dictionary (loads automatically if None)
        job_name: Job name override
    
    Returns:
        bool: True if successful
    """
    try:
        if config is None:
            config = load_config(dest_dir=os.path.dirname(dst_path))
        
        with open(src_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        if src_path.endswith('.sh') and config is not None:
            if job_name is None:
                job_name = os.path.basename(src_path).replace('.sh', '')
            content = process_shell_script(content, config, job_name)
        
        with open(dst_path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(content)
        
        return True
        
    except Exception as e:
        print(f"Error copying {src_path}: {e}")
        return False


def get_config_info(dest_dir="."):
    """
    Get information about which config file will be used
    
    Args:
        dest_dir: Destination directory
    
    Returns:
        tuple: (config_path, config_dict) or (None, None)
    """
    config_path = find_config_file(dest_dir)
    if config_path:
        config = load_config(config_path)
        return config_path, config
    return None, None
