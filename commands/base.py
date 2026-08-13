"""
GEWUM Base Command Module
Provides shared functionality for all workflow commands
"""
import os
import shutil
from abc import ABC, abstractmethod
from .template_utils import load_config, copy_with_template, get_config_info


class BaseWorkflowCommand(ABC):
    """
    Abstract base class for GEWUM workflow commands.
    Provides common file copying and execution logic.
    """
    
    @property
    @abstractmethod
    def repository_path(self) -> str:
        """Return the path to workflow-specific script repository"""
        pass
    
    @property
    @abstractmethod
    def file_map(self) -> dict:
        """Return mode-to-files mapping for workflow-specific files"""
        pass
    
    @property
    def workflow_name(self) -> str:
        """Return the workflow name for display purposes"""
        return self.__class__.__name__.replace('Command', '')
    
    def get_common_files(self, mode: str) -> list:
        """
        Return list of common files needed for the given mode.
        These files are shared across RD and SUB workflows.
        
        Args:
            mode: The operation mode (relax, post, Ehull, sym, ph, etc.)
        
        Returns:
            List of relative paths to common files
        """
        from ..config import COMMON_FILE_MAP
        return COMMON_FILE_MAP.get(mode, [])
    
    def _copy_file(self, filename: str, dest: str, copied: list, missing: list, 
                   is_common: bool = False, config=None) -> None:
        """
        Copy a single file from source to destination.
        For .sh files, apply template processing if config is available.
        
        Args:
            filename: Name of the file to copy
            dest: Destination directory
            copied: List to append successfully copied files
            missing: List to append missing files
            is_common: If True, copy from common repository
            config: Template configuration dictionary
        """
        from ..config import COMMON_REPOSITORY
        
        if is_common:
            src_path = os.path.join(COMMON_REPOSITORY, filename)
            dst_filename = os.path.basename(filename)
        else:
            src_path = os.path.join(self.repository_path, filename)
            dst_filename = filename
        
        dst_path = os.path.join(dest, dst_filename)
        
        if not os.path.exists(src_path):
            missing.append(filename)
            return
        
        dst_dir = os.path.dirname(dst_path)
        if dst_dir and not os.path.exists(dst_dir):
            os.makedirs(dst_dir, exist_ok=True)
        
        if filename.endswith('.sh') and config is not None:
            copy_with_template(src_path, dst_path, config)
        else:
            shutil.copy2(src_path, dst_path)
        
        copied.append(dst_filename)
    
    def execute(self, args) -> None:
        """
        Execute the workflow command - copy required scripts to target directory.
        
        Args:
            args: Parsed command line arguments with 'mode' and 'dest' attributes
        """
        copied = []
        missing = []
        
        os.makedirs(args.dest, exist_ok=True)
        
        config_path, config = get_config_info(args.dest)
        
        print(f" Preparing GEWUM {self.workflow_name} {args.mode} scripts...")
        print(f" Workflow repository: {self.repository_path}")
        
        if config_path:
            print(f" Using SLURM config: {config_path}")
        else:
            print(" No slurm_config.yaml found, using default values in scripts")
        
        if not os.path.exists(self.repository_path):
            print(f"[ERROR] Repository not found: {self.repository_path}")
            return
        
        common_files = self.get_common_files(args.mode)
        if common_files:
            print(f" Copying common files for mode '{args.mode}'...")
            for filename in common_files:
                self._copy_file(filename, args.dest, copied, missing, 
                               is_common=True, config=config)
        
        specific_files = self.file_map.get(args.mode, [])
        if specific_files:
            print(f" Copying {self.workflow_name}-specific files...")
            for filename in specific_files:
                self._copy_file(filename, args.dest, copied, missing, 
                               is_common=False, config=config)
        
        self._print_results(args.mode, copied, missing, args.dest)
    
    def _print_results(self, mode: str, copied: list, missing: list, dest: str) -> None:
        """Print the file copy operation results"""
        print("\n" + "=" * 60)
        print(f" GEWUM {self.workflow_name} '{mode}' scripts copy results")
        print("-" * 60)
        
        if copied:
            print(" Copied files:")
            for f in copied:
                print(f"    {f}")
        
        if missing:
            print("\n Missing files (not found in repository):")
            for f in missing:
                print(f"    {f}")
        
        print(f"\n Target location: {os.path.abspath(dest)}")
        print("=" * 60)


def setup_base_args(parser, file_map: dict, workflow_desc: str = "") -> None:
    """
    Setup common command line arguments for workflow commands.
    
    """
    mode_choices = list(file_map.keys())
    
    parser.add_argument(
        "--mode",
        required=True,
        choices=mode_choices,
        help=f"GEWUM mode to copy:\n{workflow_desc}"
    )
    parser.add_argument(
        "--dest",
        default=".",
        help="Target directory (default: current directory)"
    )
