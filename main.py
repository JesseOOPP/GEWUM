import argparse
import importlib
import subprocess
import sys
import os
import colorsys

def _rgb_to_ansi(r: int, g: int, b: int) -> str:
    """Convert RGB to ANSI 24-bit color escape sequence"""
    return f"\033[38;2;{r};{g};{b}m"

def _generate_gradient_colors(text: str, start_hue: float = 0.6, end_hue: float = 0.0) -> str:
    """Generate rainbow gradient colored text.
    
    Args:
        text: The text to colorize
        start_hue: Starting hue (0.0-1.0), default 0.6 (blue)
        end_hue: Ending hue (0.0-1.0), default 0.0 (red)
    """
    result = []
    chars = [c for c in text if c != '\n']
    total_chars = len(chars)
    
    char_idx = 0
    for line in text.split('\n'):
        colored_line = []
        for c in line:
            if c == ' ':
                colored_line.append(c)
            else:
                t = char_idx / max(total_chars - 1, 1)
                hue = start_hue + (end_hue - start_hue) * t
                if hue < 0:
                    hue += 1.0
                elif hue > 1:
                    hue -= 1.0
                
                r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 0.95)
                ansi_color = _rgb_to_ansi(int(r * 255), int(g * 255), int(b * 255))
                colored_line.append(f"{ansi_color}{c}")
                char_idx += 1
        result.append(''.join(colored_line))
    
    return '\n'.join(result) + '\033[0m'

BORDER_WIDTH = 60

HEADER_RAW = """
 ===========================================================
 =                                                         =
 =      ██████╗ ███████╗██╗    ██╗██╗   ██╗███╗   ███╗     =
 =     ██╔════╝ ██╔════╝██║    ██║██║   ██║████╗ ████║     =
 =     ██║  ███╗█████╗  ██║ █╗ ██║██║   ██║██╔████╔██║     =
 =     ██║   ██║██╔══╝  ██║███╗██║██║   ██║██║╚██╔╝██║     =
 =     ╚██████╔╝███████╗╚███╔███╔╝╚██████╔╝██║ ╚═╝ ██║     =
 =      ╚═════╝ ╚══════╝ ╚══╝╚══╝  ╚═════╝ ╚═╝     ╚═╝     =
 =                                                         =
 ===========================================================
 =        General Exploration Workflow for the             =
 =                 Utopia of Materials                     =
 ===========================================================
 =            Citation: arXiv:2604.21401                   =
 ===========================================================
"""

DEVELOPER = "Developer: Jiexi Song"
CONTACT = "Contact: songjx@szlab.ac.cn"
VERSION = "Version 1.0.0"

def print_header():
    colored_header = _generate_gradient_colors(HEADER_RAW, start_hue=0.55, end_hue=0.05)
    print(colored_header)
    print(DEVELOPER)
    print(CONTACT)
    print(VERSION)
    print("=" * BORDER_WIDTH)
    print()

def main():
    print_header()
    
    parser = argparse.ArgumentParser(
        prog="gewum",
        description="GEWUM: General Exploration Workflow for the Utopia of Materials",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="="*BORDER_WIDTH + "\n"
               "Maintainer & Lead Developer: Jiexi Song\n"
               "Contributors: Diwei Shi, Aixian She, Yanqing Qin, Zhenyu Liu\n"
               "Acknowledgments: Fengyuan Xuan and Chongde Cao\n"
               + "="*BORDER_WIDTH

    )
    
    subparsers = parser.add_subparsers(
        title="Available Commands",
        dest="command",
        metavar="<command>",
        help="Use 'gewum <command> -h' for detailed help"
    )
    
    command_mapping = {
        'cifgen': 'cifgen',
        'RD': 'RD',  
        'PT': 'PT',
        'ELA': 'ELA',
        'QHA': 'QHA',
        'TC': 'TC',
        'MD': 'MD',
        'DB': 'DB',
    }

    command_descriptions = {
        'cifgen': 'Generate INPUT for CIF Generation',
        'RD': 'Selective Random Structure Search Workflow',
        'PT': 'Perturbation Structure Search Workflow',
        'ELA': 'Elastic Constants',
        'QHA': 'Quasi-Harmonic Approximation',
        'TC': 'Thermal Conductivity',
        'MD': 'Molecular Dynamics (test)',
        'DB': 'Read-only DB inspection and batch CIF export',
    }
    
    command_modules = {}
    
    for display_name, module_name in command_mapping.items():
        try:
            full_module_name = f"gewum.commands.{module_name}"
            module = importlib.import_module(full_module_name)
            
            cmd_desc = command_descriptions.get(display_name, f"GEWUM {display_name} module")
            cmd_parser = subparsers.add_parser(
                display_name,
                help=cmd_desc,
                description=f"GEWUM {display_name} - {cmd_desc}",
                formatter_class=argparse.RawTextHelpFormatter
            )
            module.setup_args(cmd_parser)
            command_modules[display_name] = module
        except ImportError as e:
            print(f"[WARNING] GEWUM module '{display_name}' not found: {e}")
            continue
    
    parser.add_argument('--version', action='version', version=f'GEWUM {VERSION}')
    
    # Pre-scan: intercept viz mode help request before argparse captures -h
    if len(sys.argv) >= 2 and sys.argv[1] == 'RD':
        argv_str = sys.argv[2:]  # args after 'RD'
        if '--mode' in argv_str:
            mode_idx = argv_str.index('--mode')
            if mode_idx + 1 < len(argv_str):
                mode_name = argv_str[mode_idx + 1]
                if mode_name in ('viz', 'viz2') and ('-h' in argv_str or '--help' in argv_str):
                    module_map = {
                        'viz': 'gewum.src.common.postprocess.visualization',
                        'viz2': 'gewum.src.common.postprocess.viz2_analysis',
                    }
                    result = subprocess.run(
                        [sys.executable, '-m', module_map[mode_name], '--help']
                    )
                    sys.exit(result.returncode)

    args, remaining = parser.parse_known_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    try:
        print(f"Starting GEWUM {args.command} module...")
        command_modules[args.command].execute(args, remaining_args=remaining)
        print(f"\nGEWUM {args.command} completed successfully!")
        print("=" * BORDER_WIDTH)
    except KeyError:
        print(f"[ERROR] Unknown GEWUM command: {args.command}")
        print("\nAvailable commands:")
        for cmd in command_modules.keys():
            print(f"  - {cmd}")
        sys.exit(1)

if __name__ == "__main__":
    main()
