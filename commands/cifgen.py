import sys
import os
from ..config import CIFGEN_MODULE

def setup_args(parser):
    parser.add_argument("--mode", 
                        choices=["all", "oxidation", "substitute", "mutate", "dp"], 
                        required=True,
                        help="Execution mode:\n"
                             "  all: Generate all compositions\n"
                             "  oxidation: Generate compositions with common oxidation states\n"
                             "  substitute: Generate structures with element substitutions\n"
                             "  mutate: Generate perturbed structures from a template CIF\n"
                             "  dp: Generate doped structures with random atom replacements")

def execute(args, remaining_args=None):
    print(f"GEWUM CIF Generation - Mode: {args.mode}")
    
    if not os.path.exists(CIFGEN_MODULE):
        print(f"[ERROR] CIFGEN module directory not found: {CIFGEN_MODULE}")
        sys.exit(1)
    
    sys.path.insert(0, os.path.dirname(CIFGEN_MODULE))
    
    try:
        print("Starting CIF generation workflow...")
        
        if args.mode == "all":
            from cifgen_input.chemi import main as chemi_main
            from cifgen_input.conv_all_cifgen import main as conv_all_main
            chemi_main()
            conv_all_main()
        elif args.mode == "oxidation":
            from cifgen_input.chemi import main as chemi_main
            from cifgen_input.oxidation import main as oxidation_main
            from cifgen_input.conv_oxidation_cifgen import main as conv_oxidation_main
            chemi_main()
            oxidation_main()
            conv_oxidation_main()
        elif args.mode == "substitute":
            from cifgen_input.substitute import main as substitute_main
            substitute_main()
        elif args.mode == "mutate":
            from cifgen_input.mutate import main as mutate_main
            mutate_main()
        elif args.mode == "dp":
            from cifgen_input.doping import main as doping_main
            doping_main()
        else:
            print(f"[ERROR] Unknown mode: {args.mode}")
            sys.exit(1)
            
    except ImportError as e:
        print(f"[ERROR] Failed to import cifgen modules: {e}")
        print(f"[DEBUG] Python path: {sys.path}")
        print(f"[DEBUG] Files in cifgen_input: {os.listdir(CIFGEN_MODULE) if os.path.exists(CIFGEN_MODULE) else 'Directory not found'}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] CIF generation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
