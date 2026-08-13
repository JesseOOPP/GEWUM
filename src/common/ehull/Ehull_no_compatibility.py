"""GEWUM Energy Hull Calculation Module (without compatibility corrections)
Thin wrapper that delegates to the unified implementation in Ehull_compatibility.
Supports both online (MP API) and offline (local JSON) modes
"""
import logging
import argparse

logging.basicConfig(
    filename='process.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def get_above_hull_and_formation_energy(file_path, mp_api_key=None, mp_data_path=None):
    """
    Calculate energy above hull without compatibility corrections.
    Delegates to the unified implementation in Ehull_compatibility.
    
    Args:
        file_path: Path to CSV file with formula and energy columns
        mp_api_key: Materials Project API key (for online mode)
        mp_data_path: Path to offline MP JSON file (for offline mode)
        
    Returns:
        Tuple of (phase_diagram_objects, formula_energy_dict, result_dataframe)
    """
    from .Ehull_compatibility import get_above_hull_and_formation_energy as _unified_calc
    result_df = _unified_calc(file_path, mp_api_key=mp_api_key, mp_data_path=mp_data_path, use_compatibility=False)
    formulas = result_df.iloc[:, 0].tolist()
    total_energies = result_df.iloc[:, 2].tolist() if result_df.shape[1] > 2 else [0]*len(formulas)
    formula_energy_dict = dict(zip(formulas, total_energies))
    return [None]*len(formulas), formula_energy_dict, result_df


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description='Calculate energy above hull (no compatibility)')
    parser.add_argument('--input', '-i', default='0_final_result_tot.txt',
                        help='Input CSV file path')
    parser.add_argument('--output', '-o', default='Hull_result.csv',
                        help='Output CSV file path')
    parser.add_argument('--api-key',
                        help='Materials Project API key (for online mode)')
    parser.add_argument('--mp-data',
                        help='Path to offline MP JSON file (for offline mode)')
    parser.add_argument('-cor', action='store_true', default=False,
                        help='Use MP2020 energy compatibility corrections')
    args = parser.parse_args()
    
    if not args.api_key and not args.mp_data:
        parser.error("Either --api-key or --mp-data must be provided")

    from .Ehull_compatibility import get_above_hull_and_formation_energy as _unified_calc
    result_data = _unified_calc(
        file_path=args.input,
        mp_api_key=args.api_key,
        mp_data_path=args.mp_data,
        use_compatibility=args.cor
    )

    result_data.to_csv(args.output, index=False)
    logging.info("Processing completed successfully")
    print(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()
