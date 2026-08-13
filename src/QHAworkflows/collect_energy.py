import os
import numpy as np

def main():
    import argparse
    parser = argparse.ArgumentParser(description='QHA Collect Energy Data')
    parser.add_argument('--output', '-o', default='CP/v-e.dat', help='Output file for V-E data')
    parser.add_argument('--scales', type=str, default='0.95,0.96,0.97,0.98,0.99,1.00,1.01,1.02,1.03,1.04,1.05',
                        help='Comma-separated scale factors')
    args = parser.parse_args()
    
    scales = [float(s) for s in args.scales.split(',')]
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    
    with open(args.output, "w") as outfile:
        for i, scale in enumerate(scales):
            dir_path = f"{scale:.2f}"
            try:
                with open(os.path.join(dir_path, "energy.dat"), "r") as f:
                    data = f.readline().split()
                    volume = float(data[0])
                    energy = float(data[1])
                    outfile.write(f"{volume} {energy}\n")
                    print(f"Collected data from {dir_path}: V={volume}, E={energy}")
            except:
                print(f"Warning: Missing data in {dir_path}")


if __name__ == "__main__":
    main()
