import os

def write_xyz_if_not_exists(filename, atoms):
    if os.path.exists(filename):
        print(f"File '{filename}' already exists. Skipping write to avoid overwriting.")
    else:
        with open(filename, 'w') as f:
            f.write(f"{len(atoms)}\n")
            f.write("Optimized geometry\n")
            for symbol, x, y, z in atoms:
                f.write(f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n")
        print(f"Saved new file: {filename}")
