#!/usr/bin/env python3
"""
MBE SPLITTER
A simple "helper script" for setting up many-body expnansion energy calculations
using MPAC or any other method.

This script reads an XYZ file containing an N-monomer cluster, and generates all 
possible n-mers from n=1 to n=N-1, skipping N which is the full 'COMPLEX'.
Similar to XYZ SPLITTER, monomers are labelled A, B, ..., 
such that we now get n-mers/A_B_C...

i.e. it follows the format specified in kappa-SPL v1.0

USAGE:
    python3 xyz_mbe_splitter.py cluster.xyz
    python3 xyz_mbe_splitter.py cluster.xyz --verbose
    --verbose option prints "tree-like" directory organization

OUTPUT:
    - A 'COMPLEX' directory containing the original cluster XYZ file as 'complex.xyz'.
    - Directories 1-mer, 2-mer, ..., (N-1)-mer containing all possible fragment A_B_... combinations for all n

NOTE:
    - The script assumes that the atom count for each monomer is provided in the 'atlist' variable.
    - The monomer labels (e.g., A, B, C) are defined in the 'monomer_labels' variable.

Author: @etiennepalos
"""
import sys, os, shutil, argparse
from itertools import combinations

atlist = [1,3,3,3,3,3,3,3,3,3]
monomer_labels = ['A','B','C','D','E','F','G','H','I','J']
monomer_charges = [1] + [0]*(len(atlist)-1)   # assumes first monomer charged, modify at will

def read_monomers(xyz_file, atlist):
    with open(xyz_file,'r') as f:
        f.readline()   # nat
        f.readline()   # comment

        monomers = []
        for nat in atlist:
            block = [f.readline() for _ in range(nat)]
            monomers.append(block)
    return monomers


def write_fragment(indices, monomers, labels, outdir, charges):
    os.makedirs(outdir, exist_ok=True)

    frag_name = "_".join(labels[i] for i in indices)
    frag_dir = os.path.join(outdir, frag_name)
    os.makedirs(frag_dir, exist_ok=True)

    # xyz
    xyz_path = os.path.join(frag_dir, "m.xyz")
    with open(xyz_path, 'w') as f:
        total_atoms = sum(len(monomers[i]) for i in indices)
        f.write(f"{total_atoms}\n")
        f.write(f"{frag_name}\n")
        for i in indices:
            for line in monomers[i]:
                f.write(line)

    # charge
    frag_charge = sum(charges[i] for i in indices)
    with open(os.path.join(frag_dir, "charge.txt"), 'w') as f:
        f.write(str(frag_charge))

def make_complex(cluster_xyz, charges):
    os.makedirs("COMPLEX", exist_ok=True)
    shutil.copy(cluster_xyz, "COMPLEX/m.xyz")
    with open("COMPLEX/charge.txt",'w') as f:
        f.write(str(sum(charges)))

def print_tree(start_path):
    for root, dirs, files in os.walk(start_path):
        level = root.replace(start_path, '').count(os.sep)
        indent = "    " * level
        dirname = os.path.basename(root) if level > 0 else root
        print(f"{indent}\033[96m{dirname}/\033[0m")
        subindent = "    " * (level + 1)
        for f in files:
            print(f"{subindent}\033[90m{f}\033[0m")

def print_unicode_combos(r, combos, labels):
    print(f"\n{r}-mer:")
    for i, combo in enumerate(combos):
        parts = " + ".join(labels[j] for j in combo)
        branch = "└─" if i == len(combos)-1 else "├─"
        print(f"  {branch} {parts}")
    print("")


# main
def main():
    parser = argparse.ArgumentParser(description="Generate MBE fragments from a cluster XYZ.")
    parser.add_argument("xyz", type=str, help="Cluster xyz file")
    parser.add_argument("--verbose", type=int, default=0,
                        help="0 = silent, 1 = print unicode combos, 2 = print directory tree, 3 = both")
    args = parser.parse_args()

    xyz_file = args.xyz
    monomers = read_monomers(xyz_file, atlist)
    N = len(monomers)

    make_complex(xyz_file, monomer_charges)

    for r in range(1, N+1):
        r_dir = f"{r}-mer"
        os.makedirs(r_dir, exist_ok=True)

        combos = list(combinations(range(N), r))

        # verbose1: unicode combo print
        if args.verbose in (1, 3):
            print_unicode_combos(r, combos, monomer_labels)

        for combo in combos:
            write_fragment(combo, monomers, monomer_labels, r_dir, monomer_charges)

    # verbose2: filesystem tree
    if args.verbose in (2, 3):
        print("\nDirectory structure generated:\n")
        print_tree(".")
        print("")

if __name__ == "__main__":
    main()

