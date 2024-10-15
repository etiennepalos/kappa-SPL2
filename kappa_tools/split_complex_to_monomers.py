'''
XYZ SPLITTER
A simple "helper script" for MPAC interaction energy calculations
for clusters of arbitrary size (N-body).

This script reads an XYZ file containing an N-monomer cluster,
splits it into individual monomers based on a provided list of atom counts,
and organizes each monomer into its own directory (e.g., A, B, C, ...).
Additionally, it creates a directory called 'COMPLEX' and copies the original
cluster XYZ file into it as 'complex.xyz'.

i.e. it follows the format specified in kappa-SPL v1.0

USAGE:
    python3 xyz_monomer_splitter.py cluster.xyz
    
OUTPUT:
    - A 'COMPLEX' directory containing the original cluster XYZ file as 'complex.xyz'.
    - A set of directories named A, B, C, ..., where each directory contains the XYZ 
      file for the corresponding monomer (e.g., A_monomer.xyz, B_monomer.xyz, etc.).

NOTE:
    - The script assumes that the atom count for each monomer is provided in the 'atlist' variable.
    - The monomer labels (e.g., A, B, C) are defined in the 'monomer_labels' variable.

Author: Etienne Palos
'''

# Required libs 
import sys
import os
import shutil

if len(sys.argv) != 2:
    print("Usage: python3 "  + sys.argv[0] + " cluster.xyz")
    sys.exit(1)

fxyz = sys.argv[1]
atlist = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]  # N-atoms per monomer (example given is e.g. water decamer)
monomer_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

# For N-monomer complex
def create_complex_directory(xyz):
    os.makedirs("COMPLEX", exist_ok=True)  
    shutil.copy(xyz, "COMPLEX/system.xyz")  

# Function to write monomers into directories
def write_monomers(xyz, atlist, labels):
    with open(xyz, 'r') as f:
        nat = f.readline().split()[0]  
        f.readline()  
        
        mons = []
        for i in range(len(atlist)):
            m = []
            for j in range(atlist[i]):
                line = f.readline()
                m.append(line)
            mons.append(m)

    # Create directories for monomers 
    for i in range(len(atlist)):
        dir_name = labels[i]  
        os.makedirs(dir_name, exist_ok=True) 
        
        fname = f"{dir_name}/system.xyz" 
        with open(fname, 'w') as ff:
            inat = atlist[i]
            ff.write(str(inat) + "\n")
            ff.write(f"{dir_name} monomer\n")
            for l in mons[i]:
                ff.write(l)

create_complex_directory(fxyz)
write_monomers(fxyz, atlist, monomer_labels)

print("\nXYZ preparation complete. The full cluster is saved in the 'COMPLEX' directory as 'system.xyz'.")
print("Monomer files are organized in respective directories.\n")

