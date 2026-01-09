"""
Utilities for basis set superposition error (BSSE) counterpoise (CP) correction.
Implements Boys-Bernardi counterpoise correction for N-fragment systems N>=2.
"""
import os

def read_xyz_file(filepath):
    """Read an XYZ file and return the coordinate lines.
    
    Args:
        filepath (str): Path to XYZ file
        
    Returns:
        list: List of coordinate lines (without natoms/comment header)
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
        # Skip first two lines (natoms and comment)
        return [line.strip() for line in lines[2:] if line.strip()]


def get_ghost_atoms_for_fragment(fragment_dirs, fragment_index, base_dir="."):
    """Get ghost atoms for a specific fragment in counterpoise correction.
    
    For fragment i, returns all atoms from fragments j≠i as ghost atoms.
    
    Args:
        fragment_dirs (list): List of fragment directory names (e.g., ['A', 'B', 'C'])
        fragment_index (int): Index of the fragment to get ghosts for (0-based)
        base_dir (str): Base directory containing fragment directories
        
    Returns:
        str: XYZ-formatted string of ghost atoms from all other fragments
    """
    ghost_coords = []
    
    for idx, frag_dir in enumerate(fragment_dirs):
        if idx != fragment_index:
            # This is a ghost fragment
            xyz_path = os.path.join(base_dir, frag_dir, "m.xyz")
            if os.path.exists(xyz_path):
                coords = read_xyz_file(xyz_path)
                ghost_coords.extend(coords)
            else:
                # Try system.xyz as fallback
                xyz_path = os.path.join(base_dir, frag_dir, "system.xyz")
                if os.path.exists(xyz_path):
                    coords = read_xyz_file(xyz_path)
                    ghost_coords.extend(coords)
    
    return '\n'.join(ghost_coords)


def get_all_ghost_atoms(fragment_dirs, exclude_index=None, base_dir="."):
    """Get all atoms from specified fragments as ghosts.
    
    Args:
        fragment_dirs (list): List of fragment directory names
        exclude_index (int, optional): Index to exclude (for fragment calculations)
        base_dir (str): Base directory containing fragment directories
        
    Returns:
        str: XYZ-formatted string of ghost atoms
    """
    ghost_coords = []
    
    for idx, frag_dir in enumerate(fragment_dirs):
        if exclude_index is not None and idx == exclude_index:
            continue
            
        xyz_path = os.path.join(base_dir, frag_dir, "m.xyz")
        if not os.path.exists(xyz_path):
            xyz_path = os.path.join(base_dir, frag_dir, "system.xyz")
            
        if os.path.exists(xyz_path):
            coords = read_xyz_file(xyz_path)
            ghost_coords.extend(coords)
    
    return '\n'.join(ghost_coords)


def calculate_bsse_correction(energies_isolated, energies_cp):
    """Calculate BSSE correction from isolated and CP-corrected energies.
    
    Args:
        energies_isolated (list): Fragment energies in their own basis [E(A), E(B), ...]
        energies_cp (list): Fragment energies in full basis [E(A@full), E(B@full), ...]
        
    Returns:
        float: Total BSSE correction (sum over all fragments)
    """
    bsse = sum([e_cp - e_iso for e_iso, e_cp in zip(energies_isolated, energies_cp)])
    return bsse


def cp_corrected_interaction_energy(e_complex, energies_cp):
    """Calculate counterpoise-corrected interaction energy.
    
    E_int^{CP} = E(COMPLEX) - sum_{FRAGMENTS = 1}^{N} E(FRAGMENT_i @ COMPLEX)
    
    Args:
        e_complex (float): Energy of full complex
        energies_cp (list): Fragment energies in full basis
        
    Returns:
        float: CP-corrected interaction energy
    """
    return e_complex - sum(energies_cp)
