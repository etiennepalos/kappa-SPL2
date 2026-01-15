"""
Calculating the interaction energy for all MPAC functionals at once
"""

#import pyscf, numpy and numba
import numpy as np
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from all_codes.numba_all import MP2_energy_split, MP2_energy_kappa_p_SS, MP2_energy_kappa_p_OS
from all_codes.numba_all import DF_MP2_energy_split, DF_MP2_energy_kappa_p_SS, DF_MP2_energy_kappa_p_OS
from all_codes.mol_all import run_pyscf
from all_codes.mpac_fun_all import MPAC_functionals
from all_codes.constants_all import *
from kappa_codes.cp_utils import get_ghost_atoms_for_fragment
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--charge", type=int, default=0, help="the charge of the system")
    parser.add_argument("--spin", type=int, default=0, help="the spin of the system")
    parser.add_argument("--basis", type=str, default="aug-cc-pvqz",help="the basisset used in the calculations")
    parser.add_argument("--cp", action="store_true", help="perform counterpoise (Boys-Bernardi) correction for BSSE")
    parser.add_argument("--use-df", action="store_true", help="use density fitting (DF-MP2) for faster integral evaluation")
    parser.add_argument("--auxbasis", type=str, default=None, help="auxiliary basis set for DF-MP2 (auto-selected if not specified)")

args = parser.parse_args()
use_df = args.use_df
auxbasis = args.auxbasis
mols=["A","B","AB"] #A and B are fragments, AB is the complex
mpacf=["spl2","f1","f1ab","mp2"]

Ex=[]
ehf=[]
Uh=[]
rho_4_3=[]
gea_4_3=[]
rho_3_2=[]
gea_7_6=[]
E_c_int=[]
E_c_SS_k=[]
E_c_SS=[]
E_c_OS=[]
E_c_OS_k=[]
E_c_mp2tot=[]

# For counterpoise correction
if args.cp:
    Ex_cp=[]
    ehf_cp=[]
    E_c_mp2tot_cp=[]
    print("Counterpoise correction ENABLED")
else:
    print("Counterpoise correction DISABLED")

for i in range(3): #run over the fragements and complex
    run_mol=mols[i]
    #add here path to frag m.xyz file
    chkfile="chkfile_"+run_mol+".chk"
    old_pwd=os.getcwd()
    datadir=old_pwd+"/"+run_mol
    os.chdir(datadir)
    ###runs HF
    py_run=run_pyscf(atom="m.xyz",charge=args.charge,spin=args.spin,basis=args.basis)
    if use_df:
        tab, B_ia = py_run.run_eris_df(chkfile_name=chkfile, chkfile_dir=datadir, auxbasis=auxbasis)
        eris = (py_run.nocc, py_run.nvirt, py_run.e, B_ia)
    else:
        tab, eris = py_run.run_eris(chkfile_name=chkfile,chkfile_dir=datadir)
    #this prints and extracts all of the ingredients except MP2
    np.savetxt("tab.csv", tab, delimiter=",", fmt='%s')
    ehf.append(tab[0])
    Uh.append(tab[1])
    Ex.append(tab[2])
    rho_4_3.append(tab[3])
    gea_4_3.append(tab[4])
    rho_3_2.append(tab[5])
    gea_7_6.append(tab[6])
    #k1 is for same spin
    #k2 is for the opposite spin
    k1ss = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
    k2ss = k1ss
    np.savetxt("k1.csv", k1ss, delimiter=",", fmt='%s')
    np.savetxt("k2.csv", k2ss, delimiter=",", fmt='%s')
    k1s=np.array(k1ss,dtype=float)
    k2s=np.array(k2ss,dtype=float)
    if use_df:
        mp2OS = DF_MP2_energy_kappa_p_OS(*eris,k2s, 1) #calculate the opposite spin integral with kappa
        E_c_OS_k.append(mp2OS)
        np.savetxt("os.csv", mp2OS, delimiter=",", fmt='%s')
        #Spin scaled \kappa's
        mp2SS = DF_MP2_energy_kappa_p_SS(*eris,k1s, 1) #calculate the same spin integral with kappa
        E_c_SS_k.append(mp2SS)
        np.savetxt("ss.csv", mp2SS, delimiter=",", fmt='%s')
        e_mp2_split = DF_MP2_energy_split(*eris,) #calculates the same an opposite mp2 integrals
        np.savetxt("mp2.csv", e_mp2_split, delimiter=",", fmt='%s')
    else:
        mp2OS = MP2_energy_kappa_p_OS(*eris,k2s, 1) #calculate the opposite spin integral with kappa
        E_c_OS_k.append(mp2OS)
        np.savetxt("os.csv", mp2OS, delimiter=",", fmt='%s')
        #Spin scaled \kappa's
        mp2SS = MP2_energy_kappa_p_SS(*eris,k1s, 1) #calculate the same spin integral with kappa
        E_c_SS_k.append(mp2SS)
        np.savetxt("ss.csv", mp2SS, delimiter=",", fmt='%s')
        e_mp2_split = MP2_energy_split(*eris,) #calculates the same an opposite mp2 integrals
        np.savetxt("mp2.csv", e_mp2_split, delimiter=",", fmt='%s')
    E_c_SS.append(e_mp2_split[0])
    E_c_OS.append(e_mp2_split[1])
    E_c_mp2tot.append(sum(e_mp2_split))
    os.chdir(old_pwd)

# Counterpoise correction: calculate fragments in full basis
if args.cp:
    print("\n=== Running Counterpoise Correction ===")
    fragment_dirs = ["A", "B"]  # Exclude complex
    
    for i in range(2):  # Only fragments, not complex
        run_mol = mols[i]
        print(f"Calculating fragment {run_mol} in full (ghost) basis...")
        
        chkfile = "chkfile_" + run_mol + "_cp.chk"
        old_pwd = os.getcwd()
        datadir = old_pwd + "/" + run_mol
        os.chdir(datadir)
        
        # Get ghost atoms from all other fragments
        ghost_atoms_str = get_ghost_atoms_for_fragment(fragment_dirs, i, base_dir=old_pwd)
        
        # Run calculation with ghost atoms
        py_run = run_pyscf(atom="m.xyz", charge=args.charge, spin=args.spin, 
                          basis=args.basis, ghost_atoms=ghost_atoms_str)
        if use_df:
            tab_cp, B_ia_cp = py_run.run_eris_df(chkfile_name=chkfile, chkfile_dir=datadir, auxbasis=auxbasis)
            eris_cp = (py_run.nocc, py_run.nvirt, py_run.e, B_ia_cp)
        else:
            tab_cp, eris_cp = py_run.run_eris(chkfile_name=chkfile, chkfile_dir=datadir)
        
        np.savetxt("tab_cp.csv", tab_cp, delimiter=",", fmt='%s')
        ehf_cp.append(tab_cp[0])
        Ex_cp.append(tab_cp[2])
        
        # Calculate MP2 for CP
        if use_df:
            e_mp2_split_cp = DF_MP2_energy_split(*eris_cp)
        else:
            e_mp2_split_cp = MP2_energy_split(*eris_cp)
        np.savetxt("mp2_cp.csv", e_mp2_split_cp, delimiter=",", fmt='%s')
        E_c_mp2tot_cp.append(sum(e_mp2_split_cp))
        
        os.chdir(old_pwd)
    
    print("=== Counterpoise correction calculations complete ===\n")

#getting the arrays into the correct shape
E_c_SS_k=np.array(E_c_SS_k).T
E_c_OS_k=np.array(E_c_OS_k).T
E_c_OS=np.array(E_c_OS)

#initializing the MPAC functionals and calculating the HF energy difference
form_frags=MPAC_functionals(Ex[0]+Ex[1],rho_4_3[0]+rho_4_3[1],gea_4_3[0]+gea_4_3[1])
form_com=MPAC_functionals(Ex[2],rho_4_3[2],gea_4_3[2]) 
ehfdiv=ehf[2]-ehf[1]-ehf[0]
funcs=["MP2","SPL2","F1","F1ab","k-MP2","k-SPL2","k-F1","k-F1ab","ksskos-MP2","ksskos-SPL2","ksskos-F1","ksskos-F1ab","coskos-MP2","coskos-SPL2","coskos-F1","coskos-F1ab","cos-MP2","cos-SPL2","cos-F1","cos-F1ab"]

#storing all the EMP2 data
EMP2vals={
    "MP2": [E_c_mp2tot]*4,
    "k-MP2": [E_c_SS_k[5] + E_c_OS_k[5],E_c_SS_k[11] + E_c_OS_k[11],E_c_SS_k[7] + E_c_OS_k[7],E_c_SS_k[9] +E_c_OS_k[9]],
    "ksskos-MP2": [E_c_SS_k[3] + E_c_OS_k[8],E_c_SS_k[5] + E_c_OS_k[11],E_c_SS_k[4] + E_c_OS_k[8],E_c_SS_k[10] + E_c_OS_k[7]],
    "coskos-MP2":[2.1*E_c_OS_k[3],2.1*E_c_OS_k[7],2.3*E_c_OS_k[5],2.5*E_c_OS_k[4]],
    "cos-MP2": [1.7*E_c_OS,1.8*E_c_OS,2.2*E_c_OS,2*E_c_OS]
}

#calculating all the 20 functionals
for name,emp2 in EMP2vals.items():
    E_c_int.append(form_com.mp2(params[name][0],emp2[0][2])-form_frags.mp2(params[name][0],emp2[0][1]+emp2[0][0]))
    E_c_int.append(form_com.spl2(params[name][1],emp2[1][2])-form_frags.spl2(params[name][1],emp2[1][1]+emp2[1][0]))
    E_c_int.append(form_com.f1(params[name][2],emp2[2][2])-form_frags.f1(params[name][2],emp2[2][1]+emp2[2][0]))
    E_c_int.append(form_com.f1(params[name][3],emp2[3][2])-form_frags.f1(params[name][3],emp2[3][1]+emp2[3][0]))

#print json file
E_c_ints=dict(zip(funcs,kcal*(ehfdiv+np.array(E_c_int))))
print(E_c_ints) #prints out the correct E_c_int
with open("E_c_all.json","w",encoding="utf-8") as f:
    json.dump(E_c_ints,f)

# Calculate and print CP-corrected interaction energies
if args.cp:
    print("\n=== Counterpoise-Corrected Interaction Energies ===")
    
    # CP-corrected uses fragments calculated in full basis
    form_frags_cp = MPAC_functionals(Ex_cp[0]+Ex_cp[1], rho_4_3[0]+rho_4_3[1], gea_4_3[0]+gea_4_3[1])
    ehfdiv_cp = ehf[2] - ehf_cp[1] - ehf_cp[0]  # Complex - fragments@full_basis
    
    # Calculate CP-corrected energies for all functionals
    E_c_int_cp = []
    EMP2vals_cp = {
        "MP2": [E_c_mp2tot_cp]*4,
        "k-MP2": [E_c_SS_k[5] + E_c_OS_k[5], E_c_SS_k[11] + E_c_OS_k[11], 
                  E_c_SS_k[7] + E_c_OS_k[7], E_c_SS_k[9] + E_c_OS_k[9]],
        "ksskos-MP2": [E_c_SS_k[3] + E_c_OS_k[8], E_c_SS_k[5] + E_c_OS_k[11], 
                       E_c_SS_k[4] + E_c_OS_k[8], E_c_SS_k[10] + E_c_OS_k[7]],
        "coskos-MP2": [2.1*E_c_OS_k[3], 2.1*E_c_OS_k[7], 2.3*E_c_OS_k[5], 2.5*E_c_OS_k[4]],
        "cos-MP2": [1.7*E_c_OS, 1.8*E_c_OS, 2.2*E_c_OS, 2*E_c_OS]
    }
    
    for name, emp2_cp in EMP2vals_cp.items():
        E_c_int_cp.append(
            form_com.mp2(params[name][0], emp2_cp[0][2])
            - form_frags_cp.mp2(params[name][0], E_c_mp2tot_cp[0] + E_c_mp2tot_cp[1])
        )
        E_c_int_cp.append(
            form_com.spl2(params[name][1], emp2_cp[1][2])
            - form_frags_cp.spl2(params[name][1], E_c_mp2tot_cp[0] + E_c_mp2tot_cp[1])
        )
        E_c_int_cp.append(
            form_com.f1(params[name][2], emp2_cp[2][2])
            - form_frags_cp.f1(params[name][2], E_c_mp2tot_cp[0] + E_c_mp2tot_cp[1])
        )
        E_c_int_cp.append(
            form_com.f1(params[name][3], emp2_cp[3][2])
            - form_frags_cp.f1(params[name][3], E_c_mp2tot_cp[0] + E_c_mp2tot_cp[1])
        )
    
    E_c_ints_cp = dict(zip(funcs, kcal*(ehfdiv_cp+np.array(E_c_int_cp))))
    print(E_c_ints_cp)
    
    with open("E_c_all_CP.json","w",encoding="utf-8") as f:
        json.dump(E_c_ints_cp, f)
    
    # Calculate and save BSSE corrections
    bsse_corrections = {func: E_c_ints[func] - E_c_ints_cp[func] for func in funcs}
    print("\n=== BSSE Corrections [kcal/mol] ===")
    print(bsse_corrections)
    
    with open("BSSE_corrections.json","w",encoding="utf-8") as f:
        json.dump(bsse_corrections, f)

