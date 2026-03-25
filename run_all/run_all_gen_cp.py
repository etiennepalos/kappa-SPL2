"""
Generalized MPAC functionals Interaction Energies
for N-molecule clusters
E_int(N) = E_N - sum_i^N E_i  (e.g., dimer, trimer, ..., decamer)
supports different charges for different fragments

Contributors:
Etienne Palos
K.J. Daas
D.P. Kooi
S. Vuckovic 
"""

# import required libs: pyscf, numpy and numba
import numpy as np
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from all_codes.numba_all import *
from all_codes.mol_all import run_pyscf
from all_codes.mpac_fun_all import MPAC_functionals
from all_codes.constants_all import *
from kappa_codes.cp_utils import get_ghost_atoms_for_fragment
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--charge", type=int, default=0, help="the charge of the system") #good for neutral
    parser.add_argument("--charges", nargs='+', type=int, help="charge per fragment i, i= 1,...,N. Last arg is complex charge.")
    parser.add_argument("--spin", type=int, default=0, help="the spin of the system")
    parser.add_argument("--basis", type=str, default="def2-qzvppd",help="the basisset used in the calculations")
    parser.add_argument("--cp", action="store_true", help="perform counterpoise (Boys-Bernardi) correction for BSSE")
    parser.add_argument("--use-df", action="store_true", help="use density fitting (DF-MP2) for faster integral evaluation")
    parser.add_argument("--auxbasis", type=str, default=None, help="auxiliary basis set for DF-MP2 (auto-selected if not specified)")
    parser.add_argument("--fragments", nargs='+', type=str, default=None, help="list of fragment names (e.g., A B C D). Auto-detected if not provided.")

    args = parser.parse_args()
    
    use_df = args.use_df
    auxbasis = args.auxbasis

    # See kappa_tools for helper scripts to prepare your workind dir. 
    # NOTE: For now, order is assumed to be fragment1, fragment2, ..., fragmentN, complex.
    if args.fragments:
        mols = args.fragments + ["COMPLEX"]
    else:
        # Auto-detect single-letter/frag directories if possible
        detected = [d for d in os.listdir('.') if os.path.isdir(d) and (len(d) <= 2 or d.startswith('frag'))]
        detected.sort()
        if "COMPLEX" in os.listdir('.'):
            mols = detected + ["COMPLEX"]
        else:
            mols = ["A", "B", "C", "D", "COMPLEX"] # Fallback

    if args.charges:
        charges = args.charges  # Expecting charges for each fragment and the complex in same order as mol
    else:
        charges = [0] * len(mols)  # Default all 0

    from kappa_codes.output_utils import print_job_header, write_mpac_job_out
    start_time = print_job_header(mols, args.basis, use_df, args.cp, charges)

    if len(charges) != len(mols):
        raise ValueError(f"ERROR: Number of charges ({len(charges)}) must match the number of systems (Fragments + COMPLEX = {len(mols)}).")

mpacf=["spl2","f1","f1ab", "mpac25","mp2"]

# Initializing lists to store MPAC ingredients 
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
    rho_4_3_cp=[]
    gea_4_3_cp=[]
    E_c_SS_k_cp=[]
    E_c_OS_k_cp=[]
    E_c_OS_cp=[]
    E_c_mp2tot_cp=[]
    print("Counterpoise correction ENABLED")
else:
    print("Counterpoise correction DISABLED")

for i in range(len(mols)): #run over the fragements and complex
    run_mol=mols[i]
    # add here path to frag m.xyz file
    chkfile="chkfile_"+run_mol+".chk"
    old_pwd=os.getcwd()
    datadir=old_pwd+"/"+run_mol
    os.chdir(datadir)
    # Runs Hartree-Fock calculation
    py_run=run_pyscf(atom="m.xyz",charge=charges[i],spin=args.spin,basis=args.basis)
    if use_df:
        tab, B_ia = py_run.run_eris_df(chkfile_name=chkfile, chkfile_dir=datadir, auxbasis=auxbasis)
        eris = (py_run.nocc, py_run.nvirt, py_run.e, B_ia)
    else:
        tab, eris = py_run.run_eris(chkfile_name=chkfile,chkfile_dir=datadir)
    # print and extracts all of the ingredients except MP2 into tab.csv
    np.savetxt("tab.csv", tab, delimiter=",", fmt='%s')
    ehf.append(tab[0])
    Uh.append(tab[1])
    Ex.append(tab[2])
    rho_4_3.append(tab[3])
    gea_4_3.append(tab[4])
    rho_3_2.append(tab[5])
    gea_7_6.append(tab[6])
    k1ss = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
    k2ss = k1ss
    np.savetxt("k1.csv", k1ss, delimiter=",", fmt='%s')
    np.savetxt("k2.csv", k2ss, delimiter=",", fmt='%s')
    k1s=np.array(k1ss,dtype=float)
    k2s=np.array(k2ss,dtype=float)
    if use_df:
        mp2OS = DF_MP2_energy_kappa_p_OS_parallel(*eris,k2s, 1) #calculate the opposite spin integral with kappa
        E_c_OS_k.append(mp2OS)
        np.savetxt("os.csv", mp2OS, delimiter=",", fmt='%s')
        # spin scaled \kappa's
        mp2SS = DF_MP2_energy_kappa_p_SS_parallel(*eris,k1s, 1) #calculate the same spin integral with kappa
        E_c_SS_k.append(mp2SS)
        np.savetxt("ss.csv", mp2SS, delimiter=",", fmt='%s')
        e_mp2_split = DF_MP2_energy_split(*eris,) #calculates the same an opposite mp2 integrals
        np.savetxt("mp2.csv", e_mp2_split, delimiter=",", fmt='%s')
    else:
        mp2OS = MP2_energy_kappa_p_OS_parallel(*eris,k2s, 1) #calculate the opposite spin integral with kappa
        E_c_OS_k.append(mp2OS)
        np.savetxt("os.csv", mp2OS, delimiter=",", fmt='%s')
        # spin scaled \kappa's
        mp2SS = MP2_energy_kappa_p_SS_parallel(*eris,k1s, 1) #calculate the same spin integral with kappa
        E_c_SS_k.append(mp2SS)
        np.savetxt("ss.csv", mp2SS, delimiter=",", fmt='%s')
        e_mp2_split = MP2_energy_split(*eris,) #calculates the same an opposite mp2 integrals
        np.savetxt("mp2.csv", e_mp2_split, delimiter=",", fmt='%s')
    E_c_SS.append(e_mp2_split[0])
    E_c_OS.append(e_mp2_split[1])
    E_c_mp2tot.append(sum(e_mp2_split))
    os.chdir(old_pwd)

# counterpoise correction: calculate fragments in full basis
if args.cp:
    print("\n Counterpoise Correction for N-fragment system")
    fragment_dirs = mols[:-1]  # All fragments, excluding complex
    num_fragments = len(fragment_dirs)
    
    for i in range(num_fragments):  # Only fragments, not complex
        run_mol = mols[i]
        print(f"Calculating fragment {run_mol} in full (ghost) basis...")
        
        chkfile = "chkfile_" + run_mol + "_cp.chk"
        old_pwd = os.getcwd()
        datadir = old_pwd + "/" + run_mol
        os.chdir(datadir)
        
        # get ghost atoms from all other fragments
        ghost_atoms_str = get_ghost_atoms_for_fragment(fragment_dirs, i, base_dir=old_pwd)
        
        # run calculation with ghost atoms
        py_run = run_pyscf(atom="m.xyz", charge=charges[i], spin=args.spin, 
                          basis=args.basis, ghost_atoms=ghost_atoms_str)
        if use_df:
            tab_cp, B_ia_cp = py_run.run_eris_df(chkfile_name=chkfile, chkfile_dir=datadir, auxbasis=auxbasis)
            eris_cp = (py_run.nocc, py_run.nvirt, py_run.e, B_ia_cp)
        else:
            tab_cp, eris_cp = py_run.run_eris(chkfile_name=chkfile, chkfile_dir=datadir)
        
        np.savetxt("tab_cp.csv", tab_cp, delimiter=",", fmt='%s')
        ehf_cp.append(tab_cp[0])
        Ex_cp.append(tab_cp[2])
        rho_4_3_cp.append(tab_cp[3])
        gea_4_3_cp.append(tab_cp[4])
        
        # calculate MP2 for CP
        if use_df:
            mp2OS_cp = DF_MP2_energy_kappa_p_OS_parallel(*eris_cp,k2s, 1) #calculate the opposite spin integral with kappa
            E_c_OS_k_cp.append(mp2OS_cp)
            mp2SS_cp = DF_MP2_energy_kappa_p_SS_parallel(*eris_cp,k1s, 1) #calculate the same spin integral with kappa
            E_c_SS_k_cp.append(mp2SS_cp)
            e_mp2_split_cp = DF_MP2_energy_split(*eris_cp,) #calculates the same an opposite mp2 integrals
        else:
            mp2OS_cp = MP2_energy_kappa_p_OS_parallel(*eris_cp,k2s, 1) #calculate the opposite spin integral with kappa
            E_c_OS_k_cp.append(mp2OS_cp)
            mp2SS_cp = MP2_energy_kappa_p_SS_parallel(*eris_cp,k1s, 1) #calculate the same spin integral with kappa
            E_c_SS_k_cp.append(mp2SS_cp)
            e_mp2_split_cp = MP2_energy_split(*eris_cp,) #calculates the same an opposite mp2 integrals
            
        np.savetxt("mp2_cp.csv", e_mp2_split_cp, delimiter=",", fmt='%s')
        E_c_OS_cp.append(e_mp2_split_cp[1])
        E_c_mp2tot_cp.append(sum(e_mp2_split_cp))
        
        os.chdir(old_pwd)
    
    # append the uncorrected complex values to the CP lists so that index [-1] corresponds to the complex
    ehf_cp.append(ehf[-1])
    Ex_cp.append(Ex[-1])
    rho_4_3_cp.append(rho_4_3[-1])
    gea_4_3_cp.append(gea_4_3[-1])
    E_c_OS_k_cp.append(E_c_OS_k[-1])
    E_c_SS_k_cp.append(E_c_SS_k[-1])
    E_c_OS_cp.append(E_c_OS[-1])
    E_c_mp2tot_cp.append(E_c_mp2tot[-1])

    print("=== Counterpoise correction calculations complete ===\n")

# gets the arrays into the correct shape
E_c_SS_k=np.array(E_c_SS_k).T
E_c_OS_k=np.array(E_c_OS_k).T
E_c_OS=np.array(E_c_OS)

if args.cp:
    E_c_SS_k_cp=np.array(E_c_SS_k_cp).T
    E_c_OS_k_cp=np.array(E_c_OS_k_cp).T
    E_c_OS_cp=np.array(E_c_OS_cp)

# initializing the MPAC functionals and calculating the HF energy difference
if args.cp:
    form_frags_cp = MPAC_functionals(sum(Ex_cp[:-1]), sum(rho_4_3_cp[:-1]), sum(gea_4_3_cp[:-1]))
    form_com_cp=MPAC_functionals(Ex_cp[-1],rho_4_3_cp[-1],gea_4_3_cp[-1]) 
form_frags=MPAC_functionals(sum(Ex[:-1]),sum(rho_4_3[:-1]),sum(gea_4_3[:-1]))
form_com=MPAC_functionals(Ex[-1],rho_4_3[-1],gea_4_3[-1]) 

# HF interaction energy
ehfdiv=ehf[-1]-(sum(ehf[:-1]))
einthf_kcal=ehfdiv*kcal
print("Hartree-Fock Eint_HF [kcal/mol] =", einthf_kcal)

funcs=["MP2","SPL2","F1","F1ab","MPAC25","k-MP2","k-SPL2","k-F1","k-F1ab","k-MPAC25","ksskos-MP2","ksskos-SPL2","ksskos-F1","ksskos-F1ab","ksskos-MPAC25","coskos-MP2","coskos-SPL2","coskos-F1","coskos-F1ab","coskos-MPAC25","cos-MP2","cos-SPL2","cos-F1","cos-F1ab","cos-MPAC25"]

# storing all the EMP2 data
EMP2vals={
    "MP2": [E_c_mp2tot]*5,
    "k-MP2": [E_c_SS_k[5] + E_c_OS_k[5],E_c_SS_k[11] + E_c_OS_k[11],E_c_SS_k[7] + E_c_OS_k[7],E_c_SS_k[9] +E_c_OS_k[9],E_c_SS_k[5] + E_c_OS_k[5]],
    "ksskos-MP2": [E_c_SS_k[3] + E_c_OS_k[8],E_c_SS_k[5] + E_c_OS_k[11],E_c_SS_k[4] + E_c_OS_k[8],E_c_SS_k[10] + E_c_OS_k[7],E_c_SS_k[3] + E_c_OS_k[8]],
    "coskos-MP2":[2.1*E_c_OS_k[3],2.1*E_c_OS_k[7],2.3*E_c_OS_k[5],2.5*E_c_OS_k[4],2.1*E_c_OS_k[3]],
    "cos-MP2": [1.7*E_c_OS,1.8*E_c_OS,2.2*E_c_OS,2*E_c_OS,1.7*E_c_OS]
}

# calculation for N-fragment systems
for name, emp2 in EMP2vals.items():
    # Compute the interaction energy for N-fragment complex
    # form_com stores energy for the N-fragment complex E(N)
    # form_frags sotres the sum of energies of fragments A, B, ...
    E_c_int.append(
        form_com.mp2(params[name][0], emp2[0][-1])  # Energy for N-fragment complex 
        - form_frags.mp2(params[name][0], sum(emp2[0][:-1]))  # Energy for A+B+... 
    )
    E_c_int.append(
        form_com.spl2(params[name][1], emp2[1][-1])   
        - form_frags.spl2(params[name][1], sum(emp2[1][:-1]))
    )
    E_c_int.append(
        form_com.f1(params[name][2], emp2[2][-1])   
        - form_frags.f1(params[name][2], sum(emp2[2][:-1]))
    )
    E_c_int.append(
        form_com.f1(params[name][3], emp2[3][-1])   
        - form_frags.f1(params[name][3], sum(emp2[3][:-1])) 
    )
    E_c_int.append(
        form_com.f1(params[name][4], emp2[4][-1])   # MPAC25 uses f1 functional
        - form_frags.f1(params[name][4], sum(emp2[4][:-1])) 
    )

# print interaction energies in kcal/mol to json file
E_c_ints = dict(zip(funcs, kcal * (ehfdiv + np.array(E_c_int))))
with open("Eint_kcalmol_all.json", "w", encoding="utf-8") as f:
    json.dump(E_c_ints, f)

if args.cp:
    # CP-corrected uses fragments calculated in full basis
    form_frags_cp = MPAC_functionals(sum(Ex_cp[:-1]), sum(rho_4_3[:-1]), sum(gea_4_3[:-1]))
    form_com_cp = MPAC_functionals(Ex_cp[-1], rho_4_3[-1], gea_4_3[-1])
    ehfdiv_cp = ehf_cp[-1] - sum(ehf_cp[:-1])  # Complex - sum(fragments@full_basis)
    E_c_int_cp = []
    EMP2vals_cp = {
        "MP2": [E_c_mp2tot_cp]*5,
        "k-MP2": [E_c_SS_k_cp[5] + E_c_OS_k_cp[5], E_c_SS_k_cp[11] + E_c_OS_k_cp[11], 
                  E_c_SS_k_cp[7] + E_c_OS_k_cp[7], E_c_SS_k_cp[9] + E_c_OS_k_cp[9], E_c_SS_k_cp[5] + E_c_OS_k_cp[5]],
        "ksskos-MP2": [E_c_SS_k_cp[3] + E_c_OS_k_cp[8], E_c_SS_k_cp[5] + E_c_OS_k_cp[11], 
                       E_c_SS_k_cp[4] + E_c_OS_k_cp[8], E_c_SS_k_cp[10] + E_c_OS_k_cp[7], E_c_SS_k_cp[3] + E_c_OS_k_cp[8]],
        "coskos-MP2": [2.1*E_c_OS_k_cp[3], 2.1*E_c_OS_k_cp[7], 2.3*E_c_OS_k_cp[5], 2.5*E_c_OS_k_cp[4], 2.1*E_c_OS_k_cp[3]],
        "cos-MP2": [1.7*E_c_OS_cp, 1.8*E_c_OS_cp, 2.2*E_c_OS_cp, 2*E_c_OS_cp, 1.7*E_c_OS_cp]
    }
    
    for name, emp2_cp in EMP2vals_cp.items():
        E_c_int_cp.append(form_com_cp.mp2(params[name][0], emp2_cp[0][-1]) - form_frags_cp.mp2(params[name][0], sum(emp2_cp[0][:-1])))
        E_c_int_cp.append(form_com_cp.spl2(params[name][1], emp2_cp[1][-1]) - form_frags_cp.spl2(params[name][1], sum(emp2_cp[1][:-1])))
        E_c_int_cp.append(form_com_cp.f1(params[name][2], emp2_cp[2][-1]) - form_frags_cp.f1(params[name][2], sum(emp2_cp[2][:-1])))
        E_c_int_cp.append(form_com_cp.f1(params[name][3], emp2_cp[3][-1]) - form_frags_cp.f1(params[name][3], sum(emp2_cp[3][:-1])))
        E_c_int_cp.append(form_com_cp.f1(params[name][4], emp2_cp[4][-1]) - form_frags_cp.f1(params[name][4], sum(emp2_cp[4][:-1])))
    
    E_c_ints_cp = dict(zip(funcs, kcal*(ehfdiv_cp+np.array(E_c_int_cp))))
    with open("Eint_kcalmol_all_CP.json","w",encoding="utf-8") as f:
        json.dump(E_c_ints_cp, f)
    
    bsse_corrections = {func: E_c_ints[func] - E_c_ints_cp[func] for func in funcs}
    with open("BSSE_corrections_kcalmol.json","w",encoding="utf-8") as f:
        json.dump(bsse_corrections, f)

# ---------------------------------------------------------
# Formatting the professional output file ("mpac_job.out")
# ---------------------------------------------------------
write_mpac_job_out(
    funcs=funcs,
    E_c_ints=E_c_ints,
    einthf_kcal=einthf_kcal,
    filename="mpac_job.out",
    cp_enabled=args.cp,
    E_c_ints_cp=E_c_ints_cp if args.cp else None,
    bsse_corrections=bsse_corrections if args.cp else None,
    ehfdiv_cp_kcal=ehfdiv_cp*kcal if args.cp else None
)