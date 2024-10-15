"""
Generalized MPAC functionals Interaction Energies
for N-molecule clusters
E_int(N) = E_N - sum_i^N E_i  (e.g., dimer, trimer, ..., decamer)
supports different charges for different fragments

Contributors:
Etienne Palos (from v2.0)
K.J. Daas
D.P. Kooi
S. Vuckovic 
"""

#import required libs: pyscf, numpy and numba
import numpy as np
import argparse
import os
from all_codes.numba_all import *
from all_codes.mol_all import run_pyscf
from all_codes.mpac_fun_all import MPAC_functionals
from all_codes.constants_all import *
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--charge", type=int, default=0, help="the charge of the system") #good for neutral
    parser.add_argument("--charges", nargs='+', type=int, help="charge per fragment i, i= 1,...,N. Last arg is complex charge.")
    parser.add_argument("--spin", type=int, default=0, help="the spin of the system")
    parser.add_argument("--basis", type=str, default="def2-qzvppd",help="the basisset used in the calculations")

    args = parser.parse_args()

    # See kappa_tools for helper scripts to prepare your workind dir. 
    # NOTE: For now, order is assumed to be fragment1, fragment2, ..., fragmentN, complex.
    mols=["A","B","C","ABC"] 
    charges = args.charges  # Expecting charges for each fragment and the complex in same order as mol

    if len(charges) != len(mols):
        raise ValueError("ERROR: Number of charges must match the number of systems: N+1.")
        raise ValueError("Charges expected as arguments in order of fragments 1,...,N, then complex.")

print("USER PROVIDED CHARGES:",charges)
mpacf=["spl2","mp2"]

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

for i in range(len(mols)): #run over the fragements and complex
    run_mol=mols[i]
    #add here path to frag m.xyz file
    chkfile="chkfile_"+run_mol+".chk"
    old_pwd=os.getcwd()
    datadir=old_pwd+"/"+run_mol
    os.chdir(datadir)
    ### Runs Hartree-Fock calculation
    py_run=run_pyscf(atom="m.xyz",charge=charges[i],spin=args.spin,basis=args.basis)
    tab, eris = py_run.run_eris(chkfile_name=chkfile,chkfile_dir=datadir)
    # This prints and extracts all of the ingredients except MP2 into tab.csv
    np.savetxt("tab.csv", tab, delimiter=",", fmt='%s')
    ehf.append(tab[0])
    Uh.append(tab[1])
    Ex.append(tab[2])
    rho_4_3.append(tab[3])
    gea_4_3.append(tab[4])
    rho_3_2.append(tab[5])
    gea_7_6.append(tab[6])
    # k1 is for same spin
    # k2 is for the opposite spin
    k1ss = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
    k2ss = k1ss
    np.savetxt("k1.csv", k1ss, delimiter=",", fmt='%s')
    np.savetxt("k2.csv", k2ss, delimiter=",", fmt='%s')
    k1s=np.array(k1ss,dtype=float)
    k2s=np.array(k2ss,dtype=float)
    mp2OS = MP2_energy_kappa_p_OS_parallel(*eris,k2s, 1) #calculate the opposite spin integral with kappa
    E_c_OS_k.append(mp2OS)
    np.savetxt("os.csv", mp2OS, delimiter=",", fmt='%s')
    # Spin scaled \kappa's
    mp2SS = MP2_energy_kappa_p_SS_parallel(*eris,k1s, 1) #calculate the same spin integral with kappa
    E_c_SS_k.append(mp2SS)
    np.savetxt("ss.csv", mp2SS, delimiter=",", fmt='%s')
    e_mp2_split = MP2_energy_split(*eris,) #calculates the same an opposite mp2 integrals
    np.savetxt("mp2.csv", e_mp2_split, delimiter=",", fmt='%s')
    E_c_SS.append(e_mp2_split[0])
    E_c_OS.append(e_mp2_split[1])
    E_c_mp2tot.append(sum(e_mp2_split))
    os.chdir(old_pwd)

# Gets the arrays into the correct shape
E_c_SS_k=np.array(E_c_SS_k).T
E_c_OS_k=np.array(E_c_OS_k).T
E_c_OS=np.array(E_c_OS)

# Initializing the MPAC functionals and calculating the HF energy difference
form_frags=MPAC_functionals(sum(Ex[:-1]),sum(rho_4_3[:-1]),sum(gea_4_3[:-1]))
form_com=MPAC_functionals(Ex[-1],rho_4_3[-1],gea_4_3[-1]) 

# HF interaction energy
ehfdiv=ehf[-1]-(sum(ehf[:-1]))
einthf_kcal=ehfdiv*kcal
print("Hartree-Fock Eint_HF [kcal/mol] =", einthf_kcal)

funcs=["MP2","SPL2","F1","F1ab","k-MP2","k-SPL2","k-F1","k-F1ab","ksskos-MP2","ksskos-SPL2","ksskos-F1","ksskos-F1ab","coskos-MP2","coskos-SPL2","coskos-F1","coskos-F1ab","cos-MP2","cos-SPL2","cos-F1","cos-F1ab"]

# Storing all the EMP2 data
EMP2vals={
    "MP2": [E_c_mp2tot]*4,
    "k-MP2": [E_c_SS_k[5] + E_c_OS_k[5],E_c_SS_k[11] + E_c_OS_k[11],E_c_SS_k[7] + E_c_OS_k[7],E_c_SS_k[9] +E_c_OS_k[9]],
    "ksskos-MP2": [E_c_SS_k[3] + E_c_OS_k[8],E_c_SS_k[5] + E_c_OS_k[11],E_c_SS_k[4] + E_c_OS_k[8],E_c_SS_k[10] + E_c_OS_k[7]],
    "coskos-MP2":[2.1*E_c_OS_k[3],2.1*E_c_OS_k[7],2.3*E_c_OS_k[5],2.5*E_c_OS_k[4]],
    "cos-MP2": [1.7*E_c_OS,1.8*E_c_OS,2.2*E_c_OS,2*E_c_OS]
}

# Calculation for N-fragment systems
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


# Print interaction energies in kkcal/mol to json file
E_c_ints=dict(zip(funcs,kcal*(ehfdiv+np.array(E_c_int))))
print(E_c_ints) #prints out the correct E_c_int
with open("Eint_kcalmol_all.json","w",encoding="utf-8") as f:
    json.dump(E_c_ints,f)

