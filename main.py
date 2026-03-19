"""
Calculating the interaction energy for a MPAC functional defined by the user
"""

#import modules and codes
import numpy as np
import argparse
import os
from kappa_codes.numba_codes import *
from kappa_codes.mol import run_pyscf
from kappa_codes.mpac_fun import MPAC_functionals
from kappa_codes.constants import *
from kappa_codes.cp_utils import get_ghost_atoms_for_fragment

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--charge", type=int, default=0, help="the charge of the system")
    parser.add_argument("--charges", nargs='+', type=int, help="charge per fragment i, i=1,...,N. Last arg is complex charge.")
    parser.add_argument("--spin", type=int, default=0, help="the spin of the system")
    parser.add_argument("--basis", type=str, default="aug-cc-pvqz",help="the basisset used in the calculations")
    parser.add_argument("--func",type=str,default="coskos-SPL2", help="the MP AC functional used including the prefix")
    parser.add_argument("--cp", action="store_true", help="perform counterpoise (Boys-Bernardi) correction for BSSE")

args = parser.parse_args()
func = args.func.lower()
#obtaining the attributes of the inputted functional:
if "kos" in func: #checking if it uses the \kappa regularizer
    kappa=True
else:
    kappa=False
if "cos" in func: #checking if it uses the spin opposite scaling
    cos=True
else:
    cos=False
if "k-" in func: #checking if it is the original \kappa-method
    ksam=True
    kappa=True
else:
    ksam=False
if "mp2" in func: #checking which base functional is used
    mpacf="mp2"
elif "spl2" in func:
    mpacf="spl2"
elif "f1ab" in func:
    mpacf="f1ab"
elif "f1" in func:
    mpacf="f1"
elif "mpac25" in func:
    mpacf="mpac25"
else:
    raise ValueError("no valid functional provided, please use mp2, spl2, f1, f1ab, or mpac25") #gives error if the wrong functional is used

if kappa==False and cos==False and func.split(mpacf)[0]!="":
    raise ValueError("Unknown prefix use coskos-, ksskos-, k- or no prefix") #gives error if the wrong prefix is used

# Auto-detect fragments from directory structure
import string
fragment_labels = []
for letter in string.ascii_uppercase:
    if os.path.isdir(letter):
        fragment_labels.append(letter)
    else:
        break

# Check for COMPLEX directory (standard for N-mer systems)
if not os.path.isdir("COMPLEX"):
    raise ValueError(f"Could not find COMPLEX directory. Found {len(fragment_labels)} fragments: {fragment_labels}")

complex_label = "COMPLEX"

mols = fragment_labels + [complex_label]  # [A, B, ..., COMPLEX/AB]
N_fragments = len(fragment_labels)
print(f"Detected {N_fragments} fragments: {fragment_labels} with complex: {complex_label}")
Ex=[]
ehf=[]
Uh=[]
rho_4_3=[]
gea_4_3=[]
rho_3_2=[]
gea_7_6=[]
E_c_mp2=[]

# For counterpoise correction
if args.cp:
    Ex_cp=[]
    ehf_cp=[]
    E_c_mp2_cp=[]
    print("Counterpoise correction ENABLED")
else:
    print("Counterpoise correction DISABLED")

para,name=params[(kappa,cos,ksam,mpacf)] #obtain parameters for the chosen MPAC functional
print(f"the functional that will be run is: {name}")
kapcoslist=[]

while len(para)>4 or (len(para)<3 and len(para)>0): #removes the non-functional specific paremeters (i.e. removes \kappa_ss, \kappa_os and c_os)
    kapcoslist.append(para.pop())

for i in range(N_fragments + 1): #run over all fragments and complex
    run_mol=mols[i]
    #add here path to frag m.xyz file
    chkfile="chkfile_"+run_mol+".chk"
    old_pwd=os.getcwd()
    datadir=old_pwd+"/"+run_mol
    os.chdir(datadir)
    ###runs HF
    py_run=run_pyscf(atom="m.xyz",charge=args.charge,spin=args.spin,basis=args.basis)
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
    if kappa==True: #if \kappa is turned on
        #k1 is for same spin
        #k2 is for the opposite spin
        k1ss = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
        k2ss = k1ss
        np.savetxt("k1.csv", k1ss, delimiter=",", fmt='%s')
        np.savetxt("k2.csv", k2ss, delimiter=",", fmt='%s')

        k1s=np.array(k1ss,dtype=float)
        k2s=np.array(k2ss,dtype=float)
        k_os=kapcoslist[0]
        mp2OS = MP2_energy_kappa_p_OS_parallel(*eris,k2s, 1) #calculate the opposite spin integral
        np.savetxt("os.csv", mp2OS, delimiter=",", fmt='%s')
        #Spin scaled \kappa's
        if cos==False:
            mp2SS = MP2_energy_kappa_p_SS_parallel(*eris,k1s, 1) #calculate the same spin integral
            np.savetxt("ss.csv", mp2SS, delimiter=",", fmt='%s')
            k_ss=kapcoslist[1] 
            E_c_kmp2_tot= mp2SS[k1ss.index(k_ss)] + mp2OS[k2ss.index(k_os)] #take only the value that corresponds to the optimal k_ss and k_os
            E_c_mp2.append(E_c_kmp2_tot)
        else:
            c_os=kapcoslist[1]
            E_c_kmp2_cos= c_os*mp2OS[k2ss.index(k_os)] #take only the value that corresponds to the optimal c_os and k_os values
            E_c_mp2.append(E_c_kmp2_cos)
    else: #run MP2 without \kappa
        ###Runs E_c^MP2(ss) and E_c^MP2(os)
        e_mp2_split = MP2_energy_split(*eris) #cal
        np.savetxt("mp2.csv", e_mp2_split, delimiter=",", fmt='%s')
        if cos==False: #run regular MP2
            E_c_mp2_tot= sum(e_mp2_split)
            E_c_mp2.append(E_c_mp2_tot)
        else: # run spin opposite scaled mp2
            c_os=kapcoslist[0]
            E_c_mp2_cos= c_os*e_mp2_split[1]
            E_c_mp2.append(E_c_mp2_cos)

    os.chdir(old_pwd)

# Counterpoise correction: calculate fragments in full basis
if args.cp:
    print("\n=== Running Counterpoise Correction ===")
    
    for i in range(N_fragments):  # Only fragments, not complex
        run_mol = mols[i]
        print(f"Calculating {run_mol} in full (ghost) basis...")
        
        chkfile = "chkfile_" + run_mol + "_cp.chk"
        old_pwd = os.getcwd()
        datadir = old_pwd + "/" + run_mol
        os.chdir(datadir)
        
        # Get ghost atoms from all other fragments
        ghost_atoms_str = get_ghost_atoms_for_fragment(fragment_labels, i, base_dir=old_pwd)
        
        # Run calculation with ghost atoms
        py_run = run_pyscf(atom="m.xyz", charge=args.charge, spin=args.spin, 
                          basis=args.basis, ghost_atoms=ghost_atoms_str)
        tab_cp, eris_cp = py_run.run_eris(chkfile_name=chkfile, chkfile_dir=datadir)
        
        np.savetxt("tab_cp.csv", tab_cp, delimiter=",", fmt='%s')
        ehf_cp.append(tab_cp[0])
        Ex_cp.append(tab_cp[2])
        
        # Calculate MP2 with same logic as regular calculation
        if kappa==True:
            k1s=np.array(k1ss,dtype=float)
            k2s=np.array(k2ss,dtype=float)
            k_os=kapcoslist[0]
            mp2OS_cp = MP2_energy_kappa_p_OS_parallel(*eris_cp, k2s, 1)
            np.savetxt("os_cp.csv", mp2OS_cp, delimiter=",", fmt='%s')
            
            if cos==False:
                mp2SS_cp = MP2_energy_kappa_p_SS_parallel(*eris_cp, k1s, 1)
                np.savetxt("ss_cp.csv", mp2SS_cp, delimiter=",", fmt='%s')
                k_ss=kapcoslist[1]
                E_c_kmp2_tot_cp = mp2SS_cp[k1ss.index(k_ss)] + mp2OS_cp[k2ss.index(k_os)]
                E_c_mp2_cp.append(E_c_kmp2_tot_cp)
            else:
                c_os=kapcoslist[1]
                E_c_kmp2_cos_cp = c_os*mp2OS_cp[k2ss.index(k_os)]
                E_c_mp2_cp.append(E_c_kmp2_cos_cp)
        else:
            e_mp2_split_cp = MP2_energy_split(*eris_cp)
            np.savetxt("mp2_cp.csv", e_mp2_split_cp, delimiter=",", fmt='%s')
            
            if cos==False:
                E_c_mp2_tot_cp = sum(e_mp2_split_cp)
                E_c_mp2_cp.append(E_c_mp2_tot_cp)
            else:
                c_os=kapcoslist[0]
                E_c_mp2_cos_cp = c_os*e_mp2_split_cp[1]
                E_c_mp2_cp.append(E_c_mp2_cos_cp)
        
        os.chdir(old_pwd)
    
    print("=== Counterpoise correction calculations complete ===\n")

# Sum all fragment contributions
Ex_frags_sum = sum(Ex[:N_fragments])
E_c_mp2_frags_sum = sum(E_c_mp2[:N_fragments])
rho_4_3_frags_sum = sum(rho_4_3[:N_fragments])
gea_4_3_frags_sum = sum(gea_4_3[:N_fragments])

form_frags=MPAC_functionals(Ex_frags_sum, E_c_mp2_frags_sum, rho_4_3_frags_sum, gea_4_3_frags_sum) #initialize the MPAC functionals
form_com=MPAC_functionals(Ex[N_fragments], E_c_mp2[N_fragments], rho_4_3[N_fragments], gea_4_3[N_fragments])
ehfdiv = ehf[N_fragments] - sum(ehf[:N_fragments])

# Calculate standard interaction energy
if mpacf == "spl2": #calculate the interaction energy of SPL2
    E_c_int=(ehfdiv+form_com.spl2(para)-form_frags.spl2(para))*kcal

elif mpacf == "f1": #calculate the interaction energy of F1
    E_c_int=(ehfdiv+form_com.f1(para)-form_frags.f1(para))*kcal

elif mpacf == "f1ab": #calculate the interaction energy of F1[\alpha,\beta]
    E_c_int=(ehfdiv+form_com.f1(para)-form_frags.f1(para))*kcal

elif mpacf == "mpac25": #calculate the interaction energy of MPAC25
    E_c_int=(ehfdiv+form_com.f1(para)-form_frags.f1(para))*kcal

elif mpacf == "mp2": #calcullate the interaction energy of MP2
    E_c_int=(ehfdiv+form_com.mp2(para)-form_frags.mp2(para))*kcal

print(f"\nThe {name} interaction energy: {E_c_int:.6f} kcal/mol") #prints out the correct E_c_int

# Calculate and print CP-corrected interaction energy
if args.cp:
    # CP-corrected uses fragments calculated in full basis
    Ex_cp_sum = sum(Ex_cp[:N_fragments])
    E_c_mp2_cp_sum = sum(E_c_mp2_cp[:N_fragments])
    rho_4_3_frags_sum = sum(rho_4_3[:N_fragments])  # Use original fragments for grid integrals
    gea_4_3_frags_sum = sum(gea_4_3[:N_fragments])
    
    form_frags_cp = MPAC_functionals(Ex_cp_sum, E_c_mp2_cp_sum, rho_4_3_frags_sum, gea_4_3_frags_sum)
    ehfdiv_cp = ehf[N_fragments] - sum(ehf_cp[:N_fragments])  # Complex - fragments@full_basis
    
    if mpacf == "spl2":
        E_c_int_cp = (ehfdiv_cp + form_com.spl2(para) - form_frags_cp.spl2(para)) * kcal
    elif mpacf == "f1":
        E_c_int_cp = (ehfdiv_cp + form_com.f1(para) - form_frags_cp.f1(para)) * kcal
    elif mpacf == "f1ab":
        E_c_int_cp = (ehfdiv_cp + form_com.f1(para) - form_frags_cp.f1(para)) * kcal
    elif mpacf == "mpac25":
        E_c_int_cp = (ehfdiv_cp + form_com.f1(para) - form_frags_cp.f1(para)) * kcal
    elif mpacf == "mp2":
        E_c_int_cp = (ehfdiv_cp + form_com.mp2(para) - form_frags_cp.mp2(para)) * kcal
    
    bsse = E_c_int - E_c_int_cp
    
    print(f"The {name} interaction energy (CP-corrected): {E_c_int_cp:.6f} kcal/mol")
    print(f"BSSE correction: {bsse:.6f} kcal/mol")
    print(f"BSSE percentage: {abs(bsse/E_c_int)*100:.2f}%")

