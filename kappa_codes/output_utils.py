import datetime

MPAC_LOGO = r"""
 ███╗   ███╗██████╗  █████╗  ██████╗
 ████╗ ████║██╔══██╗██╔══██╗██╔════╝
 ██╔████╔██║██████╔╝███████║██║     
 ██║╚██╔╝██║██╔═══╝ ██╔══██║██║     
 ██║ ╚═╝ ██║██║     ██║  ██║╚██████╗
 ╚═╝     ╚═╝╚═╝     ╚═╝  ╚═╝ ╚═════╝
 Moller-Plesset Adiabatic Connection Functionals (MPAC)
 using PySCF
 Gernealized for NCI interactions of N-body systems

 Contributors:
 E. Palos
 K.J. Daas
 D.P. Kooi
 S. Vuckovic 
"""

def print_job_header(mols, basis, use_df, cp_enabled, charges):
    """Prints the standardized MPAC job initialization header."""
    start_time = datetime.datetime.now()
    
    print(MPAC_LOGO)
    print("="*73)
    print("               MPAC INTERACTION ENERGY CALCULATION")
    print("="*73)
    print(f"Date/Time        : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Fragments Found  : {len(mols)-1} ({', '.join(mols[:-1])})")
    print(f"Basis Set        : {basis}")
    print(f"Density Fitting  : {'ENABLED' if use_df else 'DISABLED'}")
    print(f"Counterpoise     : {'ENABLED' if cp_enabled else 'DISABLED'}")
    print(f"Charges Utilized : {charges}")
    print("="*73)
    
    return start_time

def write_mpac_job_out(funcs, E_c_ints, einthf_kcal, filename="mpac_job.out", 
                       cp_enabled=False, E_c_ints_cp=None, bsse_corrections=None, ehfdiv_cp_kcal=None):
    """Writes the standardized mpac_job.out summary table."""
    with open(filename, "w") as out:
        def write_output(text, to_console=True):
            if to_console:
                print(text)
            out.write(text + "\n")
            
        write_output(MPAC_LOGO, to_console=False)
        write_output("\n" + "="*73)
        write_output("                  FINAL MPAC INTERACTION ENERGIES")
        write_output("="*73)
        
        if cp_enabled:
            write_output(f"{'Functional':<18} | {'Uncorrected':<15} | {'CP-Corrected':<15} | {'BSSE (kcal/mol)':<15}")
            write_output("-" * 73)
            write_output(f"{'Hartree-Fock':<18} | {einthf_kcal:<15.4f} | {ehfdiv_cp_kcal:<15.4f} | {(einthf_kcal - ehfdiv_cp_kcal):<15.4f}")
            write_output("-" * 73)
            for func in funcs:
                write_output(f"{func:<18} | {E_c_ints[func]:<15.4f} | {E_c_ints_cp[func]:<15.4f} | {bsse_corrections[func]:<15.4f}")
        else:
            write_output(f"{'Functional':<18} | {'Interaction Energy (kcal/mol)':<25}")
            write_output("-" * 73)
            write_output(f"{'Hartree-Fock':<18} | {einthf_kcal:<15.4f}")
            write_output("-" * 73)
            for func in funcs:
                write_output(f"{func:<18} | {E_c_ints[func]:<15.4f}")
                
        write_output("="*73)
        write_output("                      >>> JOB COMPLETE <<<")
        write_output("="*73 + "\n")