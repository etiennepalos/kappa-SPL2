"""
A file containing the parallelized MP2 codes sped up by numba.
"""
#import numpy and numba
import numpy as np
import numba
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"

@numba.njit(parallel=True)
def MP2_energy_split(nocc, e, eri):
    """Calculating the split MP2 correlation energy in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        eri (_type_): 4d array containing the four orbital integrals in the MO basis.

    Returns:
        float, float: the same spin and opposite spin MP2 correlation energies.
    """
    energy_SS = 0. #- EMP2 SS
    energy_OS = 0. #- MP2 OS
    norb = e.shape[0] #total number of orbitals
    e_occ = e[:nocc] #occupied energies
    e_virt = e[nocc:] #virtual energies
    nvirt = norb-nocc #number of virtuals
    #loop over ijab, numba automatically figures out how to parallelize, probably only does so over i, but it can reorder
    for i in numba.prange(nocc):
        for j in numba.prange(nocc):
            for a in numba.prange(nvirt):
                for b in numba.prange(nvirt):
                    energy_SS += eri[i,a,j,b]*(eri[i,a,j,b]-eri[i,b,j,a])/(e_virt[a]+e_virt[b]-e_occ[i]-e_occ[j])
                    energy_OS += eri[i,a,j,b]*(eri[i,a,j,b])/(e_virt[a]+e_virt[b]-e_occ[i]-e_occ[j])
    return - energy_SS, - energy_OS


@numba.njit(parallel=True)
def MP2_energy_kappa_p_SS(nocc, e, eri, kappa_SS, p_SS):
    """Computes the SS MP2 correlation energies with \kappa and p in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        eri (_type_): 4d array containing the four orbital integrals in the MO basis.
        kappa_SS (float): the same spin \kappa regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097
        p_SS (float): the same spin p regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097

    Returns:
        float: the same spin correlation MP2 correlation energy for a specific \kappa and p.
    """
    energy_SS = 0. #- EMP2 SS
    #energy_OS = 0. #- MP2 OS
    norb = e.shape[0] #total number of orbitals
    e_occ = e[:nocc] #occupied energies
    e_virt = e[nocc:] #virtual energies
    nvirt = norb-nocc #number of virtuals
    for i in numba.prange(nocc):
        for j in numba.prange(nocc):
            for a in numba.prange(nvirt):
                for b in numba.prange(nvirt):
                    energy_SS += eri[i,a,j,b]*(eri[i,a,j,b]-eri[i,b,j,a])/(e_virt[a]+e_virt[b]-e_occ[i]-e_occ[j])*(1-np.exp(-kappa_SS*(e_virt[a]+e_virt[b]-e_occ[i]-e_occ[j])**p_SS))**2
    return - energy_SS 


@numba.njit(parallel=True)
def MP2_energy_kappa_p_SS_parallel(nocc, e, eri, kappa_SS, p_SS):
    """Computes the SS MP2 correlation energies with p and for different \kappa in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        eri (_type_): 4d array containing the four orbital integrals in the MO basis.
        kappa_SS (array): 1d array containing the same spin \kappa regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097
        p_SS (float): the same spin p regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097

    Returns:
        ndarray: 1d array containing the same spin MP2 correlation energies for a specific p and a range of different \kappa's.
    """
    energy_SS = np.zeros(kappa_SS.shape[0])
    for k in numba.prange(kappa_SS.shape[0]):
        energy_SS[k] = MP2_energy_kappa_p_SS(nocc, e, eri, kappa_SS[k], p_SS)
    return energy_SS


@numba.njit(parallel=True)
def MP2_energy_kappa_p_OS(nocc, e, eri, kappa_OS, p_OS):
    """Computes the OS MP2 correlation energies with \kappa and p in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        eri (_type_): 4d array containing the four orbital integrals in the MO basis.
        kappa_OS (float): the opposite spin \kappa regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097
        p_OS (float): the opposite spin p regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097

    Returns:
        float: the opposite spin correlation MP2 correlation energy for a specific \kappa and p.
    """
    # energy_SS = 0. #- EMP2 SS
    energy_OS = 0.
    norb = e.shape[0]
    e_occ = e[:nocc] #occupied energies
    e_virt = e[nocc:] #virtual energies
    nvirt = norb-nocc #number of virtuals
    for i in numba.prange(nocc):
        for j in numba.prange(nocc):
            for a in numba.prange(nvirt):
                for b in numba.prange(nvirt):
                    energy_OS += eri[i,a,j,b]*(eri[i,a,j,b])/(e_virt[a]+e_virt[b]-e_occ[i]-e_occ[j])*(1-np.exp(-kappa_OS*(e_virt[a]+e_virt[b]-e_occ[i]-e_occ[j])**p_OS))**2
    return - energy_OS

@numba.njit(parallel=True)
def MP2_energy_kappa_p_OS_parallel(nocc, e, eri, kappa_OS, p_OS):
    """Computes the OS MP2 correlation energies with p and for different \kappa in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        eri (_type_): 4d array containing the four orbital integrals in the MO basis.
        kappa_SS (array): 1d array containing the opposite spin \kappa regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097
        p_SS (float): the opposite spin p regularizer derived from Shee, J.; Loipersberger, M.; Rettig, A.; Lee, J.; Head-Gordon, M. JPCL 2021, 12, 12084–12097

    Returns:
        ndarray: 1d array containing the opposite spin MP2 correlation energies for a specific p and a range of different \kappa's.
    """
    energy_OS = np.zeros(kappa_OS.shape[0])
    for k in numba.prange(kappa_OS.shape[0]):
        energy_OS[k] = MP2_energy_kappa_p_OS(nocc, e, eri, kappa_OS[k], p_OS)
    return energy_OS

# Density Fitting (DF) / Resolution of Identity (RI) MP2 support 
@numba.njit(parallel=True)
def DF_MP2_energy_split(nocc, e, B_ia):
    """Calculating the split DF-MP2 correlation energy in parallel using 3-index integrals.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        B_ia (ndarray): 3d array of shape (nocc, nvirt, naux) containing the 3-index DF integrals.

    Returns:
        float, float: the same spin and opposite spin MP2 correlation energies.
    """
    energy_SS = 0. 
    energy_OS = 0.
    norb = e.shape[0]
    e_occ = e[:nocc]
    e_virt = e[nocc:]
    nvirt = norb - nocc
    naux = B_ia.shape[2]
    
    # Loop over occupied and virtual orbitals
    for i in numba.prange(nocc):
        for j in numba.prange(nocc):
            for a in numba.prange(nvirt):
                for b in numba.prange(nvirt):
                    # Compute (ia|jb) by contracting over auxiliary index
                    eri_iajb = 0.0
                    for P in range(naux):
                        eri_iajb += B_ia[i, a, P] * B_ia[j, b, P]
                    
                    # Compute (ib|ja) for exchange
                    eri_ibja = 0.0
                    for P in range(naux):
                        eri_ibja += B_ia[i, b, P] * B_ia[j, a, P]
                    
                    denom = e_virt[a] + e_virt[b] - e_occ[i] - e_occ[j]
                    energy_SS += eri_iajb * (eri_iajb - eri_ibja) / denom
                    energy_OS += eri_iajb * eri_iajb / denom
    
    return -energy_SS, -energy_OS


@numba.njit(parallel=True)
def DF_MP2_energy_kappa_p_SS(nocc, e, B_ia, kappa_SS, p_SS):
    """Computes the SS DF-MP2 correlation energies with κ and p in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        B_ia (ndarray): 3d array of shape (nocc, nvirt, naux) containing the 3-index DF integrals.
        kappa_SS (float): the same spin κ regularizer.
        p_SS (float): the same spin p regularizer.

    Returns:
        float: the same spin DF-MP2 correlation energy for a specific κ and p.
    """
    energy_SS = 0.
    norb = e.shape[0]
    e_occ = e[:nocc]
    e_virt = e[nocc:]
    nvirt = norb - nocc
    naux = B_ia.shape[2]
    
    for i in numba.prange(nocc):
        for j in numba.prange(nocc):
            for a in numba.prange(nvirt):
                for b in numba.prange(nvirt):
                    # Compute (ia|jb)
                    eri_iajb = 0.0
                    for P in range(naux):
                        eri_iajb += B_ia[i, a, P] * B_ia[j, b, P]
                    
                    # Compute (ib|ja)
                    eri_ibja = 0.0
                    for P in range(naux):
                        eri_ibja += B_ia[i, b, P] * B_ia[j, a, P]
                    
                    denom = e_virt[a] + e_virt[b] - e_occ[i] - e_occ[j]
                    damping = (1.0 - np.exp(-kappa_SS * denom**p_SS))**2
                    energy_SS += eri_iajb * (eri_iajb - eri_ibja) / denom * damping
    
    return -energy_SS


@numba.njit(parallel=True)
def DF_MP2_energy_kappa_p_SS_parallel(nocc, e, B_ia, kappa_SS, p_SS):
    """Computes the SS DF-MP2 correlation energies with p and for different κ in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        B_ia (ndarray): 3d array of shape (nocc, nvirt, naux) containing the 3-index DF integrals.
        kappa_SS (array): 1d array containing the same spin κ regularizers.
        p_SS (float): the same spin p regularizer.

    Returns:
        ndarray: 1d array containing the same spin DF-MP2 correlation energies for different κ's.
    """
    energy_SS = np.zeros(kappa_SS.shape[0])
    for k in numba.prange(kappa_SS.shape[0]):
        energy_SS[k] = DF_MP2_energy_kappa_p_SS(nocc, e, B_ia, kappa_SS[k], p_SS)
    return energy_SS


@numba.njit(parallel=True)
def DF_MP2_energy_kappa_p_OS(nocc, e, B_ia, kappa_OS, p_OS):
    """Computes the OS DF-MP2 correlation energies with κ and p in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        B_ia (ndarray): 3d array of shape (nocc, nvirt, naux) containing the 3-index DF integrals.
        kappa_OS (float): the opposite spin κ regularizer.
        p_OS (float): the opposite spin p regularizer.

    Returns:
        float: the opposite spin DF-MP2 correlation energy for a specific κ and p.
    """
    energy_OS = 0.
    norb = e.shape[0]
    e_occ = e[:nocc]
    e_virt = e[nocc:]
    nvirt = norb - nocc
    naux = B_ia.shape[2]
    
    for i in numba.prange(nocc):
        for j in numba.prange(nocc):
            for a in numba.prange(nvirt):
                for b in numba.prange(nvirt):
                    # Compute (ia|jb)
                    eri_iajb = 0.0
                    for P in range(naux):
                        eri_iajb += B_ia[i, a, P] * B_ia[j, b, P]
                    
                    denom = e_virt[a] + e_virt[b] - e_occ[i] - e_occ[j]
                    damping = (1.0 - np.exp(-kappa_OS * denom**p_OS))**2
                    energy_OS += eri_iajb * eri_iajb / denom * damping
    
    return -energy_OS


@numba.njit(parallel=True)
def DF_MP2_energy_kappa_p_OS_parallel(nocc, e, B_ia, kappa_OS, p_OS):
    """Computes the OS DF-MP2 correlation energies with p and for different κ in parallel.

    Args:
        nocc (integer): number of occupied orbitals.
        e (ndarray): 1d array containing the orbital energies of both the occupied and virtual orbitals.
        B_ia (ndarray): 3d array of shape (nocc, nvirt, naux) containing the 3-index DF integrals.
        kappa_OS (array): 1d array containing the opposite spin κ regularizers.
        p_OS (float): the opposite spin p regularizer.

    Returns:
        ndarray: 1d array containing the opposite spin DF-MP2 correlation energies for different κ's.
    """
    energy_OS = np.zeros(kappa_OS.shape[0])
    for k in numba.prange(kappa_OS.shape[0]):
        energy_OS[k] = DF_MP2_energy_kappa_p_OS(nocc, e, B_ia, kappa_OS[k], p_OS)
    return energy_OS
