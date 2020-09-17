
import numpy as np
from itertools import product


def bond_length_estimation(atom_list, elem, factor=1.0):
    '''
    Estimate the maximum distance between any two types of atoms in  
    atom_list using crystal radius of the atoms

    Args:
        atom_list: a list containing atomic symbols, like ['Na','Cl']
        elem: atom properties from the json file
        factor: scale factor 
    Return:
        the value of estimated maximum bond length
    '''
    pair_list = list(product(atom_list, repeat=2))
    bondlen_list = []
    for pair in pair_list:
        bondlen = 0
        for atom in pair:
            bondlen += elem[atom]['crystal_rad']
        bondlen_list.append(bondlen)
    return np.array(bondlen_list).max() * factor


def total_atomic_energy(atomic_formula, elem):
    '''
    Calculated the energy sum of the atoms, by summing all the atomic energies

    Args:
        atomic_formula: list.item() iterable pymatgen formula, typically from
                        mongodb entries
        elem: atom properties from the json file
    Return:
        total atomic energy
    '''
    total_atoms_ene = 0.0
    for atom, num in atomic_formula:
        total_atoms_ene += (elem[atom]['atom_energy']*num)
    return total_atoms_ene


def atoms_rad_statis(atomic_formula, elem, anion_num):
    '''
    Give atomic raduis statistics, i.e. standard deviation and (max-min) in a structure. 
    Currently it is designed for high-entropy ceramics materails, the statistics is on cations
    And the anion atom is found using the anion_num value

    Args:
        atomic_formula: list.item() iterable pymatgen formula
        elem: atom properties from the json file
        anion_num: number of anion atoms in a high-entropy ceramics compound
    Return:
        cation and anion atom list
        standard deviation and (max-min) of atomic and crystal radius 
    '''
    cation = []
    atom_rads = []
    cryst_rads = []
    for atom, num in atomic_formula:
        if num == anion_num :
            anion = atom
        else :
            cation.append(atom)
            atom_rads.append(elem[atom]['atomic_rad'])
            cryst_rads.append(elem[atom]['crystal_rad'])
    return anion, cation, \
            np.array(atom_rads).std(), np.array(cryst_rads).std(), \
            np.array(atom_rads).ptp(), np.array(cryst_rads).ptp()


