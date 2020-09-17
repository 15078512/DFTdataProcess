
import os
import numpy as np
import math
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar

from CrystalToolkit.composition import solid_solution_random
from CrystalToolkit.structure import poscar_fix_top_bottom, stacking_fault_generation


# note that slip is only along a1vect_uvw direction
slip_systems = [
    {
        'hkl': [1, 1, 0], 
        'a1vect_uvw': [1, -1, 0], 
        'a2vect_uvw': [0, 0, 1]
    },
    {
        'hkl': [1, 1, 1], 
        'a1vect_uvw': [1, -1, 0], 
        'a2vect_uvw': [1, 1, -2]
    },
    {
        'hkl': [1, 1, 1], 
        'a1vect_uvw': [1, 1, -2], 
        'a2vect_uvw': [1, -1, 0]
    },
]

compo_list = [
    'Ti', 'Zr', 'Nb', 'Ta',
    'Ti-Zr-Nb-Ta'
]
a2 = 0.0

for compo in compo_list:  
    stru_file = compo + 'N.cif'
    compo_name = compo + '-N'

    for slip_sys in slip_systems:
        hkl, a1vect_uvw, a2vect_uvw = slip_sys.values()

        elem_list = []
        a1_max, nstep = 0.25, 11
        nstep = 11

        if hkl == [1,1,1]:
            sizemults = [1, 1, 4]
            if a1vect_uvw == [1,1,-2]:
                a1_max, nstep = 0.5, 21
        elif hkl == [1,1,0]:
            sizemults = [1, 1, 6]

        if len(compo) > 2:
            elem_frac_site = 'Ti'
            elem_list = compo.split('-')
            sizemults = [x1 * x2 for x1, x2 in zip(sizemults, [2,2,1])]
            a1_max, nstep = 2, 21

        for a1 in np.linspace(0, a1_max, nstep):
            # every time stru should be read, otherwise it seems to be changed by the subroutine
            stru = Structure.from_file(stru_file)
            stru_fault = stacking_fault_generation(structure=stru, sizemults=sizemults, a1=a1, a2=a2, 
                                                    hkl=hkl, a1vect_uvw=a1vect_uvw, a2vect_uvw=a2vect_uvw)
            if elem_list:
                stru_fault_solidsolu = solid_solution_random(structure=stru_fault, 
                                                    elem_frac_site=elem_frac_site, elem_list=elem_list)
            else:
                stru_fault_solidsolu = stru_fault

            selective_dynamics = poscar_fix_top_bottom(stru_fault_solidsolu)
            fault_sys_poscar = Poscar(structure=stru_fault_solidsolu, selective_dynamics=selective_dynamics)
            fault_sys_poscar.sort()

            if math.isclose(a1, 0.0):
                dir = compo_name \
                        + '_'.join([''.join(map(str, hkl)), 'plane', ''.join(map(str, a1vect_uvw))]) \
                        + '_' + format(a1, '.2f')
            else:
                dir = compo_name \
                        + '_'.join([''.join(map(str, hkl)), 'plane'])

            if not os.path.isdir(dir):
                os.mkdir(dir)
                fault_sys_poscar.write_file(dir + '/POSCAR', ) 


