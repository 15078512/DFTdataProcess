
'''
2020.09.05
Tried:
    pymatgen PWOutput, only returns energy
    pymatgen PWInput.from_file and .from_string, both report error
Therefore, use ase

When use ase, got the following error
>>> a = read_espresso_out('vcrl.out')
>>> for i in a:
...     print(type(i))
... 
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "python3.8/site-packages/ase/io/espresso.py", line 169, in read_espresso_out
    for image_index in image_indexes:
TypeError: 'int' object is not iterable

dig the code, found that image_indexes should be int value, as index = -1
        image_indexes = all_config_indexes[index]

so change line 169 as the following, i.e. add [] to image_indexes to make it a one item
list, and it works well
    for image_index in [image_indexes]:
'''

import os
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.espresso import read_espresso_out
from CrystalToolkit.geometry import planar_structure_normalization, \
                                    change_lattice_constants

calc_dir = '/data/scratch/exw735/2dnets_row3'
not_planar_list = []
not_run_list = []

for dir in os.listdir(calc_dir):
    pw_out = calc_dir + '/' + dir + '/vcrl.out'
    print(pw_out, flush=True)
    strus_ase = read_espresso_out(pw_out)
    if strus_ase is not None:
        for stru_ase in strus_ase:
            stru = AseAtomsAdaptor.get_structure(stru_ase)
            stru.sort()

            is_planar, stru = planar_structure_normalization(structure=stru)
            if is_planar:
                newlatt = [[None,None,None], [None,None,None], [None,None,15.0]]
                stru = change_lattice_constants(structure=stru, lattice=newlatt)

                stru.to(fmt='poscar', filename=dir+'.vasp')
            else:
                not_planar_list.append(dir)
    else:
        not_run_list.append(dir)

with open ('/data/scratch/exw735/not_planar', 'a+') as f:
    for item in not_planar_list:
        f.write("%s\n" % item)

with open ('/data/scratch/exw735/not_run', 'a+') as f:
    for item in not_run_list:
        f.write("%s\n" % item)