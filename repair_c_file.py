#! /home/uabqut17/.conda/envs/py_env/bin/python3

chemsh_PATH = '/home/uabqut17/soft/chemshell-3.7/scripts/chemsh'

import sys, os

# argv[1] --> coord_c
# argv[2] --> atoms_pdb


if sys.argv[1] in ('-h', '-H', '--help') or len(sys.argv) != 3:
    print('Usage:')
    print('\trepair_c_file.py coords.c atoms.pdb')
    print('-h flag prints this message')
    sys.exit()

atoms_c_name = '.'.join(sys.argv[2].split('.')[:-1]) + '.c'

with open('chemsh.tmp', 'w') as chemsh_script:
    chemsh_script.write(f"read_pdb file={sys.argv[2]} coords={atoms_c_name}")

os.system(f"{chemsh_PATH} chemsh.tmp")


coords_c = open(sys.argv[1]).readlines()
atoms_c  = open(atoms_c_name).readlines()

out_c_name = '.'.join(list(sys.argv[1].split('.')[:-1])) + '_atoms.c'
out_c = open(out_c_name, 'w')

for l in range(len(atoms_c)):

    try :
        if coords_c[l].find('(null)') == 0:
            element = atoms_c[l][:2]
            out_c.write(coords_c[l].replace('(null)', element))

        else:
            out_c.write(atoms_c[l])

    except IndexError:
        out_c.write(atoms_c[l])

out_c.close()

with open('chemsh.tmp', 'w') as chemsh_script:
    chemsh_script.write(f"read_pdb file={sys.argv[2]} coords=dummy.coords\n")
    chemsh_script.write(f"write_pdb file={out_c_name[:-2]}.pdb coords={out_c_name}")


os.system(f"{chemsh_PATH} chemsh.tmp")
os.remove('chemsh.tmp')
