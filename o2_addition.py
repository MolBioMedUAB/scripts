#from _typeshed import Self
#from o2_addition import create_o2_universe, merge_protein_o2
import MDAnalysis as mda
from MDAnalysis.lib import distances as mdadist
#from MDAnalysis.core.groups import Residue
import numpy as np
import warnings
from math import dist, sqrt
import os
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='o2_addition',
#                usage
                description='Arguments for using o2_addition as CLI program')

    parser.add_argument(
        '-p',
        '--pdbfilename',
        type=str,
        action='store',
        required=True,
        help='Input the filename of the input structre (as pdb).',
        )

    parser.add_argument(
        '-i',
        '--at_num_index',
        type=int,
        action='store',
        required=True,
        help='Index number of the central atom for placing the oxygen molecule.',
        )

    parser.add_argument(
        '-d',
        '--distance',
        type=float,
        action='store',
        default=3.0,
        help='Distance for placing the oxygen molecule around the central atom. Default is 3.0 ang.',
        )

    parser.add_argument(
        '-r',
        '--resolution',
        type=str,
        action='store',
        choices=['high', 'medium', 'low'],
        default='high',
        help='Options for selecting the resolution when placing the oxygens. \'High\' generates 26 places, \'medium\' generates 18 places and \'low\' generates 6 places.',
        )

    parser.add_argument(
        '-s',
        '--save_output_pdb',
        type=str,
        action='store',
        choices=['y', 'n'],#, 'Y', 'N', 'yes', 'no', 'YES', 'NO'],
        default='y',
        help='Trigger for the (de)activation of the option of saving the output pdbs. Default is yes',
        )

    parser.add_argument(
        '-o',
        '--output_prefix',
        type=str,
        action='store',
        default='',
        help='Prefix for naming the output pdbs.',
        )

    parser.add_argument(
        '-w',
        '--hide_saving_warnings',
        type=str,
        action='store',
        choices=['y', 'n'],#, 'Y', 'N', 'yes', 'no', 'YES', 'NO'],
        default='y',
        help='Trigger for the (de)activation of the warnings when saving the pdbs',
        )


    args = vars(parser.parse_args())


    pdbfilename = args['pdbfilename']
    at_num      = args['at_num_index']
    distance    = args['distance']

    if args['resolution'] == None:
        resolution = 'high'
    else :
        resolution = args['resolution']

    prefix = args['output_prefix']

    if args['save_output_pdb'] == 'y':
        savepdb = True
    elif args['save_output_pdb'] == 'n':
        savepdb = False

    if args['hide_saving_warnings'] == 'y':
        hide_saving_warnings = True
    elif args['hide_saving_warnings'] == 'n':
        hide_saving_warnings = False




class AddOxygen:

    def __init__(self, pdbfilename, at_num, distance=3.0, resolution='high', savepdb=True, prefix='', hide_saving_warnings=True):
        """
        DESCRIPTION:
            __init__ method. Loads protein, generates coc and saves distance and resolutions parameters.
        """

        if hide_saving_warnings == True:
                warnings.filterwarnings('ignore')

        self.protein       = mda.Universe(pdbfilename)

        if hide_saving_warnings == True:
                warnings.filterwarnings('default')

        self.coc           = self.protein.select_atoms('bynum %s' % at_num).positions[0]
        self.distance      = distance
        self.resolution    = resolution
        self.savepdb       = savepdb
        self.prefix        = prefix
        self.hide_warnings = hide_saving_warnings

        self.positions = {
            'high' : np.array([
                [1,0,0],
                [-1,0,0],
                [0,1,0],
                [0,-1,0],
                [0,0,1],
                [0,0,-1],

                [(1/sqrt(2)),(1/sqrt(2)),0],
                [(1/sqrt(2)),0,(1/sqrt(2))],
                [0,(1/sqrt(2)),(1/sqrt(2))],
                [-(1/sqrt(2)),(1/sqrt(2)),0],
                [-(1/sqrt(2)),0,(1/sqrt(2))],
                [0,-(1/sqrt(2)),(1/sqrt(2))],
                [(1/sqrt(2)),-(1/sqrt(2)),0],
                [(1/sqrt(2)),0,-(1/sqrt(2))],
                [0,(1/sqrt(2)),-(1/sqrt(2))],
                [-(1/sqrt(2)),-(1/sqrt(2)),0],
                [-(1/sqrt(2)),0,-(1/sqrt(2))],
                [0,-(1/sqrt(2)),-(1/sqrt(2))],

                [(1/sqrt(3)),(1/sqrt(3)),(1/sqrt(3))],
                [-(1/sqrt(3)),(1/sqrt(3)),(1/sqrt(3))],
                [(1/sqrt(3)),-(1/sqrt(3)),(1/sqrt(3))],
                [(1/sqrt(3)),(1/sqrt(3)),-(1/sqrt(3))],
                [-(1/sqrt(3)),-(1/sqrt(3)),(1/sqrt(3))],
                [-(1/sqrt(3)),(1/sqrt(3)),-(1/sqrt(3))],
                [(1/sqrt(3)),-(1/sqrt(3)),-(1/sqrt(3))],
                [-(1/sqrt(3)),-(1/sqrt(3)),-(1/sqrt(3))]
            ]),

            'medium' : np.array([
                [1,0,0],
                [-1,0,0],
                [0,1,0],
                [0,-1,0],
                [0,0,1],
                [0,0,-1],

                [(1/sqrt(2)),(1/sqrt(2)),0],
                [(1/sqrt(2)),0,(1/sqrt(2))],
                [0,(1/sqrt(2)),(1/sqrt(2))],
                [-(1/sqrt(2)),(1/sqrt(2)),0],
                [-(1/sqrt(2)),0,(1/sqrt(2))],
                [0,-(1/sqrt(2)),(1/sqrt(2))],
                [(1/sqrt(2)),-(1/sqrt(2)),0],
                [(1/sqrt(2)),0,-(1/sqrt(2))],
                [0,(1/sqrt(2)),-(1/sqrt(2))],
                [-(1/sqrt(2)),-(1/sqrt(2)),0],
                [-(1/sqrt(2)),0,-(1/sqrt(2))],
                [0,-(1/sqrt(2)),-(1/sqrt(2))],
            ]),

            'low' : np.array([
                [1,0,0],
                [-1,0,0],
                [0,1,0],
                [0,-1,0],
                [0,0,1],
                [0,0,-1],
            ])
        }


    def create_o2_universe(self):
        n_residues = 1
        n_atoms = n_residues * 2

        resindices = np.repeat(range(n_residues), 2)

        segindices = [0] * n_residues

        oxy = mda.Universe.empty(n_atoms,
                                n_residues=n_residues,
                                atom_resindex=resindices,
                                residue_segindex=segindices,
                                trajectory=True)

        oxy.add_TopologyAttr('name', ['O1', 'O2']*n_residues)
        oxy.add_TopologyAttr('type', ['O', 'O']*n_residues)
        oxy.add_TopologyAttr('resname', ['OXY']*n_residues)
        oxy.add_TopologyAttr('resid', [1])


        return oxy


    def build_positions(self, oxy):
        """
        DESCRIPTION:

        """

        o2_interatomic = 1.209152596

        if str(self.resolution).lower() in ('high', 'h'):
            positions = self.positions['high']
        elif str(self.resolution).lower() in ('medium', 'm'):
            positions = self.positions['medium']
        elif str(self.resolution).lower() in ('low', 'l'):
            positions = self.positions['low']

        oxys = []

        for p in range(len(positions)):
            #o = self.create_o2_universe()
            oxys.append(oxy.copy())


            oxys[p].atoms.positions = np.array(
                [positions[p]*(self.distance)+self.coc,
                positions[p]*(self.distance+o2_interatomic)+self.coc])

        return oxys



    def merge_protein_oxy(self, oxys, savepdb=True, pdbfilename=pdbfilename, prefix=''):

        if prefix != '':
            self.prefix = prefix

        if self.prefix == '':

            if 'protein_oxy' not in os.listdir():
                os.mkdir('protein_oxy')

            self.prefix = 'protein_oxy/' + pdbfilename.split('.')[0]


        if savepdb != True:
            self.savepdb = savepdb



        # Finding first solvent molec (if any)
        solvent = self.protein.select_atoms('resname WAT or resname Na+ or resname NA+ or resname Cl- or resname CL-')

        if len(solvent.residues) != 0:
            first_solvent = int(str(solvent.residues[0])[(str(solvent.residues[0]).find(', ')+2):str(solvent.residues[0]).find('>')])
            last_solvent  = int(str(solvent.residues[-1])[(str(solvent.residues[-1]).find(', ')+2):str(solvent.residues[-1]).find('>')])

        else :
            first_solvent, last_solvent = None



        # Generating each universe
        #if merge_all == True:
        #    protein_all = self.protein.copy()
        #    protein_all.residues[first_solvent:last_solvent].resids = list(range(first_solvent + len(oxys), last_solvent+len(oxys)))

        # Modifying resids indexes
        if first_solvent != None and last_solvent != None:
            self.protein.residues[first_solvent:last_solvent].resids = list(range(first_solvent +1, last_solvent+1))

        proteins = []



        if isinstance(oxys, list):
            if self.hide_warnings == True:
                warnings.filterwarnings('ignore')
            for oxy in range(len(oxys)):
                if first_solvent != None and last_solvent != None:
                    oxys[oxy].residues.resids = first_solvent + 2

                    proteins.append(
                        mda.Merge(self.protein.residues[:first_solvent+1].atoms, oxys[oxy].residues[:].atoms, self.protein.residues[first_solvent+1:].atoms)
                    )

                    resnum = 1
                    for r in range(0, proteins[oxy].residues.n_residues):
                        proteins[oxy].residues[r].resid = resnum
                        resnum += 1

                    if self.savepdb == True:
                        proteins[oxy].atoms.write('%s_OXY_%s.pdb' % (self.prefix,oxy+1))

            if self.hide_warnings == True:
                warnings.filterwarnings('default')

            return proteins

        else :
            if self.hide_warnings == True:
                warnings.filterwarnings('ignore')

            if first_solvent != None and last_solvent != None:
                oxys.residues.resids = first_solvent

            protein = mda.Merge(self.protein.residues[:first_solvent-1].atoms, oxys.residues[:].atoms, self.protein.residues[first_solvent-1:].atoms)


            if self.savepdb == True:
                protein.atoms.write('%s_OXY.pdb' % (self.prefix))

            if self.hide_warnings == True:
                warnings.filterwarnings('default')

            return protein


    def merge_protein_all_oxys(self, oxys, savepdb=True, pdbfilename=pdbfilename, prefix='', sufix='_ALL_OXY'):

        if prefix != '':
            self.prefix = prefix

        if self.prefix == '':

            if 'protein_oxy' not in os.listdir():
                os.mkdir('protein_oxy')

            self.prefix = 'protein_oxy/' + pdbfilename.split('.')[0]


        if savepdb != True:
            self.savepdb = savepdb



        # Finding first solvent molec (if any)
        solvent = self.protein.select_atoms('resname WAT or resname Na+ or resname NA+ or resname Cl- or resname CL-')

        if len(solvent.residues) != 0:
            first_solvent = int(str(solvent.residues[0])[(str(solvent.residues[0]).find(', ')+2):str(solvent.residues[0]).find('>')])
            last_solvent  = int(str(solvent.residues[-1])[(str(solvent.residues[-1]).find(', ')+2):str(solvent.residues[-1]).find('>')])

        else :
            first_solvent, last_solvent = None


        protein_all = self.protein.copy()
        #protein_all.residues[first_solvent:last_solvent].resids = list(range(first_solvent + len(oxys), last_solvent + len(oxys)))

        # Modifying resids indexes
        if first_solvent != None and last_solvent != None:
            protein_all.residues[first_solvent-1:last_solvent-1].resids = list(range(first_solvent + len(oxys), last_solvent + len(oxys)))

        if self.hide_warnings == True:
            warnings.filterwarnings('ignore')

        for oxy in range(len(oxys)):
            if first_solvent != None and last_solvent != None:
                oxys[oxy].residues.resids = first_solvent + oxy

                if oxy == 0:
                    protein = mda.Merge(protein_all.residues[:first_solvent-1].atoms, oxys[oxy].residues.atoms)
                elif oxy == (len(oxys) - 1):
                    protein = mda.Merge(protein.residues[:].atoms, oxys[oxy].residues.atoms, protein_all.residues[first_solvent-1:].atoms)
                else :
                    protein = mda.Merge(protein.residues[:].atoms, oxys[oxy].residues.atoms)

        if self.savepdb == True:
            protein.atoms.write('%s_%s.pdb' % (self.prefix, sufix))

        if self.hide_warnings == True:
            warnings.filterwarnings('default')

        return protein


    def check_clashes(self, proteins, prefix=''):
        """
        DESCRIPTION:

        """

        if prefix != '':
            self.prefix = prefix



        oxy_indices = []

        fc  = open(self.prefix + '_clashes.txt', 'w')
        fnc = open(self.prefix + '_noclashes.txt', 'w')

        fc.write('The following positions have a clash with an atom of the protein or the substrate.\n\n')
        fnc.write('The following positions do not have a clash with an atom of the protein or the substrate.\n\n')
        #print('The following positions have a clash with an atom of the protein or the substrate.\n\n')
        print('The following positions do not have a clash with an atom of the protein or the substrate.\n\n')

        counter = 0
        for protein in range(len(proteins)):

            sel_oxyenv  = proteins[protein].select_atoms('(around 4 resname OXY) and not resname OXY').positions
            sel_oxy     = proteins[protein].select_atoms('resname OXY').positions

            oxyenv_dist = np.min(mdadist.distance_array(sel_oxy, sel_oxyenv))

            if oxyenv_dist < 1.4:
                fc.write('\tPosition number %s has a clash.\n' % (protein + 1))
                #print('\tPosition number %s has a clash.' % (protein + 1))
            else :
                fnc.write('\tPosition number %s does not have a clash.\n' % (protein +1))
                print('\tPosition number %s does not have a clash.' % (protein +1))
                oxy_indices.append(counter)

            counter +=1

        fc.close()
        fnc.close()

        return oxy_indices


    def run(self, saveall=True):

        #u             = AddOxygen()
        oxys                  = self.build_positions(self.create_o2_universe())#, self.distance)#, self.resolution)
        proteins              = self.merge_protein_oxy(oxys)
        not_clashing_oxys_ids = self.check_clashes(proteins)
        not_clashing_oxys = []

        for id in not_clashing_oxys_ids:
            not_clashing_oxys.append(oxys[id])

        self.merge_protein_all_oxys(not_clashing_oxys, sufix='NOCLASH_OXY')

        if saveall == True:
            self.merge_protein_all_oxys(oxys)



def main(pdbfilename, at_num, distance, resolution, savepdb, prefix, hide_saving_warnings):

    system = AddOxygen(pdbfilename, at_num, distance, resolution, savepdb, prefix, hide_saving_warnings)
    system.run()


if __name__ == '__main__':

    main(pdbfilename, at_num, distance, resolution, savepdb, prefix, hide_saving_warnings)











