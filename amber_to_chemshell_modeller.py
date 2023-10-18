# -*- coding: utf-8 -*-

# Import packages
from parmed import load_file
from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError
import sys, os
from argparse import ArgumentParser

def parser():
    parser = ArgumentParser(description="Script to adapt AMBER topology and coordinates to be used in ChemShell.")

    parser.add_argument('-p', '--parameters',
                        help='AMBER prmtop file containing topology and parameters of the system.',
                        type=str,
                        required=True
                        )

    parser.add_argument('-c', '--coordinates',
                        help='AMBER inpcrd file containing the coordinates of the system',
                        required=True,
                        type=str
                        )


    parser.add_argument('-cs', '--crop_system',
                        help="Trigger for activating the cropping of the system's solvent in a radius around a residue. If used, '-csc' has to be specified and '-csr' is optional with a default value of 17 Å.",
                        action='store_true',
                        default=False,
                        required=False
                        )
    parser.add_argument('-csc', '--crop_system_center',
                        help='Number of the residue that will be used for creating the drop of solvent',
                        required=False,
                        default=0
                        )
    parser.add_argument('-csr', '--crop_system_radius',
                    help="Radius (in Å) around the central residue specified with the '-csc' flap. Default value is 17 Å.",
                    required=False,
                    default=17
                    )

    parser.add_argument('-al', '--active_atoms_list',
                        help="Trigger for creating a tcl list containing the list of atoms around an specific atom. This is required for ChemShell QM/MM simulations.",
                        required=False,
                        action='store_true',
                        default=False
                        )
    parser.add_argument('-alc', '--active_atoms_list_center',
                        help='Number of the residue that will be used as centre for creating the active atoms list',
                        required=False,
                        default=0,
                        )
    parser.add_argument('-alr', '--active_atoms_list_radius',
                        help="Radius (in Å) around the central atom specified with the '-alc' flap. Default value is 15 Å.",
                        required=False,
                        default=15
                        )

    parser.add_argument('-ra', '--rename_atoms',
                        help="Trigger for changing conflicting atom names/types used in AMBER to atoms names understandable by ChemShell. Atom names from MCPB.py are also renamed if the original names are given with the '-ram' flag.",
                        required=False,
                        action='store_true',
                        default=False
                        )
    parser.add_argument('-ram', '--rename_atoms_MCPB',
                        help="List of atom names for changing MCPB names (like Y1, Y2, M1, etc.)",
                        nargs='+',
                        type=str,
                        required=False
                        )
    parser.add_argument('-rac', '--rename_atoms_MCPB_checker',
                        help="Trigger for checking the existance of any atom name unknown by ChemShell. It prints the list of unknown atoms and exists the program.",
                        required=False,
                        default=False,
                        action='store_true'
                        )


    args = parser.parse_args()

    if args.crop_system == True and args.crop_system_center == 0:
        sys.exit("If '-cs' is activated, the center residue for cropping the solvent has to be specified with the '-csc' flag.")

    if args.active_atoms_list == True and args.active_atoms_list_center == 0:
        sys.exit("If '-al' is activated, the center atom for creating the atoms list has to be specified with the '-alc' flag.")


def main():
    parser()

if __name__ == '__main__':
    main()