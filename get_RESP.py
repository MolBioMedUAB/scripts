import veloxchem as vlx
import os
import sys
import argparse

'''
USAGE:
    Optimisation has not been implemented yet. Use another QM software to optimise the geometry of the molecule if needed and then use this script to get the RESP charges.
    Number of cores can be set using the environment variable OMP_NUM_THREADS (export OMP_NUM_THREADS=n).

DEPENDENCIES INSTALLATION:
    Follow instructions here: https://kthpanor.github.io/echem/docs/install.html
    

'''


def parser():
    xyz_in = sys.argv[1]
    chrg   = sys.argv[2]
    mult   = sys.argv[3]


    return xyz_in, chrg, mult

def argparser():

    parser = argparse.ArgumentParser(description='Get RESP charges from a molecule')


    parser.add_argument('input', type=str, help='Input XYZ or PDB file')
    parser.add_argument('-oo', '--opt_output', type=str, default=None, help='Optimisation geometry output file name')
    parser.add_argument('-c', '--charge', default=0, type=int, help='System\'s charge [0]')
    parser.add_argument('-m', '--multiplicity', default=1, type=int, help='System\'s multiplicity [1]')
    parser.add_argument('-m_opt',  '--method_opt',  type=str, default='B3LYP', choices=['HF'] + list(vlx.available_functionals()), help='QM method for optimisation [B3LYP]')
    parser.add_argument('-m_resp', '--method_resp', type=str, default='B3LYP', choices=['HF'] + list(vlx.available_functionals()), help='QM method for RESP calculation [B3LYP]')
    parser.add_argument('-b_opt',  '--basis_opt',  type=str, default='6-31G*', help='QM basis for optimisation [6-31G*]')
    parser.add_argument('-b_resp', '--basis_resp',  type=str, default='6-31G*', help='QM basis for RESP calculation [6-31G*]')

    parser.add_argument('-solv', '--implicit_solvent', action='store_true', default=False, help='Activates the usage of implicit solvent for optimisations')

    parser.add_argument('-no_opt', '--no_opt', default=False, action='store_true', help='Deactivate geometry optimisation')
    parser.add_argument('--plot_opt', default=False, action='store_true', help='Plot SCF energy along the optimisation')
    parser.add_argument('-r', '--restart', default=False, action='store_true', help='Restart SCF calculation from checkpoint file')

    if 'OMP_NUM_THREADS' in list(os.environ):
        parser.add_argument('-n', '--num_cores', type=int, default=int(os.environ['OMP_NUM_THREADS']), help=f'Number of cores to use for the calculation (OMP_NUM_THREADS has been set to {os.environ["OMP_NUM_THREADS"]})')
    else :
        parser.add_argument('-n', '--num_cores', type=int, default=vlx.environment.cpu_count(), help=f'Number of cores to use for the calculation (default is all available cores: {vlx.environment.cpu_count()})')

    #not implemented yet
    #parser.add_argument('--non_equivalent_atoms', type=bool, default=False, action='store_false', help='Deactivate automatic search of equivalent atoms for RESP fitting')

    args = parser.parse_args()

    if args.opt_output == None:
        args.opt_output = args.input.split('.')[0] + '_opt.' + args.input.split('.')[1]
    elif '.' in args.opt_output and not args.opt_output.endswith('xyz'):
        print("Only xyz extension is allowed for optimised structure. It will be automatically changed.")
        args.opt_output = '.'.join(args.opt_output.split('.')[:-1]) + '.xyz'
    elif '.' not in args.opt_output:
        args.opt_output = args.opt_output + '.xyz'

    return args



def load_molecule(input, charge, mult):
    # molecule = vlx.Molecule.read_xyz_file(xyz_in)

    if input.endswith('.xyz'):
        molecule = vlx.Molecule.read_xyz_file(input)

    elif input.endswith('.pdb'):
        molecule = vlx.Molecule.read_pdb_file(input)

    molecule.set_charge(charge)
    molecule.set_multiplicity(mult)

    return molecule


def optimize(output, molecule, xcfun='B3LYP', basis='6-31G*', max_iter=200, restart=False, solvent=False):

    basis = vlx.MolecularBasis.read(molecule, basis, ostream=None)
    
    if molecule.get_multiplicity() == 1:
        scf_drv = vlx.ScfRestrictedDriver()
    else :
        scf_drv = vlx.ScfUnrestrictedDriver()

    if xcfun in vlx.available_functionals():
        scf_drv.xcfun = xcfun

    scf_drv.max_iter = max_iter

    if solvent:
        # Enable the CPCM solvent model for water
        cpcm_drv = vlx.CpcmDriver()
        cpcm_drv.solver.method = "direct"
        cpcm_drv.solvent.dielectric_constant = 78.39
        scf_drv.solvent_model = cpcm_drv
    
    scf_results = scf_drv.compute(molecule, basis)

    opt_drv = vlx.OptimizationDriver(scf_drv)
    opt_drv.max_iter = max_iter
    opt_drv.checkpoint_file = 'opt.chkp'

    opt_results = opt_drv.compute(molecule, basis, scf_results)

    with open(output, 'w') as optf:
        optf.write(opt_results['opt_geometries'][-1])

    return opt_results


def plot_scf_along_opt(opt_results):
    from matplotlib import pyplot as plt

    e_min_in_au = min(opt_results["opt_energies"])

    energies_in_kcalpermol = [
        (e - e_min_in_au) * vlx.hartree_in_kcalpermol() for e in opt_results["opt_energies"]
    ]

    fig, ax = plt.subplots(figsize=(6, 3))

    ax.plot(energies_in_kcalpermol, "o--")

    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"Energy relaxation (kcal/mol)")

    plt.show()

def get_scf(molecule, xcfun='B3LYP', basis='6-31G*', max_iter=200, restart=False):
    basis = vlx.MolecularBasis.read(molecule, basis, ostream=None)

    if molecule.get_multiplicity() == 1:
        scf_drv = vlx.ScfRestrictedDriver()
    else :
        scf_drv = vlx.ScfUnrestrictedDriver()

    if xcfun in vlx.available_functionals():
        scf_drv.xcfun = xcfun

    scf_drv.max_iter = max_iter
    scf_drv.checkpoint_file = 'scf.chkp'

    if 'scf.chkp' in os.listdir() and restart:
        scf_drv.restart = True

    print('Starting SCF computation...')
    scf_results = scf_drv.compute(molecule, basis)

    print('DFT finished.')

    return scf_results



def get_RESP(molecule, charge, mult, scf_results, basis='6-31G*', xcfun='B3LYP'):
    print("Starting RESP fitting...")
    
    basis = vlx.MolecularBasis.read(molecule, basis, ostream=None)

    resp_drv = vlx.RespChargesDriver()
    resp_drv.net_charge = charge
    resp_drv.multiplicity = mult

    if xcfun in vlx.available_functionals() and xcfun != 'HF':
        resp_drv.xcfun = xcfun

    #resp_drv.update_settings({"equal_charges": "1 = 3, 1 = 4"})

    resp_charges = resp_drv.compute(molecule, basis, scf_results, "resp")

    respf =  open('RESP_charges.txt', 'w')

    print("Atom num\tAtom     RESP charge")
    print(20 * "-")

    respf.write("Atom num\tAtom     RESP charge"+'\n')
    respf.write(20 * "-"+'\n')

    for n, (label, resp_charge) in enumerate(zip(molecule.get_labels(), resp_charges)):
        print(f"{n+1}\t\t{label :s} {resp_charge : 18.6f}")
        respf.write(f"{n+1}\t\t{label :s} {resp_charge : 18.6f}\n")

    print(20 * "-")
    respf.write(20 * "-"+'\n')

    print(f"Total: {resp_charges.sum() : 13.6f}")
    print(f"Total: {resp_charges.sum() : 13.6f}\n")


def main():

    #xyz_in, charge, mult = parser()
    args = argparser()
    input = args.input
    charge = args.charge
    mult   = args.multiplicity
    m_opt  = args.method_opt
    b_opt  = args.basis_opt
    m_resp = args.method_resp
    b_resp = args.basis_resp
    #non_equiv = args.non_equivalent_atomsll

    #vlx.environment.set_omp_num_threads(str(args.num_cores))

    if 'OMP_NUM_THREADS' in list(os.environ):
        if args.num_cores != int(os.environ['OMP_NUM_THREADS']):
            print(f"Warning: OMP_NUM_THREADS is set to {os.environ['OMP_NUM_THREADS']}, but you are trying to set it to {args.num_cores}. {args.num_cores} will be used instead.")
            os.system("export OMP_NUM_THREADS={args.num_cores}")
    else :
        if args.num_cores != vlx.environment.cpu_count():
            os.system("export OMP_NUM_THREADS={args.num_cores}")

    print(args.opt_output)

    molecule = load_molecule(input, charge, mult)

    if not args.no_opt:
        opt_results = optimize(args.opt_output, molecule, m_opt, b_opt, restart=args.restart)
        if args.plot_opt:
            plot_scf_along_opt(opt_results)

        molecule = load_molecule(args.opt_output, charge, mult)

    scf_results = get_scf(molecule, m_resp, b_resp)

    get_RESP(molecule, charge, mult, scf_results)


if __name__ == '__main__':
    main()
