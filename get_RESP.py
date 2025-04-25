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
    parser.add_argument('-c', '--charge', default=0, type=int, help='System\'s charge')
    parser.add_argument('-m', '--multiplicity', default=1, type=int, help='System\'s multiplicity')
    parser.add_argument('-m_opt',  '--method_opt',  type=str, default='B3LYP', choices=['HF'] + list(vlx.available_functionals()), help='QM method for optimisation')
    parser.add_argument('-m_resp', '--method_resp', type=str, default='B3LYP', choices=['HF'] + list(vlx.available_functionals()), help='QM method for RESP calculation')
    parser.add_argument('-b_opt',  '--basis_opt',  type=str, default='6-31G*', help='QM basis for optimisation')
    parser.add_argument('-b_resp', '--basis_resp',  type=str, default='6-31G*', help='QM basis for RESP calculation')

    parser.add_argument('-no_opt', '--no_opt', default=False, action='store_true', help='Deactivate geometry optimisation')
    parser.add_argument('--plot_opt', default=False, action='store_true', help='Plot SCF energy along the optimisation')
    parser.add_argument('-r', '--restart', default=False, action='store_true', help='Restart SCF calculation from checkpoint file')

    parser.add_argument('-n', '--num_cores', type=int, default=1, help='Number of cores to use for the calculation')

    #not implemented yet
    #parser.add_argument('--non_equivalent_atoms', type=bool, default=False, action='store_false', help='Deactivate automatic search of equivalent atoms for RESP fitting')

    args = parser.parse_args()

    if args.opt_output == None:
        args.opt_output = args.input.split('.')[0] + '_opt.' + args.input.split('.')[1]

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


def optimize(output, molecule, xcfun='B3LYP', basis='6-31G*', max_iter=200, restart=False):

    basis = vlx.MolecularBasis.read(molecule, basis, ostream=None)
    
    if molecule.get_multiplicity() == 1:
        scf_drv = vlx.ScfRestrictedDriver()
    else :
        scf_drv = vlx.ScfUnrestrictedDriver()

    if xcfun in vlx.available_functionals():
        scf_drv.xcfun = xcfun

    scf_drv.max_iter = max_iter
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

    print("Atom     RESP charge")
    print(20 * "-")

    for label, resp_charge in zip(molecule.get_labels(), resp_charges):
        print(f"{label :s} {resp_charge : 18.6f}")

    print(20 * "-")

    print(f"Total: {resp_charges.sum() : 13.6f}")


def main():

    #xyz_in, charge, mult = parser()
    args = argparser()
    input = args.input
    charge = args.charge
    mult   = args.multiplicity
    m_opt  = args.method_opt
    b_opt  = args.basis_opt
    #m_resp = args.method_resp
    #b_resp = args.basis_resp
    #non_equiv = args.non_equivalent_atomsll

    #vlx.environment.set_omp_num_threads(str(args.num_cores))

    if args.num_cores != vlx.environment.cpu_count():
        os.system("export OMP_NUM_THREADS={args.num_cores}")

    print(args.opt_output)

    molecule = load_molecule(input, charge, mult)

    if not args.no_opt:
        opt_results = optimize(args.opt_output, molecule, m_opt, b_opt, restart=args.restart)
        if args.plot_opt:
            plot_scf_along_opt(opt_results)

        molecule = load_molecule(args.opt_output, charge, mult)

    scf_results = get_scf(molecule, m_opt, b_opt)

    get_RESP(molecule, charge, mult, scf_results)


if __name__ == '__main__':
    main()
