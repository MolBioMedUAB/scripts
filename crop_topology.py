# -*- coding: utf-8 -*-

# Import packages
from parmed import load_file
from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError
import sys, os
from argparse import ArgumentParser

#def parser_func():
parser = ArgumentParser(description="Simple script based on parmed and MDAnalysis for cropping AMBER topologies for ChemShell QM/MM calculations.")

parser.add_argument(
    '-p', '--parameters',
    help='File containing the parameters and topology of the system in the AMBER format. Accepted file formats: prmtop.\nIt needs the coordinates file.',
    required=False
    )
parser.add_argument(
    '-c', '--coordinates',
    help='File containing the coordinates of the system in AMBER format. Accepted file formats: inpcrd, rst7.\nIt needs the paramters file.',
    required=False
    )
parser.add_argument(
    '-pdb',
    help='File containing the coordinates of the system in the PDB format. If parameters are not given, ionly the generation of the active atoms list will be available.',
    required=False
    )
parser.add_argument(
    '-o', '--output',
    help='Name for the output file (the cropped parameters and coordinates).',
    required=False,
    )

argsdict = vars(parser.parse_args())

if argsdict['coordinates'] != None and argsdict['pdb'] != None:
    print('Both coordinates and pdb have been inputed. Only coordinates from the pdb will be used.')
    argsdict['coordinates'] = None

 #   return argsdict


def options_parser(argsdict=argsdict):
    '''
    This function takes a the dictionary created using the argparse module and
    returns a dictionary with the available options. The argsdict is also output
    because it may has suffered some modifications.
    The input dictionary has to include the following keys: 'parameters',
    'coordinates', 'pdb' and 'output'.
    The output dictionary has to include the following keys: 'crop topology',
    'active list', 'parameters adaption'.
    '''

    options      = {'crop parameters' : False, 'active list' : False, 'parameters adaption' : False, 'crop pdb' : False}
    options_done = {'crop parameters' : False, 'active list' : False, 'parameters adaption' : False, 'crop pdb' : False}

    if argsdict['parameters'] != None and (argsdict['coordinates'] != None or argsdict['pdb'] != None):
        options['crop parameters']     = True
        options['active list']         = True
        options['parameters adaption'] = True

    elif argsdict['parameters'] != None and argsdict['coordinates'] == None and argsdict['pdb'] == None:
        options['parameters adaption'] = True

    elif (argsdict['parameters'] == None and argsdict['coordinates'] == None) and argsdict['pdb'] != None:
        options['active list']         = True
        options['crop pdb']            = True

    elif argsdict['parameters'] == None and argsdict['coordinates'] == None and argsdict['pdb'] == None:
        print('No filename has been input, so there is nothing to do.')
        sys.exit(0)


    return options, options_done


def crop_top(argsdict=argsdict):
    '''
    This functions takes the input parameters and coordinates/pdb and crops both
    taking a residue as the centre and specifying a radius around it.
    It saves the cropped parameters, coordinates (in AMBER format) and pdb files
    and names it with the output string or with the filename of the parameters,
    adding always the '.cropped' string.
    It returns the value of the selected radius, so it can be used as threshold
    for the active_atoms_list function.
    '''

    if argsdict['parameters'] != None:
        parameters = argsdict['parameters']

    if argsdict['coordinates'] != None:
        coordinates = argsdict['coordinates']
    if argsdict['pdb'] != None:
        coordinates = argsdict['pdb']

    if argsdict['output'] != None:
        output = argsdict['output']
    elif argsdict['output'] == None:
        output = parameters[:-7]

    if str(output + '.cropped.prmtop') in os.listdir() or str(output + '.cropped.inpcrd') in os.listdir() or str(output + '.cropped.pdb') in os.listdir():
        while True:
            quest = input('The files exist. Do you want to overwrite it ([y]/n)?')

            if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
                if str(output + '.cropped.prmtop') in os.listdir():
                    os.remove(str(output + '.cropped.prmtop'))
                if str(output + '.cropped.inpcrd') in os.listdir():
                    os.remove(str(output + '.cropped.inpcrd'))
                if str(output + '.cropped.pdb') in os.listdir():
                    os.remove(str(output + '.cropped.pdb'))
                break
            elif quest in ('no', 'NO', 'No', 'nO', 'n', 'N', 0):
                print('The cropped parameters cannot be saved. Rerun the script specifiying an output name or let overwrite the file.')
                quit()
                break
            else :
                print('Please, answer \'yes\' or \'no\'.')
                continue



    u_top = Universe(parameters)

    ligand = input('Type the number or the three letters code (only if it is a non-standard residue) of the central residue: ')
    try:
#        if type(ligand) is int:
        sel = str(u_top.select_atoms('resid %s or resname %s' % (ligand, ligand)).residues)
#        elif type(ligand) is str:
#            sel = str(u_top.select_atoms('resname %s' % ligand).residues)

        print('The selected residue is %s' % sel[24:-3])
        err = False

    except SelectionError:
        print('The selected residue does not exist. Please, reselect it.')
        err = True

    while True:
        if err == True:
            quest = 'no'
        else :
            quest = input('Is the residue correct ([y]/n)? ')

        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
            break
        elif quest in ('no', 'NO', 'No', 'nO', 'n', 'N', 0):
            if err != True:
                print('Please, reselect the ligand then.')

            err = False
            ligand = input('Type the number or the three letters code (only if it is a non-standard residue) of the central residue: ')
            try:
                if type(ligand) is int:
                    sel = str(u_top.select_atoms('resid %s' % ligand).residues)
                if type(ligand) is str:
                    sel = str(u_top.select_atoms('resname %s' % ligand).residues)
                    print('The selected residue is %s' % sel[24:-3])

            except SelectionError:
                print('The selected residue does not exist. Please, reselect it.')
                err = True

            continue
        else :
            print('Type only yes or no')
            continue


    err = False
    radius = input('Type the desired radius around the selected ligand for the water drop (in ang): ')
    try:
        float(radius)
        print('The selected radius is %s ang' % radius)

    except ValueError:
        print('Type a number')
        err = True

    while True:
        if err == True:
            quest = 'no'
        else :
            quest = input('Is this correct ([y]/n)? ')

        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
            break
        elif ('no', 'NO', 'No', 'nO', 0):
            print('Type the number again.')

            radius = input('Type the desired radius around the selected ligand for the water drop (in ang): ')
            try:
                float(radius)
                print('The selected radius is %s ang' % radius)
                err = False

            except ValueError:
                print('Type a number')
                err = True
            continue
        else :
            print('I did not understand the answer. Please, answer again.')
            continue

    del u_top, sel, quest, err

    topology = load_file(parameters, coordinates)

    topology.box = None
    topology.strip(':WAT&!:%s<:%s' % (ligand, radius))
    topology.strip(':Na+&!:%s<:%s' % (ligand, radius))
    topology.strip(':Cl-&!:%s<:%s' % (ligand, radius))

    topology.write_parm(output + '.cropped.prmtop')
    topology.save(output + '.cropped.inpcrd')
    #for res in topology.residues:
    #   res.ter = False
    topology.write_pdb(output + '.cropped.pdb')

    pdb = open(output + '.cropped.pdb', 'r').readlines()
    pdb_o = open(output + '.cropped.pdb.tmp', 'w')

    atom = int(pdb[0][6:11])
    for line in pdb:
        if 'TER' in line:
            pass
        else :
            new_line = line[0:6] + '{:5}' + line[11:]
            pdb_o.write(new_line.format(atom))
            atom += 1

    pdb_o.close()
    os.remove(output + '.cropped.pdb')
    os.rename(output + '.cropped.pdb.tmp', output + '.cropped.pdb')

    print('Cropped topology and coordinates have been saved as \'%s\', \'%s\' and \'%s\'' % (output + '.cropped.prmtop', output + '.cropped.inpcrd', output + '.cropped.pdb'))

    return radius


def crop_pdb(argsdict=argsdict):
    '''
    This functions takes the input pdb and crops it
    taking a residue as the centre and specifying a radius around it.
    It saves the cropped pdb file and names it with the output string
    or with the filename of the parameters, adding always the '.cropped' string.
    It returns the value of the selected radius, so it can be used as threshold
    for the active_atoms_list function.
    '''

    if argsdict['pdb'] != None:
        pdb = argsdict['pdb']


    if argsdict['output'] != None:
        output = argsdict['output']
    elif argsdict['output'] == None:
        output = pdb[:-4]

    if str(output + '.cropped.pdb') in os.listdir():
        while True:
            quest = input('The file exist. Do you want to overwrite it ([y]/n)?')

            if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
                if str(output + '.cropped.pdb') in os.listdir():
                    os.remove(str(output + '.cropped.pdb'))
                break
            elif quest in ('no', 'NO', 'No', 'nO', 'n', 'N', 0):
                print('The cropped parameters cannot be saved. Rerun the script specifiying an output name or let overwrite the file.')
                quit()
                break
            else :
                print('Please, answer \'yes\' or \'no\'.')
                continue

    u_pdb = Universe(pdb)

    ligand = input('Type the number or the three letters code (only if it is a non-standard residue) of the central residue: ')
    try:
        sel = str(u_pdb.select_atoms('resid %s or resname %s' % (ligand, ligand)).residues)

        print('The selected residue is %s' % sel[24:-3])
        err = False

    except SelectionError:
        print('The selected residue does not exist. Please, reselect it.')
        err = True

    while True:
        if err == True:
            quest = 'no'
        else :
            quest = input('Is the residue correct ([y]/n)? ')

        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
            break
        elif quest in ('no', 'NO', 'No', 'nO', 'n', 'N', 0):
            if err != True:
                print('Please, reselect the ligand then.')

            err = False
            ligand = input('Type the number or the three letters code (only if it is a non-standard residue) of the central residue: ')
            try:
                if type(ligand) is int:
                    sel = str(u_pdb.select_atoms('resid %s' % ligand).residues)
                if type(ligand) is str:
                    sel = str(u_pdb.select_atoms('resname %s' % ligand).residues)
                    print('The selected residue is %s' % sel[24:-3])

            except SelectionError:
                print('The selected residue does not exist. Please, reselect it.')
                err = True

            continue
        else :
            print('Type only yes or no')
            continue


    err = False
    radius = input('Type the desired radius around the selected ligand for the water drop (in ang): ')
    try:
        float(radius)
        print('The selected radius is %s ang' % radius)

    except ValueError:
        print('Type a number')
        err = True

    while True:
        if err == True:
            quest = 'no'
        else :
            quest = input('Is this correct ([y]/n)? ')

        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
            break
        elif ('no', 'NO', 'No', 'nO', 0):
            print('Type the number again.')

            radius = input('Type the desired radius around the selected ligand for the water drop (in ang): ')
            try:
                float(radius)
                print('The selected radius is %s ang' % radius)
                err = False

            except ValueError:
                print('Type a number')
                err = True
            continue
        else :
            print('I did not understand the answer. Please, answer again.')
            continue

    import sys; sys.setrecursionlimit(15000)
    cropped_ats = u_pdb.select_atoms('protein or (around %s (resid %s or resname %s)) or (resid %s or resname %s)' % (radius, ligand, ligand, ligand, ligand))
    #print(cropped_ats)

    select_resids = 'resid ' + ' or resid '.join([str(r) for r in list(set(cropped_ats.resids))])
    #print(select_resids)
    cropped = u_pdb.select_atoms(select_resids)

    cropped.write(str(output + '.cropped.pdb'))

    return radius

radius = 10000
def active_atoms_list(argsdict=argsdict,radius=radius):
    '''
    This function takes the input pdb or parameters+coordinates and saves a list
    in tcl format (it includes the 'set act' definition of the list name and the
    list itself). This list will be then useful for ChemShell calculations.
    This list contains the index of a selection of atoms made specifying a central
    atom (by its number index) and the radius around it. This radius is compared
    to the radius selected in the crop_top function if it exists or to a very
    big radius (10000) if it doesn't, so any radius canb be accepted.
    It prints the total number of selectred atoms.
    '''

    if argsdict['pdb'] != None:
        u_set_act = Universe(argsdict['pdb'])
    elif argsdict['parameters'] != None and argsdict['coordinates'] != None:
        u_set_act = Universe(argsdict['parameters'], argsdict['coordinates'])

    while True:
        try :
            carbon = input("Type the number of the central atom: ")
            sel_carbon = str(u_set_act.select_atoms("bynum %s" % carbon))
            locA = sel_carbon.find('[<') + 2
            locB = sel_carbon.find(' and segid')
            break
        except SelectionError:
            print('The selection does not exists. Please, type an atom that exists.')
            continue


    while True:
       quest = input('You selected this atom: %s. Is it correct ([y]/n)? ' % sel_carbon[locA:locB])

       quest = str(quest).lower()

       if quest not in ('', 'y', 'yes', '1', 'n', 'no', '0'):
           print("Type just \'yes\' or \'no\'.")
           continue

       if quest in ('n', 'no', '0'):
           while True:
               try:
                   carbon = input("Type the number of the central carbon: ")
                   sel_carbon = str(u_set_act.select_atoms("bynum %s" % carbon))
                   locA = sel_carbon.find('[<') + 2
                   locB = sel_carbon.find(' and segid')
                   break
               except SelectionError:
                   continue
           continue

       elif quest in ('', 'y', 'yes', '1'):
           break

#        except ValueError:
#            print("Type just \'yes\' or \'no\'.")
#            continue

    while True:
        radius_set_act = input("Which is the desired radius (in ang)? ")

        try :
            radius_set_act = float(radius_set_act)

        except ValueError:
            print('Type a number, please.')
            continue

        if radius_set_act < float(radius):
            break

        elif radius_set_act >= float(radius):
            print('The selected radius for the set_act list is smaller than the radius of solvent. Please, choose a new radius for the active atoms.')
            continue


    if argsdict['output'] != None:
        output = str(argsdict['output']) + '_' +  str(radius_set_act) + '.act_list'
    elif argsdict['output'] == None:
        if argsdict['parameters'] != None:
            output = str(argsdict['parameters'])[:-7] + '.act_list'
        elif argsdict['pdb'] != None:
            output = str(argsdict['pdb'])[:-4] + '.act_list'

    selection = u_set_act.select_atoms(str('byres around %s bynum %s' % (radius_set_act, carbon)))

    txt = open(output, 'w')
    txt.write('set act { ')
    for i in range(0, len(selection)):
        locA = str(selection[i]).find('<Atom ') + 6
        locB = str(selection[i]).find(': ')
        txt.write(str(selection[i])[locA:locB] + " ")
    txt.write("} ")
    txt.close()

    print('%s atoms have been set as active for the QM/MM calculation using ChemShell.' % (len(selection)))


def topology_adapter(argsdict=argsdict):
    '''
    This function takes the input parameters and change some conflictive atom
    types on ChemShell. Moreover it changes the atom types created using MCPB
    which describe non-common atom types as well as metals and asks the user
    for the name of the atom.
    It loops over the atom types section of the file and replaces all the
    thought-to-be-conflictive atoms:
        - hc, ha, h1
        - 2C, 3C
        - CO, CX, c, c2, c3, cx, ce, cf
        - op, os, o
        - Na+, Cl-
    Into the non-conflictive-for-ChemShell types.
    '''

    top    = open(argsdict['parameters'], 'r').readlines()
    top_out = open(str(argsdict['parameters'])[:-7] + '.mod.prmtop', 'w')

    initial = 0
    heterotypes_in = []
    metals_in      = []
    for i in range(0, len(top)):
        if '%FLAG AMBER_ATOM_TYPE' in top[i]:
            initial = i +1
        if '%FLAG TREE_CHAIN_CLASSIFICATION' in top[i]:
            final   = i
        if ('  Y' in top[i] or top[i].find('Y') == 0) and '%' not in top[i]:
            index = 0
            while index < len(top[i]):
                index = top[i].find('Y', index)
                if index == -1:
                    break
                elif index != -1:
                    if top[i].find('Y') == 0:
                        loc = 0
                        heterotypes_in.append(str(top[i])[loc:loc+3])
                    else :
                        loc = top[i].find(' Y', index-1) + 1
                        #print(loc)
                        heterotypes_in.append(str(top[i])[loc:loc+3])

                    index += 1

        if ' M' in top[i] and '%' not in top[i]:
            loc = top[i].find(' M') + 1
            try :
                int(str(top[i])[loc+1:loc+3])
                metals_in.append(str(top[i])[loc:loc+3])
            except ValueError:
                pass

    heterotypes_out = []
    for i in range(len(heterotypes_in)):
        ht_ = input("Which atom is %s: " % heterotypes_in[i])
        heterotypes_out.append('{:<3}'.format(ht_))

    metals_out = []
    for i in range(len(metals_in)):
        m_ = input('Which atom is %s: ' % metals_in[i])
        metals_out.append('{:<3}'.format(m_))

    for l in range(0, initial):
        top_out.write(top[l])

    for l in range(initial, final):
        l_ = top[l].replace('hc', 'H ')
        l_ = l_.replace('ha', 'H ')
        l_ = l_.replace('h1', 'H ')

        l_ = l_.replace('2C', 'C2')
        l_ = l_.replace('3C', 'C3')

        l_ = l_.replace('CO', 'C ')
        l_ = l_.replace('CX', 'C ')
        l_ = l_.replace('c ', 'C ')
        l_ = l_.replace('c2', 'C ')
        l_ = l_.replace('c3', 'C ')
        l_ = l_.replace('cx', 'C ')
        l_ = l_.replace('ce', 'C ')
        l_ = l_.replace('cf', 'C ')

        l_ = l_.replace('op', 'O ')
        l_ = l_.replace('os', 'O ')
        l_ = l_.replace('o ', 'O ')

        l_ = l_.replace('Na+', 'NA+')
        l_ = l_.replace('Cl-', 'CL-')

        for j in range(len(heterotypes_out)):
            l_ = l_.replace(heterotypes_in[j], heterotypes_out[j])

        for j in range(len(metals_out)):
            l_ = l_.replace(metals_in[j], metals_out[j])

        top_out.write(l_)

    for l in range(final, len(top)):
        top_out.write(top[l])

    top_out.close()
    #top.close()


def main(argsdict=argsdict):
    '''
    This function is the main function of the program. It takes the input arguments, call the
    available functions (specified using the options_parser function), passes the proper argsdict
    (it may need to be modified in order to use the new files generated in the crop_top func.).
    This function will be only run if the program is execuded as it is.
    In the future the functions present here will be converted into functions of a package.
    '''

    options, options_done = options_parser(argsdict)

    if options['crop parameters'] == True:

        while True:
            quest = input('Do you want to crop the system ([y]/n)? ')

            if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                radius = crop_top()
                options_done['crop parameters'] = True
                break
            elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                break
            else :
                print('Type just \'yes\' or \'no\'.')
                continue

    if options['crop pdb'] == True:

        while True:
            quest = input('Do you want to crop the system ([y]/n)? ')

            if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                radius = crop_pdb()
                options_done['crop pdb'] = True
                break
            elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                break
            else :
                print('Type just \'yes\' or \'no\'.')
                continue

    if options['active list'] == True:
        while True:
            quest = input('Do you want to create the tcl list of the active atoms ([y]/n)? ')

            if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                radius = 100000
                active_atoms_list(argsdict, radius)
                options_done['active list'] = True
                break
            elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                break
            else :
                print('Type just \'yes\' or \'no\'.')
                continue


    if options['parameters adaption'] == True:
        while True:
            quest = input('Do you want to adapt the atom types of the parameters file to ChemShell ([y]/n)? ')

            if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                if options_done['crop parameters'] == False:
                    topology_adapter()

                elif options_done['crop parameters'] == True:
                    while options_done['crop parameters'] == True:
                        quest2 = input('Do you want to use the new topology (1) or the source one (2) ([1]/2)? ')

                        if quest2 in ('', '1'):
                            if argsdict['output'] == None:
                                argsdict['parameters'] = str(argsdict['parameters'])[:-7] + '.cropped.prmtop'
                            elif argsdict['output'] != None:
                                argsdict['parameters'] = argsdict['output'] + '.cropped.prmtop'

                            topology_adapter()
                            options_done['crop parameters'] == False
                            break
                        elif quest2 == '2':
                            topology_adapter()
                            options_done['crop parameters'] == False
                            break
                        else :
                            print('Type just \'1\' or \'2\'.')
                            options_done['crop parameters'] == True
                            continue
                break

            elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                break
            else :
                print('Type just \'yes\' or \'no\'.')
                continue

if __name__ == '__main__':
    main()
