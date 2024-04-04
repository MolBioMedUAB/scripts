
from argparse import ArgumentParser


def parser():
    parser = ArgumentParser(description="Script to fix AMBER and MCP atom types in topology file to be used in ChemShell.")

    parser.add_argument('-p', '--parameters',
                        help='AMBER prmtop file containing topology and parameters of the system.',
                        type=str,
                        required=True
                        )
  
    parser.add_argument('-ram', '--rename_atoms_MCPB',
                        help="List of atom names for changing MCPB names (like Y1, Y2, M1, etc.)",
                        nargs='+',
                        type=str,
                        default=None,
                        required=False
                        )

    parser.add_argument('-o', '--output',
                        help='Name of the output files.',
                        type=str,
                        default=None
                        )

    args = parser.parse_args()
  
    if args.rename_atoms_MCPB == None:
        print('MCPB atom types will be automatically adapted to ChemShell atom types.')

    if args.output == None:
        args.output = args.parameters[:-7]

    return args


def topology_adapter(args):

    new_types = []
    u = Universe(args.parameters)

    top    = open(args.parameters, 'r').readlines()
    top_out = open(str(args.parameters)[:-7] + '.mod.prmtop', 'w')

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
        
    #print(heterotypes_in)
    #print(metals_in)

    if args.rename_atoms_MCPB == None and (len(heterotypes_in) > 0 or len(metals_in) > 0):

        translator = {
            'NE'  : 'N ',
            'NZ'  : 'N ',
            'NE1' : 'N  ',
            'NE2' : 'N  ',
            'ND1' : 'N  ',
            'ND2' : 'N  ',
            'O'   : 'O',
            'OW'  : 'O ',
            'OD1' : 'O  ',
            'OD2' : 'O  ',
            'OE1' : 'O  ',
            'OE2' : 'O  ',
            'OXT' : 'O  ',
            'SG'  : 'S ',
            'SD'  : 'S ',
        }

        args.rename_atoms_MCPB = []
        for t in heterotypes_in:
     #       print(u.select_atoms(f"type {t}"))
            args.rename_atoms_MCPB.append(translator[u.select_atoms(f"type {t}").names[0]])

        for t in metals_in:
            args.rename_atoms_MCPB.append(u.select_atoms(f"type {t}").names[0])


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
        l_ = l_.replace('oh', 'O ')

        l_ = l_.replace('Na+', 'NA+')
        l_ = l_.replace('Cl-', 'CL-')
        l_ = l_.replace('K+', 'K+')
        l_ = l_.replace('Ca+', 'Ca+')

        for j in range(len(args.rename_atoms_MCPB)):
            l_ = l_.replace((heterotypes_in + metals_in)[j], args.rename_atoms_MCPB[j])

        top_out.write(l_)

    for l in range(final, len(top)):
        top_out.write(top[l])

    top_out.close()

    if args.rename_atoms_MCPB != None:
        print(f"{len(args.rename_atoms_MCPB)} atom types from MCPB have been changed. ")

def main():
  args = parser()

  topology_adapter(args)
  print('prmtop fixed!')
  
main()
