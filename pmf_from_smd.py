import numpy as np
import matplotlib.pyplot as plt

import argparse

plt.style.use('ggplot')


def parser():

    parser = argparse.ArgumentParser(description='Script for obtaining PMF profiles for multi-frame steered MDs simulated with AMBER.')

    """
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--option1", help="Description for option 1")
    group.add_argument("--option2", help="Description for option 2")
    """

    parser.add_argument('-i', '--input',
                        help='Input .dat files from Amber SMD simulations.',
                        required=True,
                        nargs='+')

    parser.add_argument('-t', '--temperature',
                        help='Working temperature for calculating the PMF. 310.15 K (37ºC) is the default value.',
                        required=False,
                        default=310.15,
                        type=float,
                        )

    parser.add_argument('-n', '--num_of_cv',
                        help='Number of used collective variables',
                        required=False,
                        default=1,
                        type=int,
                        )
    
    parser.add_argument('-st', '--step',
                        help='Step for avoiding points in input files',
                        required=False,
                        default=1,
                        type=int,
                        )
    
    parser.add_argument('-cv', '--collective_variable',
                        help='Collective variable to use for Force calculation',
                        required=False,
                        default=1,
                        type=int,
                        )
    
    parser.add_argument('-pw', '--print_total_works',
                        help='Reads the last lines of the .dat file and recovers the total work done.',
                        default='highest',
                        choices=['all', 'highest', 'no'],
                        required=False
                        )
    
    parser.add_argument('-rt', '--reset_time',
                        help='Resets the time printed in the input file to 0',
                        default=True,
                        type=bool,
                        required=False
                        )

    parser.add_argument('-o', '--output',
                        help='Output files name. Default is \'pmf\'.',
                        default='pmf')

    parser.add_argument('-s', '--show',
                        help='Trigger for previsualise the generated plots',
                        default=False,
                        action='store_true')

    return parser.parse_args()

def read_input_files(input):

    smds = []
    for file in input:

        f = open(file).readlines()[4:-4] # lines 1-3 are the headers and -3-last is the total work done

        for l in range(len(f)):
            f[l] = f[l].split()

        smds.append(f)

    return smds

def read_total_works(input, print_total_works):

    if print_total_works in ('all', 'highest'):
        total_works = {}
        for file in input:
            f = open(file).readlines()[-2]
            total_works[file] = float(f.split(' ')[-1])

        total_works = sorted(total_works.items(), key=lambda x:x[1], reverse=True)
        
        print(f'\nThe profile {total_works[0][0]} has the highest work, {round(total_works[0][1], 2)} kcal/mol')

        if print_total_works in ('all'):
            print('The rest of the total works, sorted by highest to lowest, are below:')
            for i in range(1, len(total_works)):
                print(f'\tThe file {total_works[i][0]} has a total work of {round(total_works[i][1], 2)} kcal/mol')

        return


    elif print_total_works in ('no'):
        return

    
def calculate_1D_2n_order_cumulant(smds, temp=310.15, n_cv=1, cv=1, reset_time=True, step=1, output='pmf'):

    beta =1/(0.001987 * temp)

    if len(smds) == 1:
        print('No cumulant can be calculated since only one calculation has been provided.')
        print('The plot will be generated with work instead of PMF')
        output_type = 'w'

    else :
        output_type='pmf'

    smd_analysed = []
    works        = []

    f_out = open(output + '_processed.csv', 'w')

    if n_cv == 1:
        cv1_col = 1
        handle1_col = 2
        k1_col      = 3
        w_col       = 4

    elif n_cv == 2:
        if cv == 1:
            cv1_col     = 1
            handle1_col = 3
            k1_col      = 5

        elif cv == 2:
            cv1_col     = 2
            handle1_col = 4
            k1_col      = 6
        
        #cv1_col = 1
        #cv2_col = 2
        #handle1_col = 3
        #handle2_col = 4
        #k1_col      = 5
        #k2_col      = 6
        
        w_col       = 7
    

    if reset_time:
        time_step = round(float(smds[0][1][0])/100,4) - round(float(smds[0][0][0])/100,4)

    for l in range(0, len(smds[0]), step):

        tot_W = 0
        tot_W_sq = 0
        tot_F = 0
        tot_CV = 0
        works.append([])
        for smd in smds:
            tot_W    += float(smd[l][w_col])
            tot_W_sq += float(smd[l][w_col])**2
            tot_F    += float(smd[l][k1_col])*(float(smd[l][cv1_col])-float(smd[l][handle1_col]))
            tot_CV  +=  float(smd[l][cv1_col])
            works[-1].append(smd[l][w_col])


        avg_W = (1/len(smds))*tot_W
        avg_W_sq = (1/len(smds))*tot_W_sq
        avg_F = (1/len(smds))*tot_F
        avg_CV = (1/len(smds))*tot_CV

        cum = avg_W - 1/2*beta*(avg_W_sq-(avg_W**2))

        if reset_time:
            smd_analysed.append(
                [time_step*l, round(avg_CV, 2), round(avg_F, 2), round(avg_W, 2), round(avg_W_sq, 2), round(cum, 2)]
            )
        else :
            smd_analysed.append(
                [round(float(smd[l][0])/100,4), round(avg_CV, 2), round(avg_F, 2), round(avg_W, 2), round(avg_W_sq, 2), round(cum, 2)]
            )

        f_out.write(', '.join([str(v) for v in smd_analysed[-1]]) + '\n')

    f_out.close()

    return smd_analysed, output_type, works

def calculate_2D_2n_order_cumulant(smds, temp=310.15, reset_time=True, step=1, output='pmf'):

    beta =1/(0.001987 * temp)

    if len(smds) == 1:
        print('No cumulant can be calculated since only one calculation has been provided.')
        print('The plot will be generated with work instead of PMF')
        output_type = 'w'

    else :
        output_type='pmf'

    smd_analysed = []
    works        = []

    f_out = open(output + '_2D_processed.csv', 'w')
        
    cv1_col = 1
    cv2_col = 2
    handle1_col = 3
    handle2_col = 4
    k1_col      = 5
    k2_col      = 6
    w_col       = 7
    

    if reset_time:
        time_step = round(float(smds[0][1][0])/100,4) - round(float(smds[0][0][0])/100,4)

    for l in range(0, len(smds[0]), step):

        tot_W = 0
        tot_W_sq = 0
        tot_F = 0
        tot_CV1 = 0
        tot_CV2 = 0
        works.append([])
        for smd in smds:
            tot_W    += float(smd[l][w_col])
            tot_W_sq += float(smd[l][w_col])**2
            tot_F    += (float(smd[l][k1_col])*(float(smd[l][cv1_col])-float(smd[l][handle1_col]))+float(smd[l][k2_col])*(float(smd[l][cv2_col])-float(smd[l][handle2_col])))/2
            tot_CV1  +=  float(smd[l][cv1_col])
            tot_CV2  +=  float(smd[l][cv2_col])
            works[-1].append(smd[l][w_col])


        avg_W    = (1/len(smds))*tot_W
        avg_W_sq = (1/len(smds))*tot_W_sq
        avg_F    = (1/len(smds))*tot_F
        avg_CV1  = (1/len(smds))*tot_CV1
        avg_CV2  = (1/len(smds))*tot_CV2

        cum = avg_W - 1/2*beta*(avg_W_sq-(avg_W**2))

        if reset_time:
            smd_analysed.append(
                [time_step*l, round(avg_CV1, 2), round(avg_CV2, 2), round(avg_F, 2), round(avg_W, 2), round(avg_W_sq, 2), round(cum, 2)]
            )
        else :
            smd_analysed.append(
                [round(float(smd[l][0])/100,4), round(avg_CV1, 2), round(avg_CV2, 2), round(avg_F, 2), round(avg_W, 2), round(avg_W_sq, 2), round(cum, 2)]
            )

        f_out.write(', '.join([str(v) for v in smd_analysed[-1]]) + '\n')

    f_out.close()

    return smd_analysed, output_type, works


def plot_1D_timewise_pmf(smd_analysed, output='pmf', show=True, output_type='pmf', works=None):

    smd_analysed = np.array(smd_analysed)

    # Plotting force
    plt.plot(smd_analysed[:, 0], smd_analysed[:, 2])

    plt.title('Time-wise Force')
    plt.ylabel('Force')
    plt.xlabel('Time (ns)')

    plt.savefig(output + '_time-wise_force_plot.png', dpi=300)

    if show:
        plt.show()
    plt.close()

    # Plotting Work/PMF
    plt.plot(smd_analysed[:, 0], smd_analysed[:, 5])

    if output_type == 'pmf':
        plt.title('Time-wise PMF')
        plt.ylabel('PMF (kcal/mol)')
        plt.xlabel('Time (ns)')

        plt.savefig(output + '_1D_time-wise_plot_free_energy.png', dpi=300)
        if show:
            plt.show()
        plt.close()

        if works != None:
            works = np.array(works)

            for smd in range(len(works[0])):
                plt.plot(smd_analysed[:, 0], works[:, smd])

            plt.title('Time-wise W')
            plt.ylabel('W (kcal/mol)')
            plt.xlabel('Time (ns)')

            plt.savefig(output + '_1D_time-wise_plot_work.png', dpi=300)
            if show:
                plt.show()
            plt.close()

    if output_type == 'w':
        plt.title('Time-wise W')
        plt.ylabel('W (kcal/mol)')
        plt.xlabel('Time (ns)')

        plt.savefig(output + '_1D_time-wise_plot_work.png', dpi=300)
        if show :
            plt.show()
        plt.close()

def plot_1D_cvwise_pmf(smd_analysed, output='pmf', show=True, output_type='pmf', works=None):

    smd_analysed = np.array(smd_analysed)

    # Plotting force
    plt.plot(smd_analysed[:, 1], smd_analysed[:, 2])

    plt.title('CV-wise Force')
    plt.ylabel('Force')
    plt.xlabel('CV (Å)')

    plt.savefig(output + '_1D_CV-wise_force_plot.png', dpi=300)

    if show:
        plt.show()
    plt.close()

    # Plotting Work/PMF
    plt.plot(smd_analysed[:, 1], smd_analysed[:, 5])

    if output_type == 'pmf':
        plt.title('CV-wise PMF')
        plt.ylabel('PMF (kcal/mol)')
        plt.xlabel('CV (Å)')

        plt.savefig(output + '_1D_CV-wise_plot_free_energy.png', dpi=300)
        if show:
            plt.show()
        plt.close()

        if works != None:
            works = np.array(works)

            for smd in range(len(works[0])):
                plt.plot(smd_analysed[:, 1], works[:, smd])

            plt.title('CV-wise W')
            plt.ylabel('W (kcal/mol)')
            plt.xlabel('CV (Å)')

            plt.savefig(output + '_1D_CV-wise_plot_work.png', dpi=300)
            if show:
                plt.show()
            plt.close()

    if output_type == 'w':
        plt.title('Time-wise W')
        plt.ylabel('W (kcal/mol)')
        plt.xlabel('CV (Å)')

        plt.savefig(output + '_1D_CV-wise_plot_work.png', dpi=300)
        if show :
            plt.show()
        plt.close()


def plot_cv1_vs_cv2(smd_analysed, output='pmf', show=True, output_type='pmf', works=None):

    smd_analysed = np.array(smd_analysed)

    # Plotting CVs
    fig, ax = plt.subplots(1,1)

    ax.scatter(smd_analysed[:, 1], smd_analysed[:, 2])
    htmp = ax.tricontourf(smd_analysed[:, 1], smd_analysed[:, 2], smd_analysed[:, 6], levels=100, cmap='gist_rainbow')
    cbar = fig.colorbar(htmp, ax=ax, )

    if output_type == 'pmf':
        cbar.set_label('PMF (kcal/mol)')
        plt.title('2D CV-wise PMF')
        plt.ylabel('CV2 (Å)')
        plt.xlabel('CV1 (Å)')

        plt.savefig(output + '_2D_CV-wise_plot_free_energy.png', dpi=300)
        if show:
            plt.show()
        plt.close()

        if works != None:
            works = np.array(works)

            for smd in range(len(works[0])):
                plt.plot(smd_analysed[:, 1], works[:, smd])

            cbar.set_label('W (kcal/mol)')
            plt.title('2D CV-wise W')
            plt.ylabel('CV2 (Å)')
            plt.xlabel('CV1 (Å)')

            plt.savefig(output + '_1D_CV-wise_plot_work.png', dpi=300)
            if show:
                plt.show()
            plt.close()

    if output_type == 'w':
        cbar.set_label('W (kcal/mol)')
        plt.title('Time-wise W')
        plt.ylabel('CV2 (Å)')
        plt.xlabel('CV1 (Å)')

        plt.savefig(output + '_2D_CV-wise_plot_work.png', dpi=300)
        if show :
            plt.show()
        plt.close()



def main():

    args = parser()

    smds = read_input_files(args.input)
    
    read_total_works(args.input, args.print_total_works)

    smd_analysed_1D, output_type_1D, works_1D = calculate_1D_2n_order_cumulant(smds, temp=args.temperature, n_cv=args.num_of_cv, cv=args.collective_variable, reset_time=args.reset_time, step=args.step, output=args.output)
    plot_1D_timewise_pmf(smd_analysed_1D, output=args.output, show=args.show, output_type=output_type_1D, works=None)
    plot_1D_cvwise_pmf(smd_analysed_1D, output=args.output, show=args.show, output_type=output_type_1D, works=None)

    if args.num_of_cv:
        smd_analysed_2D, output_type_2D, works_1D = calculate_2D_2n_order_cumulant(smds, temp=args.temperature, reset_time=args.reset_time, step=args.step, output=args.output)
        plot_cv1_vs_cv2(smd_analysed_2D, output=args.output, show=args.show, output_type=output_type_1D, works=None)



    smd = np.array(smd_analysed_1D)

    if output_type_1D == 'pmf':
        print('Maximum PMF (kcal/mol):', np.max(smd[:,5]))
    elif output_type_1D == 'w':
        print('Maximum W (kcal/mol):', np.max(smd[:,5]))


if __name__ == '__main__':
    main()

















