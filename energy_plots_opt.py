import matplotlib.pyplot as plt
from os import listdir, environ
from sys import argv,exit
import plotext.plot as plx
from termcolor import colored

def input_file():
    while True:
        try :
            if argv[1] in listdir():
                file_name = argv[1]
                break
            elif argv[1] not in listdir():
                print("This file doesn't exist.")
                file_name = input('Which is the name of the .out file?')
                break
        except IndexError: #len(argv[1]) == 0 or len(argv[2]) == 0:
            #file_name = input('Which is the name of the .out file?')
            #break
            print('.out file is not specified. Please, specify it as an argument. (> energy_plotter file_name.out)')
            exit(0)

    return file_name

def parse_file(file_name):
    file = open(file_name, 'r')

    cyc = 0
    energy = [[],[], [], []]

    for line in file:
        if line.find('Contribution to energy from') != -1:
            if list(line.split())[4] == 'turbomole:' or list(line.split())[4] == 'gaussian:':
                energy[2].append(float(list(line.split())[-2]))

            elif list(line.split())[4] == 'dl_poly:':
                energy[3].append(float(list(line.split())[-2]))

        elif line.find('QM/MM Energy:') != -1:
            cyc += 1
            energy[0].append(cyc)
            energy[1].append(float(list(line.split())[-2]))

    while True:
        try :
            energy[1][len(energy[0])-1]
            break
        except IndexError:
            del energy[0][-1]
            break

    return energy


def write_csv(energy):
    csv = open('%s_energies.csv' % file_name, 'w+')
    csv.write("Cycle, QM/MM energy, QM energy, MM energy\n")
    for n in range(0,len(energy[0])):
        csv.write(str(energy[0][n])+','+str(energy[1][n])+','+str(energy[2][n])+','+str(energy[3][n])+'\n')
    csv.close()
    print("CSV file saved.")


def plot_energy_mtl(energy):
    plt.plot(energy[0], energy[1])
    plt.xlabel('Number of cycles')
    plt.ylabel('QM/MM Potential Energy (Hartree)')
    plt.legend(["QM/MM\nPotential\nEnergy"], bbox_to_anchor=(1.02,1), loc = 2)
    plt.grid(True)

    plt.title('QM/MM energies along the optimisation', y=1.05, loc='center')
    plt.savefig('QMMM_energies.png', transparent=False, dpi=300, bbox_inches='tight')
    plt.show()

    print("Plot for QM/MM energies saved")

    fig, ax1 = plt.subplots()
    ax1.plot(energy[0], energy[2], color='blue')#, marker='o')
    ax1.set_xlabel('Number of cycles')
    ax1.set_ylabel('QM Potential energy (Hartree)', color='blue')
    ax1.tick_params(axis='y', colors='blue')
    #ax1.yaxis.grid(linestyle='-', linewidth='0.5', color='blue')
    ax1.legend(["QM\nenergy"], bbox_to_anchor=(1.20,1), loc = 2)

    ax2=ax1.twinx()
    #ax2.plot(csv_ar[:,2], csv_ar[:,-1], color='red')
    ax2.plot(energy[0], energy[3], color='red')#, marker='^')
    #ax2.set_xlabel('O-H distance (Å)')
    #for k in range(1,i+1):
    #    ax2.plot(range(globals()['len_%s' % (k-1)]+1,globals()['len_%s' % k]+1), globals()['table_%s' % k][3], linestyle='--')
    ax2.set_ylabel('MM Potential energy (Hartree)', color='red')
    #ax2.set_yticks(mm_ticks)
    ax2.tick_params(axis='y', colors='red')
    ax2.legend(["MM\nenergy"], bbox_to_anchor=(1.20,0.90), loc = 2)
    #ax2.yaxis.grid(linestyle='-.', linewidth='0.5', color='red')

    plt.xticks()
    plt.title('QM and MM energies along the optimisation', y=1.05, loc='center')
    plt.savefig('QM_vs_MM_energies.png', transparent=False, dpi=300, bbox_inches='tight')
    plt.show()


def plot_energy_plx(energy):
    plx.scatter(energy[0], energy[1])
    plx.show()
    print(colored(str(len(energy[0])) +  ' points have been calculated so far', 'red'))

def main():
    file_name = input_file()
    energy = parse_file(file_name)
    try :
        if argv[2] == 'csv':
            write_csv(energy)

    except IndexError :
        pass

    try :
        environ['DISPLAY']
        plot_energy_mtl(energy)
    except KeyError:
        plot_energy_plx(energy)



if __name__ == "__main__":
    main()





