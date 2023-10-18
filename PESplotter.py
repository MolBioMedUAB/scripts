import matplotlib.pyplot as plt
from os import listdir, environ
from sys import argv, exit
import plotext.plot as plx
from termcolor import colored
import numpy as np
import pandas as pd

def input_file():
    """
    DESCRIPTION
        Function for setting the file_name. It can be the default (PES.plt) or an input file. It is checked if the file is in the directory and the program is quitted if it isn't.
    """

    print(argv)
    if len(argv) == 1:
        if 'PES.plt' in listdir():
            file_name = 'PES.plt'

        else :
            print('No PES.plt file in the directory. Please, specify the filename as an argument')
            exit()

    elif len(argv) == 2:
        if argv[1] in listdir():
            file_name = argv[1]

        else :
            print(file_name, 'does not exist in the directory.')
            exit()

    return file_name


def parse_file(file_name):
    df = pd.DataFrame(pd.read_csv(file_name))

    x_key = []
    y_key = []
    for key in df.keys():
        if key.find('rc') != -1:
            x_key.append(key)
        elif key.find('deltaE') != -1:
            y_key.append(key)
        else :
            pass

    if len(x_key) > 1:
        print('Detected possible coordinates:')
        for key in range(len(x_key)):
            print(key + 1, '. ', x_key[key])

        while True:
            try :
                quest = int(input('Which coordinate do you want? '))

                if quest in range(1, len(x_key) + 1):
                    x_key = x_key[quest-1]
                    break

                else :
                    print('The selected coordinate does not exist.')
                    continue

            except ValueError:
                print('Type a number')
                continue

    elif len(x_key) == 1:
        x_key = x_key[0]

    elif len(x_key) == 0:
        print('No reaction coordinate has been detected. Please, select one of the following: ')
        for k in range(len(df.keys())):
            print(k+1, ': ', df.keys()[k])

        while True:
            try :
                quest = int(input('Select the number of the column containing the reaction coordinate: '))

                if quest in range(1, len(df.keys()) + 1):
                    x_key = x_key[quest-1]
                    break

                else :
                    print('The selected column does not exist.')
                    continue

            except ValueError:
                print('Type a number')
                continue

    x = list(df[x_key])


    if len(y_key) > 1:
        print('Detected possible relative energy:')
        for key in range(len(x_key)):
            print(key + 1, '. ', y_key[key])

        while True:
            try :
                quest = int(input('Which energy do you want? '))

                if quest in range(1, len(y_key) + 1):
                    y_key = y_key[quest-1]
                    break

                else :
                    print('The selected coordinate does not exist.')
                    continue

            except ValueError:
                print('Type a number')
                continue

    elif len(y_key) == 1:
        y_key = y_key[0]

    elif len(y_key) == 0:
        print('No relative energy has been detected. Please, select one of the following: ')
        for k in range(len(df.keys())):
            print(k+1, ': ', df.keys()[k])

        while True:
            try :
                quest = int(input('Select the number of the column containing the relative energy: '))

                if quest in range(1, len(y_key) + 1):
                    y_key = y_key[quest-1]
                    break

                else :
                    print('The selected column does not exist.')
                    continue

            except ValueError:
                print('Type a number')
                continue

    y = list(df[y_key])

    return x, y, x_key, y_key

def plot_energy_plx(x,y):
    plx.plot(x,y)
    plx.show()
    print(colored(str(len(x)) + ' points of the scan have been calculated so far.', 'red'))

def plot_energy_mtl(x, y, x_key, y_key):

    plt.plot(x, y, marker='.')

    plt.xlabel(x_key)
    plt.ylabel(y_key)
    plt.grid()

    plt.show()

    print(colored(str(len(x)) + ' points of the scan have been calculated so far.', 'red'))

def main():

    file_name = input_file()
    x, y, x_key, y_key = parse_file(file_name)

    try :
        environ['DISPLAY']
        #print('MATPLOTLIB')
        plot_energy_mtl(x, y, x_key, y_key)

    except KeyError:
        plot_energy_plx(x,y)


if __name__ == '__main__':
    main()


