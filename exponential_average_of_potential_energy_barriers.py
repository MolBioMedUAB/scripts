#! /Users/mikicanyelles/.pyenv/shims/python
from scipy import constants
import  argparse
from math import exp, log

R = constants.R/4184

def sum(values):

    sum = 0
    for v in values:
        sum += float(v)

    return sum


def calculate_exponential_average(temp, barriers, R=R):

    exp_factors = []
    for barrier in barriers:
        exp_factors.append(
                exp(
                    (-float(barrier)/
                        (R*temp))
                    )
                )

    exp_avg = -R*temp*log((1/len(barriers))*sum(exp_factors))

    return exp_avg

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='exponential_average',
        description='Python script for calculating exponential averages of barriers'
    )

    parser.add_argument(
        '-T', '--temperature',
        default=298,
        type=float,
        help='Option to input the temperature (in K) for obtaining the exponential average energy barrier. By default the value is 298 K.'
    )

    parser.add_argument(
        #'-b', '--barriers',
        'barriers',
        metavar='barriers',
        action='store',
        type=str,
        nargs='+',
        help='List of single energy barriers (in kcal/mol)',
        #required=True
    )


    args = vars(parser.parse_args())
    if len(args['barriers']) == 1 and ' ' in args['barriers']:
        args['barriers'] = args['barriers'][0].split()


    print('The exponential barrier at %s K is:' % args['temperature'], round(calculate_exponential_average(args['temperature'], args['barriers']), 2), 'kcal/mol')
