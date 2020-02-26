"""Comparison of the values estimated by the different methods in the Error
class

@author: Rhamel Roomes
Template: Provided by CSC148 (fall 2018) timequeues.py at the University of
Toronto"""


from typing import List, Tuple

from Errors import Error
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from helpers_mc import average_e, rms_spec
import random
###############################################################################
# Task 3: Running timing experiments
###############################################################################
str = 0.1
sto = 2.6
ste = 0.1


path_ = '~/Documents/PHY/phy372/metropolis/1d_ma_in/'

def value_error_lists() -> Tuple[List[float], List[float], List[float]]:
    """Run timing experiments for Error.blocking, Error.Jackknife and
    Error.bootstrap.

    Return lists storing the results of the experiments.
    """

    data_ = list(map(lambda x:
                     list(pd.read_csv(path_ + '1d_ma_in{0:.1f}.csv'.format(x), index_col=0).loc
                          ['indep_energy'])[1:201], np.arange(str, sto, ste)))

    # list of times
    value_bl = []  # blocking
    value_jk = []  # jackknife
    value_bs = []  # bootstrap

    errors = [Error(rms_spec, data, True, None) for data in data_]
    for error in errors:
        value_jk.append(error.jackknife_method())
        value_bs.append(error.bootstrap_method())
        # value_bl.append(error.blocking_method(10))

    tuple_ = (value_bl, value_jk, value_bs)
    return tuple_


if __name__ == '__main__':
    g = value_error_lists()
    print(g)
    methods = ['Blocking', 'Jackknife', 'Bootstrap']
    temp = np.arange(str, sto, ste)
    with plt.style.context('dark_background'):
        # plt.plot(temp, g[0], 'ro', label='blocking')
        plt.plot(temp, g[1], 'm^',  label='jackknife')
        plt.plot(temp, g[2], 'ws', label='bootstrap')
        plt.legend()
        plt.xlabel('Temperature')
        plt.ylabel('Error Value')
        plt.title('Error of Internal energy, 200 data points')
        plt.savefig('value_comp200.png')

