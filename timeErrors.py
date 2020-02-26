"""Perform run time experiment on methods in error class for comparison purposes

Based off of CSC158 timequeue file from fall 2018
"""

from timeit import timeit
from typing import List, Tuple

# from myqueue import Queue
from Errors import Error
import matplotlib.pyplot as plt
import numpy as np
import random as r

###############################################################################
# Task 3: Running timing experiments
###############################################################################


def mean_sq(vari):
    if isinstance(vari, list):
        arr = np.array(vari.copy())
    else:
        arr = vari.copy()
    return np.mean(arr ** 2)


func = mean_sq


def _setup_errors(qsize: int, n: int) -> List[Error]:
    """Return a list of <n> Errors, each of the given size."""
    # Experiment preparation: make a list containing <n> Errors,
    # each containng a random data set of size <qsize>.
    #
    #     queue._items = list(range(qsize))

    _errors = []
    for i in range(n):
        data_ = [r.uniform(1, 400) for _ in range(qsize)]
        error = Error(func, data_, True, None)
        _errors.append(error)
    return _errors


def time_error_lists() -> Tuple[List[int], List[float], List[float],
                                List[float]]:
    """Run timing experiments for Error.blocking, Error.Jackknife and
    Error.bootstrap.

    Return lists storing the results of the experiments.
    """
    # data_sizes = [10000, 20000, 40000, 80000, 160000]
    # data_sizes = [100, 200, 400, 800, 1600]
    data_sizes = [10, 20, 40, 80, 160]
    # The number of times to call a single enqueue or dequeue operation.
    trials = 200

    # list of times
    time_bl = []  # blocking
    values_bl = []
    time_jk = []  # jackknife
    values_jk = []
    time_bs = []  # bootstrap
    values_bs = []

    # This loop runs the timing experiment.
    for data_size in data_sizes:
        errors = _setup_errors(data_size, trials)
        time = 0
        for error in errors:
            time += timeit('error.blocking_method(5)', number=1, globals=locals())
        time_bl.append(time)

    for data_size in data_sizes:
        errors = _setup_errors(data_size, trials)
        time = 0
        for error in errors:
            time += timeit('error.jackknife_method()', number=1, globals=locals())
        time_jk.append(time)

    for data_size in data_sizes:
        errors = _setup_errors(data_size, trials)
        time = 0
        for error in errors:
            time += timeit('error.bootstrap_method()', number=1, globals=locals())
        time_bs.append(time)

    tuple_ = (data_sizes, time_bl, time_jk, time_bs)
    return tuple_


if __name__ == '__main__':
    # time_queue()
    g = time_error_lists()
    plt.plot(g[0], g[1], label='Blocking Method')
    plt.plot(g[0], g[2], label='Jackknife Method')
    plt.plot(g[0], g[3], label='Bootstrap Method')
    plt.legend()
    plt.title('Comparison of error run times on Random Data Sets of varying '
              'size for a mean squared function, running 200 trials for each'
              ' size')
    # plt.show()
    plt.savefig('error_anal.png')
    plt.figure()


    # bargrapgh of the mean errors of the plot

