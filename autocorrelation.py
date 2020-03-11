import numpy as np
from typing import List


def autocorrelation(maxtime: int, time: int, quanity: np.ndarray):
    """ time displaced autocorrelation.

    === Paramaters ===
    maxtime: number of rounds the simulation is ran for
    time: the time that the autocorrelation function will be evaluated at
    magnetizaion: the meausred quanities that are being evaluated.
    """
    assert time != maxtime
    
    diff = maxtime - time

    coeff = 1 / diff
    assert quanity[:diff].shape == quanity[time:maxtime].shape
    sum1 = np.sum(np.dot(quanity[:diff], quanity[time:maxtime]))
    #sum2 = np.sum(quanity[:diff])
    #sum3 = np.sum(quanity[time:maxtime])

    #return sum1 * ( coeff - coeff ** 2)
    return sum1 * (coeff - coeff ** 2 * sum1) # * sum2 * sum3


def ic_time(maxtime, quanity):
    """ calculate the integrated correlation time.

    This time will be used to find independent values of measured quatities

    === Paramaters ===
    maxtime: number of rounds the simulation is ran for
    magnetizaion: the meausred quanities that are being evaluated.
    """
    zero_ = autocorrelation(maxtime, 0, quanity)

    a, b = 0, maxtime-1

    h = 1

    i_ = h * (1/2 * autocorrelation(maxtime, a, quanity) + 1/2 *
            autocorrelation(maxtime, b, quanity)) / zero_

    sum_ = np.sum(autocorrelation(maxtime, k, quanity) for k in range(a, b))

    i_ += h * sum_ / zero_


    return i_


