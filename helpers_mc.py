"""Python file containing the commonly used helper functions when implementing
Monte Carlo Algorithms.

Used for PHY371, Monte Carlo Simulations in Statistical Physics

@author: Rhamel Roomes-Delpeache
last updated Nov 10, 2018"""

from typing import Union, List, Any
import numpy as np
import random as r

# constants
J = 1

# Functions
def average_e(lst_: List, N) -> Union[float, int]:
    """Returns the avg of the list <lst_>"""
    return sum(lst_) / (len(lst_) * J * N)


def average(lst_: List) -> Union[float, int]:
    """Returns the avg of the list <lst_>"""
    return sum(lst_) / len(lst_) / 100


def lin_function(x, a, b):
    """model linear function for curve fitting
    """
    return a * x + b


def flip_coin(event1: Any, event2: Any) -> Any:
    """Return either <event1> or <event2>

    Makes a random decision between 2 choices with equally weighted
    probabilities.
    """
    i = r.randint(0, 1)

    if i == 0:
        return event1
    elif i == 1:
        return event2


def mean_sq(vari):
    if isinstance(vari, list):
        arr = np.array(vari.copy())
    else:
        arr = vari.copy()
    return np.mean(arr ** 2)


def rms(input_: Union[List, np.array]):
    """Returns te root mean square for the given informations.

    Uses fromula rms = sqrt{<Q^2> - <Q>^2}
    """
    if isinstance(input_, list):
        arr = np.array(input_.copy())
    else:
        arr = input_.copy()

    arr_sq = arr.copy() ** 2
    return np.sqrt(np.mean(arr_sq) - np.mean(arr) ** 2)


def two_pnt(input_: np.array):
    """Returns the two point correlation function (standard deviation
    of a set of data.

    used for specfic heat C only

    to find the specific heat per lattice site divide by N, number of lattice
    sites

    Uses fromula rms = (<Q^2> - <Q>^2) * k * beta ** 2
    """
    arr = input_.copy()

    arr_sq = arr.copy() ** 2

    return (np.mean(arr_sq) - np.mean(arr) ** 2) 


def c(beta):
    """ specific heat per spin, 1d lattice
    """
    return beta**2 * 1 / np.cosh(beta) ** 2


def sus(beta):
    """ susceptibility per spin, 1d lattice
    """
    return beta * np.exp(2 * beta)


def u(beta):
    """ energy per spin, 1d lattice
    """
    return -np.tanh(beta)
