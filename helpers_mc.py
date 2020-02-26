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
# N = 100


def autocorrelation(maxtime: int, time: int, magnetization: List) -> float:
    """Calculates the magnetic auto correlation function (3.21)

    @:param maxtime: the total amount of time my system is in equilibrium
    @:param time: time interested in
    @:param: magnetization

    ===Representation Invariants ====
    0 <= time < maxtime: a value of 0 repersents the time at equilibrium
    """
    assert time != maxtime
    maxtime_ = int(maxtime)
    time_ = int(time)
    diff = maxtime_ - time_
    coeff = 1 / diff
    times = np.arange(diff)
    sum1 = sum(magnetization[t_prime] * magnetization[t_prime + time_] for
               t_prime in times)
    sum2 = sum(magnetization[t_prime] for t_prime in times)
    sum3 = sum(magnetization[t_prime + time_] for t_prime in times)
    return coeff * sum1 - (sum2 * sum3) * coeff ** 2


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
    if np.mean(arr_sq) >= np.mean(arr) ** 2:
        return np.sqrt(np.mean(arr_sq) - np.mean(arr) ** 2)


def two_pnt_spec(input_: Union[List, np.array], beta):
    """Returns the two point correlation function (standard deviation
    of a set of data.

    used for specfic heat C only

    to find the specific heat per lattice site divide by N, number of lattice
    sites

    Uses fromula rms = (<Q^2> - <Q>^2) * k * beta ** 2
    """
    if isinstance(input_, list):
        arr = np.array(input_.copy())
    else:
        arr = input_.copy()

    arr_sq = arr.copy() ** 2

    return (np.mean(arr_sq) - np.mean(arr) ** 2) * (beta ** 2)


def two_pnt_sus(input_: Union[List, np.array], beta):
    """Returns the two point correlation function (standard deviation
    of a set of data.

    used for sus χ only

    to find the sus per lattice site:
     if using M divide by N, number of lattice sites
     if using m multiply by N
    Uses fromula rms = (<Q^2> - <Q>^2) * k * beta ** 2
    """
    if isinstance(input_, list):
        arr = np.array(input_.copy())
    else:
        arr = input_.copy()

    arr_sq = arr.copy() ** 2

    return (np.mean(arr_sq) - np.mean(arr) ** 2) * beta


def u(beta):
    beta_ = np.array(beta)
    return (-np.tanh(beta_ * J)) / (np.tanh(1 / 0.4 * J))
