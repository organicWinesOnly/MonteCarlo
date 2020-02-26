import numpy as np
from typing import List


def autocorrelation(maxtime: int, time: int, magnetization: List):
    assert time != maxtime
    diff = maxtime - time
    coeff = 1 / diff
    times = np.arange(diff)
    sum1 = sum(magnetization[t_prime] * magnetization[t_prime + time] for
               t_prime in times)
    sum2 = sum(magnetization[t_prime] for t_prime in times)
    sum3 = sum(magnetization[t_prime + time] for t_prime in times)
    return coeff * sum1 - (sum2 * sum3) * coeff ** 2


def correlation_time(maxtime, magnetization):
    zero_ = autocorrelation(maxtime, 0, magnetization)
    times = np.arange(maxtime - 1)
    integrand_ = []
    for t in times:
        integrand_.append(autocorrelation(maxtime, t, magnetization))

    return sum(integrand_)


