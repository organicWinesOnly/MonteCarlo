""" errors.py

Contains error class that computes the error of a measurment using bootstrap,
jackknife, or blocking methods
"""

import csv
import numpy as np
import random as r
from typing import TypeVar, Union, List
from helpers_mc import *


class Error:
    """Represents the error in a measurment

    ===== Attributes ====
    function: function that calculates the measurement
    variable: a list of data entries that are intrested in calculating the
               errors of. 
    independence: are the data points independent? This determines
                  how the error should be calculated if possible
    param: Parameters of the function
    """

    function: TypeVar
    variable: np.ndarray
    parameter: List
    indepedence: bool

    def __init__(self, f: TypeVar, info: np.ndarray, indepent: bool,
                 param=None: Union[list, float]) -> None:
        """Initialize Error class.
        """
        self.independence = indepent

        self.variable = info

        self.function = f
        self.param = param

    def jackknife_method(self) -> Union[float, list]:
        """Calculate the error using jackknife method
        """
        if self.independence:
            if self.param is None:
                q = self.function(self.variable)
            elif isinstance(self.param, float):
                q = self.function(self.variable, self.param)
            else:
                q = self.function(self.variable, *self.param)

            q_ = []

            data_ = self.variable.copy()

            p = len(data_)

            for i in range(p):
                item = data_[i]
                data_ = np.delete(data_, i)

                if self.param is None:
                    q_.append(self.function(data_))

                elif not isinstance(self.param, list):
                    q_ = self.function(self.variable, self.param)

                else:
                    q_.append(self.function(data_, *self.param))
                data_ = np.insert(data_, i, item)

            q_ = np.array(q_)
            sum_ = np.sum((q_ - q) ** 2)
            return np.sqrt(sum_)

    def bootstrap_method(self, num_indep: int):
        """Calculate the error using bootstrap method

        Note: The lenght of the random subset should be equal to the number of
        independent measurments made. For MC this is max time dividing by 2 *
        (correlation time)

        ==== paramater ====
        num_indep: Number of independent values in the data set
        """
        created_mesurments = np.zeros(100)

        for i in range(100):
            random_subset = [self.variable[r.randint(0, len(self.variable) - 1)]
                             for _ in range(num_indep)]
            if self.param is None:
                created_mesurments[i] = self.function(random_subset)
            else:
                created_mesurments[i] = self.function(random_subset,
                                                        *self.param)
        return rms(created_mesurments)
