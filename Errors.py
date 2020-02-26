"""Contains error class that computes the error of a measurment using bootstrap,
jackknife, or blocking methods

@author: Rhamel
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
    variable: a list of data entries that I need to find the error of.
    independence: are the data points independent? This determines
                  how the error should be calculated if possible
    param: Parameters of the function
    """

    function: TypeVar
    variable: List
    parameter: List
    indepedence: bool

    def __init__(self, f: TypeVar, info: Union[str, list], indepent: bool,
                 param: Union[None, list, float]) -> None:
        """Initialize Error class.
        """
        self.independence = indepent
        if isinstance(info, str):
            data_ = []
            with open(info) as csvfile:
                reader_ = csv.reader(csvfile)
                for line_ in reader_:
                    for i_ in range(len(line_)):
                        data_.append(float(line_[i_]))
            self.variable = np.array(data_)
        elif isinstance(info, list):
            data_ = info.copy()
            self.variable = np.array(data_)
        else:
            self.variable = info

        self.function = f
        self.param = param

    def jackknife_method(self) -> Union[float, list]:
        """Calculate the error using jackknife method
        """
        if self.independence:
            if self.param is None:
                q = self.function(self.variable)
            elif not isinstance(self.param, list):
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

    def bootstrap_method(self):
        """Calculate the error using bootstrap method

        Note: The lenght of the random subset should be equal to the number of
        indepnedent measurments made. For MC this is max time dividing by 2 *
        (correlation time)
        """
        created_mesurments = []
        for __ in range(100):
            random_subset = [self.variable[r.randint(0, len(self.variable) - 1)]
                             for _ in range(len(self.variable))]
            if self.param is None:
                created_mesurments.append(self.function(random_subset))
            else:
                created_mesurments.append(self.function(random_subset,
                                                        *self.param))
        return rms(created_mesurments)

    def blocking_method(self, number_of_blocks: int):
        """Calculate the error using blocking method
        """
        real_num_of_blocks = number_of_blocks
        while len(self.variable) % real_num_of_blocks != 0:
            real_num_of_blocks -= 1

        copy_data = self.variable.copy()
        size_of_block = int(len(self.variable) / real_num_of_blocks)
        blocks = [copy_data[(u * size_of_block): ((u+1) * size_of_block) - 1] for u in
                  range(real_num_of_blocks)]
        for k in range(len(blocks)):
            if self.param is None:
                blocks[k] = self.function(np.array(blocks[k]))
            else:
                blocks[k] = self.function(np.array(blocks[k]), *self.param)

        return np.sqrt(1 / real_num_of_blocks) * rms(blocks)


if __name__ == '__main__':

    e_ = Error(mean_sq, 'question3_2.csv', True, None)
    # e_.data = np.multiply(e_.data, e_.data)
    print(e_.jackknife_method())


    # # other part of the question
    # filename = 'question3_2.csv'
    # data = []
    # with open(filename) as csvfile:
    #     reader = csv.reader(csvfile)
    #     for line in reader:
    #         for i in range(len(line)):
    #             data.append(float(line[i]))
    #
    # data1 = np.array(data)  # first term of my expression under the radical
    # data2 = np.array(data)  # second term ------- " ----------
    #
    # data1 = data1 ** 4
    # mean_data1 = np.mean(data1)
    #
    # mean_data2 = np.mean(data2 ** 2)
    # mean_sq_d2 = mean_data2 ** 2
    #
    # var_ = np.sqrt( 1 / 11 * (mean_data1 - mean_sq_d2))
    # print(var_)
    # expect 5.05
    # blocking gives 4.98
    # bootstrap gives 2.98 -> for some reason everytime i run this the error
    #  changes highest i saw is 5.2, this ,kes semse since im taking a random
    # subset on a small set of data
    # jackknife gives 5.27
