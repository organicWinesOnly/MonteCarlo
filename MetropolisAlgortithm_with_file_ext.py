"""Implementation of the metropolis algorithm
"""


import matplotlib.pyplot as mp
from lattice import *
# from Errors import Error
from helpers_mc import *
# import scipy.optimize as sp
# import csv
import pandas as pd
from ast import literal_eval as l_e


# ===========================================================================
# Constants
# ===========================================================================
B = 0
J = 1
# k = 1
m = -3.05765101e-03
y_int = 1.12324380e+01

# ===========================================================================
# Simulation
# ===========================================================================


class Simulation:
    """Class to run the algorithm

    ==== Attributes ===
    system: All the latices to be simulated at varying temperatures
    temp: List containg the temperatures the system is run at
    magnetization:
    eq_energy:List containing List of the calculated energies of the system
                during the simulation
    eq_magnetization: List containing List of the calculated magnetization
                        of the system during the simulation
    ind_temp: List of list containing the independent energy of the simulation
    ind_magnetization: List of list containing the independent magnetization of
                        the simulation
    energy: list of energies throughout the simulation
    num_rounds: number of times to run the simulation
    time: time of simulation ran in monte carlo steps (number should == num
            rounds)
    correlation_times: list containing the measured correlation times
    n: number of spins in the system
    """
    system: List[Union[Lattice, Lattice2D, Plotts]]
    temp: List[Union[int, float]]
    energy: List[List[float]]
    magnetization: List[List[float]]
    time: int
    error_type: int
    eq_energy: List[List[float]]
    eq_magnetization: List[List[float]]
    correlation_times: List[float]
    ind_energy: List[List[float]]
    ind_magnetization: List[List[float]]
    num_rounds: int
    n: int

    def __init__(self, config):
        """Initialize class.

        range_of_temps: 2 item list. [start temperature, highest temp]
        step: temperature step
        """
        iterations = np.arange(config['range'][0], config['range'][1],
                               config['step'])
        if isinstance(config['lattice_size'], tuple):
            self.system = [Lattice2D(config['lattice_size'], config['temp'], p)
                           for p in iterations]
        else:
            self.system = [Lattice(config['lattice_size'], config['temp'], p)
                           for p in iterations]

        self.stats = {}
        self.keys = []
        for y in range(len(self.system)):
            q = {}
            q['eq_energy'] = []
            q['eq_magnetization'] = []
            q['correlation_times'] = []
            k = self.stats
            k['{0:.1f}'.format(iterations[y])] = q
            self.keys.append('{0:.1f}'.format(iterations[y]))

        self.temp = list(map(lambda x: float(x), list(self.stats.keys())))  # temperature of the system
        self.stats1 = {}
        for o in range(len(self.system)):
            u = {}
            u['indep_energy'] = []
            u['indep_mag'] = []
            h = self.stats1
            h['{0:.1f}'.format(iterations[o])] = u
            self.keys.append('{0:.1f}'.format(iterations[o]))
        self.time = 0
        self.num_rounds = config['num_rounds']
        self.correlation_times = []
        self.n = config['lattice_size']

    def run(self):
        """Run the simulation"""
        raise NotImplementedError

    def plot_system(self):
        for i in range(len(self.system)):
            lattice = self.stats[self.keys[i]]
            mp.figure()
            times_ = np.arange(self.time)
            mp.plot(times_, lattice['eq_energy'], 'b',
                    label='Internal Energy')
            mp.plot(times_, lattice['eq_magnetization'], 'r',
                    label='Magnetization')
            mp.title('Energy and Magnetization vs Time for {}K'.format(
                self.temp[i]))
            mp.legend()
            mp.savefig("startt0_{}.png".format(self.keys[i]))

    def indep_measurments(self, equil_time: Union[int, List]) -> None:
        """Caclulates the independent energy and magnetization for each lattice in the system the
        system"""
        index = 0

        for temp in self.stats:
            # latt in self.eq_energy:
            num = int(self.correlation_times[index])
            i = 0
            indep_eng = []
            indep_mag = []
            if isinstance(equil_time, int):
                while i * num < len(self.stats[temp]['eq_energy']) - equil_time:
                    indep_eng.append(self.stats[temp]['eq_magnetization'][equil_time + i
                                                              * num])
                    i += 1
            else:
                while i * num < len([self.stats[temp]['eq_energy']]) - equil_time[index]:
                    indep_eng.append(self.stats[temp]['eq_energy'][equil_time[index] + i
                                                       * num])
                    print(i * num)
                    indep_mag.append(
                        self.stats[temp]['eq_magnetization'][equil_time[index] + i * num])
                    i += 1
                index += 1
                self.stats1[temp]['indep_energy'] = indep_eng
                self.stats1[temp]['indep_mag'] = indep_mag

    def internal_energy(self) -> List[float]:
        """Calculate the internal energy"""

        lst_int_e = []

        for temp in self.stats:
            sum_ = (sum(temp['indep_energy']) / len(temp['indep_energy']) /
                    self.n)
            lst_int_e.append(sum_)
        return lst_int_e

    def tot_magnetization(self) -> List[float]:
        """Calculate the internal energy"""

        lst_int_m = []

        for temp in self.stats:
            sum_ = (sum(temp['indep_magnetization']) /
                    len(temp['indep_magnetization'])) / self.n
            lst_int_m.append(sum_)
        return lst_int_m

    def specific_heat(self) -> np.array:
        rms_ = []
        i = 0
        for temp in self.stats:
            rms_.append(rms(temp['indep_energy']) * self.system[i].beta ** 2 /
                        self.n)
            i += 1

        return rms_

    def suscepitbility(self) -> np.array:
        rms_ = []
        xo = np.zeros(self.system[0].size_tot)
        xo.fill(1)
        for temp in self.stats:
            rms_.append(rms(temp['indep_magnetization']) /
                        self.system[0].size_tot)

        return rms_

    def plot_property(self):
        """Plot property"""
        pass

    def correlation_time(self, equil_time: Union[int, List]):
        """Calculates the correlation time for each lattice.
        The time retuned is rounded to the nearest int"""
        if isinstance(equil_time, int):

            for i in range(len(self.system)):
                lattice = self.stats[self.keys[i]]
                zero_ = autocorrelation(self.time-equil_time, 0,
                                        lattice['eq_magnetization'][equil_time:])
                if zero_ == 0:
                    zero_ = 0.1
                times = np.arange(self.time - 1 - equil_time)
                integrand_ = []
                integrand_csv = []
                for t in times:
                    to_add = autocorrelation(self.time-equil_time, t,
                                             lattice['eq_magnetization']
                                             [equil_time:]) / zero_
                    to_add_ = to_add * zero_
                    integrand_.append(to_add)
                    integrand_csv.append(to_add_)
                lattice['auto_correlation_times'] = integrand_csv

                x = abs(round(sum(integrand_)))
                if x == 0:
                    x = 1
                lattice['correlation_time'] = x
                self.correlation_times.append(x)


        if isinstance(equil_time, list):
            for i in range(len(self.system)):
                lattice = self.stats[self.keys[i]]
                print(i)
                zero_ = autocorrelation(self.time-equil_time[i], 0,
                                        lattice['eq_magnetization']
                                        [equil_time[i]:])
                if zero_ == 0:
                    zero_ = 0.1
                times = np.arange(self.time - 1 - equil_time[i])
                integrand_ = []
                integrand_csv = []
                for t in times:
                    to_add = autocorrelation(self.time-equil_time[i], t,
                                             lattice['eq_magnetization']
                                             [equil_time[i]:]) / zero_
                    integrand_.append(to_add)
                    to_add_ = to_add * zero_
                    integrand_csv.append(to_add_)
                x = abs(round(sum(integrand_)))
                lattice['auto_correlation_times'] = integrand_csv

                if x == 0:
                    x = 1
                lattice['correlation_time'] = x
                self.correlation_times.append(x)


class MetropolisAlgorithm(Simulation):
    """Metropolis Algorithm"""

    def run(self):
        i = 0
        keys = list(self.stats.keys())

        for lattice_ in self.system:
            magnetization_per_lattice = []
            energy_per_lattice = []
            for _ in range(self.num_rounds):
                # generate a new lattice
                state_v = self.generate(lattice_)

                # run the algorithm and return the appropiate lattice
                self.update_lattice(lattice_, state_v)
                magnetization_per_lattice.append(lattice_.magnetization)
                energy_per_lattice.append(lattice_.energy)
                self.time += 1
            self.stats[keys[i]]['eq_magnetization'].extend(magnetization_per_lattice)
            self.stats[keys[i]]['eq_energy'].extend(energy_per_lattice)
            i += 1
        self.time /= len(self.system)
        self.time = int(self.time)
        # self.write_file()

    def write_f(self):
        for temp in self.stats:
            u = pd.DataFrame.from_dict(self.stats1[temp])
            u.T.to_csv('1d_ma_in{}.csv'.format(temp))

    def write_file(self):
        for temp in self.stats:
            u = pd.DataFrame.from_dict(self.stats[temp])
            u.T.to_csv('1d_ma{}.csv'.format(temp))

    def generate(self, lattice: Lattice):
        """generate a new state v by randomly flipping 1 spin"""
        k = r.randint(0, lattice.size - 1)  # index of spin to be flipped
        state_v = lattice.__copy__()
        spin_k = state_v.spins[k] * -1
        state_v.spins[k] = spin_k
        state_v.magnetization = sum(state_v.spins)

        state_v.energy = lattice.energy + 2 * J * spin_k * sum(
            lattice._neighbours(k))
        return state_v

    def _acceptance(self, original: Lattice, other: Lattice) -> float:
        """Return the acceptance probability according to the metropolis
        algorithm w/ single spin dynamic."""
        if other.energy - original.energy > 0:
            return np.exp(-1 * original.beta * (other.energy - original.energy))
        else:
            return 1

    def _change_state(self, original: Lattice, other: Lattice) -> None:
        """change the state to that of other
        """
        original.spins = other.spins
        original.magnetization = other.magnetization
        original.size = other.size
        original.energy = other.energy

    def update_lattice(self, original: Lattice, other: Lattice) -> None:
        """update the lattice according to the acceptance probability"""
        number = r.uniform(0, 1)
        while number == 1:
            number = r.uniform(0, 1)

        accept = self._acceptance(original, other)

        if number < accept:
            self._change_state(original, other)

# ===========================================================================
# Implementation
# ===========================================================================
# after 6000 points


if __name__ == '__main__':
    config_ = {'range': [10e-3, 3.02], 'step': 0.1, 'lattice_size': 355,
               'error_method': 12, 'temp': 1, 'num_rounds': 10000}

    sim = MetropolisAlgorithm(config_)
    sim.run()
    sim.plot_system()
    # equil_time = [2000, 2000, 5000, 7000, 6000, 8000, 8000, 6000, 4000, 5000,
    #               8000, 8000, 6000, 8500, 6000, 6000, 6000, 8000, 6002, 6000, 6000,
    #               6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000]
    # sim.correlation_time(equil_time)
    # sim.indep_measurments(equil_time)
    # sim.write_f()
