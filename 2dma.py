"""run metropolis algorithm for 2d lattice
"""

import matplotlib.pyplot as mp
from lattice import *
from helpers_mc import *
import csv

N = 64
J = 1


class Simulation:
    """Class to run the algorithm

    ==== Attributes ===
    system: All the latices to be simulated at varaying temperatures
    energy: list of energies throughout the simulation
    magnetizationL list of magnetization values throughout the simulation
    time: int representing the number of cycles in the simulation that have
    occured
    error_type:
    """
    system: List[Union[Lattice, Lattice2D, Plotts]]
    energy: List[List[float]]
    magnetization: List[List[float]]
    time: int
    error_type: int
    eq_energy: List[List[float]]
    eq_magnetization: List[List[float]]

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
        self.temp = iterations
        self.magnetization = [lattice.magnetization for lattice in self.system]
        self.energy = [lattice.energy for lattice in self.system]
        self.time = 0
        self.eq_magnetization = []
        self.eq_energy = []
        self.num_rounds = config['num_rounds']
        self.correlation_times = []
        self.ind_energy = []
        self.ind_mag = []

    def run(self):
        """Run the simulation"""
        raise NotImplementedError

    def plot_system(self):
        for lat_num in range(len(self.system)):
            mp.figure()
            times_ = np.arange(self.time)
            mp.plot(times_, self.eq_energy[lat_num], 'b',
                    label='Internal Energy')
            mp.plot(times_, self.eq_magnetization[lat_num], 'r',
                    label='Magnetization')
            # auto_corr = self.auto_correlation(times_)
            # mp.semilogy(times_, auto_corr)
            mp.title('Energy and Magnetization vs Time for {}K'.format(
                self.temp[lat_num]))
            mp.legend()
            mp.show()
            # mp.savefig('temp {}.png'.format(self.temp[lat_num]))

    def indep_energy(self, equil_time: Union[int, List]) -> None:
        """Calulates the independent energy for each lattice in the system the
        system"""
        index = 0
        indep_energy = []

        for latt in self.eq_energy:
            num = int(self.correlation_times[index])
            i = 0
            indep_measurements = []
            if isinstance(equil_time, int):
                while i * num < len(latt) - equil_time:
                    indep_measurements.append(latt[equil_time + i * num])
                    i += 1
                indep_energy.append(indep_measurements)
                index += 1
            else:
                while i * num < len(latt) - equil_time[i]:
                    indep_measurements.append(latt[equil_time[i] + i * num])
                    i += 1
                indep_energy.append(indep_measurements)
                index += 1
        self.ind_energy = indep_energy

    def indep_magnetization(self, equil_time: Union[int, List]):
        """Calculates indep magnetization"""
        index = 0
        indep_mag = []

        for latt in self.eq_magnetization:
            num = int(self.correlation_times[index])
            i = 0
            indep_measurements = []
            if isinstance(equil_time, int):
                while i * num < len(latt) - equil_time:
                    indep_measurements.append(latt[equil_time + i * num])
                    i += 1
                indep_mag.append(indep_measurements)
                index += 1
            else:
                while i * num < len(latt) - equil_time[i]:
                    indep_measurements.append(latt[equil_time[i] + i * num])
                    i += 1
                indep_mag.append(indep_measurements)
                index += 1
        self.ind_mag = indep_mag

    def internal_energy(self) -> List[float]:
        """Calculate the internal energy"""

        indep_energy = self.ind_energy.copy()
        lst_int_e = []

        for i in range(len(indep_energy)):
            sum_ = (sum(indep_energy[i]) / len(indep_energy[i])) / N
            lst_int_e.append(sum_)
        return lst_int_e

    def tot_magnetization(self) -> List[float]:
        """Calculate the internal energy"""

        indep_mag = self.ind_mag.copy()
        lst_int_m = []

        for i in range(len(indep_mag)):
            sum_ = (sum(indep_mag[i]) / len(indep_mag[i])) / N
            lst_int_m.append(sum_)
        return lst_int_m

    # def magnetization_(self, start: int) -> float:
    #     """Calculate the total magnetization"""
    #     return sum(self.energy[start:start+1000]) / (1000 * N)

    def specific_heat(self) -> np.array:
        rms_ = []
        for i in range(len(self.ind_energy)):
            rms_.append(two_pnt_spec(self.ind_energy[i], self.system[i].beta) / N)

        return rms_

    def specific_heat1(self) -> np.array:
        rms_ = []
        for i in range(len(self.ind_energy)):
            rms_.append(two_pnt_sus(self.ind_energy[i], self.system[i].beta) * N)

        return rms_
    #
    # def suscepitbility(self) -> np.array:
    #     rms_ = []
    #     xo = np.zeros(self.system[0].size_tot)
    #     xo.fill(1)
    #     for i in range(len(self.ind_mag)):
    #         rms_.append(rms(self.ind_mag[i]) / self.system[0].size_tot)
    #
    #     return rms_

    def correlation_time(self, equil_time: Union[int, List]):
        """Calculates the correlation time for each lattice.
        The time retuned is rounded to the nearest int"""
        if isinstance(equil_time, int):
            for i in range(len(self.system)):
                zero_ = autocorrelation(self.time-equil_time, 0,
                                        self.eq_magnetization[i][equil_time:])
                if zero_ == 0:
                    zero_ = 0.1
                times = np.arange(self.time - 1 - equil_time)
                integrand_ = []
                for t in times:
                    to_add = autocorrelation(self.time-equil_time, t,
                                             self.eq_magnetization[i]
                                             [equil_time:]) / zero_
                    integrand_.append(to_add)
                x = abs(round(sum(integrand_)))
                if x == 0:
                    x = 1
                self.correlation_times.append(x)

class MetropolisAlgorithm2d(Simulation):
    """Metropolis Algorithm"""

    def run(self):
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
            self.eq_magnetization.append(magnetization_per_lattice)
            self.eq_energy.append(energy_per_lattice)
        self.time /= len(self.system)
        self.time = int(self.time)
        # self.correlation_time()

    def generate(self, lattice: Lattice2D):
        """generate a new state v by randomly flipping 1 spin"""
        k = r.randint(0, lattice.size[0] - 1)  # index of spin to be flipped
        q = r.randint(0, lattice.size[1] - 1)
        state_v = lattice.__copy__()
        spin_k = state_v.spins[k][q] * -1
        state_v.spins[k][q] = spin_k
        state_v.magnetization = sum(sum(state_v.spins))

        state_v.energy = lattice.energy + 2 * J * spin_k * sum(
            lattice._neighbours((k, q)))
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

    def correlation_time(self, equil_time: Union[int, List]):
        """Calculates the correlation time for each lattice.
        The time retuned is rounded to the nearest int"""
        if isinstance(equil_time, int):
            for i in range(len(self.system)):
                iter = 0
                print("lattice: {}".format(self.system[i].temp))
                zero_ = autocorrelation(self.time-equil_time, 0,
                                        self.eq_magnetization[i][equil_time:])
                if zero_ == 0:
                    zero_ = 0.1
                times = np.arange(self.time - 1 - equil_time)
                integrand_ = []
                integrand_csv = []
                for t in times:
                    print("lattice: {}".format(self.system[i].temp), iter)
                    iter += 1
                    to_add = autocorrelation(self.time-equil_time, t,
                                             self.eq_magnetization[i]
                                             [equil_time:]) / zero_
                    to_add_ = to_add * zero_
                    integrand_.append(to_add)
                    integrand_csv.append(to_add_)
                myfile = open('corr{0:.2f}.csv'.format(self.temp[i]), 'w+')

                writer = csv.writer(myfile)
                writer.writerows(map(lambda x: [x], integrand_csv))

                x = abs(round(sum(integrand_)))
                if x == 0:
                    x = 1
                self.correlation_times.append(x)

        if isinstance(equil_time, list):
            for i in range(len(self.system)):
                zero_ = autocorrelation(self.time-equil_time[i], 0,
                                        self.eq_magnetization[i]
                                        [equil_time[i]:])
                if zero_ == 0:
                    zero_ = 0.1
                times = np.arange(self.time - 1 - equil_time[i])
                integrand_ = []
                integrand_csv = []
                for t in times:
                    to_add = autocorrelation(self.time-equil_time[i], t,
                                             self.eq_magnetization[i]
                                             [equil_time[i]:]) / zero_
                    integrand_.append(to_add)
                    to_add_ = to_add * zero_
                    integrand_csv.append(to_add_)
                x = abs(round(sum(integrand_)))

                myfile = open('corr{}.csv'.format(self.temp[i]), 'w+')

                writer = csv.writer(myfile)
                writer.writerows(map(lambda x: [x], integrand_csv))

                if x == 0:
                    x = 1
                self.correlation_times.append(x)


config_ = {'size': N, 'range': [1.0, 2.8], 'step': 0.1,
           'lattice_size': (8, 8), 'error_method': 12,
           'temp': 0, 'num_rounds': 10000, 'num_states': 2}
#
# value = []
# timess = np.arange(4000, 5999)
# with open('corr2.40.csv') as f:
#     reader = csv.reader(f)
#     for line in reader:
#         value.extend(line)
#
# val_prime = list(map(lambda x: float(x), value))
#
# mp.plot(val_prime, '.')
# mp.show()

# sim = MetropolisAlgorithm2d(config_)
# sim.run()
# sim.correlation_time(10000)
