"""Implementation of the metropolis algorithm
"""


import matplotlib.pyplot as mp
from lattice import *
from Errors import Error
from helpers_mc import *
# import scipy.optimize as sp


# ===========================================================================
# Constants
# ===========================================================================
B = 0
J = 1
# k = 1
N = 50  # number of spins
m = -3.05765101e-03
y_int = 1.12324380e+01

# ===========================================================================
# Simulation
# ===========================================================================

class Simulation:
    """Class to run the algorithm

    ==== Attributes ===
    system: All the latices to be simulated at varaying temperatures
    energy: list of energies throughout the simulation
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

    def plot_property(self, start):  #(self, property_: TypeVar, data, name: str):
        """Plot property"""
        time_ = np.array(self.time)
        # mp.plot(self.temp, self.internal_energy(start), 'g.') #, label=name)
        error = Error(self.internal_energy, start, True, None)
        return error.jackknife_method()
        # mp.errorbar(self.temp, self.internal_energy(start), yerr=error)
        # mp.show()

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

        if isinstance(equil_time, list):
            for i in range(len(self.system)):
                zero_ = autocorrelation(self.time-equil_time[i], 0,
                                        self.eq_magnetization[i][equil_time[i]:])
                if zero_ == 0:
                    zero_ = 0.1
                times = np.arange(self.time - 1 - equil_time[i])
                integrand_ = []
                for t in times:
                    to_add = autocorrelation(self.time-equil_time[i], t,
                                             self.eq_magnetization[i]
                                             [equil_time[i]:]) / zero_
                    integrand_.append(to_add)
                x = abs(round(sum(integrand_)))
                if x == 0:
                    x = 1
                self.correlation_times.append(x)


class MetropolisAlgorithm(Simulation):
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

    def generate(self, lattice: Lattice):
        """generate a new state v by randomly flipping 1 spin"""
        k = r.randint(0, lattice.size - 1)  # index of spin to be flipped
        state_v = lattice.__copy__()
        spin_k = state_v.spins[k] * -1
        state_v.spins[k] = spin_k
        state_v.magnetization = lattice.magnetization + 2 * spin_k

        state_v.energy = lattice.energy + 2 * J * spin_k * sum(
            lattice._neighbours(k))
        return state_v

    def _acceptance(self, original: Lattice, other: Lattice) -> float:
        """Return the acceptance probability according to the metropolis algorithm
        w/ single spin dynamic."""
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
    config_ = {'range': [0.11, 3.02], 'step': 0.1, 'lattice_size': N,
               'error_method': 12, 'temp': 1, 'num_rounds': 5000}
    sim = MetropolisAlgorithm(config_)
    sim.run()
    sim.correlation_time(3000)

    ##############
    # indpendent e measuments
    ##############

    # plotting
    # temps_ = np.arange(0.01, 4, 0.01)
    # np.insert(temps_, 0, 0.01)
    # mp.plot(temps_, u(1 / temps_), 'm', label='Ananlytical Solution')

    # i_energy = sim.internal_energy()
    # indep_energy_ = sim.ind_energy.copy()
    # for energy in indep_energy_:
    #     error.append(Error(rms, energy, True, None))
    # mp.errorbar(sim.temp, i_energy, yerr=[err.jackknife_method() for err
    #                                       in error], fmt='o')

    # spc = sim.specific_heat()
    # sus = sim.suscepitbility()
    # error = []
    #
    # indep_energy_ = sim.ind_energy.copy()
    # indep_mag_ = sim.ind_mag.copy()
    # for energy in indep_energy_:
    #     error.append(Error(rms, energy, True, None))
    # mp.errorbar(sim.temp, spc, yerr=[err.jackknife_method() for err
    #                                  in error], fmt='o')

    # for energy in indep_mag_:
    #     error.append(Error(rms, energy, True, None))
    # mp.errorbar(sim.temp, sus, yerr=[err.jackknife_method() for err
    #                                  in error], fmt='o')

    # mp.xlabel('Temperature 1/T (units K/J)')
    # mp.ylabel('Energy <U> (units J)')
    # mp.legend()
    # mp.title('Inverse Temperature vs Energy')
    # # mp.show()
    # mp.savefig('corr.png')

    # i_mag = sim.magnetization_(3000)
    # error_m = []
    # for mag in sim.eq_magnetization:
    #     error_m.append(Error(average, list(mag[3000:4000]), True, None))
    # mp.errorbar(sim.temp, i_mag, yerr=[err_.jackknife_method() for err_
    #                                    in error_m], fmt='o')
    #
    # mp.xlabel('Temperature 1/T (units K/J)')
    # mp.ylabel('Sus χ (units J)')
    # mp.legend()
    # mp.title('Inverse Temperature vs Sus')
    # mp.show()
    # mp.savefig('suscepy.png')

