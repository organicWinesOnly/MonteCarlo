"""Implementation of the metropolis algorithm
"""

from simulation import *

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
