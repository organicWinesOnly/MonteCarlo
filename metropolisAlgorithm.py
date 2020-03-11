"""Implementation of the metropolis algorithm with single-spin-flip dynamics.
"""


from simulation import *
# import scipy.optimize as sp


# ===========================================================================
# Constants
# ===========================================================================
B = 0
J = 1
# k = 1
N = 300  # number of spins
m = -3.05765101e-03
y_int = 1.12324380e+01

# ===========================================================================
# Simulation
# ===========================================================================

class metropolisAlgorithm(Simulation):
    """Metropolis Algorithm applied to a 1d isling model
    """
    def __init__(self, config):
        super(metropolisAlgorithm, self).__init__(config)
        self.exp = np.exp(-2 * self.beta)  #calculate exponent before hand
        self.system = [Lattice(config['lattice_size'], config['temp'], p)
                           for p in self.temp]

    def run(self):
        for j, lattice_ in enumerate(self.system):
            magnetization_mu = np.zeros(self.num_rounds)
            energy_mu = np.zeros(self.num_rounds)
            for i in range(self.num_rounds):
                # generate a new lattice

                # index of spin to be flip
                for _ in range(lattice_.size):
                    k = r.randint(0, lattice_.size - 1)  

                    # run the algorithm and return the appropiate lattice
                    self.update_lattice(lattice_, k, j)

                # magnetization for every state of the lattice in the
                # simulation
                magnetization_mu[i] = lattice_.m
                # total energy for every state of the lattice in the simulation
                energy_mu[i] = lattice_.total_energy
                self.time += 1

            # measurements for all the lattices in the system
            self.final_magnetization[j] = np.array(magnetization_mu)
            self.final_energy[j] = np.array(energy_mu)

    def _acceptance(self, original: Lattice, site: int, index) -> float:
        """Return the acceptance probability according to the metropolis algorithm
        w/ single spin dynamic."""
        delta_energy = 2 * J * original.spins[site] *\
                        np.sum(original._neighbours(site))

        if delta_energy < 0:
            return self.exp[index], delta_energy
        else:
            return 1.0, delta_energy

    def _change_state(self, original: Lattice, site: int, delta_energy: float) -> None:
        """change the state to that of other
        """
        original.spins[site] = -1 * original.spins[site]
        delta_m = 2 * original.spins[site]
        original.update(delta_energy, delta_m)

    def update_lattice(self, original: Lattice, site: int, index) -> None:
        """update the lattice according to the acceptance probability"""
        number = r.uniform(0, 1)

        while number == 1:
            number = r.uniform(0, 1)

        accept = self._acceptance(original, site, index)

        if number < accept[0]:
            self._change_state(original, site, accept[1])

# ===========================================================================
# Implementation
# ===========================================================================
# after 6000 points


#if __name__ == '__main__':
#    config_ = {'range': [10e-3, 5.02], 'step': 0.1, 'lattice_size': 300,
#               'error_method': 12, 'temp': 1, 'num_rounds': 3000}
#    sim = MetropolisAlgorithm(config_)
#    sim.run()
#    sim.correlation_time(4000)

