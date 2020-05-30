"""run metropolis algorithm for 2d lattice.

@author: organicWinesOnly
"""

from simulation import *

<<<<<<< HEAD
# ===========================================================================
# Simulation
# ===========================================================================
=======
J = 1


>>>>>>> metropolis
class metropolisAlgorithm2d(Simulation):
    """Metropolis Algorithm.
    """
    def __init__(self, config):
        """ Initialize Class.
        """
        super(metropolisAlgorithm2d, self).__init__(config)
        self.system = [Lattice2D(config['lattice_size'], config['temp'], p)
                           for p in self.temp]
        self.exp2 = np.exp(-4 * self.beta)
        self.exp1 = np.exp(-2 * self.beta)
        self.size = config['lattice_size'][0] * config['lattice_size'][1] 

    def run(self):
        for j, lattice_ in enumerate(self.system):
            magnetization_mu = np.zeros(self.num_rounds)
            energy_mu = np.zeros(self.num_rounds)
            magnetization_mu[0] = lattice_.m
            energy_mu[0] = lattice_.total_energy
            for i in range(1, self.num_rounds):
                for _ in range(lattice_.size):
                    # generate the spin flip site at random
                    kx = r.randint(0, lattice_.spins.shape[0]-1)
                    ky = r.randint(0, lattice_.spins.shape[1]-1)

                    # run the algorithm and update the lattice
                    self.update_lattice(lattice_, (kx, ky), j)
                magnetization_mu[i] = lattice_.m
                energy_mu[i] =lattice_.total_energy

            self.final_magnetization[j] = magnetization_mu
            self.final_energy[j] = energy_mu

    def _acceptance(self, original: Lattice2D, site: tuple, index):
        """Return the acceptance probability according to the metropolis
        algorithm w/ single spin dynamic."""
        delta_energy = 2 * J * original.spins[site] * \
                        np.sum(original._neighbours(site))

        if delta_energy > 1:
            return self.exp2[index], delta_energy

        elif delta_energy > 0:
            return self.exp1[index], delta_energy
        else:
            return 1.0, delta_energy

    def _change_state(self, original: Lattice, site: int, delta_energy: float) -> None:
        """change the state to that of other
        """
        delta_m = -2 * original.spins[site] 
        original.spins[site] = -1 * original.spins[site]
        original.update(delta_energy, delta_m)


    def update_lattice(self, original: Lattice, site: tuple, index: int) -> None:
        """update the lattice according to the acceptance probability"""
        number = r.uniform(0, 1)

        while number == 1:
            number = r.uniform(0, 1)

        accept = self._acceptance(original, site, index)

        if number < accept[0]:
            self._change_state(original, site, accept[1])
