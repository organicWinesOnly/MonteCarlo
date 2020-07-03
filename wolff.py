"""Wolff cluster algoritim for q-state potts spin system

last updated: May 30, 2020
"""

from metropolisAlgortithm import Simulation

J = 1
N = 100

class Wolff(Simulation):
    """Run a wolff algorithim on a 2D lattice.

    === Attributes ===
    system: set of lattices
    length_c:
    spin_value_tup: 
    """
    system: List[Plotts]
    lenght_c: List
    spin_value_tup: List

    def __init__(self, config):
        """Extend original initalizer.
        """
        super().__init__(config)
        # the range of equilibrium temperatures 
        iterations = np.arange(config['range'][0], config['range'][1],
                               config['step'])

        # a tuple implies the lattice is two dimensional
        if isinstance(config['lattice_size'], tuple):
            self.system = [Plotts(config['lattice_size'], config['temp'],
                                  config['num_states'], p) for p in iterations]
        self.length_c = []
        self.spin_value_tup = []

    def run(self):
        """Pass"""
        for p_lattice_ in self.system:
            e_energy = []
            e_mag = []
            k = 0
            for _ in range(self.num_rounds):
                copy_ = p_lattice_.__copy__()
                cluster_to_flip = self._generate_cluster(p_lattice_)
                self._flip_cluster(cluster_to_flip)
                assert copy_ != p_lattice_
                p_lattice_.update()
                e_energy.append(p_lattice_.energy)
                e_mag.append(abs(p_lattice_.magnetization))
                self.time += 1
                if not len(cluster_to_flip) / 100 >= 0.5:
                    k += 1
                if p_lattice_ is self.system[4]:
                    t2 = self.time % 100 == 2
                    t1 = self.time % 100 == 1
                    t3 = self.time % 100 == 3
                    t4 = self.time % 100 == 4
                    t5 = self.time % 100 == 5
                    t6 = self.time % 100 == 6
                    if t1 or t2 or t3 or t4 or t5 or t6:
                        copy_of_spins = p_lattice_.spins.copy()
                        spin_value = []
                        for col in copy_of_spins:
                            col1 = []
                            for spin in col:
                                col1.append(spin.spin)
                            spin_value.append(col1)
                        self.spin_value_tup.append((p_lattice_.temp, self.time,
                                                    spin_value))
            self.final_energy.append(e_energy)
            self.final_magnetization.append(e_mag)

        self.time /= len(self.system)

    def _spin_values(self, lattice: Plotts):
        """creates list of spin values needed to plot a visual repersentation of
         the cluster flipping"""
        t2 = self.time % 100 == 2
        t1 = self.time % 100 == 1
        t3 = self.time % 100 == 3
        t4 = self.time % 100 == 4
        t5 = self.time % 100 == 5
        t6 = self.time % 100 == 6
        if t1 or t2 or t3 or t4 or t5 or t6:
            copy_of_spins = lattice.spins.copy()
            spin_value = []
            for col in copy_of_spins:
                col1 = []
                for spin in col:
                    col1.append(spin.spin)
                spin_value.append(col1)
        return

    def _generate_cluster(self, p_lattice: Plotts) -> List[Spin]:
        """generate cluster"""

        y = r.randint(0, p_lattice.size[0] - 1)
        x = r.randint(0, p_lattice.size[1] - 1)

        cluster = [p_lattice.spins[y][x]]
        for spin in cluster:
            neighbours = p_lattice.neighbours(spin)
            for neighbour in neighbours:
                if not neighbour.cluster:
                    if neighbour.spin == p_lattice.spins[y][x].spin:
                        if self.p_add(p_lattice):
                            neighbour.cluster = True
                            cluster.append(neighbour)

        for spin in cluster:
            spin.considered = False
            spin.cluster = False

        self.length_c.append(len(cluster))
        return cluster

    def _flip_cluster(self, cluster) -> None:
        """flip a cluster in the lattice
        """
        new_spin = cluster[0].spin * -1
        for spin in cluster:
            spin.spin = new_spin
            spin.cluster = False

    def p_add(self, p_lattice: Plotts):
        value = np.random.random_sample()
        expont = np.e ** (-2 * p_lattice.beta * J)
        if value < 1-expont:
            return True
        return False
