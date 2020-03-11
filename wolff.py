"""Wolff cluster algoritim for q-state potts spin system
"""

import matplotlib.pyplot as mp
from metropolisAlgortithm import Simulation
from lattice import Plotts, Spin, Plotts3
from Errors import *
from typing import List

J = 1
N = 100

class Wolff(Simulation):
    """Run a wolff algorithim

    === Attributes ===
    system: 
    q_states
    length_c
    spin_value_tup: 
    """
    system: List[Plotts]
    q_states: int
    lenght_c: List
    spin_value_tup: List

    def __init__(self, config):
        """
        Extend original initalizer
        """
        super().__init__(config)
        iterations = np.arange(config['range'][0], config['range'][1],
                               config['step'])
        if isinstance(config['lattice_size'], tuple):
            self.system = [Plotts(config['lattice_size'], config['temp'],
                                  config['num_states'], p) for p in iterations]
        self.q_states = config['num_states']
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
                    # fig, ax = mp.subplots()
                    # ax.pcolor(spin_value)
                    # mp.savefig(f'{p_lattice_.temp}_{self.time}.png')
            print(k)

            self.eq_energy.append(e_energy)
            self.eq_magnetization.append(e_mag)

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

                    # if neighbour.spin == p_lattice.spins[y][x].spin and not \
                    #         neighbour.considered:
                    #     neighbour.considered = True
                    #     if self.p_add(p_lattice):
                    #         neighbour.cluster = True
                    #         sim.flip_avg += 1
                    #         cluster.append(neighbour)

        # for spin in cluster:
        #     neighbours = p_lattice.neighbours(spin)
        #     for neighbour in neighbours:
        #         if not neighbour.cluster and neighbour.spin == \
        #                 p_lattice.spins[y][x].spin:
        #             if self.p_add(p_lattice):
        #                 neighbour.cluster = True
        #                 cluster.append(neighbour)

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

    def plot_system(self):
        for lattice_num in range(len(self.system)):
            mp.figure()
            times_ = np.arange(self.time)[:1000]
            # mp.plot(times_, self.eq_energy[lattice_num][:1000], 'b.',
            #         label='Internal Energy')
            mp.plot(times_, self.eq_magnetization[lattice_num][:1000], 'r.',
                    label='Magnetization')
            # auto_corr = self.auto_correlation(times_)
            # mp.semilogy(times_, auto_corr)
            mp.title('Energy and Magnetization vs Time for {}K'.format(
                self.temp[lattice_num]))
            mp.legend()
            mp.show()


class WolffPotts(Simulation):
    """Run a wolff for 3 state potts algorithim
    === Attributes ===
    system: 
    q_states
    length_c
    spin_value_tup: 
    """
    system: List[Plotts3]
    q_states: int
    lenght_c: List
    spin_value_tup: List

    def __init__(self, config):
        """
        Extend original initalizer
        """
        super().__init__(config)
        iterations = np.arange(config['range'][0], config['range'][1],
                               config['step'])
        if isinstance(config['lattice_size'], tuple):
            self.system = [Plotts3(config['lattice_size'], config['temp'],
                                   config['num_states'], p) for p in iterations]
        self.q_states = config['num_states']
        self.lenght_c = []

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
            print(k)

            self.eq_energy.append(e_energy)
            self.eq_magnetization.append(e_mag)

        self.time /= len(self.system)

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

        self.lenght_c.append(len(cluster))
        return cluster

    def _flip_cluster(self, cluster) -> None:
        """flip a cluster in the lattice
        """
        new_spin = r.randint(1, self.q_states)
        while new_spin == cluster[0].spin:
            new_spin = r.randint(1, self.q_states)
        for spin in cluster:
            spin.spin = new_spin
            spin.cluster = False

    def p_add(self, p_lattice: Plotts):
        value = np.random.random_sample()
        expont = np.e ** (-2 * p_lattice.beta * J)
        if value < 1-expont:
            return True
        return False

    def plot_system(self):
        for lattice_num in range(len(self.system)):
            mp.figure()
            times_ = np.arange(self.time)[:1000]
            # mp.plot(times_, self.eq_energy[lattice_num][:1000], 'b.',
            #         label='Internal Energy')
            mp.plot(times_, self.eq_magnetization[lattice_num][:1000], 'r.',
                    label='Magnetization')
            # auto_corr = self.auto_correlation(times_)
            # mp.semilogy(times_, auto_corr)
            mp.title('Energy and Magnetization vs Time for {}K'.format(
                self.temp[lattice_num]))
            mp.legend()
            mp.show()


if __name__ == '__main__':
    config_ = {'size': 100, 'range': [1.2, 4.2], 'step': 0.2,
               'lattice_size': (10, 10), 'error_method': 12,
               'temp': 0, 'num_rounds': 2000, 'num_states': 2}

    sim = Wolff(config_)
    sim.run()

    sim.correlation_time(1000)
    sim.indep_magnetization(1000)
    sim.indep_energy(1000)

    spc = sim.specific_heat()
    sus = sim.suscepitbility()
    energy = sim.internal_energy()
    mag = sim.tot_magnetization()
    print(sus)

    indep_energy_ = sim.ind_energy.copy()
    indep_mag_ = sim.ind_mag.copy()

    error_spc = []
    error_sus = []
    error_mag = []
    error_en = []

    for energy_ in indep_energy_:
        error_en.append(Error(average, energy_, True, None))
    assert len(error_en) == len(energy)
    assert len(energy) == len(sim.temp)
    mp.errorbar(sim.temp, energy, yerr=[err.jackknife_method() for err
                                        in error_en], fmt='o')
    mp.title('10 x 10 lattice')
    mp.xlabel('Temperature, T (units K/J)')
    mp.ylabel('Energy, <U> units(J)')
    # mp.savefig('i2d_int_engery_wolff_with_error.png')
    mp.show()

    beta_ = [1 / (N*temp) for temp in sim.temp]
    mp.figure()
    opo = 0
    for energy_ in indep_energy_:
        bet = [beta_[opo]]
        error_spc.append(Error(two_pnt_spec, energy_, True, param=bet))
        opo += 1
    mp.errorbar(sim.temp, spc, yerr=[err.jackknife_method() for err
                                     in error_spc], fmt='o')
    mp.title('10 x 10 lattice')
    mp.xlabel('Temperature, T (units K/J)')
    mp.ylabel('Specific Heat, c')
    # mp.savefig('i2d_specific_heat_wolff_with_error.png')
    mp.show()

    mp.figure()
    for mag_ in indep_mag_:
        error_mag.append(Error(average, mag_, True, None))
    mp.errorbar(sim.temp, mag, yerr=[err.jackknife_method() for err
                                     in error_mag], fmt='o')
    mp.title('10 x 10 lattice')
    mp.xlabel('Temperature, T (units K/J)')
    mp.ylabel('Magnetization, m')
    # mp.savefig('i2d_tot_mag_wolff_with_error.png')
    mp.show()

    mp.figure()

    opo = 0
    for mag_ in indep_mag_:
        bet = [beta_[opo]]
        error_sus.append(Error(two_pnt_sus, mag_, True, None))
        opo += 1
    mp.errorbar(sim.temp, sus, yerr=[err.jackknife_method() for err
                                     in error_sus], fmt='o')
    mp.title('10 x 10 lattice')
    mp.xlabel('Temperature, T (units K/J)')
    mp.ylabel('Susceptibility, χ')
    # mp.savefig('i2d_suscep_wolff_with_error.png')
    mp.show()
