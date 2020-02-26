"""Compare (internal energy, magnetism, susceptibility, and specific heat) of
a 2d isling model with size 10x10 w/ 3000 steps, size 32x32 w/ 3000, steps, and
size 10x10 w/ 6000 steps.

Each property is plotted against the same figure"""

from wolff import Wolff
from Errors import Error
from helpers_mc import *
import matplotlib.pyplot as mp

# constants
N = [100, 32 ** 2, 100]
# configurations
config_ = {'size': 10, 'range': [1.2, 4.2], 'step': 0.2,
           'lattice_size': (10, 10), 'error_method': 12,
           'temp': 0, 'num_rounds': 3000, 'num_states': 2}
config_1 = {'size': 10, 'range': [1.2, 4.2], 'step': 0.2,
            'lattice_size': (32, 32), 'error_method': 12,
            'temp': 0, 'num_rounds': 3000, 'num_states': 2}
config_2 = {'size': 10, 'range': [1.2, 4.2], 'step': 0.2,
            'lattice_size': (10, 10), 'error_method': 12,
            'temp': 0, 'num_rounds': 6000, 'num_states': 2}

# run simulations
sim = Wolff(config_)
sim.run()

sim1 = Wolff(config_1)
sim1.run()

sim2 = Wolff(config_2)
sim2.run()

sim.correlation_time(1000)
sim.indep_magnetization(1000)
sim.indep_energy(1000)

sim1.correlation_time(1000)
sim1.indep_magnetization(1000)
sim1.indep_energy(1000)

sim2.correlation_time(2000)
sim2.indep_magnetization(2000)
sim2.indep_energy(2000)

spc = [sim.specific_heat(), sim1.specific_heat(), sim2.specific_heat()]
sus = [sim.suscepitbility(), sim1.suscepitbility(), sim2.suscepitbility()]
energy = [sim.internal_energy(), sim1.internal_energy(), sim2.internal_energy()]
magnet = [sim.tot_magnetization(), sim1.tot_magnetization(),
          sim2.tot_magnetization()]

indep_energy_ = [sim.ind_energy.copy(), sim1.ind_energy.copy(),
                 sim2.ind_energy.copy()]

indep_mag_ = [sim.ind_mag.copy(), sim1.ind_mag.copy(), sim2.ind_mag.copy()]

error_spc = []
error_sus = []
error_mag = []
error_en = []

beta_ = [1 / temp for temp in sim.temp]
opo = 0

for in_energy_ in indep_energy_:
    lst_ = []
    lst_0 = []
    for energy_ in in_energy_:
        bet = beta_[opo] / N[opo]
        lst_0.append(Error(rms_spec, energy_, True, param=bet))
        lst_.append(Error(average, energy_, True, None))
    error_en.append(lst_)
    error_spc.append(lst_0)
    opo += 1

for in_mag_ in indep_mag_:
    lst_ = []
    lst_0 = []
    for mag_ in in_mag_:
        lst_.append(Error(average, mag_, True, None))
        lst_0.append(Error(rms2, mag_, True, None))
    error_mag.append(lst_)
    error_sus.append(lst_0)


# energy
mp.errorbar(sim.temp, energy[0], yerr=[err.jackknife_method() for err
                                        in error_en[0]], fmt='o',
            label='10x10 w/ 3000 steps')
mp.errorbar(sim1.temp, energy[1], yerr=[err.jackknife_method() for err
                                         in error_en[1]], fmt='o',
            label='32x32 w/ 3000 steps')
mp.errorbar(sim2.temp, energy[2], yerr=[err.jackknife_method() for err
                                         in error_en[2]], fmt='o',
            label='10x10 w/ 3000 steps')
mp.legend()
mp.title('Energy Comparison')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Internal Energy <u>')
mp.savefig('2d_wolff_energy_comparison.png')


# mag
mp.figure()

mp.errorbar(sim.temp, magnet[0], yerr=[err.jackknife_method() for err
                                       in error_mag[0]], fmt='ro',
            label='10x10 w/ 3000 steps')
mp.errorbar(sim1.temp, magnet[1], yerr=[err.jackknife_method() for err
                                        in error_mag[1]], fmt='go',
            label='32x32 w/ 3000 steps')
mp.errorbar(sim2.temp, magnet[2], yerr=[err.jackknife_method() for err
                                        in error_mag[2]], fmt='bo',
            label='10x10 w/ 6000 steps')

mp.title('Magnetization Comparison')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Magnetization, m')
mp.legend()
mp.savefig('2d_wolff_mag_comparison.png')


# specific heat
mp.figure()

mp.errorbar(sim.temp, spc[0], yerr=[err.jackknife_method() for err
                                     in error_spc[0]], fmt='bo',
            label='10x10 w/ 3000 steps')
mp.errorbar(sim1.temp, spc[1], yerr=[err.jackknife_method() for err
                                      in error_spc[1]], fmt='go',
            label='32x32 w/ 3000 steps')
mp.errorbar(sim2.temp, spc[2], yerr=[err.jackknife_method() for err
                                      in error_spc[2]], fmt='ro',
            label='10x10 w/ 6000 steps')
mp.legend()
mp.title('Specific Heat Comparison')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Specific Heat, c')
mp.savefig('2d_wolff_spc_comparison.png')


# sus
mp.figure()

mp.errorbar(sim.temp, sus[0], yerr=[err.jackknife_method() for err
                                     in error_sus[0]], fmt='bo',
            label='10x10 w/ 3000 steps')
mp.errorbar(sim1.temp, sus[1], yerr=[err.jackknife_method() for err
                                      in error_sus[1]], fmt='go',
            label='32x32 w/ 3000 steps')
mp.errorbar(sim1.temp, sus[2], yerr=[err.jackknife_method() for err
                                      in error_sus[2]], fmt='ro',
            label='10x10 w/ 6000 steps')
mp.legend()
mp.title('Susceptibility Comparison')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Susceptibility, χ')
mp.savefig('2d_wolff_sus_comparison.png')
