"""simulate a 3 state potts model using wollf algorithm"""

from wolff import WolffPotts
from Errors import Error
from helpers_mc import *
import matplotlib.pyplot as mp

config_ = {'size': 10, 'range': [0.2, 4.2], 'step': 0.05,
           'lattice_size': (10, 10), 'error_method': 12,
            'temp': 0, 'num_rounds': 3000, 'num_states': 3}

sim = WolffPotts(config_)
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
                                    in error_en], fmt='go')
mp.title('10 x 10 lattice')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Energy, <U> units(J)')
mp.savefig('_int_engery_wolff_with_error.png')

beta_ = [1 / temp for temp in sim.temp]
mp.figure()
opo = 0
for energy_ in indep_energy_:
    bet = beta_[opo] / N
    error_spc.append(Error(rms_spec, energy_, True, param=bet))
    opo += 1
mp.errorbar(sim.temp, spc, yerr=[err.jackknife_method() for err
                                 in error_spc], fmt='bo')
mp.title('10 x 10 lattice')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Specific Heat, c')
mp.savefig('3state_specific_heat_wolff_with_error.png')

mp.figure()
for mag_ in indep_mag_:
    error_mag.append(Error(average, mag_, True, None))
mp.errorbar(sim.temp, mag, yerr=[err.jackknife_method() for err
                                 in error_mag], fmt='ro')
mp.title('10 x 10 lattice')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Magnetization, m')
mp.savefig('3state_tot_mag_wolff_with_error.png')

mp.figure()

for mag_ in indep_mag_:
    error_sus.append(Error(rms2, mag_, True, None))
mp.errorbar(sim.temp, sus, yerr=[err.jackknife_method() for err
                                 in error_sus], fmt='mo')
mp.title('10 x 10 lattice')
mp.xlabel('Temperature, T (units K/J)')
mp.ylabel('Susceptibility, χ')
mp.savefig('3state_suscep_wolff_with_error.png')
