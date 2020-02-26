"""Simulating 1d isling model at equilibrium using Metropolis Algorithm
"""
import matplotlib.pyplot as mp
from helpers_mc import *
import pandas as pd

# import csv
# constants
B = 0
J = 1
# k = 1
N = 175  # number of spins


def analytic_spec_heat(temp):
    return 1 / temp ** 2 * (1 / np.cosh(1 / temp)) ** 2

# config_ = {'range': [0.1, 4.0], 'step': 0.1, 'lattice_size': N,
#            'error_method': 12, 'temp': 1, 'num_rounds': 7000}
# sim = MetropolisAlgorithm(config_)
# sim.run()
# sim.correlation_time(equil_time=4000)
# sim.indep_measurments(4000)
#
# equil_times = [4000 for _ in range(14)]
# equil_times.insert(0, 6000)
#
value = []
#path_ = '~/Documents/PHY/phy372/metropolis/coor/'
timess = np.arange(4000, 6999)
temps = np.arange(0.001, 3.0, 0.001)
dataframes = [pd.read_csv('1d_ma_in/1d_ma_in{0:.1f}.csv'.format(x)) for x
              in np.arange(0.1, 3.0, 0.1)]

list_of_mag = []
list_of_energy = []
for item in dataframes:
    list_of_mag.append(list(item.loc[1][1:]))
    list_of_energy.append(list(item.loc[0][1:]))

beta_ = [1/q for q in np.arange(0.1, 3.0, 0.1)]
spec_heat = [rms_spec(list_of_mag[i], beta_[i]) for i in
             range(len(np.arange(0.1, 3.0, 0.1)))]
int_energy = [average_e(list_of_mag[i]) for i in
             range(len(np.arange(0.1, 3.0, 0.1)))]
eng_temp = []


# for i in range(1, len(list_of_energy)):
    # for j in range(len(list_of_energy[i])):
    #     list_of_energy[i][j] /= 10
    # eng_temp.append(rms_spec(list_of_energy[i], beta_[i]))

mp.plot(np.arange(0.1, 3.0, 0.1), int_energy, 'ro')

# mp.plot(np.arange(0.001, 3.0, 0.001), analytic_spec_heat(temps))

# mp.plot(np.arange(0.1, 3.0, 0.1), spec_heat, 'ro')



# with open('1d_ma_in/ind.61.csv') as f:
#     reader = csv.reader(f)
#     for line in reader:
#         value.extend(line)
#
# val_prime = list(map(lambda x: float(x), value))
#
# mp.plot(val_prime, '.')
#
# mp.savefig('t1.61.png')
# # #
# sim = MetropolisAlgorithm(config_)
# sim.run()
# sim.correlation_time(4000)
# sim.indep_measurments(4000)
# k = sim.eq_magnetization
# spc = sim.specific_heat()
# sus = sim.suscepitbility()
# print(sus)
# energy = sim.internal_energy()
# mag = sim.tot_magnetization()
#
# indep_energy_ = sim.ind_energy.copy()
# indep_mag_ = sim.ind_mag.copy()
#
# error_spc = []
# error_sus = []
# error_mag = []
# error_en = []
#
# beta_ = [1 / temp for temp in sim.temp]
# opo = 0
# for energy_ in indep_energy_:
#     error_en.append(Error(average, energy_, True, None))
#     bet = beta_[opo] / N
#     error_spc.append(Error(rms_spec, energy_, True, param=bet))
#     opo += 1

# ##################
# # Energy
# ##################
# temps_ = np.arange(0.01, 4, 0.01)
# np.insert(temps_, 0, 0.01)
# mp.plot(temps_, u(1 / temps_), 'm', label='Ananlytical Solution')
# mp.errorbar(sim.temp, energy, yerr=[err.jackknife_method() for err
#                                     in error_en], fmt='bo')
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Energy, <U> units(J)')
# mp.savefig('1d_int_engery_ma_with_error.png')
#
#
##################
# Specific Heat
##################
# path_ = '~/Documents/PHY/phy372/metroplis/1d_ma/'

#
# mp.figure()
# mp.errorbar(sim.temp, spc, yerr=[err.jackknife_method() for err
#                                  in error_spc], fmt='o')
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Specific Heat, c')
# mp.savefig('1d_specific_heat_ma_with_error.png')
#
#
# for mag_ in indep_mag_:
#     error_mag.append(Error(average, mag_, True, None))
#     error_sus.append(Error(rms2, mag_, True, None))
# mp.errorbar(sim.temp, mag, yerr=[err.jackknife_method() for err
#                                  in error_mag], fmt='o')

#
# ##################
# # Magnetization
# ##################
# mp.figure()
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Magnetization, m')
# mp.savefig('1d_tot_mag_ma_with_error.png')
#
#
# ##################
# # Susceptibility
# ##################
#
# mp.figure()
# mp.errorbar(sim.temp, sus, yerr=[err.jackknife_method() for err
#                                  in error_sus], fmt='o')
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Susceptibility, χ')
# mp.savefig('1d_suscep_ma_with_error.png')
