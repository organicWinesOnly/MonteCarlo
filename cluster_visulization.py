from wolff import *

config_ = {'size': 32**2, 'range': [1.2, 3.2], 'step': 0.2,
           'lattice_size': (32, 32), 'error_method': 12,
           'temp': 0, 'num_rounds': 4000, 'num_states': 2}

sim = Wolff(config_)
sim.run()

# sim.correlation_time(1000)
# sim.indep_magnetization(1000)
# sim.indep_energy(1000)
#
# spc = sim.specific_heat()
# sus = sim.suscepitbility()
# energy = sim.internal_energy()
# mag = sim.tot_magnetization()
# print(sus)
#
# indep_energy_ = sim.ind_energy.copy()
# indep_mag_ = sim.ind_mag.copy()
#
# error_spc = []
# error_sus = []
# error_mag = []
# error_en = []
#
# for energy_ in indep_energy_:
#     error_en.append(Error(average, energy_, True, None))
# assert len(error_en) == len(energy)
# assert len(energy) == len(sim.temp)
# mp.errorbar(sim.temp, energy, yerr=[err.jackknife_method() for err
#                                     in error_en], fmt='o')
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Energy, <U> units(J)')
# # mp.savefig('i2d_int_engery_wolff_with_error.png')
# mp.show()
#
# beta_ = [1 / (N*temp) for temp in sim.temp]
# mp.figure()
# opo = 0
# for energy_ in indep_energy_:
#     bet = [beta_[opo]]
#     error_spc.append(Error(two_pnt_spec, energy_, True, param=bet))
#     opo += 1
# mp.errorbar(sim.temp, spc, yerr=[err.jackknife_method() for err
#                                  in error_spc], fmt='o')
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Specific Heat, c')
# # mp.savefig('i2d_specific_heat_wolff_with_error.png')
# mp.show()
#
# mp.figure()
# for mag_ in indep_mag_:
#     error_mag.append(Error(average, mag_, True, None))
# mp.errorbar(sim.temp, mag, yerr=[err.jackknife_method() for err
#                                  in error_mag], fmt='o')
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Magnetization, m')
# # mp.savefig('i2d_tot_mag_wolff_with_error.png')
# mp.show()
#
# mp.figure()
#
# opo = 0
# for mag_ in indep_mag_:
#     bet = [beta_[opo]]
#     error_sus.append(Error(two_pnt_sus, mag_, True, bet))
#     opo += 1
# mp.errorbar(sim.temp, sus, yerr=[err.jackknife_method() for err
#                                  in error_sus], fmt='o')
# mp.title('10 x 10 lattice')
# mp.xlabel('Temperature, T (units K/J)')
# mp.ylabel('Susceptibility, χ')
# # mp.savefig('i2d_suscep_wolff_with_error.png')
# mp.show()
#
# mp.figure()

# lst = [sim.spin_value_tup[i] for i in range(300, 400)]
# print(lst)
for tupp in sim.spin_value_tup:
    mp.figure()
    mp.pcolor(tupp[2])
    mp.title('32 x 32 lattice')
    mp.colorbar()
    mp.savefig('{:.2}_{}.png'.format(tupp[0], tupp[1]))
