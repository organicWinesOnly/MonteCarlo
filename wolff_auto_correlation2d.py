from wolff import *
import csv
import matplotlib.pyplot as plt


config_ = {'size': 10, 'range': [0.2, 4.2], 'step': 0.2,
           'lattice_size': (10, 10), 'error_method': 12,
            'temp': 0, 'num_rounds': 3000, 'num_states': 3}

# sim = Wolff(config_)
# sim.run()
# equil_times = [1000 for _ in range(20)]
# sim.correlation_time(equil_times)
# sim.write_file()
str = 0.2
sto = 4.2
ste = 0.2
path_ = '~/Documents/PHY/phy372/wolff/'

data_ = pd.read_csv(path_ + '2d_wolff0.4.csv')
# data_ = list(map(lambda x:
#                      list(pd.read_csv(path_ + '2d_wolff{0:.1f}.csv'.format(x),
#                                       index_col=0)), np.arange(str, sto, ste)))
# for i in range(5):
#     with plt.style.context('dark_background'):
#         # plt.plot(temp, g[0], 'ro', label='blocking')
#         plt.plot(data_[i], 'm^', label='jackknife')
#         plt.legend()
#         plt.xlabel('timestep')
#         plt.ylabel('Autocorrelation Value')
#         plt.title('Autocorrelation value, 200 data points')
#         plt.savefig('auto_corr_2dwolff.png')

# value = []
# times = np.arange(4000, 6999)
# with open('') as f:
#     reader = csv.reader(f)
#     for line in reader:
#         value.extend(line)

# val_prime = list(map(lambda x: float(x), value))

# mp.plot(val_prime, '.')
