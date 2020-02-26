"""testing the results for the ma_confused file; the previous metropolis
algorithm simulation
"""

from ma_confused import *

N = 750

config_ = {'range': [0.11, 3.02], 'step': 0.1, 'lattice_size': N,
               'error_method': 12, 'temp': 0, 'num_rounds': 6500}

sim = MetropolisAlgorithm(config_)
sim.run()
sim.correlation_time(3000)
