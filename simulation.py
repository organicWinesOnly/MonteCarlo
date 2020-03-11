""" simulation.py

File contains abstract class, Simulation, which will be used throughout the
directory.
"""

from lattice import *
from Errors import *
from helpers_mc import *
from autocorrelation import *

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
class Simulation:
    """Class to run an algorithm.

    ==== Attributes ===
    system: All the latices to be simulated at varying temperatures
    temp: List of temperatures the system is run at
    mag_spin: the magnetization per spin for every lattize in the system
    tot_energy: total energy for every lattice is the system
    final_energy:List containing List of the calculated energies of the system
                during the simulation
    final_magnetization: List of Lists of the calculated magnetization
                        of the system during the simulation
    ind_temp: List of list containing the independent energy of the simulation
    ind_magnetization: List of list containing the independent magnetization of
                        the simulation
    mag_spin: list of magnetization per spin for each lattice
    energy: list of energies throughout the simulation
    num_rounds: number of times to run the simulation
    time: time of simulation ran in monte carlo steps (number should == num
            rounds)
    correlation_times: list containing the measured correlation times
    n: number of spins in the system
    equilib_time: List containing the time in (monte carlo steps) for a lattice
    in the system to reach equilibrium during the simulation. Value must be
    hard coded.
    """

    system: List[Union[Lattice, Lattice2D, Plotts]]
    temp: List[Union[int, float]]
    tot_energy: np.ndarray
    mag_spin: np.ndarray
    time: int
    error_type: int
    final_energy: List[np.ndarray]
    final_magnetization: List[np.ndarray]
    correlation_times: List[float]
    ind_energy: List[np.ndarray]
    ind_mag: List[np.ndarray]
    num_rounds: int
    n: int
    equilib_time: List

    #config_ = {'range': [10e-3, 5.02], 'step': 0.1, 'lattice_size': 300,
    #           'error_method': 12, 'temp': 1, 'num_rounds': 3000}
    def __init__(self, config):
        """Initialize class.

        range_of_temps: 2 item list. [start temperature, highest temp]
        step: temperature step
        """
        self.system = [] # is intialized by the individual algorithms

        self.temp = np.arange(config['range'][0], config['range'][1],
                               config['step'])
        self.mag_spin = np.zeros(len(self.temp))
        self.energy_spin = np.zeros(len(self.temp))
        self.num_rounds = config['num_rounds']
        self.time = 0
        self.final_magnetization = np.zeros((self.temp.shape[0],
            self.num_rounds))
        self.final_energy = np.zeros((self.temp.shape[0],self.num_rounds))
        self.correlation_times = np.zeros(self.temp.shape)
        self.ind_energy = [0 for _ in range(self.temp.shape[0])]
        self.ind_mag = [0 for _ in range(self.temp.shape[0])]
        self.size = config['lattice_size']
        self.beta = np.array([1/t for t in self.temp])
        self.equilib_time = []
        self.error = config['error']

    def run(self):
        """Run the simulation"""
        raise NotImplementedError


    def indep_meausurements(self): #, equil_time: Union[int, List]):
        """ Finds the independent values of magnetization and energy for each 
        lattice in the system.
        """
        index = 0
        indep_mag = []
        indep_e = []

        for j in range(len(self.system)):
            num = abs(int(self.correlation_times[j]))

            i = 0
            indep_e_measurements = []
            indep_m_measurements = []
            
            while i * num < self.num_rounds:  # - equil_time:
                if i * num < self.equilib_time[j]:
                    i += 1
                else:
                    indep_m_measurements.append(self.final_magnetization[j][i * num])
                    indep_e_measurements.append(self.final_energy[j][i * num])
                    i += 1
            self.ind_mag[j] = np.array(indep_m_measurements)
            self.ind_energy[j] = np.array(indep_e_measurements)


    def energy_per_spin(self, independent=False) -> None:
        """Calculate the internal energy"""

        if independent == False:
            arr = self.final_energy
            for i, energy_sim in enumerate(arr):
                self.energy_spin[i] = np.mean(np.abs(energy_sim[self.equilib_time[i]:]))
        else:
            arr = self.ind_energy
            for i, energy_energy in enumerate(arr):
                self.energy_spin[i] = np.mean(np.abs(energy_sim))

    def magnetization_per_site(self, independent=False) -> None:
        """Calculate the magnetizaion per site averged over all the states in
        the simulation.

        The output is a 1d arraya containing all the values for each
        configuration.
        """
        if independent == False:
            arr = self.final_magnetization
            for i, mag_sim in enumerate(arr):
                self.mag_spin[i] = np.mean(np.abs(mag_sim[self.equilib_time[i]:]))
        else:
            arr = self.ind_mag
            for i, mag_sim in enumerate(arr):
                self.mag_spin[i] = np.mean(np.abs(mag_sim))

    def specific_heat(self) -> np.array:
        """ calculate the specific heat

        correceteed march 2020
        """
        c = np.zeros(len(self.system))
        for i in range(len(self.system)):
            c[i] = self.beta[i] ** 2 * two_pnt(self.ind_energy[i])  / self.size
        return c

    def suscepitbility(self) -> np.array:
        """ Calculate the susceptibility.
        """
        sus = np.zeros(len(self.system))
        for i in range(len(self.system)):
            sus[i] = self.beta[i] * two_pnt(np.abs(self.ind_mag[i]))  * self.size
        return sus

    def correlation_time(self) -> None:
        """Calculates the correlation time for each lattice.

        The time retuned is rounded to the nearest int.

        This function finds the correlation time of each function then
        mutiplies it by 2"""

        for i, mag_sim in enumerate(self.final_magnetization):
            self.correlation_times[i] =\
                            2 * ic_time(self.num_rounds-self.equilib_time[i],
                            np.abs(mag_sim[self.equilib_time[i]:])) 
