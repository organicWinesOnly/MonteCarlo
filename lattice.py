""" lattice.py

File containg lattices of different dimenions for different purposes.
Each lattice tpe is represented by a class.
The different lattices contained are: 
    Ising: represents a 1D isling lattice
    Ising2D: represents a 2D isling lattice
    Potts: represents a 1D q-pott lattice
    Potts2D: represents a 2D q-pott lattice

@author: organicWinesOnly
"""
import numpy as np
from typing import *
import random as r
from helpers_mc import flip_coin


class Spin:
    """
    Repersent a spin in the lattice. Used for wolff algorithm.

    ==== attributes ===
    spin: spin value. Value can only be a one or zero.
    cluster: True if spin is in a cluster
    location: (row, column) location of spin in lattice
    considered: True if spin was already considered to join a latttice
    """
    spin: int
    cluster: bool
    location: tuple
    considered: bool

    def __init__(self, q: int, location_: tuple, set_=0):
        # set_ = 0 implies a value of the spin has no been assigned
        if set_ == 0:
            self.spin = flip_coin(-1, 1)
        else:
            self.spin = set_
        self.cluster = False
        self.location = location_
        self.considered = False


class Lattice:
    """1D lattice

    ==== Attributes ===
    spins: collection of Spins representing the lattice
    size: size of lattice
    total_energy: total energy per spin
    m: mean magnetizaion per spin
    """
    spins: np.ndarray #List[Spin]
    size: int
    m: float
    beta: float
    temp: int
    total_energy: int

    def __init__(self, n: int, temp: int, run_at: float) -> None:
        """Initialize Lattice class

        === parameters ===
        n: size of the lattice
        temp: the temp you want the sytem to start at, 0 or -1
              A zero value implies that temperature is 0 k. -1 implies the
              temperature is starting at positive infinity
        run_at: temperture you want the system to come to

        Representation Invariants:
        temp = 0 or 1
        """
        self.size = n  # needed?

        # build lattice
        # spin class is not needed for 1d metropolis
        if temp == 0:
            # k = flip_coin(-1, 1)
            self.spins = np.array([1 for _ in range(n)])
        else:
            self.spins = np.array([flip_coin(-1, 1) for _ in range(n)])

        self.m = np.sum(self.spins) / n
        self.beta = 1 / run_at

        interaction_ = np.zeros(n)

        for k in range(n):
            spin_j = self._neighbours(k)
            interaction_[k] = np.sum(spin_j) * self.spins[k]

        self.total_energy = -1 * np.sum(interaction_) /self.size # J = 1

    def _neighbours(self, site: int) -> np.ndarray:
        """Return neighbouring spins [<left>, <right>].

        Periodic Boundary conditions are applied.
        """
        neighbours = np.zeros(2)

        if site + 1 == self.spins.shape[0]:
            right_ = self.spins[0]
        else:
            right_ = self.spins[site + 1]

        if site - 1 < 0:
            left_ = self.spins[-1]
        else:
            left_ = self.spins[site - 1]

        neighbours[0], neighbours[1] = left_, right_

        return neighbours

    def update(self, delta_energy, delta_m) -> None:
        """ Recalculate attributes of the lattice that depend on spin
        """
        self.m += delta_m / self.size # not true?
        self.total_energy += delta_energy / self.size

class Lattice2D:
    """Create a 2d lattice

     === Attributes ===
        size: amount of spins in the lattice
        energy: total energy in the lattice
        m: magnetization per spin
        temp: starting temperature
        beta: inverse of temp we want our system to come to at equilibrium
        spins: ndarray repersenting the spins in our lattice

        Repersentation Invariants:
        temp = 0 or 1
        """
    size: Tuple
    total_energy: float
    m: float
    temp: int
    spins: np.ndarray

    def __init__(self, n: tuple, initial_temp: int, run_at: float) -> None:
        """Initailize a 2d lattice

        === parameters ===
        n: size of the lattice
        temp: the temp you want the sytem to start at, 0 or -1
        run_at: temperture you want the system to come to
        """
        self.spins = np.zeros(n)
        self.size = self.spins.size
        self.temp = initial_temp
        self.beta = 1 / run_at

        if initial_temp == 0:
            k = flip_coin(-1, 1)
            self.spins.fill(k)
        else:
            for i in range(self.spins.shape[0]):
                for j in range(self.spins.shape[1]):
                    self.spins[i, j] = flip_coin(-1, 1)

        self.m = np.sum(self.spins) / self.size

        interaction_ = np.zeros(self.spins.shape)
        for i in range(self.spins.shape[0]):
            for j in range(self.spins.shape[1]):
                nn_spins = self._neighbours((i, j))
                interaction_[i, j] = np.sum(nn_spins) * self.spins[i, j]

        self.total_energy = -1 * np.sum(interaction_) / self.size

    def _neighbours(self, site: tuple) -> List[int]:
        """Return neighbouring spins [<left>, <right>, <above> , <below>].

        Periodic Boundary conditions are applied
        """
        x = site[0]  # corridnates of lattice
        y = site[1]

        # left right
        if x + 1 == self.spins.shape[0]:
            right_ = self.spins[0][y]
        else:
            right_ = self.spins[x + 1][y]

        if x - 1 < 0:
            left_ = self.spins[-1][y]
        else:
            left_ = self.spins[x - 1][y]

        # up, down
        if y + 1 == self.spins.shape[1]:
            up_ = self.spins[x][0]
        else:
            up_ = self.spins[x][y + 1]

        if y - 1 < 0:
            down_ = self.spins[x][-1]
        else:
            down_ = self.spins[x][y - 1]

        return [left_, right_, up_, down_]

    def update(self, delta_energy, delta_m) -> None:
        """ Recalculate attributes of the lattice that depend on spin
        """
        self.m += delta_m / self.size # not true?
        self.total_energy += delta_energy / self.size 


class Plotts:
    """Create a 2d lattice that can spin in q directions

     === Attributes ===
        size: 2 dim tupple repersenting (num of rows, num of coloums)
        energy: total energy of lattice
        magnetization: total magnetization of the lattice
        temp: we want our sytem to come to equilibrium at
        beta: inverse of temp we want are ystem to come to equilibrium at
        spins: ndarray repersenting the spins in our lattice


        Repersentation Invariants:
        temp = 0 or 1
        """
    size: Tuple
    energy: float
    magnetization: int
    temp: int
    spins: np.ndarray

    def __init__(self, n: tuple, q: int, initial_temp: int, run_at: float)\
            -> None:
        """Initailize a 2d lattice

        === parameters ===
        n: size of the lattice
        temp: the temp you want the sytem to start at, 0 or -1
        run_at: temperture you want the system to come to
        """
        self.size = n
        columns = []
        self.intial_temp = initial_temp
        mag = []
        self.state = q
        self.size_tot = n[1] * n[0]

        k = 0
        for i in range(n[0]):
            if initial_temp == 0:
                if k == 0:
                    k = flip_coin(1, -1 )# r.randint(1, q)
                u = [Spin(q, (i, i_), set_=k) for i_ in range(n[1])]
                mag.append(sum(item.spin for item in u))
            else:
                u = [Spin(q, (i, i_)) for i_ in range(n[1])]
                mag.append(sum(item.spin for item in u))
            columns.append(u)

        self.spins = columns
        # assert self.spins.shape == self.size

        self.temp = run_at
        self.magnetization = sum(mag)

        interaction_ = []
        for column in self.spins:
            for spin_loc in column:
                spin_j_right = self.neighbours(spin_loc)[1].spin
                interaction_.append(spin_j_right *
                                    spin_loc.spin)
                spin_j_down = self.neighbours(spin_loc)[3].spin
                interaction_.append(spin_j_down *
                                    spin_loc.spin)
                spin_j_left = self.neighbours(spin_loc)[2].spin
                interaction_.append(spin_j_left *
                                    spin_loc.spin)
                spin_j_up = self.neighbours(spin_loc)[0].spin
                interaction_.append(spin_j_up *
                                    spin_loc.spin)

        self.beta = 1 / run_at
        self.energy = -1 * sum(interaction_) * self.beta

    def __copy__(self):
        """Return a copy of the lattice"""
        other = Plotts(self.size, self.state, self.intial_temp, self.temp)
        other.spins = self.spins.copy()
        other.magnetization = self.magnetization
        other.energy = self.energy
        # assert other.magnetization == self.magnetization
        # assert other.spins is not self.spins
        return other

    def neighbours(self, spin: Spin) -> List[Spin]:
        """Return neighbouring spins [<left>, <right>, <above> , <below>].

        Periodic Boundary conditions are applied
        """
        x = spin.location[0]  # corridnates of lattice
        y = spin.location[1]

        # left right
        if x + 1 == self.size[0]:
            right_ = self.spins[0][y]
        else:
            right_ = self.spins[x + 1][y]

        if x - 1 < 0:
            left_ = self.spins[-1][y]
        else:
            left_ = self.spins[x - 1][y]

        # up, down
        if y + 1 == self.size[1]:
            up_ = self.spins[x][0]
        else:
            up_ = self.spins[x][y + 1]

        if y - 1 < 0:
            down_ = self.spins[x][-1]
        else:
            down_ = self.spins[x][y - 1]

        return [left_, right_, up_, down_]

    def update(self):
        interaction_ = []
        for column in self.spins:
            for spin_loc in column:
                spin_j_right = self.neighbours(spin_loc)[1].spin
                interaction_.append(spin_j_right *
                                    spin_loc.spin)
                spin_j_down = self.neighbours(spin_loc)[-1].spin
                interaction_.append(spin_j_down *
                                    spin_loc.spin)
                spin_j_left = self.neighbours(spin_loc)[2].spin
                interaction_.append(spin_j_left *
                                    spin_loc.spin)
                spin_j_up = self.neighbours(spin_loc)[0].spin
                interaction_.append(spin_j_up *
                                    spin_loc.spin)

        self.energy = -1 * sum(interaction_) * self.beta
        magnet = []
        for i in range(self.size[0]):
            magnet.append(sum(self.spins[i][p].spin for p in range(self.size[1])))

        self.magnetization = sum(magnet)


class Plotts2:
    """Create a 2d lattice that can spin in q directions

     === Attributes ===
        size: 2 dim tupple repersenting (num of rows, num of coloums)
        energy: total energy of lattice
        magnetization: total magnetization of the lattice
        temp: we want our sytem to come to equilibrium at
        beta: inverse of temp we want are ystem to come to equilibrium at
        spins: ndarray repersenting the spins in our lattice


        Repersentation Invariants:
        temp = 0 or 1
        """
    size: Tuple
    energy: float
    magnetization: int
    temp: int
    spins: np.ndarray

    def __init__(self, n: tuple, q: int, initial_temp: int, run_at: float)\
            -> None:
        """Initailize a 2d lattice

        === parameters ===
        n: size of the lattice
        temp: the temp you want the sytem to start at, 0 or -1
        run_at: temperture you want the system to come to
        """
        self.size = n
        columns = []
        self.intial_temp = initial_temp
        mag = []
        self.state = q
        self.size_tot = n[1] * n[0]

        for i in range(n[0]):
            if initial_temp == 0:
                k = r.randint(1, q)
                u = [Spin(q, (i, i_), set_=k) for i_ in range(n[1])]
                mag.append(sum(item.spin for item in u))
            else:
                u = [Spin(q, (i, i_)) for i_ in range(n[1])]
                mag.append(sum(item.spin for item in u))
            columns.append(u)

        self.spins = columns
        # assert self.spins.shape == self.size

        self.temp = run_at
        self.magnetization = sum(mag)

        interaction_ = []
        for column in self.spins:
            for spin_loc in column:
                spin_j_right = self.neighbours(spin_loc)[1].spin
                interaction_.append(spin_j_right *
                                    spin_loc.spin)
                spin_j_down = self.neighbours(spin_loc)[3].spin
                interaction_.append(spin_j_down *
                                    spin_loc.spin)
                spin_j_left = self.neighbours(spin_loc)[2].spin
                interaction_.append(spin_j_left *
                                    spin_loc.spin)
                spin_j_up = self.neighbours(spin_loc)[0].spin
                interaction_.append(spin_j_up *
                                    spin_loc.spin)

        self.beta = 1 / run_at
        self.energy = -1 * sum(interaction_) * self.beta

    def __copy__(self):
        """Return a copy of the lattice"""
        other = Plotts(self.size, self.state, self.intial_temp, self.temp)
        other.spins = self.spins.copy()
        other.magnetization = self.magnetization
        other.energy = self.energy
        # assert other.magnetization == self.magnetization
        # assert other.spins is not self.spins
        return other

    def neighbours(self, spin: Spin) -> List[Spin]:
        """Return neighbouring spins [<left>, <right>, <above> , <below>].

        Periodic Boundary conditions are applied
        """
        x = spin.location[0]  # corridnates of lattice
        y = spin.location[1]

        # left right
        if x + 1 == self.size[0]:
            right_ = self.spins[0][y]
        else:
            right_ = self.spins[x + 1][y]

        if x - 1 < 0:
            left_ = self.spins[-1][y]
        else:
            left_ = self.spins[x - 1][y]

        # up, down
        if y + 1 == self.size[1]:
            up_ = self.spins[x][0]
        else:
            up_ = self.spins[x][y + 1]

        if y - 1 < 0:
            down_ = self.spins[x][-1]
        else:
            down_ = self.spins[x][y - 1]

        return [left_, right_, up_, down_]

    def update(self):
        interaction_ = []
        for column in self.spins:
            for spin_loc in column:
                spin_j_right = self.neighbours(spin_loc)[1].spin
                interaction_.append(spin_j_right *
                                    spin_loc.spin)
                spin_j_down = self.neighbours(spin_loc)[-1].spin
                interaction_.append(spin_j_down *
                                    spin_loc.spin)
                spin_j_left = self.neighbours(spin_loc)[2].spin
                interaction_.append(spin_j_left *
                                    spin_loc.spin)
                spin_j_up = self.neighbours(spin_loc)[0].spin
                interaction_.append(spin_j_up *
                                    spin_loc.spin)

        self.energy = -1 * sum(interaction_) * self.beta
        magnet = []
        for i in range(self.size[0]):
            magnet.append(sum(self.spins[i][p].spin for p in range(self.size[1])))

        self.magnetization = sum(magnet)


if __name__ == '__main__':
    pass
