import numpy as np
import matplotlib.pyplot as plt
import pandas as ps
from scipy.constants import k, eV, h, atomic_mass
from scipy.integrate import trapezoid

# setting up constant
k_B = k/eV

# function for internal energy
def internal_energy(Z, T):
    Beta = 1/(k_B* T)
    U = -np.gradient(np.log(Z), Beta)
    return U

# function for heat capacity
def heat_capa(U, T):
    Cv = np.gradient(U, T)
    return Cv