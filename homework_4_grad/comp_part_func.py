import numpy as np
import pandas as pd
from scipy.constants import k, eV, h, atomic_mass
from scipy.integrate import trapezoid

#constants/givens
k_B = k/eV
epsilon = 0.0103
sigma = 3.4
T = np.linspace(10, 1000, 100)
V = 1000
Beta = 1/(k_B* T)
m = 39.948 * atomic_mass
lam = np.sqrt((Beta*h**2)/(2*np.pi*m))
r = np.linspace(3, 6, 100) 

# setting up lennard jones potential (see optimize_Argon_dimer)
def lennard_Jones_potential(r, epsilon, sigma): # default values for Argon
    V_LJ = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) # Lennard-Jones equation
    return V_LJ

# setting up partition_function
def partition_function(T):
    Z = []
    V_LJ = lennard_Jones_potential(r, epsilon, sigma)

    # generating a loop to store a value for each temepurature
    for temp in T:
        # re setting the constants to go over all Temps
        Beta = 1 / (k_B*temp)
        lam = np.sqrt((Beta*h**2) / (2*np.pi*m))

        #C = (1 / h**6) * (1 / lam**6)  # Pre-factor from momenta integration was breaking my code computre could not handle this had to estimate it as 1
        # surface area of sphereical elements
        Vr = 4 * np.pi*r**2
        # was trying to put restrictions on my exponets
        exponent = -Beta*V_LJ
        exponent = np.clip(exponent, -700, 700)
        integrand = Vr*np.exp(exponent)*.5

        # run the intergal with respect to r
        integral = trapezoid(integrand, r)
        Z.append(integral)
    Z = np.array(Z)
    return Z
