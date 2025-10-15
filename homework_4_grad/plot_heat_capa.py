import numpy as np
import matplotlib.pyplot as plt
import pandas as ps
from scipy.constants import k, eV, h, atomic_mass
from scipy.integrate import trapezoid

#constants
k_B = k/eV
epsilon = 0.0103
sigma = 3.4
T = np.linspace(10, 1000, 100)
V = 1000
Beta = 1/(k_B* T)
m = 39.948 * atomic_mass
lam = np.sqrt((Beta*h**2)/(2*np.pi*m))
r = np.linspace(3, 6, 100) 

def lennard_Jones_potential(r, epsilon, sigma): # default values for Argon
    V_LJ = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) # Lennard-Jones equation
    return V_LJ

def partition_function(T):
    Z = []
    V_LJ = lennard_Jones_potential(r, epsilon, sigma)
    for temp in T:
        Beta = 1 / (k_B * temp)
        lam = np.sqrt((Beta * h**2) / (2 * np.pi * m))
        C = (1 / h**6) * (1 / lam**6)  # Pre-factor from momenta integration
        Vr = 4 * np.pi * r**2
        exponent = -Beta * V_LJ
        exponent = np.clip(exponent, -700, 700)
        integrand = Vr * np.exp(exponent)
        integral = trapezoid(integrand, r)
        Z_p = C * V * integral
        Z.append(Z_p)
    Z = np.array(Z)
    Z /= np.max(Z)
    return Z



def entropy(Z, T):
    Beta = 1/(k_B* T)
    dBeta = Beta[1] - Beta[0]
    U = -np.gradient(np.log(Z), dBeta)
    return U

def heat_capa(U, T):
    Cv = np.gradient(U, T)
    return Cv

Z = partition_function(T)
U = entropy(Z, T)
Cv = heat_capa(U, T)

plt.plot(T, Cv, label = 'heat capacity', color = 'blue')
plt.xlabel('kelvin')
plt.ylabel('Heat capacity (eV/K)')
plt.title('Heat Capacity vs Temperature')
plt.grid(True)
plt.legend()
plt.show()

#more debugging
print(np.isnan(Cv).sum(), np.isinf(Cv).sum())
print(Cv[:10])
Z = partition_function(T)
print("Z min, max:", np.min(Z), np.max(Z))
print("Any zeros or negatives in Z?", np.any(Z <= 0))
print("NaNs in Z:", np.isnan(Z).sum())