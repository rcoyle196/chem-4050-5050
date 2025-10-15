import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import k, eV, h, atomic_mass, c
from scipy.integrate import trapezoid

# find Cv max
# dimer of LJ in cubic box
# LJ is given by U = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) 
# sum of U is the total potential energies of the system
# Z = (1/h**6)*(1/lamba**6) intergal dx,dy,dz,dx2,dy2,dz2 e^-(-beta*U(x,y,z,x2,y2,z2))
# switch to sphereical coordinate (V/2*lambda) intergal (a -> b) 4*pi*r^2 e^(-Beta*U)
# lambda = h/squrt(2*pi*m*k_b*T)

#constants/givens and adjust units
k_B = k/eV #eV/K
h_ev = h/eV   #eV*s
epsilon = 0.0103 #ev
sigma = 3.4 #ang
T = np.linspace(10, 1000, 100) #k
V = 1000 #ang^3
Beta = 1/(k_B* T) #1/eV
m = (39.95* atomic_mass*c**2)/eV #mass to eV
lam = h_ev/np.sqrt(2*np.pi*m*k_B*T) #1/eV
r = np.linspace(3, 6, 100) #ang


# setting up lennard jones potential (see optimize_Argon_dimer)
def lennard_Jones_potential(r, epsilon, sigma): # default values for Argon
    U = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) # Lennard-Jones equation
    return U

# setting up partition_function using spericah coordinates
def partition_function(T):
    Z = []
    U = lennard_Jones_potential(r, epsilon, sigma)
    for dt in T:
        Beta = 1/(k_B*dt)
        lam = h_ev/np.sqrt(2*np.pi*m*k_B*dt)
        Cnum = V*4*np.pi
        Cdom = 2*lam**6
        C = Cnum/Cdom
        intergal = np.exp(-Beta*U)*r**2
        dz = trapezoid(intergal, r) * C
        Z.append(dz)
    return np.array(Z) 

