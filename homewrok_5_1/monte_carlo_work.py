import numpy as np
import matplotlib.pyplot as plt
import pandas as ps
from scipy import stats as st

# psi(r,0) = 1/(4*squrt(pi)*a^3/2)*(r/a)*cos(0)e^-r/(2a)
# S(R) = intergal, intergal, intergal, psi^*(x, y, z + R/2)psi(x,y,x - R/2)dxdydz
# constants given

R = 1 #bohr radi
x = 0 
y = 0 
z_0 = -R/2
z_1 = R/2
z = R
r = np.sqrt(x**2 + y**2 + z**2)
psi_z = (1/(4*np.sqrt(np.pi)*2)*z*np.exp(-r/2))

# function for 2p orbital on z axis
def psi_2p_z(x, y, z):
    r = np.sqrt(x**2+y**2+z**2)
    psi_z = (1/(4*np.sqrt(2*np.pi))*z*np.exp(-r/2))
    return psi_z

# random distrubition
np.random.seed(42)
a = 0
b = 20 
R = 2
N = [10,100,1000,10000,100000, 1000000, 10000000]
ave = []
S_R = []
V = 2*(b)**3
for n in N:
    x = np.random.uniform(a, b, n)
    y = np.random.uniform(a, b, n)
    z = np.random.uniform(a, b, n)
    intergrand= psi_2p_z(x, y, z + R/2) * psi_2p_z(x, y, z -R/2)
    intergral = np.average(intergrand)*V
    S_R.append(intergral)
    

# importance sampling (generate and find which function and shape to use)
x = 0
y = 0
z = np.linspace(-10,10,1000) # two identical, so only worry about one side for importatance sampling
intergrand= psi_2p_z(x, y, z + R/2) * psi_2p_z(x, y, z -R/2)*100 #scale it up to show up on plot
importance_sampling = st.gamma.pdf(z,3,1)


# generate the S_R value with importance sampling
S_R2 = []
for n in N:
    u = 3.5
    x = st.gamma.rvs(u, scale = 1, size = n)
    y = st.gamma.rvs(u, scale = 1, size = n)
    z = st.gamma.rvs(u, scale = 1, size = n)
    num = psi_2p_z(x, y, z + R/2) * psi_2p_z(x, y, z -R/2)
    den = st.gamma.pdf(x, u, scale = 1)* st.gamma.pdf(y, u, scale = 1)* st.gamma.pdf(z, u, scale = 1)
    intergrand2 = (num/den)
    intergral2 = np.average(intergrand2)
    S_R2.append(intergral2)


# generate S_R as function of seperation distanaces
N_2 = 100000
r = np.linspace(0.5, 20, 39)
S_R3 = []
for dr in r:
    u = 3.5
    x = st.gamma.rvs(u, scale = 1, size = N_2)
    y = st.gamma.rvs(u, scale = 1, size = N_2)
    z = st.gamma.rvs(u, scale = 1, size = N_2)
    num = psi_2p_z(x, y, z + dr/2) * psi_2p_z(x, y, z -dr/2)
    den = st.gamma.pdf(x, u, scale = 1)* st.gamma.pdf(y, u, scale = 1)* st.gamma.pdf(z, u, scale = 1)
    intergrand3 = (num/den)
    intergral3 = np.average(intergrand3)
    S_R3.append(intergral3)

