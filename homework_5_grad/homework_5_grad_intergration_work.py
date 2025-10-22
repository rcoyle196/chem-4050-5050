import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st


R = 1.4 #bhor radius

def psi_1s(x, y, z, Z=1, ao=1):
    r = np.sqrt(x**2 + y**2 + z**2)
    #psi = (1/np.sqrt(np.pi*ao**3)) * np.exp(-r/ao)
    psi_Z = (1/np.sqrt(np.pi*ao**3)) * Z**(3/2) * np.exp(-Z*r/ao)
    return psi_Z

def laplacian_psi_1s(x, y, z, Z=1, ao=1):
    r = np.sqrt(x**2 + y**2 + z**2)
    psi = psi_1s(x, y, z, Z, ao)
    #laplacian = (1/ao**2)*psi + (2/r)*(psi*1/ao)
    laplacian_Z = (Z**2/ao**2)*psi + (2/r)*(-Z/ao)*psi
    return laplacian_Z

np.random.seed(42)
a = -7
b = 7 
R = 2
N = [10,100,1000,10000,100000, 1000000, 10000000]
kii = []
V = (2*(b-a))**2
for n in N:
    x = np.random.uniform(a, b, n)
    y = np.random.uniform(a, b, n)
    z = np.random.uniform(a, b, n)
    intergrand= -.5 * psi_1s(x, y, z,) * laplacian_psi_1s(x, y, z)
    intergral = np.average(intergrand)*V
    kii.append(intergral)

# find importance sampling function
# two identical, so only worry about one side for importatance sampling
x = np.linspace(0,2,1000)
y = 0
z = 0
intergrand= -.5 * psi_1s(x, y, z,) * laplacian_psi_1s(x, y, z) /10
importance_sampling = st.gamma.pdf(x, .5, scale = .5)

# generate the S_R value with importance sampling
kii_2 = []
for n in N:
    u = 1
    x = st.gamma.rvs(u, scale = .5, size = n)
    y = st.gamma.rvs(u, scale = .5, size = n)
    z = st.gamma.rvs(u, scale = .5, size = n)
    num = -.5*psi_1s(x, y, z) * laplacian_psi_1s(x,y,z)
    den = st.gamma.pdf(x, u, scale = .5)* st.gamma.pdf(y, u, scale = .5)* st.gamma.pdf(z, u, scale = .5)
    intergrand2 = (num/den)
    intergral2 = np.average(intergrand2)
    kii_2.append(intergral2)



np.random.seed(42)
a = -7
b = 7 
R = 2
N = [10,100,1000,10000,100000, 1000000, 10000000]
kij = []
V = (2*(b-a))**2
for n in N:
    x = np.random.uniform(a, b, n)
    y = np.random.uniform(a, b, n)
    z = np.random.uniform(a, b, n)
    intergrand= -.5*psi_1s(x, y, z + R/2) *laplacian_psi_1s(x, y, z -R/2)
    intergral = np.average(intergrand)*V
    kij.append(intergral)

# find importance sampling function
x = 0
y = 0
z = np.linspace(0,1,1000) # identical left and right side of intergral
intergrand_kij = -.5*psi_1s(x, y, z) * laplacian_psi_1s(x, y, z)/10 #scale it up to show up on plot
importance_sampling_kij = st.gamma.pdf(z,.05,scale = 1)

# importance sampling to find kij
kij_2 = []
for n in N:
    u_kij = .05
    x = st.gamma.rvs(u_kij, scale = .5, size = n)
    y = st.gamma.rvs(u_kij, scale = .5, size = n)
    z = st.gamma.rvs(u_kij, scale = .5, size = n)
    num = -.5*psi_1s(x, y, z + R/2) * laplacian_psi_1s(x, y, z -R/2)
    den = st.gamma.pdf(x, u, scale = .5)* st.gamma.pdf(y, u, scale = .5)* st.gamma.pdf(z, u, scale = .5)
    intergrand_kij_2 = (num/den)
    intergral_kij_2 = np.average(intergrand_kij_2)
    kij_2.append(intergral_kij_2)