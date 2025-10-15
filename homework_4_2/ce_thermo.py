import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.constants import k, eV

# Constants
k_B = k / eV  # Boltzmann constant in eV/K
T = np.linspace(300, 2000, 1700)
Beta = 1/(k_B * T)

#Ce3+ 
#isolated gives
g_iso = np.array([14])
E_iso = np.array([0])

# spin orbit coupling (SOC) 
# 2 energy states
g_soc = np.array([6,8])
E_soc = np.array([0,.28]) #eV difference between Low and high

# spin orbit and crystal feild splitting (SOC & CFS)
# further split into 5 energy states
g_cfs = np.array([4,2,2,4,2])
E_cfs = np.array([0, .12, .35, .41, .42])

# function to find and store values of Z
def partition(g, E, T):
    Z = []
    for dt in T:
        Beta = 1/(k_B*dt)
        dz =  np.sum(g * np.exp(-Beta* E)) 
        Z.append(dz)
    return np.array(Z)

# function to find internal energy
def internal_energy(Z, T):
    U = []
    Beta = 1/(k_B*T)
    U = -np.gradient(np.log(Z), Beta) 
    return U

#function to find free energy
def free_energy(Z, T):
    F = -k_B * T * np.log(Z)
    return F

# function to find entropy
def entropy(F, T):
    S = -np.gradient(F, T)
    return S

# storing all the data points fround from the fucntions in different systems for CSV and graph
Z_iso = partition(g_iso, E_iso, T)
U_iso = internal_energy(Z_iso, T)
F_iso = free_energy(Z_iso, T)
S_iso = entropy(F_iso, T)

Z_soc = partition(g_soc, E_soc, T)
U_soc = internal_energy(Z_soc, T)
F_soc = free_energy(Z_soc, T)
S_soc = entropy(F_soc, T)

Z_cfs = partition(g_cfs, E_cfs, T)
U_cfs = internal_energy(Z_cfs, T)
F_cfs = free_energy(Z_cfs, T)
S_cfs = entropy(F_cfs, T)

#plot interal vs free energy of each system not meaningful that is why they are not saved graph just wanted to see how they reacted
plt.plot(T, U_iso, label = 'internal energy', color = 'blue', linestyle = '--')
plt.plot(T, F_iso, label = 'free energy', color = 'red', linestyle = '--')
plt.xlabel('K')
plt.ylabel('eV')
plt.title('therodynamic processes of Ce in isolated system vs temperature')
plt.grid(True)
plt.legend()
plt.show()

plt.plot(T, U_soc, label = 'internal energy', color = 'blue', linestyle = '--')
plt.plot(T, F_soc, label = 'free energy', color = 'red', linestyle = '--')
plt.xlabel('K')
plt.ylabel('eV')
plt.title('therodynamic processes of Ce in spin orbit cuppling system vs temperature')
plt.grid(True)
plt.legend()
plt.show()

plt.plot(T, U_cfs, label = 'internal energy', color = 'blue', linestyle = '--')
plt.plot(T, F_cfs, label = 'free energy', color = 'red', linestyle = '--')
plt.xlabel('K')
plt.ylabel('eV')
plt.title('therodynamic processes of Ce in spin orbit cuppling and crystal feild splitting system vs temperature')
plt.grid(True)
plt.legend()
plt.show()

#plot comparsion of entropy between each system
plt.plot(T, S_iso, label = 'entropy isolated', color = 'red')
plt.plot(T, S_soc, label = 'entropy spin orbit cuppling', color = 'blue')
plt.plot(T, S_cfs, label = 'entropy spin orbit and crystal feild', color = 'green')
plt.xlabel('temperature in K')
plt.ylabel('entropy in eV/k')
plt.title('entropy of Ce in three systems')
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\entropy.png', dpi=300)
plt.show()

#plot comparision of interal energy between each system
plt.plot(T, U_iso, label = 'interal energy isolated', color = 'red')
plt.plot(T, U_soc, label = 'interal energy spin orbit cuppling', color = 'blue')
plt.plot(T, U_cfs, label = 'interal energy spin orbit and crystal feild', color = 'green')
plt.xlabel('temperature in K')
plt.ylabel('energy in eV')
plt.title('interal energy of Ce in three systems')
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\internal_energy.png', dpi=300)
plt.show()

#plot free energy comparision between each system
plt.plot(T, F_iso, label = 'free energy isolated', color = 'red')
plt.plot(T, F_soc, label = 'free energy spin orbit cuppling', color = 'blue')
plt.plot(T, F_cfs, label = 'free energy spin orbit and crystal feild', color = 'green')
plt.xlabel('Temperature in K')
plt.ylabel('energy in eV')
plt.title('free energy of Ce in three systems')
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\free_energy.png', dpi=300)
plt.show()

data = {
    "Temperature (K)": T,
    "U_iso (eV)": U_iso,
    "U_soc (eV)": U_soc,
    "U_cfs (eV)": U_cfs,
    "F_iso (eV)": F_iso,
    "F_soc (eV)": F_soc,
    "F_cfs (eV)": F_cfs,
    "S_iso (eV/K)": S_iso,
    "S_soc (eV/K)": S_soc,
    "S_cfs (eV/K)": S_cfs
}

df = pd.DataFrame(data)
df.to_csv('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\Ce_thermo.csv', index=False)
