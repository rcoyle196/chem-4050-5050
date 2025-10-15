import numpy as np
import matplotlib.pyplot as plt
import pandas as ps
from scipy.constants import k, eV

# Constants
k_B = k / eV  # Boltzmann constant in eV/K
T = np.linspace(300, 2000, 1700)
beta = 1/(k_B * T)

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


# isolated system calulation
def isolated_system(g, E, T, beta):
    Z =  np.sum(g * np.exp(-np.outer(beta, E)), axis=1) 
    U = -np.gradient(np.log(Z), beta) #interal energy
    F = -k_B * T * np.log(Z) #free energy
    S = -np.gradient(F,T) #entropy
    return Z, U, F, S

#SOC
def SOC_system(g, E, T, beta):
    Z = np.sum(g * np.exp(-np.outer(beta, E)), axis=1)
    U = -np.gradient(np.log(Z), beta) 
    F = -k_B * T * np.log(Z)
    S = -np.gradient(F,T)
    return Z, U, F, S

# CFS & SOC
def CFS_system(g, E, T, beta):
    Z = np.sum(g * np.exp(-np.outer(beta, E)), axis=1) 
    U = -np.gradient(np.log(Z), beta) 
    F = -k_B * T * np.log(Z)
    S = -np.gradient(F,T)
    return Z, U, F, S

#unpack each function and their values (I only needed one function)
Z_iso, U_iso, F_iso, S_iso = isolated_system(g_iso, E_iso, T, beta)
Z_soc, U_soc, F_soc, S_soc = SOC_system(g_soc, E_soc, T, beta)
Z_cfs, U_cfs, F_cfs, S_cfs = CFS_system(g_cfs, E_cfs, T, beta)

#plot interal vs free energy of each system
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
plt.xlabel('K')
plt.ylabel('eV/k')
plt.title('entropy of Ce')
plt.grid(True)
plt.legend()
#plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\entropy.png', dpi=300)
plt.show()

#plot comparision of interal energy between each system
plt.plot(T, U_iso, label = 'interal energy isolated', color = 'red')
plt.plot(T, U_soc, label = 'interal energy spin orbit cuppling', color = 'blue')
plt.plot(T, U_cfs, label = 'interal energy spin orbit and crystal feild', color = 'green')
plt.xlabel('K')
plt.ylabel('eV')
plt.title('interal energy of Ce')
plt.grid(True)
plt.legend()
#plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\internal_energy.png', dpi=300)
plt.show()

#plot free energy comparision between each system
plt.plot(T, F_iso, label = 'free energy isolated', color = 'red')
plt.plot(T, F_soc, label = 'free energy spin orbit cuppling', color = 'blue')
plt.plot(T, F_cfs, label = 'free energy spin orbit and crystal feild', color = 'green')
plt.xlabel('K')
plt.ylabel('eV')
plt.title('free energy of Ce')
plt.grid(True)
plt.legend()
#plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\free_energy.png', dpi=300)
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

df = ps.DataFrame(data)
df.to_csv('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_2\\Ce_thermo.csv', index=False)