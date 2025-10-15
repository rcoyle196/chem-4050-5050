import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.constants import k, eV, h, atomic_mass
from scipy.integrate import trapezoid
from comp_ther_prop import internal_energy, heat_capa
from comp_part_func import partition_function

#constants
T = np.linspace(10, 1000, 100)

# unwraping the function and assigining values to each 
Z = partition_function(T)
U = internal_energy(Z, T)
Cv = heat_capa(U, T)

#plot heat capacity vs temp
plt.plot(T, Cv, label = 'heat capacity', color = 'blue')
plt.xlabel('temperature in K')
plt.ylabel('energy over temperuature eV/K')
plt.title('Heat Capacity vs Temperature')
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_grad\\heat capacity vs temp', dpi=300)
plt.show()

#plot internal energy vs temp
plt.plot(T, U, label = 'internal energy', color = 'red')
plt.xlabel('temperuature in K')
plt.ylabel('energy in eV')
plt.title('internal energy vs temeperature')
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_grad\\interal energy vs temp', dpi=300)
plt.show()

#generate CSV files
pd.DataFrame({'Temperature (K)': T, 'Partition Function (Z)': Z}).to_csv('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_grad\\partition_function.csv', index=False)
pd.DataFrame({'Temperature (K)': T, 'Heat Capacity(eV/K)': Cv, 'internal energy(eV)': U}).to_csv('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_grad\\inetreal energy and heat capacity.csv', index=False)