import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
import pandas as pd
from compute_work_isothermal import ideal_gas
from compute_work_adiabatic import adiabatic_gas
#gives for ideal gas
n = 1.0
R = 8.314
T = 300.0
Vi = 0.1
gamma = 1.4

# unpack work and final volume for the graph and CSV file
work_iso, dv_iso = ideal_gas(Vi, R, T, n)
work_adi, dv_adi = adiabatic_gas(Vi, R, gamma, n, T)

# print to make sure things are working
print("comparision")
print(f"first few line for iso work {work_iso[:10]}\n")
print(f"firt few lines for adi work {work_adi[:10]}\n")

# plot each type of work done as a messure of work vs volume
plt.plot(dv_iso, work_iso, label = 'isothermic', color = 'blue')
plt.plot(dv_adi, work_adi, label = 'adiabatic', color = 'red')
plt.title("work done adiabatic vs isotermic")
plt.xlabel("volume meters cubed")
plt.ylabel("work in joules")
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_1\\gas_law', dpi=300)
plt.show()


df_iso = pd.DataFrame({"dv meter cubed": dv_iso, "Work_isothermal jouls": work_iso})
df_adi = pd.DataFrame({"dv meter cubed": dv_adi, "Work_adiabatic jouls": work_adi})

df_iso.to_csv('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_1\\work_isothermal.csv', index=False)
df_adi.to_csv('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_4_1\\work_adiabatic.csv', index=False)