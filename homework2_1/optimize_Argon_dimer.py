import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import minimize

# Part 1 
# Define the Lennard-Jones potential function
def lennard_Jones_potential(r, epsilon = 0.01, sigma = 3.4): # default values for Argon
    V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) # Lennard-Jones equation
    return V

# use optimization to find min or equilibrium
r = np.linspace(3, 6, 100) # range of distances in Angstroms
result = minimize(
    fun=lennard_Jones_potential, # function to minimize
    x0 = 4, # initial guess
    method="Nelder-Mead", # optimization method
    tol=1e-6 # tolerance for termination
)
print (f" this is equilibrium distance {result["x"][0]}") # equilibrium in Angstroms 
print ("The angle between the two atoms is 180 degrees since it is a dimer")

# Plot the Lennard-Jones potential
plt.plot(r, lennard_Jones_potential(r), linestyle="--") # plot the potential
plt.scatter(result["x"], lennard_Jones_potential(result["x"]), color="red", label="equilibrium point") # plot the equilibrium point
plt.xlabel("r (Angstroms)") # x-axis label
plt.ylabel("V (potential energy)") # y-axis label
plt.title("Lennard-Jones Potential") # title
plt.legend() # legend
plt.grid(True) # grid
plt.show() # show the plot

# xyz file
import os
from writexyz import write_xyz_if_not_exists
r_Ar = result["x"][0] # equilibrium distance
atoms = [
    ("Ar",0.000, 0.000, 0.000), # position of atom 1
    ("Ar",r_Ar, 0.000, 0.000)  # position of atom 2
    ]
write_xyz_if_not_exists("homework_2_1\\argon_dimer_optimized.xyz", atoms)

