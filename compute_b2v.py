import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from optimize_Argon_dimer import lennard_Jones_potential
import pandas as pd
from writecsv import write_dict_to_csv

#PV = RT(1 + B2v/V + B3v^2/V^2 + ...)
#B2v(T) = -2πNA ∫[0 to ∞] (exp(-u(r)/kT) - 1) r^2 dr
#u(r) = potential energy function
#k = boltzmann constant 0.0019872041 kcal/(mol*K)

# Define potential energy functions
# Hard Sphere
def hardsphere(r, sigma = 3.4):
    if r < sigma:
        V = 1000
    else:
        V = 0
    return V

# Square Well
def squaredwell(r, sigma = 3.4, epsilon = 0.01, lambda_ = 1.5):
    if r < sigma:
        V = 1000
    elif r >= sigma and r < lambda_ * sigma:
        V = -epsilon
    else:
        V = 0
    return V
# Lennard-Jones (load from homework 2_1)

# Calculate B2v using numerical integration
def B2v_function(potential, T):
    k = 0.0019872041 # Boltzmann constant in kcal/(mol*K)
    N = 6.02214076e23 # Avogadro's number
    r = np.linspace(0.0001, 5*3.4, 1000) # Angstroms linespace
    dr = 5*3.4 / 1000 # step size
    u = np.array([potential(rp) for rp in r]) #set up array for each r value in order to caluclate u(r)
    integrand = (np.exp(-u / (k * T)) - 1) * r**2 #given equation
    B2v = -2 * np.pi * N * trapezoid(integrand, r, dr, axis=0) #trapezoidal rule for integration over setp size
    print(B2v)
    return B2v

# Test the B2v function with different potentials at 100K
print("Testing B2v values at 100K") 
print("Lennard-Jones")
B2v_function(lennard_Jones_potential, 100) 
print("Square Well")
B2v_function(squaredwell, 100)
print("Hard Sphere")
B2v_function(hardsphere, 100)

# Plot B2v vs Temperature for each potential
temperatures = np.linspace(0,800,100) # Temperature range from 0K to 800K
potentials = [hardsphere, squaredwell, lennard_Jones_potential] # List of potential functions
potential_names = ["Hard Sphere", "Square Well", "Lennard-Jones"] # Corresponding names for the potentials
B2v_values = {name: [] for name in potential_names} # Dictionary to store B2v values for each potential
# Calculate B2v for each potential at each temperature
for T in temperatures:
    for potential, name in zip(potentials, potential_names): #zip to iterate over both lists simultaneously
        B2v = B2v_function(potential, T) #calculate B2v
        B2v_values[name].append(B2v) #store the result in the dictionary
write_dict_to_csv(temperatures, B2v_values, x_label="Temperature (K)", filename="homework_2_2.csv")

# Plotting
for name in potential_names: 
    plt.plot(temperatures, B2v_values[name], marker='o', label=name)
plt.xlabel("Temperature (K)")
plt.ylabel("B2v")
plt.title("Second Virial Coefficient vs Temperature")
plt.legend()
plt.grid(True)
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_2_2\\B2v_plot.png', dpi=300)
plt.show()

