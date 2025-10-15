import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from optimize_Argon_dimer import lennard_Jones_potential #importing the Lennard-Jones potential function from the dimer file

# Part 2
# Define the potential energy function for a three-atom system
def potential_energy(cord1):
    [r1, x, y] = cord1 # r1 is distance between atom 1 and 2, (x,y) are coordinates of atom 3
    pos1 = np.array([0,0]) # position of atom 1
    pos2 = np.array([r1,0]) # position of atom 2
    pos3 = np.array([x,y]) # position of atom 3
    d1 = np.linalg.norm(pos3 - pos1,2) # distance between atom 1 and 3
    d2 = np.linalg.norm(pos3 - pos2,2) # distance between atom 2 and 3
    d3 = np.linalg.norm(pos2 - pos1,2) # distance between atom 1 and 2
    # use the Lennard-Jones potential to calculate the total potential energy of the system
    Vt = lennard_Jones_potential(d1) + lennard_Jones_potential(d2) + lennard_Jones_potential(d3)
    return Vt

#minimize function to find the optimal geometry
result2 = minimize( 
    fun=potential_energy, #function to minimize
    x0 = [4,1,2], #initial guess for [r1, x, y]
    method="Nelder-Mead", #optimization method
    tol=1e-20, #tolerance for termination
    )

# print(result2["x"[0]])
# make a function to take the x-values found at min point (r1, x, y) and calculate the distances and angle
r1 = result2["x"][0] # distance between atom 1 and 2
x = result2["x"][1] # x-coordinate of atom 3
y = result2["x"][2] # y-coordinate of atom 3
def geometry(result2):
    pos1 = np.array([0,0]) # position of atom 1
    pos2 = np.array([r1,0]) # position of atom 2
    pos3 = np.array([x,y]) # position of atom 3
    d1 = np.linalg.norm(pos3 - pos1,2) # distance between atom 1 and 3
    d2 = np.linalg.norm(pos3 - pos2,2) # distance between atom 2 and 3
    d3 = np.linalg.norm(pos2 - pos1,2) # distance between atom 1 and 2

    # calculate the angle between the three atoms using the dot product

    dot = np.dot((pos1 - pos2),(pos1 - pos3)) # dot product between vectors 1-2 and 1-3
    angle = np.arccos(dot/(d1*d2)) * (180/np.pi) # angle in degrees
     # print the results
    print(f"the distance between atom 1 and 2 is {d3}")
    print(f"the distance between atom 1 and 3 is {d1}")
    print(f"the distance between atom 2 and 3 is {d2}")
    print(f"the angle between atom 1, 3, and 2 is {angle}")
    print(f"the potential energy of the system is {potential_energy(result2['x'])}")
    return d1, d2, d3, angle

# call the geometry function and print the results    
print(geometry(result2))
print("Argon dimer is a near perfect equilateral triangle within my pramenters")

#xyz file
import os
from writexyz import write_xyz_if_not_exists


atoms = [
    ("Ar", 0.000, 0.000, 0.000), # position of atom 1
    ("Ar", r1, 0.000, 0.000), # position of atom 2
    ("Ar", x, y, 0.000) # position of atom 3
]
write_xyz_if_not_exists("homework_2_1\\Argon_trimer_optimized.xyz", atoms)
