import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from optimize_Argon_dimer import lennard_Jones_potential #importing the Lennard-Jones potential function from the dimer file
from optimize_Argon_trimer import geometry
#importing from homework 1_2
from homework_1_2 import bond_angle_function, bond_length_function

#grabbing and sort the bond lengths from homework 1_2
list_bondlength_H2, list_bondlength_H2O, list_bondlength_benz = bond_length_function()
uniq_list_H2 = sorted(set(list_bondlength_H2)) 
uniq_list_H2O = sorted(set(list_bondlength_H2O))
uniq_list_benze = sorted(set(list_bondlength_benz))
print("\n this is a unique array for bondlength found from itteration")
print(f"H2{uniq_list_H2}")
print(f"H2O{uniq_list_H2O}")
print(f"benzene{uniq_list_benze}")

# grabbing and sort the bond angles from homework 1_2
list_bond_angle_H2O, list_bond_angle_benze = bond_angle_function()
unique_angle_list_H2O = sorted(set(list_bond_angle_H2O))
unique_angle_list_benze = sorted(set(list_bond_angle_benze))
print(f"\nangle for water excluding replicate bonds for water \n{unique_angle_list_H2O}")
print(f"\nangle for benzene excluding replicate bonds for benzene \n{unique_angle_list_benze}")

#using the lennard jones potential on the unique bond lengths found from homework 1_2
def lennard_jones_geometry():
    water_length = []
    benzene_length = []
    for i in uniq_list_benze:
        lennard_Jones_potential(i)
        benzene_length.append(i)
    for j in list_bondlength_H2O:
        lennard_Jones_potential(j)
        water_length.append(j)
    print (f"this is water {water_length}")
    print (f"this is benzene {benzene_length}")
    # setting up arrays for benzene to find angle
    pos1_b = np.array([0,0])
    pos2_b = np.array([benzene_length[0],0])
    pos3_b = np.array([benzene_length[1],benzene_length[2]])
    d1_b = np.linalg.norm(pos3_b - pos1_b,2)
    d2_b = np.linalg.norm(pos3_b - pos2_b,2)
    d3_b = np.linalg.norm(pos2_b - pos1_b,2)
    dot = np.dot((pos1_b - pos2_b),(pos1_b - pos3_b)) # dot product between vectors 1-2 and 1-3
    angle = np.arccos(dot/(d1_b*d2_b)) * (180/np.pi) # angle in degrees
    print(f"in benzene")
    print(f"the distance between atom 1 and 2 is {d3_b}")
    print(f"the distance between atom 1 and 3 is {d1_b}")
    print(f"the distance between atom 2 and 3 is {d2_b}")
    print(f"the angle between atom 1, 3, and 2 is {angle}")
    # setting up arrays in water to get angle 
    pos1_w = np.array([0,0])
    pos2_w = np.array([water_length[0],0])
    pos3_w = np.array([water_length[1],water_length[2]])
    d1_w = np.linalg.norm(pos3_w - pos1_w,2)
    d2_w = np.linalg.norm(pos3_w - pos2_w,2)
    d3_w = np.linalg.norm(pos2_w - pos1_w,2)
    dotw = np.dot((pos1_w - pos2_w),(pos1_w - pos3_w)) # dot product between vectors 1-2 and 1-3
    anglew = np.arccos(dotw/(d1_b*d2_b)) * (180/np.pi) # angle in degrees
    print(f"in water")
    print(f"the distance between atom 1 and 2 is {d3_w}")
    print(f"the distance between atom 1 and 3 is {d1_w}")
    print(f"the distance between atom 2 and 3 is {d2_w}")
    print(f"the angle between atom 1, 3, and 2 is {anglew}")

print(lennard_jones_geometry())






#I am not sure if this is what you meant by using the lennard jones potential
#did it mean just from part 1?
#would not make sense to use for bonded compounds so I don't fully understand the question
#potential_energy(list_bondlength_H2O[0:3])
#result3 = minimize(
 #   fun=potential_energy,
 #   x0 = [2,1,2],
 #   method="Nelder-Mead",
 #   tol=1e-20,
 #   )
#geometry(result3)
#print(geometry(result3))

