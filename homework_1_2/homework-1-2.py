import numpy as np

#import quard
molecule_H2 = {
    "H1": np.array([0.0000, 0.0000, 0.0000]), #unit angstrums
    "H2": np.array([0.0000, 0.0000, 0.7414]),
}

molecule_H2O = {
    "O1": np.array([0.0000, 0.0000, 0.1173]),
    "H2": np.array([0.0000, 0.7572, -0.4692]),
    "H3": np.array([0.0000, -0.7572, -0.4692]),
}
molecule_benz = {
    "C1": np.array([0.0000,	1.3970,	0.0000]),
    "C2": np.array([1.2098,	0.6985,	0.0000]),
    "C3": np.array([1.2098,	-0.6985, 0.0000]),
    "C4": np.array([0.0000,	-1.3970, 0.0000]),
    "C5": np.array([-1.2098, -0.6985, 0.0000]),
    "C6": np.array([-1.2098, 0.6985,	0.0000]),
    "H7": np.array([0.0000,	2.4810,	0.0000]),
    "H8": np.array([2.1486,	1.2405,	0.0000]),
    "H9": np.array([2.1486,	-1.2405, 0.0000]),
    "H10": np.array([0.0000, -2.4810, 0.0000]),
    "H11": np.array([-2.1486, -1.2405, 0.0000]),
    "H12": np.array([-2.1486, 1.2405, 0.0000]),
}
print (type(molecule_benz))
#print (f"H2 {molecule_H2} \nH2O{molecule_H2O} \nbenzene {molecule_benz}")

#manual bond length import
bond_length_H2 = np.linalg.norm(molecule_H2["H1"] - molecule_H2["H2"],2)
bond_length_H2O = np.linalg.norm(molecule_H2O["O1"] - molecule_H2O["H2"],2)
bond_length_H20_2 = np.linalg.norm(molecule_H2O["O1"] - molecule_H2O["H3"],2)
bond_length_benz = np.linalg.norm(molecule_benz["C1"] - molecule_benz["C2"],2)

print (f"d_H2 (H1 & H2) {bond_length_H2:.3f}\nd_H2O (H2 & O1) {bond_length_H2O:.3f} \nd_benzene (C1 and C2) {bond_length_benz:.3f}")

#bond lenghts fit
print("testing fit of bonds")
if bond_length_benz < 2:
    print(f"d_benzene {bond_length_benz:.3f}")
else:
    print("warrning bond to long")

if bond_length_H2 < 2:
    print(f"d_H2 {bond_length_H2:.3f}")
else:
    print("warring too long")

if bond_length_H2O < 2:
    print(f"d_H2O {bond_length_H2O:.3f}")
else:
    print("warrning ugh too long")

print()
print()


#setup for water example
set_up_vector_1 = molecule_H2O["O1"] - molecule_H2O["H2"]
set_up_vector_2 = molecule_H2O["O1"] - molecule_H2O["H3"]
dot_H2O = np.dot(set_up_vector_1, set_up_vector_2)
angle_H2O = np.arccos (dot_H2O/(bond_length_H2O * bond_length_H20_2)) * 180/np.pi

print(f"\nangle H2O {angle_H2O:.3f}")

if angle_H2O < 90:
    print ("acute")
elif angle_H2O == 90:
    print("right")
else:
    print("obtuse")

def bond_length_function():
    list_bondlength_H2 = []
    list_bondlength_H2O = []
    list_bondlength_benz = []
    
    print ("\nlength for H2")
    for (index, value) in (molecule_H2.items()): 
        for (index2, value2) in (molecule_H2.items()):
            if index == index2:
                continue
            bond_length_H = np.linalg.norm(value - value2,2) 
            if bond_length_H > 0.1:
                list_bondlength_H2.append(bond_length_H)
            print(f"{index} - {index2} length {bond_length_H}")
        else:
            break
            #print("thats not a bond")

    print (f"\n length for water")      
    for (index, value) in (molecule_H2O.items()): 
        for (index2, value2) in (molecule_H2O.items()):
            if index == index2:
                continue
            bond_length_H2O = np.linalg.norm(value - value2,2)
            
            if bond_length_H2O >= 0.1 and bond_length_H2O <1.5:
                print(f"{index}-{index2} length {bond_length_H2O}")
                list_bondlength_H2O.append(bond_length_H2O)
            else:
                break
               # print(f"{index} and {index2} dont bind")
                

    print (f"\n length for benzene")
    for (index, value) in (molecule_benz.items()): 
        for (index2, value2) in (molecule_benz.items()):
                if index == index2:
                    continue
                bond_length_benz = np.linalg.norm(value - value2,2)
                if bond_length_benz >= 0.1 and bond_length_benz <2:
                    print(f"{index}-{index2} length {bond_length_benz}")
                    list_bondlength_benz.append(bond_length_benz)
                    #print(f"{index} and {index2} dont bind")
    return list_bondlength_H2, list_bondlength_H2O, list_bondlength_benz

list_bondlength_H2, list_bondlength_H2O, list_bondlength_benz = bond_length_function()
uniq_list_H2 = sorted(set(list_bondlength_H2))
uniq_list_H2O = sorted(set(list_bondlength_H2O))
uniq_list_benze = sorted(set(list_bondlength_benz))

print("\n this is a unique array for bondlength found from itteration")
print(f"H2{uniq_list_H2}")
print(f"H2O{uniq_list_H2O}")
print(f"benzene{uniq_list_benze}")
#print(list_bondlength_benz[0])

def bond_angle_function():
    list_bond_angle_H2O = []
    list_bond_angle_benze = []
    print (f"\n angles for water")      
    for (index, value) in (molecule_H2O.items()): 
        for (index2, value2) in (molecule_H2O.items()):
            for(index3, value3) in (molecule_H2O.items()):
                if len({index, index2, index3}) <3:
                    continue
                bond_length_H2O = np.linalg.norm(value - value2,2)
                bond_length_H20_2 = np.linalg.norm(value - value3,2)

                if not (0.1<bond_length_H2O<1 and 0.1<bond_length_H20_2<1):
                    continue
                dot_product_setup = value2 - value
                dot_product_setup_2 = value3 - value
                dot_product = np.dot(dot_product_setup, dot_product_setup_2)
                angle = np.arccos(dot_product/(bond_length_H2O*bond_length_H20_2))*180/np.pi
                print(f"{index}-{index2}-{index3} have angle of {angle}")
                list_bond_angle_H2O.append(angle)
                
    print(f"\n angles for benzene")
    for (index, value) in (molecule_benz.items()): 
        for (index2, value2) in (molecule_benz.items()):
            for(index3, value3) in (molecule_benz.items()):
                if len({index, index2, index3}) <3:
                    continue

                bond_length_benz = np.linalg.norm(value - value2,2)
                bond_length_benz_2 = np.linalg.norm(value - value3,2)

                if not(0.1<bond_length_benz<2 and 0.1<bond_length_benz_2<2):
                    continue
                dot_product_setup = value - value2
                dot_product_setup_2 = value - value3
                dot_product = np.dot(dot_product_setup, dot_product_setup_2)
                angle = np.arccos(dot_product/(bond_length_benz*bond_length_benz_2))*180/np.pi
                print(f"{index}-{index2}-{index3} have angle of {angle}")
                list_bond_angle_benze.append(angle)

    return list_bond_angle_H2O, list_bond_angle_benze

list_bond_angle_H2O, list_bond_angle_benze = bond_angle_function()
unique_angle_list_H2O = sorted(set(list_bond_angle_H2O))
unique_angle_list_benze = sorted(set(list_bond_angle_benze))

print(f"\nangle for water excluding replicate bonds for water \n{unique_angle_list_H2O}")
print(f"\nangle for benzene excluding replicate bonds for benzene \n{unique_angle_list_benze}")
            
