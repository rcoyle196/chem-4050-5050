import matplotlib.pyplot as plt
import numpy as np

#test of five 
#n_5 = 5
#I_5 = np.eye(n_5)
#I_off_diag_5 = np.eye(n_5, k=1) + np.eye(n_5, k=-1)
#lap_5 = (-2*I_5+I_off_diag_5)
#print(lap_5)

#setting up varriables
hbar_r = 1 
mass_e_r = 1 #mass of electron
L_r = 1 # (Bohr radii)
neg_L_2 = -.5
L_2 = .5
n = 2000

#space grid and dx
space_grid = np.linspace(neg_L_2, L_2, n)
dx = L_r/(n-1)

# identity matrix
I = np.eye(n)
I_off_diag = np.eye(n, k=1) + np.eye(n, k=-1)
lap_m = (-2*I+I_off_diag)
lap = (lap_m)/(dx**2)
#print(lap)

#hamilition
H = -1*(1/2)*lap
#print(H)

#eigenvalues/vectors
eigenvalues, eigenvectors = np.linalg.eig(H)
eigenvalues = np.argsort(eigenvalues)
eigenvectors = eigenvectors[:, eigenvalues]

#extract first 7 energy levels
def waves():
    for i in range(8):
        print(eigenvalues[:i])

#plot 5 lowest energy state wave functions#for i in range():
    for j in range(6):
        plt.plot(space_grid, eigenvectors[:,j], color = 'green', linewidth = '2', linestyle = ':' )
        plt.title (f"energy state {j}")
        plt.xlabel ('dx \nangstroms or bhor radii')
        plt.ylabel ('psi')
        plt.grid(True)
        plt.legend()
        plt.show()
        
    for k in range(6):
        plt.plot(space_grid, eigenvectors[:,k]**2*1000, color = 'blue', linewidth = '2', linestyle = '-' )
        plt.title (f"normilized energy density {k}")
        plt.xlabel ('dx \n angstroms or bhor radii')
        plt.ylabel ('psi')
        plt.grid(True)
        plt.legend()
        plt.show()

waves()