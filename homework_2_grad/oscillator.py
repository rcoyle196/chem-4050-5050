# V = 1/2mw^2x^2
# V = Vx = D[1-exp(-a(x-re))]^2

import numpy as np
import matplotlib.pyplot as plt

h_bar = 1
m = 1
w = 1
D = 10
B = np.sqrt(1/(2*D))
L = 40
L_linspace = np.linspace(-L/2, L/2, 2000)
dx = L/(2000-1)
V_harmonic = 1/2 * m * w**2 * L_linspace**2
V_morse = D * (1 - np.exp(-B * L_linspace))**2
#construct matrix and hamiltonian
I = np.eye(2000)
I_off_diag = np.eye(2000, k=1) + np.eye(2000, k=-1)
lap = 1/(dx**2) * (-2*I + I_off_diag)
H_harmonic = -h_bar**2/(2*m) * lap + np.diag(V_harmonic)
H_morse = -h_bar**2/(2*m) * lap + np.diag(V_morse)
#print(H_harmonic)
# eigenvalues/vectors harmonic
store_eigenvalues_harmonic = []
store_eigenvectors_harmonic = []
eigenvalues_harmonic, eigenvectors_harmonic = np.linalg.eig(H_harmonic)
eigenvalues = np.argsort(eigenvalues_harmonic)
eigenvalues_harmonic = eigenvalues_harmonic[(eigenvalues)]
eigenvectors_harmonic = eigenvectors_harmonic[:,(eigenvalues)]
store_eigenvalues_harmonic.append(eigenvalues_harmonic)
store_eigenvectors_harmonic.append(eigenvectors_harmonic)

# eigenvalues/vectors morse
store_eigenvalues_morse = []
store_eigenvectors_morse = []
eigenvalues_morse, eigenvectors_morse = np.linalg.eig(H_morse)
eigenvalues = np.argsort(eigenvalues_morse)
eigenvalues_morse = eigenvalues_morse[(eigenvalues)]
eigenvectors_morse = eigenvectors_morse[:,(eigenvalues)]
store_eigenvalues_morse.append(eigenvalues_morse)
store_eigenvectors_morse.append(eigenvectors_morse)

def plot_harmonic():
    for i in range(10):
        plt.plot(L_linspace, eigenvectors_harmonic[:,i], label=f"n={i}")
    plt.title("Harmonic Oscillator Wavefunctions")
    plt.xlabel("x (arbitrary units)")
    plt.ylabel("Wavefunction (ψ)")
    plt.grid(True)
    plt.legend()
    plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_2_grad\\har1.png', dpi=300)
    plt.show()
plot_harmonic()

def plot_morse():
    for j in range(10):
        plt.plot(L_linspace, eigenvectors_morse[:,j], label=f"n={j}")
    plt.title("Morse Oscillator Wavefunctions")
    plt.xlabel("x (arbitrary units)")
    plt.ylabel("Wavefunction (ψ)")
    plt.grid(True)
    plt.legend()
    plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_2_grad\\anh_1.png', dpi=300)
    plt.show()
plot_morse()

def plot_harmonic_2():
    for k in range(10):
        plt.plot(L_linspace, eigenvectors_harmonic[:,k]**2, label=f"n={k}")
    plt.title("Harmonic Oscillator Wavefunctions")
    plt.xlabel("x (arbitrary units)")
    plt.ylabel("Wavefunction (ψ**2)")
    plt.grid(True)
    plt.legend()
    plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_2_grad\\har2.png', dpi=300)
    plt.show()
plot_harmonic_2()

def plot_morse_2():
    for l in range(10):
        plt.plot(L_linspace, eigenvectors_morse[:,l]**2, label=f"n={l}")
    plt.title("Morse Oscillator Wavefunctions")
    plt.xlabel("x (arbitrary units)")
    plt.ylabel("Wavefunction (ψ**2)")
    plt.grid(True)
    plt.legend()
    plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_2_grad\\anh2.png', dpi=300)
    plt.show()
plot_morse_2()
