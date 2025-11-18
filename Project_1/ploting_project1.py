import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import N_A, k, eV
from monte_carlo_project_1 import plot_lattice, run_simulation

# Constants
#ideal  
En_I = -0.1
Eh_I = -0.1 
Enn_I = 0
Ehh_I = 0
Enh_I = 0

# repulsive
En_r = -0.1
Eh_r = -0.1
Enn_r = 0.05
Ehh_r = 0.05
Enh_r = 0.05

# attractive
En_a = -0.1
Eh_a = -0.1
Enn_a = -0.05
Ehh_a = -0.05
Enh_a = -0.05

# immiscible
En_im = -0.1
Eh_im = -0.1
Enn_im = -0.05
Ehh_im = -0.05
Enh_im = 0.05

# like dissolves unlikes
En_l = -0.1
Eh_l = -0.1
Enn_l = 0.05
Ehh_l = 0.05
Enh_l = -0.05

# Parameters for simulations
size = 4
n_steps = 1000
mu_B = -0.1
mus_A = np.linspace(-0.2, 0, 7)
Ts = np.linspace(1, 1199, 7)
params_ideal = []
params_repulsive = []
params_attractive = []
params_immiscible = []
params_like_dissolves_unlikes = []


# generate changes in energy and temperature parameters and run simulations for each condition
# store parameter sets for phase diagram plotting
# plot example lattice for each condition (final state)
# ideal condition
for ideal in [(En_I, Eh_I, Enn_I, Ehh_I, Enh_I)]:
    for mu_A in mus_A:
        for T in Ts:
            param_set = {
                'mu_A': mu_A,
                'mu_B': mu_B,
                'T': T,
                'epsilon_A': En_I,
                'epsilon_B': Eh_I,
                'epsilon_AA': Enn_I,
                'epsilon_BB': Ehh_I,
                'epsilon_AB': Enh_I
                }
            params_ideal.append(param_set)


# repulsive condition
for repulsive in [(En_r, Eh_r, Enn_r, Ehh_r, Enh_r)]:
    for mu_A in mus_A:
        for T in Ts:
            param_set = {
                'mu_A': mu_A,
                'mu_B': mu_B,
                'T': T,
                'epsilon_A': En_r,
                'epsilon_B': Eh_r,
                'epsilon_AA': Enn_r,
                'epsilon_BB': Ehh_r,
                'epsilon_AB': Enh_r
                }
            params_repulsive.append(param_set)

# attractive condition
for attractive in [(En_a, Eh_a, Enn_a, Ehh_a, Enh_a)]:
    for mu_A in mus_A:
        for T in Ts:
            param_set = {
                'mu_A': mu_A,
                'mu_B': mu_B,
                'T': T,
                'epsilon_A': En_a,
                'epsilon_B': Eh_a,
                'epsilon_AA': Enn_a,
                'epsilon_BB': Ehh_a,
                'epsilon_AB': Enh_a
                }
            params_attractive.append(param_set)

# immiscible condition
for immiscible in [(En_im, Eh_im, Enn_im, Ehh_im, Enh_im)]:
    for mu_A in mus_A:
        for T in Ts:
            param_set = {
                'mu_A': mu_A,
                'mu_B': mu_B,
                'T': T,
                'epsilon_A': En_im,
                'epsilon_B': Eh_im,
                'epsilon_AA': Enn_im,
                'epsilon_BB': Ehh_im,
                'epsilon_AB': Enh_im
                }
            params_immiscible.append(param_set)

# like dissolves unlikes condition
for like_dissolves_unlikes in [(En_l, Eh_l, Enn_l, Ehh_l, Enh_l)]:
    for mu_A in mus_A:
        for T in Ts:
            param_set = {
                'mu_A': mu_A,
                'mu_B': mu_B,
                'T': T,
                'epsilon_A': En_l,
                'epsilon_B': Eh_l,
                'epsilon_AA': Enn_l,
                'epsilon_BB': Ehh_l,
                'epsilon_AB': Enh_l
                }
            params_like_dissolves_unlikes.append(param_set)

# generate function for plotting phase diagrams
def plot_phase_diagram(mus_A, Ts, coverage_A, coverage_B, title):
    fig, ax = plt.subplots(2, 3, figsize=(13, 8))

    # plot coverage for A (hydrogen)
    ax[0,0].pcolormesh(mus_A, Ts, coverage_A.T, shading='auto', cmap='Reds')
    ax[0,0].set_title(f'Coverage Hydrogen')
    ax[0,0].set_xlabel('Chemical Potential Hydrogen')
    ax[0,0].set_ylabel('Temperature (K)')
    ax[0,0].invert_yaxis()
    colorbar_A = fig.colorbar(plt.cm.ScalarMappable(cmap='Reds'), ax=ax[0,0])
    colorbar_A.set_label('Hydrogen')

    # plot coverage for B (nitrogen)
    ax[0,1].pcolormesh(mus_A, Ts, coverage_B.T, shading='auto', cmap='Blues')
    ax[0,1].set_title(f'Coverage Nitrogen')
    ax[0,1].set_xlabel('Chemical Potential Hydrogen')
    ax[0,1].set_ylabel('Temperature (K)')
    ax[0,1].invert_yaxis()
    colorbar_B = fig.colorbar(plt.cm.ScalarMappable(cmap='Blues'), ax=ax[0,1])
    colorbar_B.set_label('Nitrogen')

    # plot total coverage
    total_coverage = coverage_A + coverage_B
    ax[0,2].pcolormesh(mus_A, Ts, total_coverage.T, shading='auto', cmap='Purples')
    ax[0,2].set_title(f'Total Coverage')
    ax[0,2].set_xlabel('Chemical Potential Hydrogen')
    ax[0,2].set_ylabel('Temperature (K)')
    ax[0,2].invert_yaxis()
    colorbar_total = fig.colorbar(plt.cm.ScalarMappable(cmap='Purples'), ax=ax[0,2])
    colorbar_total.set_label('Total Coverage')

    # plot example lattices

    ax[1,0] = plot_lattice(final_lattice[0, 3], ax[1,0], f"mu_A = -0.2 eV, T = 1 K")
    ax[1,1] = plot_lattice(final_lattice[3, 3], ax[1,1], f"mu_A = -0.1 eV, T = 600 K")
    ax[1,2] = plot_lattice(final_lattice[6, 3], ax[1,2], f"mu_A = -0 eV, T = 1,999 K")

    # generate visual 
    plt.suptitle(f'Phase Diagram - {title}', fontsize=16)
    plt.tight_layout()
    return ax


# Example of plotting phase diagram for ideal case
np.random.seed(42)
final_lattice = np.zeros((len(mus_A), len(Ts), size, size), dtype=int)
coverage_A_ideal = np.zeros((len(mus_A), len(Ts)))
coverage_B_ideal = np.zeros((len(mus_A), len(Ts)))

# loop through ideal parameters
for i, param_set in enumerate(params_ideal):
    i_mu = i // len(Ts)
    i_T = i % len(Ts)
    lattice, cov_A, cov_B = run_simulation(size, n_steps, param_set)
    final_lattice[i_mu, i_T] = lattice
    mean_cov_A = np.mean(cov_A[int(n_steps/2):])
    mean_cov_B = np.mean(cov_B[int(n_steps/2):])
    coverage_A_ideal[i_mu, i_T] = mean_cov_A
    coverage_B_ideal[i_mu, i_T] = mean_cov_B
plot_phase_diagram(mus_A, Ts, coverage_A_ideal, coverage_B_ideal, 'Ideal Mixture')
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_1\\ideal_mixture', dpi=300)
plt.show()

# Example of plotting phase diagram for repulsive case
coverage_A_repulsive = np.zeros((len(mus_A), len(Ts)))
coverage_B_repulsive = np.zeros((len(mus_A), len(Ts)))

# loop through repulsive parameters
for i, param_set in enumerate(params_repulsive):
    i_mu = i // len(Ts)
    i_T = i % len(Ts)
    lattice, cov_A, cov_B = run_simulation(size, n_steps, param_set)
    final_lattice[i_mu, i_T] = lattice
    mean_cov_A = np.mean(cov_A[int(n_steps/2):])
    mean_cov_B = np.mean(cov_B[int(n_steps/2):])
    coverage_A_repulsive[i_mu, i_T] = mean_cov_A
    coverage_B_repulsive[i_mu, i_T] = mean_cov_B
plot_phase_diagram(mus_A, Ts, coverage_A_repulsive, coverage_B_repulsive, 'Repulsive Interactions')
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_1\\repulsive_mixture', dpi=300)
plt.show()

# Example of plotting phase diagram for attractive case
coverage_A_attractive = np.zeros((len(mus_A), len(Ts)))
coverage_B_attractive = np.zeros((len(mus_A), len(Ts)))

# loop through attractive parameters
for i, param_set in enumerate(params_attractive):
    i_mu = i // len(Ts)
    i_T = i % len(Ts)
    lattice, cov_A, cov_B = run_simulation(size, n_steps, param_set)
    final_lattice[i_mu, i_T] = lattice
    mean_cov_A = np.mean(cov_A[int(n_steps/2):])
    mean_cov_B = np.mean(cov_B[int(n_steps/2):])
    coverage_A_attractive[i_mu, i_T] = mean_cov_A
    coverage_B_attractive[i_mu, i_T] = mean_cov_B
plot_phase_diagram(mus_A, Ts, coverage_A_attractive, coverage_B_attractive, 'Attractive Interactions')
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_1\\attraction_mixture', dpi=300)
plt.show()

# Example of plotting phase diagram for immiscible case
coverage_A_immiscible = np.zeros((len(mus_A), len(Ts)))
coverage_B_immiscible = np.zeros((len(mus_A), len(Ts)))

# loop through immiscible parameters
for i, param_set in enumerate(params_immiscible):
    i_mu = i // len(Ts)
    i_T = i % len(Ts)
    lattice, cov_A, cov_B = run_simulation(size, n_steps, param_set)
    final_lattice[i_mu, i_T] = lattice
    mean_cov_A = np.mean(cov_A[int(n_steps/2):])
    mean_cov_B = np.mean(cov_B[int(n_steps/2):])
    coverage_A_immiscible[i_mu, i_T] = mean_cov_A
    coverage_B_immiscible[i_mu, i_T] = mean_cov_B
plot_phase_diagram(mus_A, Ts, coverage_A_immiscible, coverage_B_immiscible, 'Immiscible Species')
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_1\\Immiscible', dpi=300)
plt.show()

# Example of plotting phase diagram for like dissolves unlikes case
coverage_A_like_dissolves_unlikes = np.zeros((len(mus_A), len(Ts)))
coverage_B_like_dissolves_unlikes = np.zeros((len(mus_A), len(Ts)))

# loop through like dissolves unlikes parameters
for i, param_set in enumerate(params_like_dissolves_unlikes):
    i_mu = i // len(Ts)
    i_T = i % len(Ts)
    lattice, cov_A, cov_B = run_simulation(size, n_steps, param_set)
    final_lattice[i_mu, i_T] = lattice
    mean_cov_A = np.mean(cov_A[int(n_steps/2):])
    mean_cov_B = np.mean(cov_B[int(n_steps/2):])
    coverage_A_like_dissolves_unlikes[i_mu, i_T] = mean_cov_A
    coverage_B_like_dissolves_unlikes[i_mu, i_T] = mean_cov_B
plot_phase_diagram(mus_A, Ts, coverage_A_like_dissolves_unlikes, coverage_B_like_dissolves_unlikes, 'Like Dissolves Unlikes')
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_1\\Like_dissolves_unlike', dpi=300)
plt.show()