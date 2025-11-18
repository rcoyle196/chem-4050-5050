import numpy as np 
import random
from scipy.constants import N_A, k, eV


#set up lattice 
def initlize_lattice(size):
    lattice = np.zeros((size, size), dtype=int)
    return lattice

# compute neighbors with periodic boundary conditions
def compute_neighbors(size):
    neighbors = {}
    for x in range(size):
        for y in range(size):
            neighbors[(x,y)] = [
                ((x-1) % size, y % size),
                ((x+1) % size, y % size),
                (x % size, (y-1) % size),
                (x % size, (y+1) % size)
            ]
    return neighbors

# compute interaction energy for a given site and particle type
def compute_interaction_energy(lattice, site, particle, neighbors, epsilon_AA, epsilon_BB, epsilon_AB):
    x, y = site
    interaction_energy = 0
    for neighbor in neighbors[(x,y)]:
        neighbor_state = lattice[neighbor]
        if particle == 1 and neighbor_state == 1:
            interaction_energy += epsilon_AA
        elif particle == 2 and neighbor_state == 2:
            interaction_energy += epsilon_BB
        elif (particle == 2 and neighbor_state ==1) or (particle ==1 and neighbor_state == 2):
            interaction_energy += epsilon_AB
    return interaction_energy

# add, remove, swap moves inside the lattice
def attempt_move(lattice, N_A, N_B, N_empty, neighbors, params):
    # extract parameters
    size = lattice.shape[0]
    N_sites = size * size
    kB = 8.62e-5
    beta = 1 / (kB * params['T'])
    mu_A = params['mu_A']
    mu_B = params['mu_B']
    epsilon_A = params['epsilon_A']
    epsilon_B = params['epsilon_B']
    epsilon_AB = params['epsilon_AB']
    epsilon_AA = params['epsilon_AA']
    epsilon_BB = params['epsilon_BB']
    extract_paramaters = [mu_A, mu_B, epsilon_A, epsilon_B, epsilon_AB, epsilon_AA, epsilon_BB]
    # randomly choose move type
    rand = random.random()
    # perform move addition and caluclate energy change (and if it is accepted)
    if rand < 1/3: #add
        if N_empty == 0:
            return N_A, N_B, N_empty
        site = random.choice(np.argwhere(lattice == 0))
        particle = random.choice([1, 2])
        if particle == 1:
            mu = mu_A
            epsilon = params['epsilon_A']
            N_S = N_A
        else:
            mu = mu_B
            epsilon = params['epsilon_B']
            N_S = N_B
        delta_Ea = epsilon + compute_interaction_energy(lattice, tuple(site), particle, neighbors, params['epsilon_AA'], params['epsilon_BB'], params['epsilon_AB'])
        acc = min(1, (N_empty) /(N_S + 1) * np.exp(-beta * (delta_Ea - mu)))
        r = np.random.rand()
        # perform addition
        if r < acc:
            lattice[site[0], site[1]] = particle
            if particle == 1:
                N_A += 1
            else:
                N_B += 1
            N_empty -= 1  
    # remove move and caluclate energy change (and if it is accepted)  
    elif rand >= 2/3: #remove
        if N_sites - N_empty == 0:
            return N_A, N_B, N_empty
        site = random.choice(np.argwhere(lattice != 0))
        particle = lattice[site[0], site[1]]
        if particle == 1:
            mu = mu_A
            epsilon = epsilon_A
            N_S = N_A
        else:
            mu = mu_B
            epsilon = epsilon_B
            N_S = N_B
        delta_Er = -(epsilon + compute_interaction_energy(lattice, tuple(site), particle, neighbors, params['epsilon_AA'], params['epsilon_BB'], params['epsilon_AB']))
        acc = min(1, N_S /(N_empty + 1) * np.exp(-beta * (delta_Er + mu)))
        r = np.random.rand()
        # perform removal
        if r < acc:
            lattice[site[0], site[1]] = 0
            if particle == 1:
                N_A -= 1
            else:
                N_B -= 1
            N_empty += 1
    else: #swap and caluclate energy change (and if it is accepted)
        site_A = np.argwhere(lattice == 1)
        site_B = np.argwhere(lattice == 2)
        if len(site_A) == 0 or len(site_B) == 0:
            return N_A, N_B, N_empty
        site_A = random.choice(site_A)   
        site_B = random.choice(site_B)
        E_Ao = compute_interaction_energy(lattice, site_A, 1, neighbors, params['epsilon_AA'], params['epsilon_BB'], params['epsilon_AB']) 
        E_Bo = compute_interaction_energy(lattice, site_B, 2, neighbors, params['epsilon_AA'], params['epsilon_BB'], params['epsilon_AB'])
        E_B = compute_interaction_energy(lattice, site_B, 1, neighbors, params['epsilon_AA'], params['epsilon_BB'], params['epsilon_AB']) 
        E_A = compute_interaction_energy(lattice, site_A, 2, neighbors, params['epsilon_AA'], params['epsilon_BB'], params['epsilon_AB'])
        delta_E = (E_A + E_B) - (E_Ao + E_Bo)
        acc = min(1, np.exp(-beta * (delta_E - (mu_B - mu_A))))
        r = np.random.rand()
        # perform swap 
        if r < acc:
            lattice[site_B[0], site_B[1]] = 1
            lattice[site_A[0], site_A[1]] = 2
    return N_A, N_B, N_empty

# run simulation for given size, steps, and parameters and update coverage
def run_simulation(size, n_steps, params):
    lattice = initlize_lattice(size)
    neighbors = compute_neighbors(size)
    N_site = size * size
    N_A = 0
    N_B = 0
    N_empty = N_site
    coverage_A = np.zeros(n_steps)
    coverage_B = np.zeros(n_steps)
    for step in range(n_steps):
        N_A, N_B, N_empty = attempt_move(lattice, N_A, N_B, N_empty, neighbors, params)
        coverage_A[step] = N_A / N_site
        coverage_B[step] = N_B / N_site
    return lattice, coverage_A, coverage_B

# plotting function
def plot_lattice(lattice, ax, title):
    size = lattice.shape[0]
    for x in range(size):
        for y in range(size):
            if lattice[x, y] == 1:
                ax.plot(x + 0.5, y + 0.5, 'ro')
            elif lattice[x, y] == 2:
                ax.plot(x + 0.5, y + 0.5, 'bo')
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticks(np.arange(0, size + 1, 1), minor=True)
    ax.set_yticks(np.arange(0, size + 1, 1), minor=True)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=2)   
    ax.set_title(title)
    return ax



