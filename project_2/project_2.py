def paramaters():
    # Simulation parameters
    dt = 0.0001  # Time step
    total_steps = 10000  # Number of steps
    box_size = 100.0  # Size of the cubic box
    k = 5  # Spring constant #k5 
    mass = 1.0  # Particle mass
    r0 = 1.0  # Equilibrium bond length
    target_temperature = 0.1  # Target temperature
    rescale_interval = 100  # Steps between velocity rescaling
    n_particles = 10  # Number of particles p=25
    epsilon_repulsive = .25  # Depth of repulsive LJ potential e.75
    epsilon_attractive =   0.5  # Depth of attractive LJ potential
    sigma = 1.0  # LJ potential parameter
    kb = 1.0  # Boltzmann constant
    return (n_particles, box_size, r0, k, sigma, epsilon_attractive, 
            epsilon_repulsive, mass, target_temperature, dt, 
            total_steps, rescale_interval, kb)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


parmas = paramaters()
(n_particles, box_size, r0, k, sigma, epsilon_attractive, 
epsilon_repulsive, mass, target_temperature, dt, 
total_steps, rescale_interval, kb) = parmas

# perodic boundry
def apply_pbc(position, box_size):
    return position % box_size

def min_image(positionI, positionJ, box_size):
    displacment = positionI - positionJ
    displacment -= box_size * np.round(displacment / box_size)
    return displacment


#intial position and Velo
def initialize_chain(n_particles, box_size, r0):
    positions = np.zeros((n_particles, 3))
    current_position = np.array([box_size / 2, box_size / 2, box_size / 2])
    positions[0] = current_position
    for i in range(1, n_particles):
        direction = np.random.normal(0, 1, 3)
        direction /= np.linalg.norm(direction)
        next_position = current_position + r0 * direction
        next_position = apply_pbc(next_position, box_size)
        positions[i] = next_position
        current_position = next_position
    return positions

def initialize_velocities(n_particles, target_temp, mass):
    velocity = np.random.normal(0, np.sqrt(kb * target_temp / mass), (n_particles, 3))
    velocity -= np.mean(velocity, axis=0)
    return velocity



# harmonic functions
def compute_harmonic_forces(positions, k, r0, box_size):
    forces = np.zeros_like(positions)
    for i in range(len(positions) - 1):
        displacment = min_image(positions[i+1], positions[i], box_size)
        distance = np.linalg.norm(displacment)
        force_magnitude = -1 *k * (distance - r0)
        force = force_magnitude * (displacment / distance)
        forces[i] -= force
        forces[i + 1] += force
    return forces

# lenard Jones
def compute_lenard_jones_forces(positions, epsilon, sigma, box_size, interaction_type):
    forces = np.zeros_like(positions)
    # avoid stacking
    cutoff_r = 2**(1/6) * sigma 
    cutoff_a = 2.5 * sigma
    # loops to determine LJ force and type
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            if interaction_type == 'repulsive' and abs(i-j) == 2:
                eps = epsilon_repulsive
                cutoff = cutoff_r
            elif interaction_type == 'attractive' and abs(i-j) > 2:
                eps = epsilon_attractive
                cutoff = cutoff_a
            else: 
                continue
            displacment = min_image(positions[j], positions[i], box_size)
            distance = np.linalg.norm(displacment)
            if distance < cutoff and distance > 1e-12:
                force_magnitude = 24 * eps * ((sigma/distance)**12 - 0.5 * (sigma/distance)**6) / distance
                force = force_magnitude * (displacment / distance)
                forces[i] -= force
                forces[j] += force
    return forces

def compute_forces(positions, k, r0, box_size, epsilon_repulsive, epsilon_attractive, sigma):

    forces_harmonic = compute_harmonic_forces(positions, k, r0, box_size)
    forces_repulsive = compute_lenard_jones_forces(positions, epsilon_repulsive, sigma, box_size, 'repulsive')
    forces_attractive = compute_lenard_jones_forces(positions, epsilon_attractive, sigma, box_size, 'attractive')
    total_forces = forces_harmonic + forces_repulsive + forces_attractive
    return total_forces

# velocity verlet integration
def velocity_verlet(positions, velocities, forces, dt, mass, box_size):
    velocities += 0.5 * forces / mass * dt
    positions += velocities * dt
    positions = apply_pbc(positions, box_size)  
    forces_new = compute_forces(positions, k, r0, box_size, epsilon_repulsive, epsilon_attractive, sigma)  
    velocities += 0.5 * forces_new / mass * dt
    return positions, velocities, forces_new

# velocity rescaling thermostat
def rescale_velocities(velocities, target_temp, mass):
    kinetic_energy = 0.5 * mass * np.sum(np.linalg.norm(velocities, axis=1)**2)
    n_particles = velocities.shape[0]
    current_temp = (2/3) * kinetic_energy / (n_particles * kb) 
    scaling_factor = np.sqrt(target_temp / current_temp)
    velocities *= scaling_factor
    return velocities

# analyze results
def compute_potential_energy(positions, k, r0, box_size, epsilon_repulsive, epsilon_attractive, sigma):
    potential_energy = 0.0
    cutoff_rep = 2**(1/6) * sigma
    cutoff_att = 2.5 * sigma
    for i in range(len(positions) - 1):
        for j in range(i + 1, len(positions)):
            displacment = min_image(positions[j], positions[i], box_size)
            distance = np.linalg.norm(displacment) 
            if distance < 1e-12:
                continue
            if abs(i-j) == 1:
                harmonic_potential = 0.5 * k * (distance - r0)**2
                potential_energy += harmonic_potential
                continue
            elif abs(i-j) == 2 and distance < cutoff_rep:
                LJ_potential = 4 * epsilon_repulsive * ((sigma / distance)**12 - (sigma/distance)**6 + 0.25)
                potential_energy += LJ_potential
                continue
            if abs(i-j) > 2 and distance < cutoff_att:
                LJ_potential = 4 * epsilon_attractive * ((sigma / distance)**12 - (sigma/distance)**6 )
                potential_energy += LJ_potential
                continue
    return potential_energy


def calculate_radius_of_gyration(positions):
    center_of_mass = positions.mean(axis=0)
    rg_squared = np.mean(np.linalg.norm(positions - center_of_mass, axis=1)**2)
    rg = np.sqrt(rg_squared)
    return rg

def calculate_end_to_end_distance(positions):
    displacment = min_image(positions[-1], positions[0], box_size)
    Ree = np.linalg.norm(displacment)
    return Ree

def find_positions_at_temperatures(temperatures, n_particles, box_size, r0, k, sigma, epsilon_repulsive, epsilon_attractive, mass, dt, total_steps, rescale_interval): 
    positions_at_temperatures = {} 
    for T in temperatures:
        positions = initialize_chain(n_particles, box_size, r0)
        velocities = initialize_velocities(n_particles, T, mass) 
        total_forces = compute_forces(positions, k, r0, box_size, epsilon_repulsive, epsilon_attractive, sigma) 
        for step in range(total_steps):
            positions, velocities, total_forces = velocity_verlet(positions, velocities, total_forces, dt, mass, box_size) 
            if step % rescale_interval == 0: 
                velocities = rescale_velocities(velocities, T, mass)
        positions_at_temperatures[T] = positions
    return positions_at_temperatures


parmas = paramaters()
(n_particles, box_size, r0, k, sigma, epsilon_attractive, 
epsilon_repulsive, mass, target_temperature, dt, 
total_steps, rescale_interval, kb) = parmas

# Initialize positions and velocities
positions = initialize_chain(n_particles, box_size, r0)
velocities = initialize_velocities(n_particles, target_temperature, mass)



positions_data = []
velocities_data = []
total_forces = compute_forces(positions, k, r0, box_size, epsilon_repulsive, epsilon_attractive, sigma)
for step in range(total_steps):
    # Integrate equations of motion
    positions, velocities, total_forces = velocity_verlet(positions, velocities, total_forces, dt, mass, box_size)

    # Apply thermostat
    if step % rescale_interval == 0:
        velocities = rescale_velocities(velocities, target_temperature, mass)
    # (Optional) Store data for analysis
    positions_data.append(positions.copy())
    velocities_data.append(velocities.copy())

# Plot results
# plot positions
plt.figure()
ax = plt.subplot(111, projection='3d')
ax.plot(positions[:, 0], positions[:, 1], positions[:, 2])
plt.title('Positions of Particles')
plt.xlabel(r'X Position ($\sigma$)')
plt.ylabel(r'Y Position ($\sigma$)')
ax.set_zlabel(r'Z Position ($\sigma$)')
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\chained_positions_p30.png', dpi=300)
plt.show()

# plot velocities
plt.figure()
plt.plot(velocities_data[0][:, 0], label='X Velocity', marker='o', linestyle='-')
plt.plot(velocities_data[0][:, 1], label='Y Velocity', marker='o', linestyle='-')
plt.plot(velocities_data[0][:, 2], label='Z Velocity', marker='o', linestyle='-')
plt.title('Initial Velocities of Particles')
plt.xlabel('Particle Index')
plt.ylabel(r'Velocity ($\sigma/\tau$)')
plt.legend()
plt.grid(True)
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\initial_velocities_p30.png', dpi=300)
plt.show()

# example of analysis
# Arrays to store properties
temperatures = np.linspace(0.1, 1.0, 10)
Rg_values = []
Ree_values = []
potential_energies = []
equilibration = total_steps // 2


for T in temperatures:
    # Set target temperature
    target_temperature = T
    # (Re-initialize positions and velocities)
    positions = initialize_chain(n_particles, box_size, r0)
    velocities = initialize_velocities(n_particles, target_temperature, mass)
    total_forces = compute_forces(positions, k, r0, box_size, epsilon_repulsive, epsilon_attractive, sigma)
    U_array = []
    Rg_array = []
    Ree_array = []
    temp_history = []

    for step in range(total_steps):
        positions, velocities, total_forces = velocity_verlet(positions, velocities, total_forces, dt, mass, box_size)
        kinetic_energy = 0.5 * mass * np.sum(np.linalg.norm(velocities, axis=1)**2)
        T_inst = (2 / 3) * kinetic_energy / (n_particles * kb)
        temp_history.append(T_inst)
        if step % rescale_interval == 0:
            velocities = rescale_velocities(velocities, target_temperature, mass)
        if step >= equilibration and step % 10 == 0:
            U_array.append(compute_potential_energy(positions, k, r0, box_size, epsilon_repulsive, epsilon_attractive, sigma))
            Rg_array.append(calculate_radius_of_gyration(positions))
            Ree_array.append(calculate_end_to_end_distance(positions))


    # Compute properties
    Rg_values.append(np.mean(Rg_array))
    Ree_values.append(np.mean(Ree_array))
    potential_energies.append(np.mean(U_array))
# find transition points (based on derivative change)
# reterive data
Rg_values = np.array(Rg_values)
Ree_values = np.array(Ree_values)
U_values = np.array(potential_energies)

# get derivatives
dRg_dT = np.gradient(Rg_values, temperatures)
dRee_dT = np.gradient(Ree_values, temperatures)
dU_dT = np.gradient(U_values, temperatures)
# find transition points
transition_point_Rg = temperatures[np.argmax(np.abs(dRg_dT))]
transition_point_Ree = temperatures[np.argmax(np.abs(dRee_dT))]
transition_point_U = temperatures[np.argmax(np.abs(dU_dT))]
# print transition points
print(f"Transition point for End to end distance at T={transition_point_Ree:.2f}")
print(f"Transition point for U at T={transition_point_U:.2f}")
print(f"Transition point for Radius of gyration at T={transition_point_Rg:.2f}")

# plot transition points
plt.figure()
plt.plot(temperatures, dRee_dT, color = 'blue', label='dRee/dT', linestyle='--')
plt.xlabel(r'Temperature ($\epsilon/k_B$)')
plt.ylabel(r'Derivative ($\sigma / (\epsilon/k_B)$)')
plt.title('Derivative of End-to-End Distance vs Temperature')
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\end_to_end_distance_derivative_p30.png', dpi=300)
plt.show()

plt.figure()
plt.plot(temperatures, dU_dT, color = 'red', label='dRg/dT', linestyle='--')
plt.xlabel(r'Temperature ($\epsilon/k_B$)')
plt.ylabel(r'Derivative ($\epsilon / (\epsilon/k_B)$)')
plt.title('Derivative of Potential Energy vs Temperature')
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\potential_energy_derivative_p30.png', dpi=300)
plt.show()

# Plotting
plt.figure()
plt.plot(temperatures, Rg_values, label='Radius of Gyration', color = 'red')
plt.xlabel(r'Temperature ($\epsilon/k_B$)')
plt.ylabel(r'Radius of Gyration ($\sigma$)')
plt.title('Radius of Gyration vs Temperature')
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\radius_of_gyration_p30.png', dpi=300)
plt.show()

plt.figure()
plt.plot(temperatures, Ree_values, label='End-to-End Distance', color = 'blue')
plt.axvline(x=transition_point_Ree, color='green', linestyle='--', label='Transition Point')
plt.xlabel(r'Temperature ($\epsilon/k_B$)')
plt.ylabel(r'End-to-End Distance ($\sigma$)')
plt.title('End-to-End Distance vs Temperature')
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\end_to_end_distance_p30.png', dpi=300)
plt.show()

plt.figure()
plt.plot(temperatures, potential_energies, label='Potential Energy', color = 'green')
plt.xlabel(r'Temperature ($\epsilon/k_B$)')
plt.ylabel(r'Potential Energy ($\epsilon$)')
plt.title('Potential Energy vs Temperature')
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\potential_energy_p30.png', dpi=300)
plt.show()

plt.figure()
plt.plot(temp_history)
plt.axhline(target_temperature, color='r', linestyle='--')
plt.xlabel(r'Time Step ($\tau$)')
plt.ylabel(r'Temperature ($\epsilon/k_B$)')
plt.title('Temperature Stability')
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\temp_vs_target_temp_png', dpi=300)
plt.show()



# Example of finding positions at 3 different temperatures
temperatures = [0.1, 0.5, 1.0]
positions_at_temperatures = find_positions_at_temperatures(temperatures, n_particles, box_size, r0, k, sigma, epsilon_repulsive, epsilon_attractive, mass, dt, total_steps, rescale_interval)
# Plotting positions at different temperatures
for T, positions in positions_at_temperatures.items():
    state = 'folded' if T < transition_point_Ree else 'unfloded'
    plt.figure()
    ax = plt.subplot(111, projection='3d')
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], label=f'Temperature = {T}', color='blue', marker='o', linestyle='-')
    plt.title(f'{state} Final Positions and state of Particles at Temperature {T}')
    plt.xlabel(r'X Position ($\sigma$)')
    plt.ylabel(r'Y Position ($\sigma$)')
    ax.set_zlabel(r'Z Position ($\sigma$)')
    plt.legend()
    plt.savefig(f'C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\project_2\\positions_at_temperature_{T}_p30.png', dpi=300)
    plt.show()


print(f"temperature of transition point is {transition_point_Ree:.2f}")

