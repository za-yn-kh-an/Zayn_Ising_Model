import numpy as np
import matplotlib.pyplot as plt
import random as random
from scipy.constants import Boltzmann

def initial_state(L):
    """
    This function creates the initial state of the lattice by creating an
    L x L lattice and randomly allocating each point +1 or -1 to correspond to spin up or spin down.

    Parameters
    ----------
    L : INT
        Defines the size of the lattice sides.

    Returns
    -------
    lattice_points : INT
        A numpy array where each point represents a different spin.
    L : INT
        Defines the size of the lattice sides.
    """
    
    """
    lattice_points = np.zeros((L,L))
    
    
    up_or_down = [-1,1]
    
    for i in range(L):
        for j in range(L):
            lattice_points[i,j] = random.choice(up_or_down)
    """
    
    lattice_points = np.ones((L,L))
    
    
    return lattice_points, L
    

def Hamiltonian(lattice_points, L, T):
    """
    Calculates the Hamiltonian of each microstate of the lattice.
    It also sets out periodic boundary conditions for the lattice edges.
    
    Parameters
    ----------
    lattice_points : Array
        An array that contains the spins at each point of the lattice.
    L : INT
        Defines the size of the lattice sides.
    T : FLOAT 64
        Defines the temperature

    Returns
    -------
    energy : FLOAT 64
        The Hamiltonian of a given microstate of the lattice.
    weight : FLOAT.64
        Statistical weighting of the lattice microstate.
    """

    energy = 0; J = 1
    
    
    for i in range(L):
        for j in range(L):
            reference_point = lattice_points[i,j]
            above = lattice_points[i,j-1]
            left = lattice_points[i-1,j]

            if i == L-1:
                right = lattice_points[0,j]
            else:
                right = lattice_points[i+1,j]
                
            if j == L-1:
                below = lattice_points[i,0]
            else: 
                below = lattice_points[i, j+1]
                
            energy += -J * reference_point * (above+below+right+left) #J is ang. mom. (high J means less magnetisation)
            energy = energy/2
                    
    return energy



def metropolis(L, flips, T, r = 0.5):
    """
    This function calculates the inital parameters of the lattice and then
    randomly changes one spin at a time. If the change is accepted then the
    magnetisation is calculated and applied to the results

    Parameters
    ----------
    L : INT
        Length of the lattice sides.
    reps : INT
        number of times the simulation is run.
    T : FLOAT64
        The temperature of the simulation.

    Returns
    -------
    sum_mag_density : FLOAT 64
        The sum of the magnetisation densities at a given temperature.
    z : INT
        The partition function of the system.
    N : INT
        Number of iterations the model has gone through.
    """
    #mag_densities = []
    #partition_function = 0
    N = 0
    #H_list = []
    #k = Boltzmann
    k = 1
    
    
    lattice_points, length = initial_state(L)
    initial_energy = Hamiltonian(lattice_points, length, T)
    initial_mag_density = np.abs(np.sum(lattice_points))/(L**2)
    plt.imshow( lattice_points , cmap = 'viridis', animated=False)
    plt.title(f"Initial state, Energy = {initial_energy}, Mag density = {initial_mag_density}")
    plt.show()

    
    #print(H)
    #print(w)
    #H_list.append(energy)
    #print(A)
    #mag_densities.append(initial_mag_density); #partition_function += 1 ; N += 1
    
    for l in range(flips):
        i = random.randint(0, L-1); j = random.randint(0, L-1)
        lattice_points[i,j] = -1 * lattice_points[i,j]
        new_energy = Hamiltonian(lattice_points, L, T)
        #new_mag_density = sum(sum(lattice_points))/(L**2) 
        #print(w/sum_w)
        
        energy_change = new_energy - initial_energy
        weight = np.exp(-energy_change/(k*T))
        r = random.random()
        
        if energy_change <= 0:
            #mag_densities.append(new_mag_density) ; #partition_function += 1 ; N +=1
            initial_energy = new_energy
            #H_list.append(energy)
        elif (energy_change > 0 and weight >= r):
            initial_energy = new_energy
        else: lattice_points[i,j] = -1 * lattice_points[i,j]; N += 1
    

        if l%(flips/10) == 0:
            plt.imshow( lattice_points , cmap = 'viridis', animated=False)
            plt.title(f"T = {T}, Flip = {l}/{flips}, E = {new_energy}")
            plt.show()

   
    
    """Check this"""
    #mean_mag_density = sum_mag_density/N
    #mean_mag = np.mean(mag_densities)
    final_mag_density = np.abs(np.sum(lattice_points))/(L**2)
    
    return final_mag_density, N, lattice_points#, partition_function

def data(T_initial, T_final, T_step, reps, L, flips, r = 0.5):
    """
    This function uses the metropolis functions to generate the magnetisation densities of the lattice at different temperatures, as well as the standard deviation for each.

    Parameters
    ----------
    T_initial : FLOAT64
        The initial temperature of the lattice.
    T_final : FLOAT64
        Final temperature of the lattice.
    T_step : FLOAT64
        The size of the increment as we increase temperature for each simulation.
    reps : INT
        Number of times the simulation is run at each temperature.
    L : INT
        Length of the lattice sides.
    flips : INT
        Number of spin flips per simulations.

    Returns
    -------
    T : NP ARRAY
        An array of all of the temperatures the sim has been run at.
    A : NP ARRAY
        An array of the magetisation density at each temperature.
    std : NP ARRAY
        The standard deviation of each mag density.

    """
    T = np.arange(T_initial, T_final , T_step)
    T_len = len(T)
    
    A = np.zeros(T_len)
    std = np.zeros(T_len)
    for i in range(T_len):
        B = np.zeros(reps)
        for j in range(reps):
            B[j] += metropolis(L, flips, T[i], r)[0]
            print(f"T = {T[i]}/{T[-1]}, Rep {j+1}/{reps}")
        A[i] = np.mean(B)
        print(A[i])
        std[i] = np.std(B)
    return T, A, std

def plot(T, A, std):
    plt.plot(T,A)
    plt.xlabel("Temperature")
    plt.ylabel("Magnetisation density")
    plt.errorbar(T, A, std)
    #plt.vlines(2.26, 0, 1, colors = "r", linestyles = "--", label = "Expected Critical Temperature")
    plt.title("Temperature vs. Magnetisation Density")

#data(0.5, 4.5, 0.5, 10, 20, 10000)