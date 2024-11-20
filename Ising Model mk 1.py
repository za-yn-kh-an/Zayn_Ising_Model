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
    
    
    lattice_points = np.zeros((L,L))
    up_or_down = [-1,1]

    
    for i in range(L):
        for j in range(L):
            lattice_points[i,j] = random.choice(up_or_down)
    
    """
    lattice_points = np.ones((L,L))
    """
    
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
    H : FLOAT.64
        The Hamiltonian of a given microstate of the lattice.
    w : FLOAT.64
        DESCRIPTION.
    """
    #k = Boltzmann
    k = 1
    delta_H = 0; J = 1
    
    
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
                
            delta_H += -J * reference_point * (above+below+right+left) #J is ang. mom. (high J means less magnetisation)
            w = np.exp(-delta_H/(k*T))
                    
    return delta_H, w



def metropolis(L,reps, T):
    """
    This function calculates the inital parameters of the lattice and then
    randomly changes one spin at a time. If the change is accepted then the
    magnetisation is calculated and applied to the results

    Parameters
    ----------
    L : INT
        DESCRIPTION.
    reps : INT
        DESCRIPTION.
    T : INT
        DESCRIPTION.

    Returns
    -------
    sum_A : FLOAT 64
        The sum of the magnetisation densities.
    z : INT
        DESCRIPTION.
    N : INT
        Number of iterations the model has gone through.
    """
    sum_A = 0
    z = 0
    N = 0
    H_list = []
    sum_w = 0
    
    lattice_points, length = initial_state(L)
    H,w = Hamiltonian(lattice_points, length, T)
    
    #print(H)
    #print(w)
    
    H_list.append(H)
    sum_w += w
    
    A = sum(sum(lattice_points))/(L**2)
    #print(A)
    sum_A += A; z +=1 ; N +=1
    
    for k in range(reps):
        i = random.randint(0, L) - 1; j = random.randint(0, L) - 1
        lattice_points[i,j] = -1 * lattice_points[i,j]
        H,w = Hamiltonian(lattice_points, L, T)
        new_A = sum(sum(lattice_points))/(L**2) 
        #print(w/sum_w)
        
        if H <= 0:
            sum_A += new_A ; z += 1 ; N +=1
            H_list.append(H)
            sum_w += w
        elif H > 0 and w >= 0.5:
            sum_A += new_A ; z += 1 ; N +=1
            H_list.append(H)
            sum_w += w
        else: lattice_points[i,j] = -1 * lattice_points[i,j]; N += 1 
    
    observed_A = sum_A/N   
    
    return observed_A, z, N

""" FIX ACCEPTANCE CONDITION """

T = np.arange(0.5,3,0.2)
T_len = len(T)


A = np.zeros(T_len)
std = np.zeros(T_len)
reps = 10
for i in range(T_len):
    B = np.zeros(reps)
    for j in range(reps):
        B[j] += (metropolis(20, 20000, T[i])[0])
    A[i] = np.mean(B)
    std[i] = np.std(B)



plt.plot(T,A)
plt.xlabel("Temperature")
plt.ylabel("Magnetisation density")
plt.errorbar(T, A, std)
plt.vlines(2.26, 0, 1, colors = "r", linestyles = "--", label = "Expected Critical Temperature")
plt.title("Temperature vs. Magnetisation Density")