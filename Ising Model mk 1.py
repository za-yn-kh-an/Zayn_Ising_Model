import numpy as np
import matplotlib.pyplot as plt
import random as random


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
    return lattice_points, L

def Hamiltonian(lattice_points, L):
    """
    Calculates the Hamiltonian of each microstate of the lattice.
    It also sets out periodic boundary conditions for the lattice edges.
    
    Parameters
    ----------
    lattice_points : Array
        An array that contains the spins at each point of the lattice.
    L : INT
        Defines the size of the lattice sides.

    Returns
    -------
    H : FLOAT.64
        The Hamiltonian of a given microstate of the lattice.
    w : FLOAT.64
        DESCRIPTION.

    """
    H = 0; J = 1
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
            
            H += J * reference_point * (above+below+right+left)
    w = np.exp(-H/8)
    return H, w

    

def metropolis(L,reps):
    sum_A = 0
    z = 0
    N = 0
    
    lattice_points, length = initial_state(L)
    H,w = Hamiltonian(lattice_points, length)
    
    z +=1 ; N +=1
    A = w/z; sum_A += A
    
    for k in range(reps):
        i = random.randint(0, L) - 1; j = random.randint(0, L) - 1
        lattice_points[i,j] = -1 * lattice_points[i,j]
        H,w = Hamiltonian(lattice_points, L)
        p = w/z
        if p >= 0.5:
            sum_A += p ; z += 1 ; N +=1
    return sum_A, z, N