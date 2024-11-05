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
    H : FLOAT.64
        The Hamiltonian of a given microstate of the lattice.
    w : FLOAT.64
        DESCRIPTION.
    """
    #k = Boltzmann
    k = 1
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
                
            H += J * reference_point * (above+below+right+left) #J is ang. mom. (high J means less magnetisation)
            w = np.exp(-H/(k*T))
                    
    return H, w



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
    sum_A += A; z +=1 ; N +=1
    
    for k in range(reps):
        i = random.randint(0, L) - 1; j = random.randint(0, L) - 1
        lattice_points[i,j] = -1 * lattice_points[i,j]
        H,w = Hamiltonian(lattice_points, L, T)
        
        #print(w/sum_w)
        
        if w/sum_w >= 0.2:
            sum_A += A ; z += 1 ; N +=1
            H_list.append(H)
            sum_w += w
        
    return sum_A, z, N



T = np.arange(0.2,3.2,0.2)
T_len = len(T)
"""
A = np.zeros(15)
reps = 1000
for i in range(15):
    for j in range(reps):
        A[i] += (metropolis(4, 1000, T[i])[0])
    A[i] = A[i]/reps
#myInt = 100
#new_A = [x / myInt for x in A]
"""

#%%

C = np.zeros(T_len)
reps = 1000
for j in range(reps):
    B = np.zeros(T_len)
    for i in range(T_len):
        B[i] += (metropolis(4, 1000, T[i])[0])
    C = C + B
A = np.divide(C,reps)



plt.plot(T,A)
plt.xlabel("Temperature")
plt.ylabel("Magnetisation Density")