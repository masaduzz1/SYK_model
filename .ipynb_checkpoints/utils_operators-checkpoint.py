#import usual libraries
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import sparse
from numpy import linalg as LA
from scipy import linalg as sLA
from scipy.sparse.linalg import eigs
from IPython.display import display, Latex
from decimal import *
import pandas as pd
import math
from qiskit.visualization import plot_histogram
from qiskit.visualization import *
from qiskit.circuit.library import iSwapGate
from scipy.linalg import expm
from scipy.stats import unitary_group
from IPython.core.interactiveshell import InteractiveShell
#from bristol.ensembles import Circular
from itertools import cycle
from functools import reduce
mycolor=cycle('krcbm')
mymarker=cycle('oxs^v')

# Paulis:  create Pauli to test commutations
def sigX():
    return np.array([[0,1],[1,0]])
def sigY(): 
    return np.array([[0,-1.j],[1.j,0]])
def sigZ():
    return np.array([[1,0],[0,-1]])
def sminus():
    return 0.5*(sigX()+1.0j*sigY() )
def splus():
    return 0.5*( sigX()-1.0j*sigY() ) 

# Pauli-in multiqubit system $\sigma^{\rm tag}_{q}$ 
#$\mathrm{tag}=x,y,z,+,-$
#$q=$ position of lattice
def sigma_creator(nq,q,tag):
    if q>=nq:
        print("error: q greater than nq")
        quit()
    iden=np.eye(2)
    
    if tag==0:
        Sa=sigX()
    elif tag==1:
        Sa=sigY()
    elif tag==2:
        Sa=sigZ()
    elif tag==3:
        Sa=splus()
    elif tag==4:
        Sa=sminus()
    
    if q==0:
        S=Sa
        for i in range(1,nq):
             S= np.kron(S,iden)
        return S
    else:
        S=iden
        for i in range(1,nq):
            if q!=i:
                S=np.kron(S,iden)
            else:
                S=np.kron(S,Sa)
        
        return S
    
    
simpl = 0 # No simplified SYK Hamiltonian 

def make_Majorana(N):

    # Make 'N' Majorana fermions. Set of N Hermitian matrices psi_i, i=1,..N
    # obeying anti-commutation relations {psi_i,psi_j} = δ_{ij}
    
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    I = np.array([[1, 0], [0, 1]])

    psi = dict()

    if simpl == 0: 
        for i in range(1, N+1):

            if (i % 2) == 1:
                matlist = [Z] * int((i-1)/2)
                matlist.append(X)
                matlist = matlist + [I] * int((N/2 - (i+1)/2))
                psi[i] = 1/np.sqrt(2)*reduce(np.kron, matlist)
            else:
                matlist = [Z] * int((i - 2) / 2)
                matlist.append(Y)
                matlist = matlist + [I] * int((N/2 - i/2))
                psi[i] = 1/np.sqrt(2)*reduce(np.kron, matlist)

    if simpl == 1: 
        for i in range(1, N+1):

            if (i % 2) == 1:
                matlist = [I] * int((i-1)/2)
                matlist.append(X)
                matlist = matlist + [I] * int((N/2 - (i+1)/2))
                psi[i] = 1/np.sqrt(2)*reduce(np.kron, matlist)
            else:
                matlist = [I] * int((i - 2) / 2)
                matlist.append(Y)
                matlist = matlist + [I] * int((N/2 - i/2))
                psi[i] = 1/np.sqrt(2)*reduce(np.kron, matlist)


    for i in range(1, N+1):
        for j in range(1, N+1):

            if simpl == 0:
            # Not checking this for simplified H yet.  

                if i != j:
                    if np.allclose(psi[i] @ psi[j], -psi[j] @ psi[i]) == False:
                        print ("Does not satisfy algebra")

                if i == j:
                    if np.allclose(psi[i] @ psi[j] + psi[j] @ psi[i], np.eye(int(2**(N/2)))) == False:
                        print ("Does not satisfy algebra for i=j")

    return psi


def SYK_Hamiltonian(psi, N, instances, J_squared):
    
    # Creates multiple realisations of the SYK Hamiltonian
    # Variance of couplings is given by 'J_squared * 3!/N^3'.

    H = 0
    J = dict()
    sigma_sq = 6.*J_squared/(N**3)
    sigma = math.sqrt(sigma_sq)
    np.random.seed(seed=213423)
    
    for i in range(1, N+1):
        for j in range(i+1, N+1):
            for k in range(j+1, N+1):
                for l in range(k+1, N+1):
                    J[i, j, k, l] = np.random.normal(loc=0, scale=sigma,size=instances)
                    M = psi[i] @ psi[j] @ psi[k] @ psi[l]
                    H = H + np.array([element * M for element in J[i, j, k, l]])

    return H

def single_SYK_Hamiltonian(psi, N, J_squared):
    
    # Creates multiple realisations of the SYK Hamiltonian
    # Variance of couplings is given by 'J_squared * 3!/N^3'.

    H = 0
    sigma_sq = 6.*J_squared/(N**3)
    sigma = math.sqrt(sigma_sq)
    
    

    # Generate all possible combinations of indices
    indices = np.array(np.meshgrid(range(1, N+1), range(1, N+1), range(1, N+1), range(1, N+1))).T.reshape(-1, 4)

    # Select unique combinations where i < j < k < l
    indices = indices[np.all(np.diff(indices, axis=1) > 0, axis=1)]

    # Generate random values for J outside the loop
    Js = np.random.normal(loc=0, scale=sigma, size=len(indices))

    # Compute M for each combination and add to H
    for idx, (i, j, k, l) in enumerate(indices):
        M = np.dot(np.dot(np.dot(psi[i], psi[j]), psi[k]), psi[l])
        H += Js[idx] * M


    return H

