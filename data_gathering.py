# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:27:59 2026

@author: eduvs
"""

# Hubbard tree automata, data gathering

import math
import numpy as np
from numpy import linalg
from fractions import Fraction
import matplotlib.pyplot as plt

import hubbard_tree_automata as hta

####

def size_of_mating_with(theta1, N):
    # Fixed an angle theta1
    # Prints the sizes of all the matings of theta1 with a/2**N

    denom = 2 ** N
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.mating_dyadics(theta1, theta2)
        size = A.size
        
        print(f'Mating {theta1} with {theta2}: size = {size}')
        
    return

####

def leading_eig(M):
    # Input is a Perron-Frobenius matrix (np array)
    # Output is the leading real eigenvalue
    
    eigenvalues = np.linalg.eigvals(M)
    
    maximum = 0
    
    for eigenvalue in eigenvalues:
        if np.isreal(eigenvalue):
            if eigenvalue.real > maximum:
                maximum = eigenvalue
    
    return maximum

####

def entropy(M):
    # Input is an adajcency matrix for a directed graph with loops
    # Output is the entropy
    
    # since the leading eigenvalue is already real, may take its real part
    # Otherwise python complains
    
    leading_eigenvalue = leading_eig(M).real
    
    return math.log(leading_eigenvalue)
    
####

def entropy_of_mating_with(theta1, N):
    # Fixed an angle theta1
    # Prints the entropies of all the matings of theta1 with a/2**N

    denom = 2 ** N
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.mating_dyadics(theta1, theta2)
        ent = round(entropy(A.matrix), 8)
        
        print(f'Mating {theta1} with {theta2}: entropy = {ent}')
        
    return

####

def plot_zero_entropy_angles(theta1, N):
    # Fixed an angle theta1
    # Outputs a plot of all dyadic angleswith denominator 2**N whose mating with theta1 has zero entropy
    
    denom = 2 ** N
    
    dots_x_red = []
    dots_y_red = []
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.mating_dyadics(theta1, theta2)
        h = entropy(A.matrix)
        
        if h == 0:
            dots_x_red.append(math.cos(theta2 * math.tau))
            dots_y_red.append(math.sin(theta2 * math.tau))
            
    plt.plot(dots_x_red, dots_y_red, ',')
    plt.gca().set_aspect("equal")
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.suptitle(f"Zero entropy on matings with {theta1}")
    plt.show()

    
    return

####

def ignore_trivial(A):
    # Input is an automaton A representing ray connections in a mating
    # Removes the self-transitions, corresponding to the trivial angles 1/2 and 0,
    # And then trims the automaton
    # Usefulness: If the only ray connections are the trivial ones and preimages, the new automaton is empty
    
    M = A.matrix
    states = A.states
    size = A.size
    
    new_M = M.copy()
    
    for j in range(0, size):
        if new_M[j][j] == 1:
            new_M[j][j] = 0
    
    new_A = hta.trim(hta.automaton(new_M, states))
    
    return new_A

####

def plot_no_non_trivial_ray_connections(theta1, N):
    
    denom = 2 ** N
    
    dots_x = []
    dots_y = []
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        B = ignore_trivial(hta.mating_dyadics(theta1, theta2))
        
        if B.size == 0:
            dots_x.append(math.cos(theta2 * math.tau))
            dots_y.append(math.sin(theta2 * math.tau))
            
    plt.plot(dots_x, dots_y, '.')
    plt.gca().set_aspect("equal")
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.suptitle(f"No non-trivial ray connections on mating with {theta1}")
    plt.show()

####

def main1():
    
    theta1 = Fraction(1,4)
    M = ''
    
    for n in range(6, 20):
        
        l = []
    
        for a in range(1, 2**n):
            
            theta2 = Fraction(a, 2**n)
            A = hta.mating_dyadics(theta1, theta2)
            
            if A.size == 4 and A.size != 2:
                l.append(a)
                
        m = max(l)
        
        M = M + ' ' + str(m) + '/' + str(2**n) + '\n'
    
        print(M)
    
    return

####

def main2():
    
    dots_x_red = []
    dots_y_red = []
    
    dots_x_blu = []
    dots_y_blu = []
    
    theta = Fraction(1,4)
    
    for a in range(1, 128):
        A = hta.mating_dyadics(theta, Fraction(a,128))
        M = A.matrix
        
        if A.size > 4 or A.size == 2:
            dots_x_blu.append(math.cos(Fraction(a, 128) * math.tau))
            dots_y_blu.append(math.sin(Fraction(a, 128) * math.tau))
        
        eigenvalues = np.linalg.eigvals(M)
        
        maximum = 0
        
        for eigenvalue in eigenvalues:
            if np.isreal(eigenvalue) and eigenvalue > maximum:
                maximum = eigenvalue
        
        entropy = round(math.log(maximum), 8)
        
        if entropy > 0:
            
            dots_x_red.append(math.cos(Fraction(a, 128) * math.tau))
            dots_y_red.append(math.sin(Fraction(a, 128) * math.tau))
    
    plt.plot(dots_x_red, dots_y_red, 'ro')
    plt.plot(dots_x_blu, dots_y_blu, 'g^')
        
    return

####

def main3():
    
    for n in range(6,20):
        
        theta = Fraction(5, 2**n)
        
        print(hta.forbidden_region_dyadic(theta))
        
    return
