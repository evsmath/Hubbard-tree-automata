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

import scipy
from scipy import sparse
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

import hubbard_tree_automata as hta

####

# Recall: we have constructed hubbard_tree_automaton(theta)
# and two the two mating constructions: product_mating and union_mating
# union_mating seems to be faster

#######################
# AUXILIARY FUNCTIONS #
#######################

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

def is_weakly_connected(M):
    
    n_components, labels = connected_components(csgraph = M, directed = True, connection = 'weak')
    
    if n_components == 1:
        return True
    
    else:
        return False
    
####

##################
# DATA GATHERING #
##################

####

def size_of_mating_with_dyadics(theta, N):
    # Fixed an angle theta
    # Prints the sizes of all the mating automata of theta with a/2**N

    denom = 2 ** N
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.union_mating(theta, theta2)
        size = A.size
        
        print(f'Mating {theta} with {theta2}: size = {size}')
        
    return

####

def entropy_of_mating_with(theta, N):
    # Fixed an angle theta
    # Prints the entropies of all the matings of theta with a/2**N

    denom = 2 ** N
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.mating_dyadics(theta, theta2)
        ent = round(entropy(A.matrix), 8)
        
        print(f'Mating {theta} with {theta2}: entropy = {ent}')
        
    return

####

def plot_zero_entropy_angles(theta1, N):
    # Fixed an angle theta1
    # Outputs a plot of all dyadic angle swith denominator 2**N whose mating with theta1 has zero entropy
    
    denom = 2 ** N
    
    plt.grid()
    
    dots_x = []
    dots_y = []
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.mating_dyadics(theta1, theta2)
        h = entropy(A.matrix)
        
        if h == 0:
            dots_x.append(math.cos(theta2 * math.tau))
            dots_y.append(math.sin(theta2 * math.tau))
            
    plt.plot(dots_x, dots_y, '.')
    plt.gca().set_aspect("equal")
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.suptitle(f"Zero entropy on matings with {theta1}")
    plt.show()
    
    return

####

def plot_zero_size(theta1, N):
    # Fixed an angle theta1
    # Outputs a plot of all dyadic angles swith denominator 2**N whose mating with theta1 has size non-zero
    
    denom = 2 ** N
    
    plt.grid()
    
    dots_x = []
    dots_y = []
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.union_mating(theta1, theta2)
        
        if A.size != 0:
            dots_x.append(math.cos(theta2 * math.tau))
            dots_y.append(math.sin(theta2 * math.tau))
            
    plt.plot(dots_x, dots_y, '.')
    plt.gca().set_aspect("equal")
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.suptitle(f"Zero size on matings with {theta1}")
    plt.show()
    
    return

####

def plot_zero_size_and_neg_F(theta1, N):
    # Fixed an angle theta1
    # Outputs a plot of all dyadic angles swith denominator 2**N whose mating with theta1 has size non-zero
    
    denom = 2 ** N
    F = hta.forbidden_region(theta1)
    neg_F = hta.negated_collection(F)
    
    plt.grid()
    
    dots_x_zero_size = []
    dots_y_zero_size = []
    
    dots_x_in_neg_F = []
    dots_y_in_neg_F = []
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.union_mating(theta1, theta2)
        
        if A.size != 0:
            dots_x_zero_size.append(math.cos(theta2 * math.tau))
            dots_y_zero_size.append(math.sin(theta2 * math.tau))
            plt.plot(math.cos(theta2 * math.tau), math.sin(theta2 * math.tau), '.')
            
        in_forbidden_arc = False
        
        f = len(F)
        i = 0
        while i < f and in_forbidden_arc == False:
            if hta.in_arc(theta2, neg_F[i]):
                dots_x_in_neg_F.append(math.cos(theta2 * math.tau))
                dots_y_in_neg_F.append(math.sin(theta2 * math.tau))
                in_forbidden_arc = True
            i += 1
              
    plt.plot(dots_x_zero_size, dots_y_zero_size, '.', color = 'blue')
    plt.plot(dots_x_in_neg_F, dots_y_in_neg_F, '.', color = 'yellow')        
    plt.gca().set_aspect("equal")
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.suptitle(f"Blue: non-zero size on matings with {theta1} \n Yellow: negated forbidden region")
    plt.show()
    
    return

####

def plot_zero_entropy_angles_height(theta1, N):
    # Fixed an angle theta1
    # Outputs a plot of all dyadic angle swith denominator 2**N whose mating with theta1 has zero entropy
    # x coordinate is the angle in [0,1], y coordinate is the height of the denominator
    
    denom = 2 ** N
    
    angles = []
    heights = []
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        A = hta.mating_dyadics(theta1, theta2)
        h = entropy(A.matrix)
        
        if h == 0:
            angles.append(theta2)
            heights.append(1/math.log2(theta2.denominator))
            
    plt.plot(angles, heights, '.')
    # plt.gca().set_aspect("equal")
    
    plt.xlim(0, 1)
    
    plt.xlabel('Angle theta = a/2^m')
    plt.ylabel('Height 1/m')
    
    plt.suptitle(f"Dyadic matings with {theta1} having zero entropy")
    plt.show()
    
    return

####

def plot_only_dyadic_ray_connections(theta1, N):
    
    denom = 2 ** N
    
    dots_x = []
    dots_y = []
    
    for a in range(1, denom):
        
        theta2 = Fraction(a, denom)
        B = ignore_trivial(hta.mating_dyadics(theta1, theta2))
        
        if B.size == 0:
            dots_x.append(math.cos(theta2 * math.tau))
            dots_y.append(math.sin(theta2 * math.tau))
            
    plt.plot(dots_x, dots_y, ',')
    plt.gca().set_aspect("equal")
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.suptitle(f"No non-dyadic ray connections on mating with {theta1}")
    plt.show()

####

def info(theta1, N):
    # Given N and theta
    # Prints all info on the mating with theta and a/2**N for all a
    
    denom = 2**N
    
    for a in range(1, denom):
        theta2 = Fraction(a, denom)
        A = hta.mating(theta1, theta2)
        M = A.matrix
        n_components, labels = connected_components(csgraph = M, directed = True, connection = 'weak')
        ent = entropy(M)
        if 5/7 < theta2 and theta2 < 6/7:
            conjugate = True
        else:
            conjugate = False
            
        B = ignore_trivial(A)

        
        print(f'Mating of {theta1} with {theta2}:\n - entropy = {ent}\n - size = {A.size}\n - size ignoring trivial = {B.size}\n - weak components = {n_components}\n conjugate limbs = {conjugate}\n')
    
    return


def info_self_matings(N):
    # Given N
    # Prints all info the mating with a/2**N wwith itself for all a
    
    denom = 2**N
    
    for a in range(1, denom):
        theta = Fraction(a, denom)
        A = hta.mating_dyadics(theta, theta)
        M = A.matrix
        n_components, labels = connected_components(csgraph = M, directed = True, connection = 'weak')
        ent = entropy(M)
            
        B = ignore_trivial(A)

        
        print(f'Mating of {theta} with {theta}:\n - entropy = {ent}\n - size = {A.size}\n - size ignoring trivial = {B.size}\n - weak components = {n_components}\n')
    
    return

####
    
def not_weakly_connected_up_to(N):
    
    denom = 2 ** N
    
    for a in range(1, denom):
        for b in range(a, denom):
            
            if a + b != denom: # Weirdness if conjugate: more than 2 components possibly
            
                theta1 = Fraction(a, denom)
                theta2 = Fraction(b, denom)
                
                A = hta.mating_dyadics(theta1, theta2)
                M = A.matrix
                
                n_components, labels = connected_components(csgraph = M, directed = True, connection = 'weak')
                
                if not is_weakly_connected(M):
                    
                    print(f'Mating {theta1} with {theta2} not weakly connected, size = {A.size}, number = {n_components}')
                
    return

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

def no_ray_connections_in_mating_with_dyadic(theta, N):
    
    power = 2**N
    
    for a in range(1, power):
        
        A = hta.union_mating(theta, Fraction(a, power))
        
        if A.size == 0:
            print(Fraction(a, power))
            
    return
    
####

def main4(N):
    for n in range(3, N+1):
        dyadic = Fraction(1, 2**(n-1))
        non = Fraction(1, 3 * (2**(n-2)))
        
        A = hta.union_mating(dyadic, non)
        
        print(A.size)
    return

####

def size_of_self_mating_up_to(N):
    
    for denom in range(2, N+1):
        for a in range(1, denom):
            
            if math.gcd(a, denom) == 1:
                
                theta = Fraction(a,denom)
                
                A = hta.union_mating(theta, theta)
                s = A.size
                
                print(f'Self-mating of {theta}: size = {s}')

    return

####

def main6(theta, N):
    
    for n in range(3, N+1):
        denom = 2 ** n
        
        first_non_zero_size_found = False
        a = 1
        
        while a < denom and first_non_zero_size_found == False:
            
            theta2 = Fraction(a, denom)
            A = hta.union_mating(theta, theta2)
            size = A.size
            
            b = hta.dyadic_to_binary(theta2)
            
            if size != 0:
                print(f'first non-zero size in mating with {theta}: is {b}')
                first_non_zero_size_found = True
                
            a += 1
        
    return

def main7(theta, N):
    # Conjecture: mating of theta with phi, where phi is in -F, where F is the forbidden region for theta
    # has empty mating automaton (except for specific tip cases)
    
    F = hta.forbidden_region(theta)
    neg_F = hta.negated_collection(F)
    f = len(F)
    
    for denom in range(2, N):
        if denom % 2 == 1:
            for a in range(1, denom):
                if math.gcd(a, denom) == 1:
                    phi = Fraction(a, denom)
                    
                    phi_in_neg_F = False
                    
                    i = 0
                    while i < f and phi_in_neg_F == False:
                        arc = neg_F[i]
                        if hta.in_arc(phi, arc):        
                            phi_in_neg_F = True
                        else:
                            i += 1
                            
                    # guaranteed phi is reduced form with odd denominator
                    if phi_in_neg_F == True:
                        A = hta.union_mating(theta, phi)
                        size = A.size
                        
                        print(f'size of mating of {theta} with {phi}: is {size}')

    return

# Conjecture has been verified for theta = 1/6 up to denominator 80, false for 1/4
# true for airplane






