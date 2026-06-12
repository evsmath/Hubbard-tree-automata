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

def main():
    
    theta = Fraction(1,4)
    
    for a in range(1, 32):
        M = hta.mating_dyadics(theta, Fraction(a,32)).matrix
        
        eigenvalues = np.linalg.eigvals(M)
        
        maximum = 0
        
        for eigenvalue in eigenvalues:
            if np.isreal(eigenvalue) and eigenvalue > maximum:
                maximum = eigenvalue
        
        entropy = round(math.log(maximum), 8)
        
        print(f'Fraction: {a}/32, entropy = {entropy}')
        
    return

####

def main1():
    
    theta = Fraction(1,4)
    
    M = ''
    
    for n in range(6, 20):
        
        l = []
    
        for a in range(1, 2**n):
            A = hta.mating_dyadics(theta, Fraction(a,2**n))
            
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

####