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

import hubbard tree automata as hta

def main3():
    
    for n in range(6,20):
        
        theta = Fraction(5, 2**n)
        
        print(hta.forbidden_region_dyadic(theta))
        
    return