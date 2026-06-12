# -*- coding: utf-8 -*-
"""
Created on Mon May 25 16:33:06 2026

@author: Eduardo
"""

import math
import numpy as np
from numpy import linalg
from fractions import Fraction
import matplotlib.pyplot as plt

####

###############
# LAMINATIONS #
###############

####

class arc:
    # Defines non-empty arcs of the circle, oriented counterclockwise
    
    def __init__(self, left_endpoint, right_endpoint, closed):
        # endpoints are in the form of Fractions, assumed to be in [0,1), 0 = 1
        # closed is a boolean; True for closed arcs, False for open
        
        self.left_endpoint = left_endpoint
        self.right_endpoint = right_endpoint
        
        if left_endpoint == right_endpoint:
            self.degenerate = True
            
        else:
            self.degenerate = False
        
        if left_endpoint <= right_endpoint:
            self.ordered = True
        else:
            self.ordered = False
            
        self.closed = closed
        
        # Arc is either open or closed, not mixed
        # Convention: open degenerate arcs are circle minus endpoint
        # closed degenerate arcs are just the endpoint
        
    def __repr__(self):
        left = str(self.left_endpoint.numerator) + '/' + str(self.left_endpoint.denominator)
        right = str(self.right_endpoint.numerator) + '/' + str(self.right_endpoint.denominator)
        
        if self.closed:
            return '[' + left + ', ' + right + ']'
        else:
            return '(' + left + ', ' + right + ')'
        
    def __str__(self):
        
        left = str(self.left_endpoint.numerator) + '/' + str(self.left_endpoint.denominator)
        right = str(self.right_endpoint.numerator) + '/' + str(self.right_endpoint.denominator)
        
        if self.closed:
            return '[' + left + ', ' + right + ']'
        else:
            return '(' + left + ', ' + right + ')'
            
####

def are_equal(arc1, arc2):
    # Decides if two arcs are equal or not
    # Must be both open or both closed to consider equal as well
    
    if arc1.left_endpoint == arc2.left_endpoint:
        if arc1.right_endpoint == arc2.right_endpoint:
            if arc1.closed == arc2.closed:
                return True
            
    return False

####

def image_arc(arc1):
    # Gives the image arc of an arc
    # Assumes that the arc covers less than half of the circle; true in the relevant applications below
    
    new_left_endpoint = (2 * arc1.left_endpoint) % 1
    new_right_endpoint = (2 * arc1.right_endpoint) % 1
    
    image = arc(new_left_endpoint, new_right_endpoint, arc1.closed)
    
    return image

####

def in_arc(angle, arc):
    # Angle is received in the form of Fraction in [0,1)
    
    left_endpoint, right_endpoint = arc.left_endpoint, arc.right_endpoint
    
    if arc.closed: # Arc is closed
        if arc.degenerate: # Applying convention for closed degenerate arcs
            if angle == left_endpoint:
                return True
            return False
        
        if arc.ordered:
            if left_endpoint <= angle and angle <= right_endpoint:
                return True
            return False
        
        else:
            
            if angle <= right_endpoint or angle >= left_endpoint:
                return True
            return False # Observe that this also contemplates the case [t, 0]
        
    else: # Arc is open
        if arc.degenerate: # Applying convention for open degenerate arcs
            if angle != left_endpoint:
                return True
            else:
                return False
        
        if arc.ordered:
            if left_endpoint < angle and angle < right_endpoint:
                return True
            else:
                return False
            
        else:
            if angle < right_endpoint or angle > left_endpoint:
                return True
            else:
                return False

####

def arc_inclusion(arc_1, arc_2):
    # Arguments are two arcs, decides if the first is a subset of the second
    A, B = arc_1.left_endpoint, arc_1.right_endpoint
    C, D = arc_2.left_endpoint, arc_2.right_endpoint
    
    # Subdivides into the possible cases of open, closed and degenerate
    
    if arc_2.closed: # arc_2 open
        if arc_1.closed: # arc_1 clsoed
            if in_arc(A, arc_2):
                arc_3 = arc(A,D, True)
                if in_arc(B, arc_3):
                    return True
                return False
            return False
        
        else: #arc_1 open
            if in_arc(A, arc_2) and not arc_1.degenerate:
                arc_3 = arc(A,D, True)
                if in_arc(B, arc_3):
                    return True
                return False
            return False
    
    else: # arc_2 open
        if arc_1.closed: # arc_1 closed
            if in_arc(A, arc_2):
                arc_3 = arc(A,D, False)
                if B == A or in_arc(B, arc_3):
                    return True
                return False
            return False
        
        else: # arc_1 open
            if A == C or in_arc(A, arc_2):
                arc_3 = arc(A, D, False)
                if B == D or in_arc(B, arc_3):
                    return True
                return False
            return False
        
####

class leaf:
    # Defines the leaf object class
    # Endpoints are Fractions
    
    def __init__(self, endpoint_A, endpoint_B):
        
        # orders the endpoints in [0,1). May omit later
        if endpoint_A > endpoint_B: 
            endpoint_A, endpoint_B = endpoint_B, endpoint_A
        
        self.endpoint_A = endpoint_A
        self.endpoint_B = endpoint_B
        
        # degenerate leaf: just the endpoint
        if endpoint_A == endpoint_B:
            self.degenerate = True
        else:
            self.degenerate = False
            
        self.open_arcs = [arc(endpoint_A, endpoint_B, False), arc(endpoint_B, endpoint_A, False)]
        self.closed_arcs = [arc(endpoint_A, endpoint_B, True), arc(endpoint_B, endpoint_A, True)]
        # arcs are open and ordered counterclockwise
        # if leaf is degenerate, open/closed arcs are the same. May change later
        
    def __repr__(self):
        
        endpoint_A_string = str(self.endpoint_A.numerator) + '/' + str(self.endpoint_A.denominator)
        endpoint_B_string = str(self.endpoint_B.numerator) + '/' + str(self.endpoint_B.denominator)
        
        return 'leaf(' + endpoint_A_string + ', ' + endpoint_B_string + ')'
        
    def __str__(self):
        
        endpoint_A_string = str(self.endpoint_A.numerator) + '/' + str(self.endpoint_A.denominator)
        endpoint_B_string = str(self.endpoint_B.numerator) + '/' + str(self.endpoint_B.denominator)
        
        return 'leaf(' + endpoint_A_string + ', ' + endpoint_B_string + ')'
    
####

# TO SELF: may need to define function of when leaves are equal, or touch at endpoints.

def leaf_in_arc(leaf_1, arc_1):
    # Decides if a leaf is fully contained in a single arc, including endpoints if arc is closed
    
    A, B = leaf_1.endpoint_A, leaf_1.endpoint_B
    if in_arc(A, arc_1) and in_arc(B, arc_1):
        return True
    else:
        return False

####

def major_leaf(theta):
    # Given a rational of type Fraction, returns the leaf connecting its two preimages
    
    endpoint_A = theta / 2 
    endpoint_B = ((1 + theta) / 2)
    
    return leaf(endpoint_A, endpoint_B)

####

####################################
# KNEADING SEQUENCES AND ADDRESSES #
####################################

####

def kneading_string(theta):
    # Given a dyadic rational of type Fraction, returns the kneading sequence
    # Parentheses indicate the period
    
    iterate = theta
    
    diameter = major_leaf(theta)
    
    open_arcs = diameter.open_arcs
    
    for arc in open_arcs:
        if in_arc(theta, arc):
            arc_1 = arc
            
    kneading = ''
    
    while iterate.denominator % 2 == 0:
        if in_arc(iterate, arc_1):
            kneading += '1'
        else:
            kneading += '0'
        iterate = (2 * iterate) % 1
        
    kneading += '('
    
    first_periodic = iterate
    
    if in_arc(iterate, arc_1):
        kneading += '1'
    else:
        kneading += '0'
        
    iterate = (2 * iterate) % 1
        
    while iterate != first_periodic:
        if in_arc(iterate, arc_1):
            kneading += '1'
        else:
            kneading += '0'
        iterate = (2 * iterate) % 1 
        
    kneading += ')'
    
    return kneading

####

def kneading_list(theta):
    # Given a dyadic rational of type Fraction, returns the kneading sequence
    # First entry in the list is preperiod
    # Second is period
    
    iterate = theta
    
    diameter = major_leaf(theta)
    
    open_arcs = diameter.open_arcs
    
    for arc in open_arcs:
        if in_arc(theta, arc):
            arc_1 = arc
            
    kneading_preperiod = ''
    
    while iterate.denominator % 2 == 0:
        if in_arc(iterate, arc_1):
            kneading_preperiod += '1'
        else:
            kneading_preperiod += '0'
        iterate = (2 * iterate) % 1
        
    kneading_period = ''
    
    first_periodic = iterate
    
    if in_arc(iterate, arc_1):
        kneading_period += '1'
    else:
        kneading_period += '0'
        
    iterate = (2 * iterate) % 1
        
    while iterate != first_periodic:
        if in_arc(iterate, arc_1):
            kneading_period += '1'
        else:
            kneading_period += '0'
        iterate = (2 * iterate) % 1 
        
    k_list = [kneading_preperiod, kneading_period]
    
    return k_list

####

def internal_address_dyadic(theta):
    # returns the internal address of a dyadic angle theta
    # After the '...' is the sequence of consecutive integers following the last one
    
    k_list = kneading_list(theta)
    
    preperiod = k_list[0]
    period = k_list[1]
    
    
    m = len(preperiod)
    
    kneading = preperiod + (m * period)
    
    S = 1
    
    int_address = [1]
    
    while S < m:
        
        j = S + 1
        difference_found = False
    
        while not difference_found:
            eta_lower = kneading[j - 1 - S]
            eta_upper = kneading[j - 1]
            
            if eta_lower == eta_upper:
                j += 1
                
            else:
                difference_found = True
                S = j
                int_address.append(S)
                
    int_address.append('...')
                
    return int_address

####

def denominators_angled_internal_address_dyadic(theta):
    # returns the denominators in the angled internal address of a dyadic theta
    # Everything after the last 2 is a continuing list of 2's
    
    k_list = kneading_list(theta)
    
    preperiod = k_list[0]
    period = k_list[1]
    
    m = len(preperiod)
    
    kneading = preperiod + (m * period)
    
    int_address = internal_address_dyadic(theta)
    
    M = len(int_address) - 2
    k = 0
    
    denominators = []
    
    while k < M:
        
        S_k = int_address[k]
        S_next = int_address[k+1]
        
        r = S_next % S_k
        
        if r == 0:
            
            q_k = S_next // S_k
            
        else:
            
            r_i = r
            
            j = r_i + 1
            difference_found = False
        
            while not difference_found:
                eta_lower = kneading[j - 1 -r_i]
                eta_upper = kneading[j - 1]
                
                if eta_lower == eta_upper:
                    j += 1
                    
                else:
                    difference_found = True
                    r_i = j
                    
            # This computes the next element in the orbit of r.
                
            S_k_in_orbit = False
            
            while (not S_k_in_orbit) and r_i <= S_k:
                
                if S_k == r_i:
                    S_k_in_orbit = True
                    
                else:
                    j = r_i + 1
                    difference_found = False
                
                    while not difference_found:
                        eta_lower = kneading[j - 1 - r_i]
                        eta_upper = kneading[j - 1]
                        
                        if eta_lower == eta_upper:
                            j += 1
                            
                        else:
                            difference_found = True
                            r_i = j
                            
                    # This iterates the orbit of r_i
                    # If r_i ever surpasses S_k, or is equal, stop
                    
            if S_k_in_orbit:
                q_k = ((S_next - r) // S_k) + 1
            else:
                q_k = ((S_next - r) // S_k) + 2
                
        denominators.append(q_k)
        k += 1
    
    denominators.append(2)
    denominators.append('...')
    
    return denominators
    
####

def numerators_angled_internal_address_dyadic(theta):
    # returns the numerators in the angled internal address of a dyadic angle theta
    # Everything after the last 1 is a continuing list of 1's
    
    int_address = internal_address_dyadic(theta)
    denominators = denominators_angled_internal_address_dyadic(theta)
    
    M = len(int_address) - 2
    
    numerators = []
    
    k = 0
    
    while k < M:
        S_k = int_address[k]
        multiplier = 2 ** (S_k)
        q_k = denominators[k]
        
        count = 0
        theta_now = theta
        
        for i in range(q_k - 1):
            
            if theta_now <= theta:
                 count += 1
                 
            theta_now = (multiplier * theta_now) % 1
        
        p_k = count
        numerators.append(p_k)
        k += 1
        
    numerators.append(1)
    numerators.append('...')
    
    return numerators

####

def dyadic_to_binary(fraction):
    # input is a dyadic fraction, returns its terminating binary expansion
    
    num = fraction.numerator
    denom = fraction.denominator
    
    m = 0
    d = denom
    
    while d != 1:
        d = d // 2 
        m += 1
        
    A = format(num, 'b')
    l = len(A)
    
    binary_string = '.' + ((m - l) * '0') + A

    return binary_string

    # Note that len() will returns m + 1, since there is a '.' at the start

####

#######################
# AUXILIARY FUNCTIONS #
#######################

####

def binary_period(alpha):
    # Input is a periodic angle theta in the form of a Fraction
    # Output is the string of bits forming its period for its binary expansion
    
    if alpha == 0:
        return '0'
    
    s = ''
    
    if alpha < Fraction(1,2):
        s += '0'
    else:
        s += '1'
    
    angle = (2 * alpha) % 1
    
    while angle != alpha:
        
        if angle < Fraction(1,2):
            s += '0'
        else:
            s += '1'
        
        angle = (2 * angle) % 1
        
    return s
    
####

def preperiodic_orbit(theta):
    # Given a preperiodic angle under doubling in the form of a fraction
    # Returns its orbit
    
    angle = theta
    
    orb = []
    
    while angle.denominator % 2 == 0:
        orb.append(angle)
        
        angle = (2 * angle) % 1
        
    first_periodic = angle
    
    orb.append(angle)
    angle = (2 * angle) % 1
    
    while angle != first_periodic:
        
        orb.append(angle)
        angle = (2 * angle) % 1
    
    return orb

def partner_angle(theta):
    # Given a periodic angle under doubling in the form of a fraction
    # Returns the partner angle
    
    if theta == 0:
        return Fraction(0,1)
    
    orb = preperiodic_orbit(theta)
    
    m = major_leaf(theta)
    closed_arcs = m.closed_arcs
    
    if m.endpoint_A == orb[-1]:
        p = m.endpoint_B
        
    else:
        p = m.endpoint_A
        
    k = -1
    
    current_theta_orb = orb[-1]
    current_pullback = p
    
    while current_theta_orb != theta:
        
        k -= 1
        
        current_theta_orb = orb[k]
        
        for arc in closed_arcs:
            if in_arc(current_theta_orb, arc):
                current_arc = arc
        
        preimages = [current_pullback / 2, (current_pullback + 1) / 2]
        
        for preimage in preimages:
            if in_arc(preimage, current_arc):
                current_pullback = preimage
    
    n = len(orb)
    
    denom = (2 ** n) - 1
    
    C = denom * current_pullback
    
    num = C.__round__()
    
    partner = Fraction(num,denom)
    
    return partner
        
####

def cardioid_angle(t):
    # Input is t = p/q rotation number in the form of a Fraction
    # Output is the angle theta_minus landing at the p/q-sublimb of the main cardioid
    # Recall that theta_plus = theta_minus + 1/(2^q -1)
    
    q = t.denominator
    
    b_0 = 2*q + 2
    
    partial_sum = 0
    
    power = 1
    
    for b in range(1, b_0 + 1):
        
        power = power * 2
        
        partial_sum += (math.ceil(b*t) - 1) / power
        
    q_power = (2 ** q) - 1
    
    a_minus = math.ceil(q_power * partial_sum)
        
    return Fraction(a_minus, q_power)

####

def cardioid_tuning(alpha, t):
    # Input is a periodic angle landing at the root of a hyperbolic component and t = p/q rotation number
    # Output is the angle theta_minus(t) landing at the root of the p/q-sublimb, from the left
    
    partner = partner_angle(alpha)
    
    if alpha < partner:
        alpha_minus = alpha
        alpha_plus = partner
    else:
        alpha_minus = partner
        alpha_plus = alpha
    
    u_minus = binary_period(alpha_minus)
    u_plus = binary_period(alpha_plus)
    
    card_minus = cardioid_angle(t)
    
    tuning_reference = binary_period(card_minus)
    
    s = ''
    j = 0
    
    for j in tuning_reference:
        if j == '0' :
            s += u_minus
        if j == '1' :
            s += u_plus
        
    num = int(s, 2)
    
    angle = Fraction(num, (2 ** len(s)) - 1)
    
    return angle

####

def periodic_branch_portraits_dyadic(theta):
    
    int_address = internal_address_dyadic(theta)
    denominators = denominators_angled_internal_address_dyadic(theta)
    numerators = numerators_angled_internal_address_dyadic(theta)
    
    M = len(int_address) - 2
    
    bin_string = dyadic_to_binary(theta)
    
    branch_satellite_orbits = []
    
    k = 0
    
    # First, case of k = 0:
        
    q_0 = denominators[0]
    
    if q_0 >= 3:
        
        p_0 = numerators[0]
        
        theta_minus = cardioid_angle(Fraction(p_0, q_0))
        
        orbit_portrait = [[]]
        orbit_angle = theta_minus
        
        for j in range(q_0):
            
            orbit_portrait[0].append(orbit_angle)
            orbit_angle = (2 * orbit_angle) % 1
    
        branch_satellite_orbits.append(orbit_portrait)
    
    k += 1
    
    while k < M:
        
        q_k = denominators[k]
               
        if q_k >= 3:
            
            S_k = int_address[k]
            p_k = numerators[k]
            
            approx_string = bin_string[1:S_k + 1]
            
            approx_num = int(approx_string, 2)
            
            angle = Fraction(approx_num, (2 ** S_k) - 1)
            
            tuned_angle = cardioid_tuning(angle, Fraction(p_k, q_k))
            
            orbit_angle = tuned_angle
            
            orbit_portrait = []
            
            for i in range(S_k):
                orbit_portrait.append([])
            
            for j in range(q_k):
                for i in range(S_k):
                    
                    orbit_portrait[i].append(orbit_angle)
                    orbit_angle = (2 * orbit_angle) % 1
                    
            branch_satellite_orbits.append(orbit_portrait)
        
        k += 1
            
    return branch_satellite_orbits

####

################################
# FINDING THE FORBIDDEN REGION #
################################

def preimage_angles_landing_at_preperiodic_point(angles):
    # Input is a collection of angles landing together at a periodic point
    # Output is the collection of angles landing together at its strictly preperiodic preimage
    
    preimage_angles = []
    
    for angle in angles:
        
        if angle.numerator % 2 == 0:
            preimage_angle = (1 + angle) / 2
        else:
            preimage_angle = angle / 2
            
        preimage_angles.append(preimage_angle)
    
    return preimage_angles

####

def next_in_cyclic_order(alpha, angles):
    # Given a list, set, or dictionarty of angles and an angle alpha within it
    # Returns the next angle in cyclic order
    
    # In practice: assume set of angles >= 3 in size, and 0 is not an element
    
    found_one_larger = False
    
    best_so_far = -1
    
    for angle in angles:
        if angle > alpha:
            if found_one_larger == False:
                found_one_larger = True
                best_so_far = angle
                
            else:
                if angle < best_so_far:
                    best_so_far = angle
    
    if found_one_larger == False:
        next_angle = min(angles)
        
    else:
        next_angle = best_so_far
    return next_angle
    
####

def contains(open_arc, P):
    # Given an open arc and a set P
    # Decides if that arc contains points of P or not
    
    for p in P:
        if in_arc(p, open_arc):
            return True
    
    return False

####

def is_possible_attaching_point(angles_landing_at_point, P):
    # Given a set of angles landing at a point and a set P of points
    # Decides if at least two of the corresponding arcs contain points of P
    # If only one arc did, then the point is not a candidate for being on the tree
    
    at_least_one_arc_with_points_of_P = False
    at_least_two_arcs_with_points_of_P = False
    
    for angle in angles_landing_at_point:
        
        next_angle = next_in_cyclic_order(angle, angles_landing_at_point)
        open_arc = arc(angle, next_angle, False)
        
        if contains(open_arc, P):
            if at_least_one_arc_with_points_of_P:
                at_least_two_arcs_with_points_of_P = True
                
            at_least_one_arc_with_points_of_P = True
            
    if at_least_two_arcs_with_points_of_P:
        return True
    
    else:
        return False
    
####
    
def is_branch_point(angles_landing_at_point, P):
    # Given a set of angles landing at a point and a set P of points
    # Decides if the point is a branch point of the corresponding tree or not
    
    at_least_one_arc_with_points_of_P = False
    at_least_two_arcs_with_points_of_P = False
    at_least_three_arcs_with_points_of_P = False
    
    for angle in angles_landing_at_point:
        
        next_angle = next_in_cyclic_order(angle, angles_landing_at_point)
        open_arc = arc(angle, next_angle, False)
        
        if contains(open_arc, P):
            if at_least_two_arcs_with_points_of_P:
                at_least_three_arcs_with_points_of_P = True
                
            if at_least_one_arc_with_points_of_P:
                at_least_two_arcs_with_points_of_P = True
                
            at_least_one_arc_with_points_of_P = True
            
    if at_least_three_arcs_with_points_of_P:
        return True
    
    else:
        return False

####

def preimages_with_labels(angles_with_labels, P, diameter):
    # Input is a dictionary of angles A_i with labels T(A_i). Then angles are assumed to be of rays landing
    # at a given (preperiodic) point. Also as an input is a set of points P and a diameter
    
    # T(A_i) is true if the arc (A_i, A_{i+1}) in the cyclic order is a preimage of a forbidden arc,
    # and so doesn't need to be counted
    
    # Output are two dictionaries for the two preimages of the point, with corresponding preimage angles
    # If T(A_i) == True, then T(A'_i) == T(A''_i) == True for the two preimage angles
    # If T(A_i) == False but (A_i, A_{i+1}) does not contain points of P, so that it is a forbidden arc,
    # then T(A'_i) == T(A''_i) == True
    # Else, T(A'_i) == T(A''_i) == False
    
    open_arcs = diameter.open_arcs
    
    reference_arc = open_arcs[0]
    preimage_0 = {}
    preimage_1 = {}
    
    for angle in angles_with_labels:
        
        preimage_angle_A = angle / 2
        preimage_angle_B = (1 + angle) / 2
        
        next_angle = next_in_cyclic_order(angle, angles_with_labels)
        open_arc = arc(angle, next_angle, False)
        containment = contains(open_arc, P)
        
        if in_arc(preimage_angle_A, reference_arc):
            
            if angles_with_labels[angle] == True:
                preimage_0[preimage_angle_A] = True
                preimage_1[preimage_angle_B] = True
                
            else:
                if containment:
                    preimage_0[preimage_angle_A] = False
                    preimage_1[preimage_angle_B] = False
                else:
                    preimage_0[preimage_angle_A] = True
                    preimage_1[preimage_angle_B] = True
            
            
        else:
            if angles_with_labels[angle] == True:
                preimage_1[preimage_angle_A] = True
                preimage_0[preimage_angle_B] = True
                
            else:
                if containment:
                    preimage_1[preimage_angle_A] = False
                    preimage_0[preimage_angle_B] = False
                else:
                    preimage_1[preimage_angle_A] = True
                    preimage_0[preimage_angle_B] = True
    
    return [preimage_0, preimage_1]

####

def forbidden_region_dyadic(theta):

    P = preperiodic_orbit(theta)
    diameter = major_leaf(theta)
    periodic_branch_portraits = periodic_branch_portraits_dyadic(theta)
    
    candidate_points = []
    
    for portrait in periodic_branch_portraits:
        for angles_landing_at_point in portrait:
            
            preperiodic_preimage = preimage_angles_landing_at_preperiodic_point(angles_landing_at_point)
            
            candidate_point = {}
            
            for angle in preperiodic_preimage:
                candidate_point[angle] = False 
                
                # This says that the corresponding arcs are not yet preimages of forbidden arcs
            
            candidate_points.append(candidate_point)
        
    # Now, we have created the first list of candidates for attaching points of the Hubbard tree,
    # And for the corresponding arcs that could be forbidden
    # Note that a forbidden region has not yet been identified, and so all arcs are labeled False
    # As they are not yet preimages of forbidden arcs

    forbidden_arcs = []
    no_more_attaching_points = False
    
    while not no_more_attaching_points:
        
        new_list_of_candidate_points = []
        
        for candidate_point in candidate_points:
                
            if is_possible_attaching_point(candidate_point, P):
                # the arcs must separate at least two points of P, otherwise the point is not a possible attaching
                # point, nor do we add its preimages
            
                for angle in candidate_point:
                    
                    next_angle = next_in_cyclic_order(angle, candidate_point)
                    open_arc = arc(angle, next_angle, False)
                    
                    if contains(open_arc, P) == False and (candidate_point[angle] == False):
                        forbidden_arcs.append(open_arc)
                        # Add the forbidden arcs, if it is not already a preimage of a forbidden arc
                        
                if is_branch_point(candidate_point, P): 
                    # if it is a branch point of H, then adds the two preimages with labels
                    # Recall that when taking the preimage points, the preimages of forbidden arcs get the label
                    # True, so that they are no counted as different forbidden arcs
                    
                    preimages = preimages_with_labels(candidate_point, P, diameter)
                    new_list_of_candidate_points += preimages
                    
        candidate_points = new_list_of_candidate_points
        
        if candidate_points == []:
            no_more_attaching_points = True
            
    # Obtained list of forbidden arcs
    # Now we sort the forbidden arcs in cyclic order within the unit circle

    arcs_w_left_endpoints = {}
    
    for forbidden_arc in forbidden_arcs:
        arcs_w_left_endpoints[forbidden_arc] = forbidden_arc.left_endpoint
            
    sorted_arcs_dict = {key: value for key, value in sorted(arcs_w_left_endpoints.items(), key = lambda item: item[1])}
    
    sorted_arcs = []
    for forbidden_arc in sorted_arcs_dict:
        sorted_arcs.append(forbidden_arc)
        
    return sorted_arcs

####

#######################
# AUTOMATA AND MATING #
#######################

####

class automaton:
    
    def __init__(self, M, states):
    # defines an automaton with transition graph given by a matrix M
    # states are numbered 0 to n-1, where M is n x n
    
    # second input is an n x 3 matrix: first column is the name of the state, second is the number in M,
    # third is the bit digit assigned

    # By convention, origins of a transition are columns, destinations are rows:
    # alings with matrix multiplication M * v on the left
        
        Matrix = np.array(M) # converts into array, if it isn't already
        
        self.matrix = Matrix
        self.size = len(M)
        self.states = states
        
    def __repr__(self):
        
        return str(self.matrix)
    
    def __str__(self):
        
        return str(self.matrix)
        
####

def hubbard_automaton(theta):
    # creates the automaton for theta based on the forbidden region
    # Need to keep track of isolated angles/consecutive forbidden arcs
    
    # Finding the isolated angles. Note that the last forbidden inte
    
    forbidden_arcs = forbidden_region_dyadic(theta)
    A = len(forbidden_arcs)
    
    # Note that the forbidden arcs are already ordered counterclockwise, starting from 0
    
    if A == 0:
        M = np.ones((2,2), dtype = int)
        
        states = [['A_0', 0, 0], ['A_1', 1, 1]]
        
        return automaton(M, states)
    
    # Janky solution: image of the arc [0,1/2] would mistakenly return the degenerate arc [0,0].
    # Only case where an arc between marked points would have length >= 1/2 is for theta = Fraction(1,2) anyway
    
    degenerate_arcs = []
    
    for i in range(A-1):
        if forbidden_arcs[i].right_endpoint == forbidden_arcs[i+1].left_endpoint:
            degenerate_arcs.append( arc(forbidden_arcs[i].right_endpoint, forbidden_arcs[i+1].left_endpoint, True) )
    
    # Finds the isolated angles / degenerate closed arcs within two consecutive forbidden regions (open arcs)
    
    marked_points = [Fraction(0,1), Fraction(1,2)]
    
    for forbidden_arc in forbidden_arcs:
        
        marked_points += preperiodic_orbit(forbidden_arc.left_endpoint)
        marked_points += preperiodic_orbit(forbidden_arc.right_endpoint)
        
    marked_points = list(set(marked_points))
    marked_points.sort()
    
    m = len(marked_points)
    
    # Finds the set of all marked points, to construct the Markov partition
    
    states = []
    dict_addresses_to_arcs = {}
    
    i = 0
    c = 0
    
    while c < m - 1:
        closed_arc = arc(marked_points[c], marked_points[c+1], True)
        
        not_a_forbidden_arc = True 
        
        for forbidden_arc in forbidden_arcs:
            
            if (closed_arc.left_endpoint == forbidden_arc.left_endpoint and closed_arc.right_endpoint == forbidden_arc.right_endpoint):
                
                not_a_forbidden_arc = False
        
        if not_a_forbidden_arc:
                
            name = 'A_' + str(i)
            address = i
                
            if closed_arc.right_endpoint <= Fraction(1,2):
                bit = 0
            else:
                bit = 1
                    
            states.append([name, address, bit])
            dict_addresses_to_arcs[i] = closed_arc
                
            i += 1
            
        c += 1
        
    
    # Perform for last arc [t,0]:
    
    closed_arc = arc(marked_points[m-1], marked_points[0], True)
    
    not_a_forbidden_arc = True
    
    for forbidden_arc in forbidden_arcs:
        
        if (closed_arc.left_endpoint == forbidden_arc.left_endpoint and closed_arc.right_endpoint == forbidden_arc.right_endpoint):
            not_a_forbidden_arc = False
            
    if not_a_forbidden_arc:
    
        name = 'A_' + str(i)
        address = i
        bit = 1 # The last arc in counterclockwise order is always labeled 1
                
        states.append([name, address, bit])
        dict_addresses_to_arcs[i] = closed_arc
            
    i += 1
    
    # Now we need to attribute names and addresses to the degenerate arcs
    
    P_degenerate_arcs = []
    
    for degenerate_arc in degenerate_arcs:
        
        angle = degenerate_arc.left_endpoint
        
        name = 'I_' + str(i)
        address = i
        
        if angle < Fraction(1,2):
            bit = 0
        else:
            bit = 1
            
            # Note that 1/2 and 0 cannot be isolated angles.
            
        states.append([name, address, bit])
        dict_addresses_to_arcs[i] = degenerate_arc
        
        i += 1
    
        orbit = preperiodic_orbit(angle)
        
        for angle_in_orbit in orbit[1:]:
            
            name = 'P_' + str(i)
            address = i
            
            if angle_in_orbit < Fraction(1,2):
                bit = 0
            else:
                bit = 1
            
            states.append([name, address, bit])
            dict_addresses_to_arcs[i] = arc(angle_in_orbit, angle_in_orbit, True)
            
            P_degenerate_arcs.append(angle_in_orbit)
            
            i += 1

            
    # Here, we have finished adding all the states with their names, addresses and bits
    # We also keep track of the isolated angles and their names
    
    # Now we must construct the transition matrix
    # For this, we see how each arc-state maps over other arc states, including degenerate arcs
    # Next, for each isolated angle, we just map it along its own separate orbit
    
    S = len(states)
    
    M = np.zeros((S, S), dtype = int)

    for column in dict_addresses_to_arcs:
        
        arc_column = dict_addresses_to_arcs[column]
        image_column = image_arc(arc_column)
        
        # Each arc will be of length less than 1/2 total length of the circle, so there is no problem in taking images
        # The case for theta = 1/2 was discussed earlier
        
        for row in dict_addresses_to_arcs:
            
            arc_row = dict_addresses_to_arcs[row]
            
            if arc_inclusion(arc_row, image_column):
                
                if arc_row.degenerate == False:
                    M[row][column] = 1
                
                else:
                    angle_row = arc_row.left_endpoint
                    
                    if (angle_row not in P_degenerate_arcs) or (image_column.degenerate == True):
                        M[row][column] = 1
                
                # For non-degenerate arcs, will always map (bijectively) over a union of arcs
                # If the destination arc is degenerate: should only be mapped to if it is an isolated angle
                # or the origin is degenerate itself
    
    # dict_addresses_to_arcs is also pretty useful!
                
    A = automaton(M, states)
                    
    return A

####

def flip(A):
    # Input is an automaton with names, output is the automaton with same names and bit assignments flipped
    
    M = A.matrix
    states = A.states
    
    new_M = M.copy()
    flipped_states = []
    
    #Here we copy both the transition matrix (which in numpy copies also the lists inside, I assume)
    #And the list of states, to flip the bits without chaning the original ones
    
    for state in states:
        new_name = '~'+ state[0]
        new_address = state[1]
        flipped_bit = 1 - state[2]
        
        flipped_states.append([new_name, new_address, flipped_bit])
    
    return automaton(new_M, flipped_states)

####

def intersect_automata(A, B):
    # Input is two automata with names A and B
    # Input is the automaton formed by all pairs of states with matching bits,
    # and all allowable transitions
    
    M = A.matrix
    N = B.matrix
    
    states_A = A.states
    states_B = B.states
    
    paired_states = []
    paired_addresses = []
    
    i = 0
    for state_A in states_A:
        for state_B in states_B:
            if state_A[2] == state_B[2]:
                state_A_name, state_B_name = state_A[0], state_B[0]
                new_state_bit = state_A[2] # because the bits are the same
                
                new_state_name = '(' + state_A_name + ', ' + state_B_name + ')'
                
                new_state = [new_state_name, i, new_state_bit]
                paired_states.append(new_state)
                
                paired_addresses.append([state_A[1], state_B[1]])
                
                i += 1
         
    # We use the list of addresses (i,j) that each allowable state corresponds to to construct the new
    # transition matrix
    # The "new" addresses of the paired states are ordered starting from 0.
    # By the construction below, they will correspond to the appropriate rows and columns of the transition matrix
    
    P = []
    
    for pair_dest_address in paired_addresses:
        row = []
        for pair_origin_address in paired_addresses:
            origin_address_A = pair_origin_address[0]
            dest_address_A = pair_dest_address[0]
            
            origin_address_B = pair_origin_address[1]
            dest_address_B = pair_dest_address[1]
            
            if M[dest_address_A][origin_address_A] == 1 and N[dest_address_B][origin_address_B] == 1:
                row.append(1)
            else:
                row.append(0)
                
        P.append(row)
        
    P = np.array(P)

    return automaton(P, paired_states)

####

def trim(A):
    # Input is an automaton representing a transition graph
    # Output is the automaton with the states which eventually have no allowed transitions removed
    # In other words: remove the basis vectors representing states that are in the generalized
    # eigenspace of 0
    
    M = A.matrix
    states = A.states
    n = A.size
    
    trimmed_M = M.copy()
    
    # no need for anything more fancy to copy the numpy matrix of states M?
    
    trimmed_states = []
    
    for state in states:
        trimmed_states.append(state.copy())
        
    # This now genuinely creates a copy of the list of states (each of which is a copy)
    
    guaranteed_fully_trimmed = False
    
    while guaranteed_fully_trimmed == False:
    
        found_zero_column = False
        
        j = 0
        
        while found_zero_column == False and j < n:
            
            found_transition_from_origin = False
            
            i = 0
            
            while found_transition_from_origin == False and i < n:
                
                origin_state = trimmed_states[j]
                origin_address = origin_state[1]
                
                destination_state = trimmed_states[i]
                destination_address = destination_state[1]
                
                # the i-th state is the destination
                # the j-th state is the origin
                # the entry [1] of the state is its address in the transition matrix M
                
                if trimmed_M[destination_address][origin_address] != 0:
                    
                    found_transition_from_origin = True
                else:
                    i += 1
                    
            if found_transition_from_origin == False:
                    
                found_zero_column = True
                
                # found that j is the index of the origin state, which represents a column of zeros
                
                zero_column_address = origin_address
                    
            else:
                j += 1
                    
        if found_zero_column == True:
            
            trimmed_states.pop(j)
            
            for state in trimmed_states:
                if state[1] > zero_column_address:
                    state[1] -= 1
            
            # Here is where we assume M is an numpy matrix (or array).
            
            trimmed_M = np.delete(trimmed_M, (zero_column_address), axis = 1)
            trimmed_M = np.delete(trimmed_M, (zero_column_address), axis = 0)
            
            n -= 1
            
        else:
            guaranteed_fully_trimmed = True
                
    return automaton(trimmed_M, trimmed_states)


####

def mating_automata(A,B):
    # Input is two automata representing rays landing on the Hubbard trees
    # Output is the automata of ray connections between the Hubbard trees in the mating
    
    B_flipped = flip(B)
    
    mating = intersect_automata(A, B_flipped)
    
    trimmed_mating = trim(mating)
    
    return trimmed_mating

####

def mating_dyadics(theta1, theta2):
    # Given two dyadic angles, whose corresponding parameter rays land at parameters of the Mandelbrot set
    # Returns the automaton corresponding to ray connections between the Hubbard trees in the mating
    # of the corresponding polynomials
    
    A1 = hubbard_automaton(theta1)
    A2 = hubbard_automaton(theta2)
    
    A = mating_automata(A1, A2)
    
    return A

####

# new idea: trim mating matrix to get rid of self.transitions and states that aren't transitioned to
# seeks to find only the periodic cycles, check how many, how many components, etc

####

##################
# GATHERING DATA #
##################

####

def main():
    
    theta = Fraction(1,4)
    
    for a in range(1, 32):
        M = mating_dyadics(theta, Fraction(a,32)).matrix
        
        eigenvalues = np.linalg.eigvals(M)
        
        maximum = 0
        
        for eigenvalue in eigenvalues:
            if np.isreal(eigenvalue) and eigenvalue > maximum:
                maximum = eigenvalue
        
        entropy = round(math.log(maximum), 8)
        
        print(f'Fraction: {a}/32, entropy = {entropy}')
        
    return

def main1():
    
    theta = Fraction(1,4)
    
    M = ''
    
    for n in range(6, 20):
        
        l = []
    
        for a in range(1, 2**n):
            A = mating_dyadics(theta, Fraction(a,2**n))
            
            if A.size == 4 and A.size != 2:
                l.append(a)
                
        m = max(l)
        
        M = M + ' ' + str(m) + '/' + str(2**n) + '\n'
    
        print(M)
    
    return

def main2():
    
    dots_x_red = []
    dots_y_red = []
    
    dots_x_blu = []
    dots_y_blu = []
    
    theta = Fraction(1,4)
    
    for a in range(1, 128):
        A = mating_dyadics(theta, Fraction(a,128))
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

def main3():
    
    for n in range(6,20):
        
        theta = Fraction(5, 2**n)
        
        print(forbidden_region_dyadic(theta))
        
    return


