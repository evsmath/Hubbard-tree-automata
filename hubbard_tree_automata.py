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

def fractions_equal(a , b):
    # Given two Fractions, checks if they are equal
    
    if a.numerator == b.numerator and a.denominator == b.denominator:
        return True
    
    else:
        return False

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
    # Assumes that the arc covers less than half of the circle
    
    # This will always be true in the applications below, but be careful
    
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

def image_leaf(l):
    # Given a leaf, outputs its image
    
    endpoint_A, endpoint_B = l.endpoint_A, l.endpoint_B
    
    image_A = (2 * endpoint_A) % 1
    image_B = (2 * endpoint_B) % 1
    
    return leaf(image_A, image_B)
    
####

###################################################
# KNEADING SEQUENCES, ADDRESSES AND COMBINATORICS #
###################################################

####

def principal_period(s):
    # Given a string s
    # Returns the principal period of s (Stack Exchange)
    
    i = (s+s).find(s, 1, -1)
    
    if i == -1:
        return s
    
    else:
        return s[:i]

####

def period(theta):
    # Given theta periodic under doubling
    # Computes its period
    
    iterate = (2 * theta) % 1
    n = 1
    
    while not fractions_equal(iterate, theta):
        iterate = (2 * iterate) % 1
        n += 1
        
    return n

####

def kneading_list(theta):
    # Given a dyadic rational of type Fraction, returns the kneading sequence
    # First entry in the list is preperiod
    # Second is period
    
    # Note: for preperiodic theta, the kneading_period of theta may be a strict divisor of the orbit period
    
    if theta == 0: # Convention
        return ['', '1']
    
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
        
        if (2 * iterate) % 1 == theta:
            kneading_period += '*'
        elif in_arc(iterate, arc_1):
            kneading_period += '1'
        else:
            kneading_period += '0'
        iterate = (2 * iterate) % 1 
        
    # Computes principal period for preperiodic:
        
    if theta.denominator % 2 == 0:
        
        kneading_period = principal_period(kneading_period)
        
    k_list = [kneading_preperiod, kneading_period]
    
    return k_list

####

def partner_angle(theta):
    # Given a periodic angle under doubling in the form of a fraction
    # Returns the partner angle
    
    # See MathOverflow post that proves this algorithm works
    # Idea: pullback the non-periodic preimage of theta along the itinerary of theta
    # The partner angle will be sufficiently close to the n-th pullback, if n is the period
    
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
    
    # Proof of this formula by Douady + tail estimates on the infinite series
    
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

def upper_kneading(theta):
    # Given an angle theta periodic under doubling
    # Computes the upper_kneading sequence of theta
    
    # The trick: find the angle of tuning with the basilica
    # Even though this lands at a parabolic point on the hyperbolic component that theta lands at the root of,
    # The kneading sequence will coincide with the upper kneading sequence of theta up to the n-th entry
    
    theta_minus_one_half = cardioid_tuning(theta, Fraction(1,2))
    
    n = period(theta)

    kneading_theta_minus_one_half = kneading_list(theta_minus_one_half)[1]

    return kneading_theta_minus_one_half[: n ]

####        

##########
# DYADIC #
##########

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

##############
# NON-DYADIC #
##############

####

def internal_address(theta):
    # Input is a rational angle theta, in the form of a Fraction
    # If theta is periodic under doubling, outputs the (terminating) internal address
    # If theta is preperiodic, outputs the internal address up to some S_k,
    # where we may guarantee that for all k' >= k, the denominator q_{k'} = 2
    
    # Biggest difficulty: finding a sufficient 'stop' condition that guarantees we won't miss denominators q >= 3
    
    k_list = kneading_list(theta)
    
    preperiod = k_list[0]
    period = k_list[1]
    
    m = len(preperiod)
    l = len(period)
    
    if m == 0: # periodic kneading sequence: know that it terminates at S_k = l
    
        kneading = period
        S = 1
        address = [1]
        
        while S < l:
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
                    address.append(S)
        return address
        
    else: # strictly preperiodic: infinite kneading sequence
    # Lemma: if S_k < m, then S_{k+1} < m + m * l
    
        kneading = preperiod + (m * period)
    
        S = 1
        address = [1]
        
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
                    address.append(S)
            
        # Guarantees that S >= m
        # Lemma: if S_k >= m, then S_{k+1} <= S_k + m + l
        
        while S < m + l:
            
            while len(kneading) < S + m + l:
                kneading += period
                
            # Guarantees that len(kneading) >= S + m + l >= S_{k+1}; difference will be found
            
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
                    address.append(S)
                    
        # Guarantees that S >= m + l
        # Hence, combined with the previous lemma,
        # in the formula for the denominators q_k, r = S_{k+1} - S_k, and (S_{k+1} - r)/S_k = 1
        
        # If nu^k is the approximating periodic block of length S_k, it will end in B_k, an initial segment of B
        # Lemma: B_k determines S_{k+1} - S_k = r = r_k, and therefore the next initial segment B_{k+1}
        # By finiteness, the collection of initial segments eventually repeats
        
        # Let {C_1, ..., C_t} be this collection, with corresponding {r_1, ..., r_t}
        # If S_k is in the rho-orbit of r_j, mark that C_j as "good" 
        # Any other S_{k'}, k'> k, corresponding to this block C_j will have denominator q = 2
        # This is because the rho-orbits of r_j and S_k merge
        
        # Once all the C_j are marked good, can stop computing the internal address
        
        # Remarks: to keep track of the initial segment B_k, only need to keep track of its length (S_k - m) % l
        # We also keep track the corresponding r = r_k, which we find by computing S_{k+1} still
        
        pairs_lengths_r = []
        
        repeated_length = False
        
        while not repeated_length:
            
            while len(kneading) < S + m + l:
                kneading += period
                
            # Guarantees that len(keading) >= S + m + l >= S_{k+1}
            
            S_k = S
            length_of_B_k = (S_k - m) % l
            
            j = S + 1
            difference_found = False
            
            while not difference_found:
                eta_lower = kneading[j - 1 - S]
                eta_upper = kneading[j - 1]
                
                if eta_lower == eta_upper:
                    j += 1
                
                else:
                    difference_found = True
            
            S_next = j
            r_k = S_next - S_k
            
            pairs_lengths_r.append( (length_of_B_k, r_k) )
            
            S = S_next
            address.append(S)
                    
            length_of_B_next = (S - m) % l
            
            # Now, check if this length has already appeared, and keep track of the first time it appeared as an index
            # There probably is a more efficient way of doing this withs dicts and proper indexing
            
            L = len(pairs_lengths_r)
            
            for i in range(L):
                if length_of_B_next == pairs_lengths_r[i][0]:
                    repeated_length = True
                    index = i
                
        # Now the current S = S_k has a "periodic" initial segment B_k,
        # whose length has appeared before in the list [(l_i, r_i)]
        
        # create list of the periodic orbit of the lengths, and corresponding r, using index

        periodic_lengths_r = pairs_lengths_r[index : ]
        p = len(periodic_lengths_r)
        
        # Now need to compute internal address until all the lengths are marked as "good"
        
        # Recall that r in {1, ..., m+l}
        # If an r = r_k repeats, and S_k was attested to be in the rho-orbit of r, don't need to compute the rho_orbit again for that r
        
        r_seen_and_S_in_orbit = {}
        for R in range(1, m+l+1):
            r_seen_and_S_in_orbit[R] = False
        
        i = 0 # cyclic counter. Note that periodic_lengths_r[0] = (l_k, r_k), for the current S = S_k
        all_good = False
        
        good_indices_dict = {}
        number_of_good_indices = 0
        
        for t in range(p):
            good_indices_dict[t] = 0 # 0 represents not good
        
        while not all_good:
            
            while len(kneading) < S + m + l:
                kneading += period
                
            # Guarantees that len(keading) >= S + m + l >= S_{k+1}
            
            S_k = S
            r_k = periodic_lengths_r[i][1]
            S_next = S_k + r_k
            
            # Here we have S_k, S_{k+1}, and the difference r_k from the list of pairs obtained
            
            # Need to check if S_k is in the rho-orbit of r. If yes, mark index i as good; if not, continue
            # Recall that the S_k are just the rho-orbit of 1, which is is the internal address
            # Hence, we only need to see if the orbits merge
            
            # For a given iterate r_orbit, starting with r = r_k, check if it is in the internal address
            
            # Keep track of those r for which S_k is in the orbit; skip these calculations for future S_{k'}
            
            if r_seen_and_S_in_orbit[r_k] == True:
                # if this r was already seen, and its orbit already merges with the internal address,
                # Then we skip calculations below and just mark the index as good. Denominator guaranteed to be 2
                    
                if good_indices_dict[i] == 0:
                    number_of_good_indices += 1
                    good_indices_dict[i] = 1
                
                    # Here, we increase the counter of the number of good indices by 1 exactly if we the index was
                    # not good, and we changed it to good. If it was already good, we don't increase the counter
                        
        
            r_orbit = r_k # change index to treat as iteration, to find rho-orbit of r_k
            
            S_k_in_the_orbit_of_r_k = False
            
            while not (S_k_in_the_orbit_of_r_k) and r_orbit <= S_k:
                
                address_index = 0
                    
                while address[address_index] <= r_orbit:
                        
                    if address[address_index] == r_orbit: # checks to see if the current r_orbit appears in the internal address
                            
                        r_seen_and_S_in_orbit[r_k] = True 
                        S_k_in_the_orbit_of_r_k = True
                            
                        if good_indices_dict[i] == 0:
                            number_of_good_indices += 1
                            good_indices_dict[i] = 1
                                
                            # Again increase the counter of the number of good indices by 1 exactly if we the index was
                            # not good, and we changed it to good. If it was already good, we don't increase the counter
                    
                    address_index += 1
                    
                if not S_k_in_the_orbit_of_r_k: 
                    # r_orbit hasn't appearead in the internal address computed so far,
                    # Find the iterate rho(r_orbit)
                
                    j = r_orbit + 1
                    difference_found = False
                
                    while not difference_found:
                        
                        # Since I haven't proven the claim ** above, for safety, include a check
                        
                        if j > len(kneading):
                            kneading += period
                        
                        eta_lower = kneading[j - 1 -r_orbit]
                        eta_upper = kneading[j - 1]
                        
                        if eta_lower == eta_upper:
                            j += 1
                            
                        else:
                            difference_found = True
                            r_orbit = j
                            
                    # This iterates the orbit of r, indexed as r_orbit
                    # If r_orbit ever surpasses S_k, or is equal, stop
                     
            # Now check if all indices are good:
                
            if number_of_good_indices == p:
                
                all_good = True
                
                S = S_next
                address.append(S)
                
                # This guarantees that, for preperiodic theta, the last entry computed has denominator 2.
                
            else:
                S = S_next
                address.append(S)
                
                i = (i+1) % p # Recall that it is a cyclic index, ranging over the period of lengths
            
        # The loop only ends when all indices are good: must be the case eventually.
        # Guarantees that for any other entry of the internal address not computed, the denominator is 2
        
        address.append('...')
        
        return address
    
    # Remark: There is considerable redundancy. Could cut out the extraneous trailing 2's if needed

####

def denominators_angled_internal_address(theta):
    # Input is a rational theta in the form of a Fraction
    # Output is the list of denominators in the angled internal address of theta
    # If theta has odd denominator, the list is finite. if not, all denominators after the listed ones are 2
    
    # Recall that, in the function internal_address, for preperiodic theta, the last S computed has denominator 2
    
    k_list = kneading_list(theta)
    
    preperiod = k_list[0]
    period = k_list[1]
    
    m = len(preperiod)
    l = len(period)
    
    address = internal_address(theta)
    
    if address[-1] == '...':
        address.pop(-1)
    
    highest_S = address[-1]
    
    D = max(math.ceil((highest_S / l) + 1), m)
    kneading = preperiod + (D * period)
    
    # This guarantees that we can still compute S_{k+1} from kneading. Possibly redundant
    
    M = len(address) - 1    
    k = 0
    
    denominators = []
    
    while k < M:
        
        S_k = address[k]
        S_next = address[k+1]
        
        r_k = S_next % S_k
        
        # As in internal_address, need to compute the rho-orbit of the r and see if it merges with the internal address
        # Can speed up computations by seeing if the iterate of r is present in the internal address before S_k
        
        if r_k == 0:
            
            q_k = S_next // S_k
            
        else:
            
            r_orbit = r_k # change index to treat as iteration, to find rho-orbit of r_k
            
            S_k_in_the_orbit_of_r_k = False
            
            while not (S_k_in_the_orbit_of_r_k) and r_orbit <= S_k:
                
                address_index = 0
                    
                while address[address_index] <= r_orbit:
                        
                    if address[address_index] == r_orbit: # checks to see if the current r_orbit appears in the internal address
                            
                        S_k_in_the_orbit_of_r_k = True
                            
                    address_index += 1
                    
                if not S_k_in_the_orbit_of_r_k: 
                    # r_orbit hasn't appearead in the internal address computed so far,
                    # Find the iterate rho(r_orbit)
                
                    j = r_orbit + 1
                    difference_found = False
                
                    while not difference_found:
                        
                        # Since I haven't proven the claim ** above, for safety, include a check
                        
                        if j > len(kneading):
                            kneading += period
                        
                        eta_lower = kneading[j - 1 -r_orbit]
                        eta_upper = kneading[j - 1]
                        
                        if eta_lower == eta_upper:
                            j += 1
                            
                        else:
                            difference_found = True
                            r_orbit = j
                            
                    # This iterates the orbit of r, indexed as r_orbit
                    # If r_orbit ever surpasses S_k, or is equal, stop
                     
            if S_k_in_the_orbit_of_r_k:
                q_k = ((S_next - r_k) // S_k) + 1
            else:
                q_k = ((S_next - r_k) // S_k) + 2
                
        denominators.append(q_k)
        k += 1
        
    if m == 0:
        denominators.append(1)
        
    else:
        denominators.append(2)
        denominators.append('...')

    return denominators

####

def numerators_angled_internal_address(theta):
    # returns the numerators in the angled internal address of an angle theta
    # Everything after the last 1 is a continuing list of 1's
    
    k_list = kneading_list(theta)
    
    preperiod = k_list[0]
    
    m = len(preperiod)
    
    address = internal_address(theta)
    denominators = denominators_angled_internal_address(theta)
    
    if m != 0:
        M = len(address) - 2
    else:
        M = len(address) - 1
    
    numerators = []
    
    k = 0
    
    while k < M:
        q_k = denominators[k]
        
        if q_k == 2:
            p_k = 1
        else:    
            S_k = address[k]
            multiplier = 2 ** (S_k)
            
            count = 0
            theta_now = theta
            
            for i in range(q_k - 1):
                
                if theta_now <= theta:
                     count += 1
                     
                theta_now = (multiplier * theta_now) % 1
            
            p_k = count
            
        numerators.append(p_k)
        k += 1
        
    if m != 0:
        numerators.append(1) # This is because we guaranteed that the last entry in the internal address had denominator 2
        numerators.append('...')
        
    else:
        numerators.append(0)
    
    return numerators

####

#######################
# AUXILIARY FUNCTIONS #
#######################

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

def preperiodic_to_binary(fraction):
    # input is a fraction
    # Output is a list, first entry is the binary preperiod, second is the binary period
    
    iterate = fraction

    binary_preperiod = ''
    
    while iterate.denominator % 2 == 0:
        if 0 <= iterate and iterate < Fraction(1,2):
            binary_preperiod += '0'
        else:
            binary_preperiod += '1'
        iterate = (2 * iterate) % 1
        
    binary_period = ''
    
    first_periodic = iterate
    
    if 0 <= iterate and iterate < Fraction(1,2):
        binary_period += '0'
    else:
        binary_period += '1'
        
    iterate = (2 * iterate) % 1
        
    while iterate != first_periodic:
        if 0 <= iterate and iterate < Fraction(1,2):
            binary_period += '0'
        else:
            binary_period += '1'
        iterate = (2 * iterate) % 1 
        
    b_list = [binary_preperiod, binary_period]
    
    return b_list

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
    # Returns its finite orbit as a list
    
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

####

def periodic_branch_portraits_dyadic(theta):
    # Input is a dyadic angle theta
    # Output is the list of all periodic orbit portraits (Milnor) for the polynomial z^2 + c_theta,
    # having rotation number p/q, q >= 3
    # These correspond to periodic branch points on the Hubbard tree and their external arguments
    
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
        
        # Lemma: For any theta, suppose that for some k we have q_k >= 3.
        # Then, if phi_k is the periodic approximation of theta of period S_k,
        # phi_k lands at the root of the hyperbolic component corresponding to S_k
        
        # This allows us, by tuning a cardioid angle, to find the rays landing at the q_k * S_k satellite bifurcation
        
        q_k = denominators[k]
               
        if q_k >= 3:
            
            S_k = int_address[k]
            p_k = numerators[k]
            
            approx_string = bin_string[1:S_k + 1]
            
            approx_num = int(approx_string, 2)
            
            angle = Fraction(approx_num, (2 ** S_k) - 1)
            
            # This works because of our estimates for dyadics:
            # If m is the preperiod, and k is such that q_k >= 3,
            # then S_k =< (m-1)/2 < m
            # bin_string has length m, the preperiod of theta, so cutting it off at S_k doesn't leave out anything
            # (even the trailing zeros which would be necessary to find the approximating periodic angle)
            
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

def periodic_branch_portraits(theta):
    # Input is a RATIONAL angle theta
    # Output is the list of all periodic orbit portraits (Milnor) for the polynomial z^2 + c_theta,
    # having rotation number p/q, q >= 3
    # These correspond to periodic branch points on the Hubbard tree and their external arguments
    
    # Remark: need to be careful when theta is an angle landing at a satellite component of bifurcation q > 2
    
    address = internal_address(theta)
    denominators = denominators_angled_internal_address(theta)
    numerators = numerators_angled_internal_address(theta)
    
    if theta.denominator % 2 == 0:
        M = len(address) - 1 # recall that for preperiodic theta, the internal address ends with '...'
    else:
        M = len(address)
        
    preperiod_bin_string = preperiodic_to_binary(theta)[0]
    period_bin_string = preperiodic_to_binary(theta)[1]
    
    bin_string = preperiod_bin_string + period_bin_string
    
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
    
    # now k = 1, iterate it to M - 1
    
    while k < M:
        
        # Lemma: For any theta, suppose that for some k we have q_k >= 3.
        # Then, if phi_k is the periodic approximation of theta of period S_k,
        # phi_k lands at the root of the hyperbolic component corresponding to S_k
        
        # This allows us, by tuning a cardioid angle, to find the rays landing at the q_k * S_k satellite bifurcation
        
        q_k = denominators[k]
               
        if q_k >= 3:
            
            S_k = address[k]
            p_k = numerators[k]
            
            while S_k > len(bin_string):
                
                bin_string += period_bin_string
                
            # This guarantees that the binary string for theta has at least S_k entries
            # Hence we may take the approximation
            
            # Differently from dyadics, the string doesn't start with a dot: .b1b2...
            
            approx_string = bin_string[:S_k]
            
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

####

def preimages(angles_at_point, diameter):
    # Input is a collection of angles landing at a given point, and a diameter (major leaf)
    # Output are the two preimage points, represented by the angles landing on them,
    # respecting the diamater partition
    
    open_arcs = diameter.open_arcs
    
    reference_arc = open_arcs[0]
    preimage_0 = []
    preimage_1 = []
    
    for angle in angles_at_point:
        
        preimage_angle_A = angle / 2
        preimage_angle_B = (1 + angle) / 2
        
        if in_arc(preimage_angle_A, reference_arc):
            
            preimage_0.append(preimage_angle_A)
            preimage_1.append(preimage_angle_B)
            
        else:
            
            preimage_1.append(preimage_angle_A)
            preimage_0.append(preimage_angle_B)
    
    return [preimage_0, preimage_1]

####

def angles_landing_at_P_preperiodic(theta):
    # Given an angle theta preperiodic under doubling
    # Finds all the rays landing on the postcritical set P, grouped into points
    
    orbit = preperiodic_orbit(theta)
    kneading = kneading_list(theta)
    diameter = major_leaf(theta)
    
    m = len(kneading[0]) # preperiod
    l = len(kneading[1]) # orbit period
    
    n = len(orbit) - m # ray period
    
    periodic_orbit = orbit[m:]
    
    r = n // l # number of rays landing at each point
    
    if r > 1: # critical value falls onto a satellite orbit
        
        angles_at_P = []
        for i in range(l):
            angles_at_P.append([])
            
        for i in range(n):
            angle = periodic_orbit[i]
            
            angles_at_P[i % l].append(angle)
            
        index = m - 1
        
        current_point = angles_at_P[0]
        
        while index >= 0:
            
            preimage_0, preimage_1 = preimages(current_point, diameter)
            
            found_angle_of_orbit_in_preimage_0 = False
            i = 0
            
            while i < r and found_angle_of_orbit_in_preimage_0 == False:
                angle = preimage_0[i]
                
                if fractions_equal(orbit[index], angle):
                
                    found_angle_of_orbit_in_preimage_0 = True
                    
                i += 1
                
            if found_angle_of_orbit_in_preimage_0 == True:
                current_point = preimage_0
                
            else:
                current_point = preimage_1
            
            angles_at_P = [current_point] + angles_at_P
                
            index -= 1
    
        return angles_at_P

    else: # r = 1, the periodic portrait of the critical orbit is trivial or primitive
    # To see if the portrait is non-trivial: finds all partners of the periodic angles
    # if theta is in the angle sector determined by some periodic angle and its partner,
    # Then a portrait is realized. It must be primitive because the satellite case is 
    # Remark: only one portrait is realized, because the characteristic leaf determines it
    
    # Note: if theta is dydic, always returns trivial portrait, as expected
    # Because in the cose below, angle = partner = 0
        
        found_nontrivial_portrait = False
        
        j = 0
        
        while found_nontrivial_portrait == False and j < n:
            
            angle = periodic_orbit[j]
            partner = partner_angle(angle)
            
            if angle < partner:
                if angle < theta and theta < partner:
                    # record that they land together
                    
                    found_nontrivial_portrait = True
                    
                    char_angle = angle
                    char_partner = partner
                    
                        
            else:
                if partner < theta and theta < angle:
                    # record that they land together
                    
                    found_nontrivial_portrait = True
                    
                    char_angle = angle
                    char_partner = partner
        
            j += 1
        
        # Because the satellite case is already considered above, we know that the
        # portrait, if found, must be primitive
        
        angles_at_P = []
        
        if found_nontrivial_portrait == False:
            
            for element in orbit:
                
                angles_at_P.append([element])
                
                return angles_at_P
                
        else: # Presence of primitive portrait
            # char_angle and char_partner for the characteristic leaf for the corresponding primitive orbit
            # j -1 is the index in the periodic orbit
            
            start_j = j - 1
            
            angles_at_P = []
            current_angle, current_partner = char_angle, char_partner
            
            while j - 1 < n:
                
                angles_at_P.append([current_angle, current_partner])
                
                j += 1
                
                current_angle = (2 * current_angle) % 1
                current_partner = (2* current_partner) % 1
                
            # Now j = n, return to start of period
            
            for i in range(0, start_j):
                
                angles_at_P = [[current_angle, current_partner]] + angles_at_P 
                
                current_angle = (2 * current_angle) % 1
                current_partner = (2* current_partner) % 1
            
            # Now need to pullback to the preperiodic ones
            
            index = m - 1
            
            current_point = angles_at_P[0]
            
            while index >= 0:
                
                preimage_0, preimage_1 = preimages(current_point, diameter)
                
                found_angle_of_orbit_in_preimage_0 = False
                i = 0
                
                while i < r and found_angle_of_orbit_in_preimage_0 == False:
                    angle = preimage_0[i]
                    
                    if fractions_equal(orbit[index], angle):
                    
                        found_angle_of_orbit_in_preimage_0 = True
                        
                    i += 1
                    
                if found_angle_of_orbit_in_preimage_0 == True:
                    current_point = preimage_0
                    
                else:
                    current_point = preimage_1
                
                angles_at_P = [current_point] + angles_at_P
                    
                index -= 1
        
            return angles_at_P

####

def angles_landing_at_roots_periodic(theta):
    # Input is a periodic theta under doubling
    # Output is the collection of angles landing at the roots of the periodic Fatou components
    
    # Less computationally intensive way of finding the period of theta:
    
    n = 1
    current_orbit = (2 * theta) % 1
    
    while not fractions_equal(current_orbit, theta):
    
        n += 1
        current_orbit = (2 * current_orbit) % 1
    
    # now n is the period
        
    
    angles_at_roots = []
    
    current_angle = theta
    current_partner = partner_angle(theta)
    
    for i in range(n):
        
        angles_at_roots.append([current_angle, current_partner])
        
        current_angle = (2 * current_angle) % 1
        current_partner = (2 * current_partner) % 1
        
    return angles_at_roots
    
####

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

def forbidden_region(theta):
    # Input is a rational angle theta
    # Output is the forbidden_region of the Hubbard tree of f_{c_{theta}}
    
    # To find the forbidden region, need to find the attaching points of the preimage Hubbard tree
    # Divide into two cases: attaching point is in P or not
    
    # Not in P: it is a branch point of the preimage tree, and maps into a periodic branch cycle
    # Need to find all periodic branch portraits, iterate preimages until leaves stop separating P
    
    # In P: Find those points in P where the Hubbard tree gains a new branch at its image
    # if f hyperbolic, only happens once, if f subhyperbolic, may happen many times
    
    # Danger: Orbit of P could fall into a periodic branch orbit (ex.: 9/56), want to treat these
    # separately. This is because that determine whether the point is a branch point,
    # or possible attaching point, don't work well if the point itself is in P
    # (this was not possible in the dyadic case)

    P = preperiodic_orbit(theta)
    diameter = major_leaf(theta)
    per_branch_portraits = periodic_branch_portraits(theta)
    
    m = len(kneading_list(theta)[0])
    
    if theta.denominator % 2 == 0:
        
        periodic = False
        angles_at_P = angles_landing_at_P_preperiodic(theta)
        
    else:
        periodic = True
        angles_at_roots = angles_landing_at_roots_periodic(theta)
    
    # We want to exclude the branch portrait corresponding to the periodic orbit of preperiodic theta
    # This only happens if theta is preperiodic and three or more rays land at each point of P
    
    if periodic == False:
        if len(angles_at_P[0]) >= 3:
        
            found_postcritical_portrait = False
            
            num_portraits = len(per_branch_portraits)
            i = 0
            
            while i < num_portraits and (found_postcritical_portrait == False):
                
                portrait = per_branch_portraits[i]
                num_points_in_portrait = len(portrait)
                j = 0
                
                while j < num_points_in_portrait and (found_postcritical_portrait == False):
                    
                    angles_landing_at_point = portrait[j]
                    num_angles_landing_at_point = len(angles_landing_at_point)
                    k = 0
                    
                    while k < num_angles_landing_at_point and (found_postcritical_portrait == False):
                        
                        angle = angles_landing_at_point[k]
                        
                        num_forward_orbit_theta = len(P)
                        l = 0
                        
                        while l < num_forward_orbit_theta and (found_postcritical_portrait == False):
                            
                            postcritical_angle = P[l]
                            
                            if fractions_equal(angle, postcritical_angle):
                                found_postcritical_portrait == True
                                index_of_postcritical_portait = i
                                
                            l += 1    
                        k += 1            
                    j += 1                
                i += 1        
                
            per_branch_portraits.pop(index_of_postcritical_portait)
    
    # Now we've excluded the postcritical portrait from those branch point
    # Treat postcritical portrait of P separately after
    
    candidate_points = []
    
    for portrait in per_branch_portraits:
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
            
    # Obtained list of forbidden arcs for the preperiodic preimages of periodic branch points
    # Now treat attaching points in P separately
    
    if periodic == False:
        if len(angles_at_P[0]) >= 2:
            
            # Need to find postcritical attaching points
            # They are exactly the postcritical points c_i for which f(c_i) has more branches than c_i
            # Still measure this with respect to leaves and separation of P, but ignore the point itself
            
            for i in range(m):
                c_i = angles_at_P[i]
                f_c_i = angles_at_P[i+1]
                
                # If a leaf of c_i doesn't separate P \ {angles at c_i},
                # and its image leaf of f_c_i separates P \ {angles at f_c_i},
                # then it corresponds to a forbidden arc
                
                # Remark: because of non-periodicity, only need to exclude, from P, the angle P[i] landing at c_i
                
                P_excluded_theta_i = P.copy()
                P_excluded_theta_i.pop(i)
                
                if i != m-1:
                    # only need to exclude more angles from f_c_i when it is the first periodic point
                    
                    P_excluded_f_c_i = P.copy()
                    P_excluded_f_c_i.pop(i+1)
                    
                else:
                    # If i == m-1, last preperiodic point in the orbit
                    # Will need to exclude more angles from P in this case: all of those who land at f_c_i
                    
                    P_excluded_f_c_i = []
                    
                    for iterate_of_theta in P:
                        for image_angle in f_c_i:
                            if not fractions_equal(iterate_of_theta, image_angle):
                                P_excluded_f_c_i.append(iterate_of_theta)
                        
                    # This removes from P all those angles landing at f_c_i, to consider separation of leaves
                
                for angle in c_i:
                    next_angle = next_in_cyclic_order(angle, c_i)
                    
                    open_arc = arc(angle, next_angle, False)
                    
                    # Check if arc doesn't contain points of P_excluded_theta_i
                    # Then if the image arc has points of P_excluded_f_c_i inside and outside
                    
                    N = len(P_excluded_theta_i)
                    counter1 = 0
                    
                    doesnt_contain_points = True
                    
                    while counter1 < N and doesnt_contain_points:
                        
                        if in_arc(P_excluded_theta_i[counter1], open_arc):
                            doesnt_contain_points = False
                        counter1 += 1
                    
                    if doesnt_contain_points:
                        
                        # Now need to check if the image arc separates P_excluded_theta_i_plus_1
                        # Remark: arc will always be < 1/2 in length
                        # If it had length >= 1/2, image would map non-injectively, hence c_0 is in the arc
                        # But because c_0 is in the tree, there must be endpoints of the tree in the arc too
                        # This implies we can consier image_arc without problem
                        
                        open_arc_image = image_arc(open_arc)
                        
                        N_image = len(P_excluded_f_c_i)
                        
                        counter2 = 0
                        points_in_open_image_arc = 0
                        points_notin_open_image_arc = 0
                        
                        while counter2 < N_image and (points_in_open_image_arc == 0 or points_notin_open_image_arc == 0):
                        
                            if in_arc(P_excluded_f_c_i[counter2], open_arc_image):
                                points_in_open_image_arc += 1
                            else:
                                points_notin_open_image_arc += 1
                            counter2 += 1
                        
                        if not (points_in_open_image_arc < 1 or points_notin_open_image_arc < 1):
                            # so that there are PCF points inside and outside of the image arc:
                            # add the leaf as a forbidden arc!
                            # either open_arc0 or open_arc1, whichever has 0 points in it
                            
                            forbidden_arcs.append(open_arc)
    
    else: # periodic = True
        # Repeat roughly same process: find the open arcs that dont contain points of P \ {...},
        # and whose image separates P \ {...}
        # Remark: careful with satellite orbits, may overcount forbidden arcs of the roots
        # don't need to consider satellite cases 
        if len(angles_at_roots) == 2:
            
            return
        
        # NEED TO FIX ISSUE WITH UPPER KNEADING SEQUENCES
    
    
    
    
    # Now we sort the forbidden arcs in cyclic order within the unit circle
    
    # Only nuance: if periodic, in the case of a rabbit, a root point of a periodic Fatou component
    # may coincide with a periodic satellite branch point

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
    # Note that the forbidden arcs are already ordered counterclockwise, starting from 0
    
    forbidden_arcs = forbidden_region_dyadic(theta)
    A = len(forbidden_arcs)
    
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
    
    i = 0 # The address
    c = 0 # Index for the arcs between the marked points
    
    while c < m - 1:
        closed_arc = arc(marked_points[c], marked_points[c+1], True)
        
        not_a_forbidden_arc = True 
        
        for forbidden_arc in forbidden_arcs:
            
            if (closed_arc.left_endpoint == forbidden_arc.left_endpoint and closed_arc.right_endpoint == forbidden_arc.right_endpoint):
                
                not_a_forbidden_arc = False
                
        # Verifies if the arc between the marked points is forbidden or not
        # If not, adds it to the list of states with its name and address
        
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
        
    # Adds last arc between marked points, of the form [t,0]:
    
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
        
    number_of_non_degenerate_arc_states = i
            
    i += 1
    
    # Now we must to attribute names and addresses to the degenerate arcs and their images.
    # The solution to understanding the degenerate arcs in the Markov partition is that each one,
    # representing and isolated angle landing on the Hubbard tree, must fall into a periodic cycle
    # Hence we add an extra state of the automaton for each element of its orbit, and this orbit is preperiodic
    
    degenerate_orbits_addresses = [] # Keeps track of the orbits of the isolated angles / degenerate arcs
    
    for degenerate_arc in degenerate_arcs:
        
        degenerate_orbit = []
        
        angle = degenerate_arc.left_endpoint
        
        name = 'I_' + str(i)
        address = i
        
        if angle < Fraction(1,2):
            bit = 0
        else:
            bit = 1
            
            # Note that 1/2 and 0 are never isolated angles / degenerate arcs, so there's no ambiguity with bits
            
        states.append([name, address, bit])
        dict_addresses_to_arcs[i] = degenerate_arc
        
        degenerate_orbit.append(i)
        
        i += 1
    
        orbit = preperiodic_orbit(angle)
        
        # Now we add the states for each element in the forward orbit of the isolated angles.
        
        for angle_in_orbit in orbit[1:]:
            
            name = 'P_' + str(i)
            address = i
            
            if angle_in_orbit < Fraction(1,2):
                bit = 0
            else:
                bit = 1
            
            states.append([name, address, bit])
            dict_addresses_to_arcs[i] = arc(angle_in_orbit, angle_in_orbit, True)
            
            degenerate_orbit.append(i)
            
            i += 1
        
        degenerate_orbits_addresses.append(degenerate_orbit)
        
    # degenerate_orbits_addresses keeps track of the addresses of the orbits of the degenerate arcs / isolated angles
            
    # We have finished adding all the states with their names, addresses and bits.
    # We also kept track of the isolated angles and their names.
    
    # Now we must construct the transition matrix.
    # For this, we see how each arc-state maps over other arc states, and over the isolated angles
    # Next, for each isolated angle, we just map it along its own separate orbit and created states
    
    # Important remark: some of the P-states may actually have degenerate arcs that coincide.
    # Ex.: Fraction(15,16)
    # In this case, it is important that the orbit of the I and P-states only follows itself, and doesn't
    # lead into these other P-states from orbits of different isolated angles
    # To do this, we abuse the fact the orbits of isolated angles are indexed with consecutive addresses
    
    S = len(states)
    
    M = np.zeros((S, S), dtype = int)
    
    # By convention, M[i][j] means the transition from address j to i
    # So that matrix multiplication of states corresponds to column vectors: M*v
    
    # For non-degenerate arcs, will always map (bijectively) over a union of arcs
    
    for j in range(number_of_non_degenerate_arc_states + 1): # Here we range only over non-degenerate origin states first
        
        arc_j = dict_addresses_to_arcs[j]
        image_j = image_arc(arc_j)
        
        for i in dict_addresses_to_arcs:
            
            arc_i = dict_addresses_to_arcs[i]
            
            if arc_i.degenerate == False: # if destination arc is not degenerate, proceed as normal
                if arc_inclusion(arc_i, image_j):
                    M[i][j] = 1
                    
            if arc_i.degenerate == True: # if destination arc is degenerate, should only be mapped to if its an isolated angle from the forbidden region
                if arc_inclusion(arc_i, image_j) and arc_i in degenerate_arcs:
                    M[i][j] = 1
                    
    # Now we need to consider transitions from degenerate origin states.
    # For this, we use the list degenerate_orbits_addresses and the orbits induced
    
    for degenerate_orbit in degenerate_orbits_addresses:
        
        length = len(degenerate_orbit)
        
        for l in range( length - 1):
            M[degenerate_orbit[l+1], degenerate_orbit[l]] = 1
            
        M[degenerate_orbit[0]][degenerate_orbit[-1]] = 1
    
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



