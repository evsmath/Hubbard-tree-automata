# hubbard_tree_automata
Given a normalized quadratic complex polynomial and its associated Hubbard tree, the set of angles whose rays land on the Hubbard tree
can be expressed in the form of a finite state automaton.

The python code in this repository computes this automaton; the relevant mathematical facts are contained in https://arxiv.org/abs/2611.04190. It also
computes the product mating automaton of two given quadratic polynomials, enconding the set of direct ray connections between their Hubbard trees.

Relevant functions and classes:
- forbidden_region(theta): computes the forbidden region of the Hubbard tree, which determines how the branches of the preimage Hubbard tree are attached to it, and what are the enclosed forbidden angle intervals.
- automaton: class containing two attributes, a numpy matrix representing the transition graph of an automaton, and a list containing the names, addresses, and bit values of the states of the automaton.
- automaton_and_addresses(theta): function that outputs the Hubbard tree automaton for the PCF polynomial given by the parameter external angle theta, and the dicitionary containing the angle intervals corresponding to the states of the automaton
- hubbard_tree_automaton(theta): function that outputs just the Hubbard tree automaton.
- mating(theta1, theta2): computes the mating automaton of two quadratic PCF polynomials respectively given by theta1 and theta2, according to the product automaton and trimming construction.

# data_gathering
Collection of function computing entropy, number of states, weakly connected components and other properties of the mating automata.


