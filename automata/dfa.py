#!/usr/bin/env python3

from automata.finite_automata import Finite_Automata
from typing import Set, Tuple, List 
import heapq

class DFA(Finite_Automata):
    def __init__(
        self,
        states : set,
        initial_state : str, 
        acceptance_states : set,
        alphabet : set,
        transitions : List[Tuple[str, str, str]] = []):

        Finite_Automata.__init__(self, states, initial_state, alphabet, transitions)

        self.construct_transition_function_hashmap()
        self.construct_reachability_hashmap()
