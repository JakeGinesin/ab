#!/usr/bin/env python3
from __future__ import annotations

from automata.finite_automata import Finite_Automata
from typing import Set, Tuple, List, Union, Dict
import heapq
from z3 import *

"""
follows from the exact formal definition of a NFA
which is a 5-tuple, (Q, Sigma, T, q0, F)
where T : Q x Sigma -> P(Q) (where P(Q) is the powerset of Q)

in this library we include epsilon transitions. e.g., we can
make a transition without consuming a character from a word.
"""
class NFA(Finite_Automata):
    def __init__(
        self,
        states : set,
        initial_state : str, 
        acceptance_states : set,
        alphabet : set,
        transitions : List[Tuple[str, str, str]] = []):

        Finite_Automata.__init__(self, states, initial_state, acceptance_states, alphabet, transitions)

        assert "" not in self.alphabet, "epsilon transitions (defined just as empty transitions, "") not allowed"
        self.alphabet.append('')

        self.construct_transition_function_hashmap()
        self.construct_reachability_hashmap()

    def test_word(self, word : str, current_state : str) -> Set[str]:
        """
        test_word: word, initial state (opt) -> (set of acceptance states reached)
        """
        assert current_state in self.states, "test_word: initial state not in NFA"

        if word == [] and current_state not in self.acceptance_states : return set()
        if word == [] and current_state in self.acceptance_states:
            ret = set()
            ret.add(current_state)
            return ret

        wordc = word.copy()
        current_char = wordc.pop(0)
        assert current_char in self.alphabet, "test_word: character not in alphabet"
        reachable = self.transition_function(current_state, current_char)
        reachable_epsilon = self.transition_function(current_state, '')
        accept_ret = set()
        for state in reachable:
            accept_ret = accept_ret.union(self.test_word(wordc, state))

        for state in reachable_epsilon:
            accept_ret = accept_ret.union(self.test_word(word.copy(), state))

        # must consume entire word before accepting
        if len(reachable) == 0 and len(reachable_epsilon) == 0 and word != []:
            return set()

        return accept_ret

    def to_single_acceptance_state(self) -> NFA:
        """
        returns a NFA that accepts the same langauge as self, 
        but the automata itself has a single accepting state.
    
        this is useful when reasoning about the language of the NFA with SMT/LP
                    
        also, this has the side effect of removing dead (unconnected) states
        """
        if not hasattr(self, "reachability_hashmap") : self.construct_reachability_hashmap()

        s = []
        s.append(self.initial_state)
        discovered=set()

        nfa_initial = self.initial_state
        nfa_alphabet = self.alphabet.copy()
        nfa_alphabet.remove('')
        nfa_acceptance_state = "s" + str(len(self.states) + 1)
        nfa_transitions = self.transitions.copy()

        end_state_transition_mapping = {}
        for state in self.states : end_state_transition_mapping[state] = set()
        for (sA, i, sB) in self.transitions : end_state_transition_mapping[sB].add((sA, i, sB))

        while s != []:
            v = s.pop()
            if v not in discovered:
                if v in self.acceptance_states:
                    for (sA, i, sB) in end_state_transition_mapping[v]:
                        nfa_transitions.append((sA, i, nfa_acceptance_state))

                # and, we don't add v to acceptance states of the created nfa
                discovered.add(v)
                reachable = self.reachability_hashmap[v]
                for state, char in reachable : s.append(state)

        nfa_states = set(discovered)
        nfa_states.add(nfa_acceptance_state)

        new_transitions = []
        for transition in nfa_transitions:
            if transition[0] in nfa_states and transition[2] in nfa_states:
                new_transitions.append(transition)

        return NFA(nfa_states, nfa_initial, [nfa_acceptance_state], nfa_alphabet, new_transitions)

    def minimal_parikh_image(self, init_state : str = None, usable_states : set = set()) -> Dict[str, int]:
        """
        a parikh image is simply a vector representing the number of characters in any given accepted string.
        finding the *minimal* parikh image for a NFA/buchi automata is trivially done using SMT/ILP. 
        see https://arxiv.org/abs/1002.1464 for details

        this is also in relation to parikh automata 
        """

        if init_state == None : init_state = self.initial_state
        # note: we only consider states relevant to the language of the automata
        if usable_states == set() : usable_states = self.get_wellconnected_states()

        # if init state is an accepting state, our parikh image is empty
        if init_state in self.acceptance_states : return {char: 0 for char in self.alphabet} 

        incoming_transitions, outgoing_transitions = {}, {}
        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()
    
        z3_transitions = []
        for t in self.transitions:
            if t[0] in usable_states and t[2] in usable_states : z3_transitions.append(t)

        # maintain a mapping between transition labels and actual transitions
        transition_mapping = {}
        transition_labels = ["t" + str(i) for i in range(len(z3_transitions))]

        count = 0
        for (sA, i, sB) in z3_transitions:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            outgoing_transitions[sA].add("t" + str(count))
            incoming_transitions[sB].add("t" + str(count))

            count = count + 1

        # plug flow constraints into z3

        solver = Optimize()
        tosolve = {}
        accepting_states_to_solve = set()
        for t in transition_labels : tosolve[t] = Int(str(t))
        for t in tosolve.values() : solver.add(t >= 0) # transitions must not have a negative value

        solver.minimize(Sum(list(tosolve.values())))

        for state in usable_states:
            if state in self.acceptance_states : accepting_states_to_solve.add(state)
            elif state == init_state:
                # want to ensure initial state has one more outgoing than incoming
                solver.add(Sum([tosolve[t] for t in outgoing_transitions[state]]) - 0 == Sum([tosolve[t] for t in incoming_transitions[state]]) + 1)
            else:
                # want to add [all outgoing] == [all incoming] to solver 
                solver.add(Sum([tosolve[t] for t in outgoing_transitions[state]]) == Sum([tosolve[t] for t in incoming_transitions[state]]))

        # if we have 2 accepting states, we only want our parikh image to consider the closer one
        solver.add(Or([Sum([tosolve[t] for t in outgoing_transitions[state]]) + 1 == Sum([tosolve[t] for t in incoming_transitions[state]]) + 0 for state in accepting_states_to_solve]))


        if solver.check() == sat:
            m = solver.model()
            solns = {var: m[var] for var in tosolve.values()}
            
            char_total = {}
            for char in self.alphabet : char_total[char] = 0

            for transition, value in solns.items():
                sA, i, sB = transition_mapping[str(transition)]
                char_total[i] += int(str(value))

            return char_total

        else:
            return {}

    def parikh_image_inclusion(self, parikh_image = Dict[str, int], init_state : str = None, usable_states : set = set()) -> bool:
        """
        determines if a given parikh image is included in the language of the buchi automata. 
        this problem can nicely be encoded as an ILP problem, and is NP hard via: https://arxiv.org/abs/1002.1464
        """

        if init_state == None : init_state = self.initial_state
        # note: we only consider states relevant to the language of the automata
        if usable_states == set() : usable_states = self.get_wellconnected_states()
        for char in self.alphabet : assert char in parikh_image, "parikh image doens't have all the required characters"

        # if init state is an accepting state, our parikh image is empty
        if init_state in self.acceptance_states : return {char: 0 for char in self.alphabet} 

        incoming_transitions, outgoing_transitions = {}, {}
        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        z3_transitions = []
        for t in self.transitions:
            if t[0] in usable_states and t[2] in usable_states : z3_transitions.append(t)

        # maintain a mapping between transition labels and actual transitions
        transition_mapping = {}
        transition_labels = ["t" + str(i) for i in range(len(z3_transitions))]

        # map characters to certain transitions
        character_mapping = {}
        for char in self.alphabet : character_mapping[char] = set()

        count = 0
        for (sA, i, sB) in z3_transitions:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            outgoing_transitions[sA].add("t" + str(count))
            incoming_transitions[sB].add("t" + str(count))
            character_mapping[i].add(transition_labels[count])

            count = count + 1

        # plug flow constraints into z3

        solver = Solver()
        tosolve = {}
        accepting_states_to_solve = set()
        for t in transition_labels : tosolve[t] = Int(str(t))
        for t in tosolve.values() : solver.add(t >= 0) # transitions must not have a negative value

        for char in self.alphabet:
            solver.add(
                Sum([tosolve[t] for t in character_mapping[char]]) == parikh_image[char]
            )

        for state in usable_states:
            if state in self.acceptance_states : accepting_states_to_solve.add(state)
            elif state == init_state:
                # want to ensure initial state has one more outgoing than incoming
                solver.add(Sum([tosolve[t] for t in outgoing_transitions[state]]) - 0 == Sum([tosolve[t] for t in incoming_transitions[state]]) + 1)
            else:
                # want to add [all outgoing] == [all incoming] to solver 
                solver.add(Sum([tosolve[t] for t in outgoing_transitions[state]]) == Sum([tosolve[t] for t in incoming_transitions[state]]))

        # if we have 2 accepting states, we only want our parikh image to consider the closer one
        solver.add(Or([Sum([tosolve[t] for t in outgoing_transitions[state]]) + 1 == Sum([tosolve[t] for t in incoming_transitions[state]]) + 0 for state in accepting_states_to_solve]))
       
        if solver.check() == sat:
            m = solver.model()
            solns = {var: m[var] for var in tosolve.values()}
            
            char_total = {}
            for char in self.alphabet : char_total[char] = 0

            for transition, value in solns.items():
                sA, i, sB = transition_mapping[str(transition)]
                char_total[i] += int(str(value))

            return char_total

        else:
            return {}

