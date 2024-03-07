#!/usr/bin/env python3
# Various algorithms to test intersection emptiness between automata

from automata.nfa import NFA
from automata.buchi import Buchi

from z3 import *
from typing import Union, List
import heapq
import itertools

def dijkstra_nfa_intersection_emptiness(nfa_array : List[NFA]) -> Union[str, None]:
    """
    the canonical implementation for testing intersection emptiness between
    nondeterministic finite automata via djikstra and on-the-fly construction

    either returns a string violating intersection emptiness, or returns None,
    which indicates intersection emptiness holds
    """
    assert len(nfa_array) > 1, "must take intersection emptiness between two or more automata"
    alphabet = nfa_array[0].alphabet
    for i in range(1, len(nfa_array)) : assert alphabet == nfa_array[i].alphabet, "intersection emptiness checking requires equivalent alphabets"

    # initialize 
    initial_state = tuple(nfa.initial_state for nfa in nfa_array)

    # we choose to maintain a dict to keep track of words while exploring states
    # alternatively, we could backtrack once we hit our desired condition 
    distance, word = {}, {}

    queue = []
    distance[initial_state] = 0
    word[initial_state] = "" 
    # push initial state onto the heap
    heapq.heappush(queue, (0, initial_state))

    while queue:
        current_distance, current_state = heapq.heappop(queue)

        # if we've already been to this state before, no need to re-explore it. 
        # note, we already push all possible explorable states onto the stack
        if current_distance > distance[current_state] : continue

        # acceptance condition: check if all NFAs being intersected are at an acceptance state
        if all([current_state[i] in nfa_array[i].acceptance_states for i in range(len(current_state))]): 
            return word[current_state]

        # go through alphabet, construct adjacent states on the fly
        reachable = []
        for char in alphabet:

            """
            when constructing the next on-the-fly state, we require for any given character
            that all NFAs make a transition. because of non-determinism, we observe for any given
            character we may be able to construct multiple possible states
            """
            new_states, break_loop = [], False
            new_state = [None]*len(nfa_array)
            for i in range(len(nfa_array)):
                new_state[i] = nfa_array[i].transition_function(current_state[i], char)

                if new_state == []:
                    break_loop = True
                    break

            if break_loop : continue

            # now, we have a new_state with (array, ..., array)
            new_state = tuple(new_state)
            for new_state in list(itertools.product(*new_state)):
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""
                reachable.append((new_state, char))

        for (state, char) in reachable:
            new_distance = distance[current_state] + 1

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state))
    
    return None

def dijkstra_buchi_intersection_emptiness(buchi_array : List[Buchi]) -> Union[str, None]:
    """
    the canonical implementation for testing intersection emptiness between
    nondeterministic buchi automata via dijkstra and on-the-fly construction

    either returns a string violating intersection emptiness (in the case of a buchi automata,
    a lasso containing an acceptance state, or returns None,
    which indicates intersection emptiness holds
    """

    assert len(buchi_array) > 1, "must take intersection emptiness between two or more automata" 
    alphabet = buchi_array[0].alphabet
    for i in range(1, len(buchi_array)) : assert alphabet == buchi_array[i].alphabet, "intersection emptiness checking requires equivalent alphabets"



    initial_state = tuple(buchi.initial_state for buchi in buchi_array)
    prev_state_mapping = {}
    prev_character_mapping = {}
    word = {}

    queue = []
    word[initial_state] = "" 
    heapq.heappush(queue, (None, initial_state))
    
    explored_states = set()
    tried_combos = set()
    explored_combos = set()
    prev_state_mapping = {}

    while queue:
        previous_state, current_state = heapq.heappop(queue)
        prev_state_mapping[current_state] = previous_state

        if current_state in explored_states and not (current_state, previous_state) in tried_combos:
            # if so, iterate backwards until we hit the same state
            hit_acceptance = False
            iterate = True
            iterate_state = current_state
            loop_word = ""
            hit_states = set()
            while iterate:
                loop_word+=prev_character_mapping[iterate_state]
                iterate_state = prev_state_mapping[iterate_state]
                if iterate_state in hit_states : break

                if all([iterate_state[i] in buchi_array[i].acceptance_states for i in range(len(iterate_state))]):
                    hit_acceptance = True
                    iterate = False
                else:
                    if iterate_state == current_state : iterate = False

                hit_states.add(iterate_state)
                if iterate_state == initial_state : break

            if hit_acceptance : return word[current_state] + loop_word[::-1]
            else : tried_combos.add((current_state, previous_state))

        if (previous_state, current_state) in explored_combos : continue
            
        explored_combos.add((previous_state, current_state))
        explored_states.add(current_state)

        # go through alphabet, construct adjacent states on the fly
        reachable = []
        for char in alphabet:

            # exact same on-the-fly construction for NFAs
            new_states, break_loop = [], False
            new_state = [None]*len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state[i], char)

                if new_state == []:
                    break_loop = True
                    break

            if break_loop : continue

            # now, we have a new_state with (array, ..., array)
            new_state = tuple(new_state)
            for new_state in list(itertools.product(*new_state)):
                if not new_state in word : word[new_state] = ""
                prev_character_mapping[new_state] = char
                reachable.append((new_state, char))

        for (state, char) in reachable:
            word[state] = word[current_state] + char
            heapq.heappush(queue, (current_state, state))

    return None

def z3_dijkstra_nfa_intersection_emptiness(nfa_array : List[NFA]) -> Union[str, None]:
    """
    we reason about the automata we're intersecting individually in order to navigate the 
    on-the-fly constructed search space (which is exponential in size)
    """

    assert len(nfa_array) > 1, "must take intersection emptiness between two or more automata"
    alphabet = nfa_array[0].alphabet
    for i in range(1, len(nfa_array)) : assert alphabet == nfa_array[i].alphabet, "intersection emptiness checking requires equivalent alphabets"
    for nfa in nfa_array : assert len(nfa.acceptance_states) >= 1, "all NFAs must have at least one acceptance state)"

    # initialize 
    initial_state = tuple(nfa.initial_state for nfa in nfa_array)

    # we choose to maintain a dict to keep track of words while exploring states
    # alternatively, we could backtrack once we hit our desired condition 
    distance, word = {}, {}

    queue = []
    distance[initial_state] = 0
    word[initial_state] = "" 
    # push initial state onto the heap
    heapq.heappush(queue, (0, initial_state))

    well_connected_states = {}
    for nfa in nfa_array : well_connected_states[nfa] = nfa.get_wellconnected_states()

    """
    use the minimal parikh image to get the distance to the acceptance state
    """
    def heuristic(state):
        ret_sum = 0
        for i in range(len(state)):
            initial_state, usable_states = state[i], well_connected_states[nfa_array[i]]
            if len(nfa_array[i].acceptance_states) >= 1 and initial_state == list(nfa_array[i].acceptance_states)[0] : continue
            parikh_img = nfa_array[i].minimal_parikh_image(initial_state, usable_states)
            if parikh_img == {} : return float('inf')
            for (char, value) in parikh_img.items() : ret_sum+=value

        return ret_sum


    while queue:
        current_distance, current_state = heapq.heappop(queue)

        # if we've already been to this state before, no need to re-explore it. 
        # note, we already push all possible explorable states onto the stack
        if current_distance > distance[current_state] : continue

        # acceptance condition: check if all NFAs being intersected are at an acceptance state
        if all([current_state[i] in nfa_array[i].acceptance_states for i in range(len(current_state))]): 
            return word[current_state]

        # go through alphabet, construct adjacent states on the fly
        reachable = []
        for char in alphabet:

            """
            when constructing the next on-the-fly state, we require for any given character
            that all NFAs make a transition. because of non-determinism, we observe for any given
            character we may be able to construct multiple possible states
            """
            new_states, break_loop = [], False
            new_state = [None]*len(nfa_array)
            for i in range(len(nfa_array)):
                new_state[i] = nfa_array[i].transition_function(current_state[i], char)

                if new_state == []:
                    break_loop = True
                    break

            if break_loop : continue

            # now, we have a new_state with (array, ..., array)
            new_state = tuple(new_state)
            for new_state in list(itertools.product(*new_state)):
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""
                reachable.append((new_state, char))

        for (state, char) in reachable:
            new_distance = distance[current_state] + heuristic(state) + 1

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state))
    
    return None

def z3_dijkstra_buchi_intersection_emptiness(buchi_array : List[Buchi]) -> Union[str, None]:
    """
    we reason about the automata we're intersecting individually in order to navigate the 
    on-the-fly constructed search space (which is exponential in size).

    uses a sophisticated heuristic to underapproximate the distance to the buchi acceptance condition
    (i.e. discovering a lasso)
    """

    assert len(buchi_array) > 1, "must take intersection emptiness between two or more automata" 
    alphabet = buchi_array[0].alphabet
    for i in range(1, len(buchi_array)) : assert alphabet == buchi_array[i].alphabet, "intersection emptiness checking requires equivalent alphabets"

    initial_state = tuple(buchi.initial_state for buchi in buchi_array)
    prev_state_mapping = {}
    prev_character_mapping = {}
    distance, word = {}, {}

    queue = []
    distance[initial_state] = 0
    word[initial_state] = "" 
    explored = set()
    heapq.heappush(queue, (0, initial_state))

    # initialize various dicts to hold all the heurisitc-required stuff
    well_connected_states, hit_acceptance_states, incoming_transitions_all, outgoing_transitions_all, transition_mapping_all = {}, {}, {}, {}, {}

    for buchi in buchi_array:
        well_connected_states[buchi] = buchi.get_wellconnected_states()
        incoming_transitions, outgoing_transitions, transition_mapping = {}, {}, {}

        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        usable_transitions = []
        for (sA, i, sB) in buchi.transitions: 
            if sA in usable_states and sB in usable_states:
                usable_transitions.append((sA, i, sB))

        transition_labels = [i for i in range(len(usable_transitions))]

        count = 0
        for (sA, i, sB) in ts:
            transition_mapping[transition_labels[count]] = (sA, i, sB)

            outgoing_transitions[sA].add(count)
            incoming_transitions[sB].add(count)

            count = count + 1

        incoming_transitions_all[buchi] = incoming_transitions
        outgoing_transitions_all[buchi] = outgoing_transitions
        transition_mapping_all[buchi] = transition_mapping



    while queue: 
        current_distance, current_state = heapq.heappop(queue)

        # termination condition: check if we already saw the current state
        if current_state in explored:
            # if so, iterate backwards until we hit the same state
            hit_acceptance = False
            iterate = True
            iterate_state = current_state
            loop_word = ""
            while iterate:
                loop_word+=prev_character_mapping[iterate_state]
                iterate_state = prev_state_mapping[iterate_state]

                if all([iterate_state[i] in buchi_array[i].acceptance_states for i in range(len(iterate_state))]):
                    hit_acceptance = True
                    iterate = False
                else:
                    if iterate_state == current_state : iterate = False

                if iterate_state == initial_state : break

            if hit_acceptance : return word[current_state] + loop_word[::-1]

        explored.add(current_state)

        # if we've already been to this state before, no need to re-explore it. 
        # note, we already push all possible explorable states onto the stack
        if current_state in explored : continue

        # go through alphabet, construct adjacent states on the fly
        reachable = []
        for char in alphabet:

            # exact same on-the-fly construction for NFAs
            new_states, break_loop = [], False
            new_state = [None]*len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state[i], char)

                if new_state == []:
                    break_loop = True
                    break

            if break_loop : continue

            # now, we have a new_state with (array, ..., array)
            new_state = tuple(new_state)
            for new_state in list(itertools.product(*new_state)):
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""
                prev_state_mapping[new_state] = current_state
                prev_character_mapping[new_state] = char
                reachable.append((new_state, char))

        for (state, char) in reachable:
            new_distance = distance[current_state] + 1

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state))

    return None

