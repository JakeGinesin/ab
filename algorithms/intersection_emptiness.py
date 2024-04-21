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


# it seems i'll need to store captured in a dict... 
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
                    if iterate_state == current_state: 
                        iterate = False

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

"""
From here: to be refined further. 
The code below is quite messy. You have been warned.
"""


def expr_gurobi_nfa_intersection(nfa_array):
    fa = nfa_array[0].alphabet
    for i in range(1, len(nfa_array)) : assert fa == nfa_array[i].alphabet

    initial_state = tuple(nfa.initial_state for nfa in nfa_array)

    distance = {}
    word = {}

    queue = []
    heapq.heappush(queue, (0, initial_state))
    distance[initial_state] = 0
    word[initial_state] = "" 

    # initialize states to be considered in heuristic by each nfa
    heuristic_wellconnected_states = {}
    nfa_array_single_state = []

    # Ax = b
    nfa_matrix_A = {}
    nfa_matrix_b = {}

    nfa_matrix_column_mapping = {}

    transition_mapping = {}
    transition_mapping_reverse = {}

    # initialize NFAs
    for nfa in nfa_array:
        nfa_s = nfa.to_single_acceptance_state() # O(V+E)... + NFA construction time
        nfa_array_single_state.append(nfa_s)

        # this isn't the fastest solution. 
        # O(V+E) * 2 
        # should optimize with one DFS (a single O(V+E) to eliminate dead states
        nfa_s_acceptance_state = list(nfa_s.acceptance_states)[0]
        usable_states = nfa_s.reachable_from(nfa_s.initial_state)
        acceptance_usable = nfa_s.reachable_to(nfa_s_acceptance_state)
        usable_states = usable_states.intersection(acceptance_usable)

        heuristic_wellconnected_states[nfa_s] = usable_states

        # initialize the matrix
        nfa_matrix_A[nfa_s] = []
        nfa_matrix_b[nfa_s] = []
        nfa_matrix_column_mapping[nfa_s] = {}

        # map the transitions
        transition_mapping[nfa_s] = {}
        transition_mapping_reverse[nfa_s] = {}
        transition_labels = ["t" + str(i) for i in range(len(nfa_s.transitions))]

        count = 0
        for (sA, i, sB) in nfa_s.transitions:
            transition_mapping[nfa_s][transition_labels[count]] = (sA, i, sB)
            transition_mapping_reverse[nfa_s][(sA, i, sB)] = transition_labels[count]
            count+=1

        index = 0
        for state in heuristic_wellconnected_states[nfa_s]: 
            # ok, we set the column for the state
            nfa_matrix_column_mapping[nfa_s][state] = index
            
            # now we need to set the row for the transitions
            r = [0] * len(nfa_s.transitions) * 2
            for t in nfa_s.end_state_transition_mapping[state]:
                label = transition_mapping_reverse[nfa_s][t]
                num = int(label[1:])
                r[num]+=1

            for t in nfa_s.beginning_state_transition_mapping[state]:
                label = transition_mapping_reverse[nfa_s][t]
                num = int(label[1:])
                r[num + len(nfa_s.transitions)]-=1

            nfa_matrix_A[nfa_s].append(r)
            
            a = [None]
            # Ax=(b)
            # if state == nfa_s.initial_state : a[0] = -1
            if state == nfa_s.initial_state : a[0] = 0
            elif state == nfa_s_acceptance_state : a[0] = 1
            else : a[0] = 0
            nfa_matrix_b[nfa_s].append(a)
            index+=1

    nfa_arr_len = len(nfa_array_single_state)

    # state is a tuple of length len(nfa_array) all with the current state of each nfa
    # so, all we need to do is plug each state and each corresponding set of usable (connected) states into the parikh image heuristic
    def heuristic(state):
        ret_sum = 0
        for i in range(nfa_arr_len):

            nfa_h = nfa_array_single_state[i]
            initial_state, well_connected_states = state[i], heuristic_wellconnected_states[nfa_array_single_state[i]]
            
            # passing case
            if initial_state == list(nfa_h.acceptance_states)[0] : continue

            # since we don't set the initial state in the matrix, we set it now
            # (we use the mapping between states and rows of the matrix to instafind)
            A_matrix = nfa_matrix_A[nfa_h].copy()
            b_matrix = nfa_matrix_b[nfa_h].copy()

            b_matrix[nfa_matrix_column_mapping[nfa_h][initial_state]][0] = -1

            A = np.matrix(A_matrix)
            b = np.matrix(b_matrix)

            model = Model('matrix1')
            num_constrs, num_vars = A.shape

            model.setParam('OutputFlag', 0)
            x = model.addVars(num_vars)
            model.update()

            for i in range(num_constrs):
                model.addConstr(quicksum(A[i, j]*x[j] for j in range(num_vars)) == b[i])

            # Specify objective: minimize the sum of x
            model.setObjective(x.sum(), GRB.MINIMIZE)

            # Optimize model
            model.optimize()

            total_dist = 0

            if model.status == GRB.INFEASIBLE:
                return float('inf')

            elif model.status == GRB.OPTIMAL:
                for v in model.getVars():
                    ii = int(v.varName[1:])
                    value=v.x

                    # incoming transitions
                    if ii < len(nfa_h.transitions):
                        total_dist+=value

                ret_sum+=total_dist

        return ret_sum

    total = 0
    while queue:
        total+=1
        current_distance, current_state = heapq.heappop(queue)
        print(queue)

        # acceptance condition: check if all states are in acceptance etc
        if all([current_state[i] in nfa_array_single_state[i].acceptance_states for i in range(len(current_state))]): 
            print(total)
            return word[current_state]

        if current_distance > distance[current_state] : continue

        reachable = []
        for char in fa:
            new_states = []
            break_loop = False

            new_state = [None]*nfa_arr_len
            for i in range(nfa_arr_len):
                new_state[i] = nfa_array_single_state[i].transition_function(current_state[i], char)

                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))
            for new_state in l:
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""
                reachable.append((new_state, char))

        for (state, char) in reachable:
            # where state is an array of the current states of the nfas 
            # in nfa_array_single_state
            # associated by index
            new_distance = distance[current_state] + heuristic(state)

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state))

    return None

def fast_gurobi_nfa_intersection(nfa_array):
    # 1. assert all alphabets are the same
    fa = nfa_array[0].alphabet
    for i in range(1, len(nfa_array)) : assert fa == nfa_array[i].alphabet
    # 2. create initial state for on-the-fly stuff
    initial_state = tuple(nfa.initial_state for nfa in nfa_array)

    # 3. run A*
    distance = {}
    word = {}

    queue = []
    heapq.heappush(queue, (0, initial_state))
    distance[initial_state] = 0
    word[initial_state] = ""

    # initialize states to be considered in heuristic by each nfa
    heuristic_wellconnected_states = {}
    nfa_array_single_state = []

    nfa_matrix_A = {}
    nfa_matrix_b = {}
    gurobi_models = {}
    nfa_matrix_x = {}

    # transition_mapping_all = {}
    # transition_mapping_reverse_all = {}

    state_mapping_all = {}
    # state_mapping_reverse_all = {}

    nfa_s_array = []

    for nfa_s in nfa_array:
        nfa = nfa_s.to_single_acceptance_state() # O(V+E) + NFA construction time (linear)
        nfa_s_array.append(nfa)

        nfa_acceptance_state = list(nfa.acceptance_states)[0]
        usable_states = nfa.reachable_from(nfa.initial_state)
        acceptance_usable = nfa.reachable_to(nfa_acceptance_state)
        usable_states = usable_states.intersection(acceptance_usable)

        heuristic_wellconnected_states[nfa] = usable_states

        # initialize matrices
        # nfa_matrix_A[nfa] = np.zeros( (len(usable_states), len(nfa.transitions)) )
        # nfa_matrix_b[nfa] = np.zeros( (len(usable_states), 1) )
        A_init = np.zeros( (len(usable_states), len(nfa.transitions)) )
        b_init = np.zeros( (len(usable_states), 1) )

        state_mapping = {}
        state_mapping_reverse = {}
        # transition_mapping = {}
        # transition_mapping_reverse = {}

        sc = 0
        for state in usable_states:
            state_mapping[state] = sc
            # state_mapping_reverse[sc] = state
            sc+=1

        tc = 0
        for (sA, i, sB) in nfa.transitions:
            if not sA in usable_states or sB in usable_states : continue
            # transition_mapping[(sA, i, sB)] = tc
            # transition_mapping_reverse[tc] = (sA, i, sB)

            A_init[state_mapping[sA]][tc]-=1
            A_init[state_mapping[sB]][tc]+=1

            tc+=1

        b_init[state_mapping[nfa_acceptance_state]] = 1

        # re-map the dicts and stuff
        nfa_matrix_A[nfa] = A_init
        nfa_matrix_b[nfa] = b_init

        # transition_mapping_all[nfa] = transition_mapping
        # transition_mapping_reverse_all[nfa] = transition_mapping_reverse

        state_mapping_all[nfa] = state_mapping

        # gurobi initialization
        model = Model('matrix-model')
        num_constrs, num_vars = A_init.shape

        model.setParam('OutputFlag', 0)
        x = model.addMVar(num_vars)
        model.setObjective(x.sum(), GRB.MINIMIZE)

        gurobi_models[nfa] = model
        nfa_matrix_x[nfa] = x

    def heuristic(state):
        ret_sum = 0
        for i in range(len(nfa_s_array)):
            nfa = nfa_s_array[i]
            initial_state, well_connected_states = state[i], heuristic_wellconnected_states[nfa_s_array[i]]
            nfa_acceptance_state = list(nfa.acceptance_states)[0]

            if initial_state == nfa_acceptance_state : continue

            model = gurobi_models[nfa]
            A = nfa_matrix_A[nfa]
            b = nfa_matrix_b[nfa].copy()
            b[state_mapping_all[nfa][initial_state]] = -1
            x = nfa_matrix_x[nfa]

            for c in model.getConstrs():
                model.remove(c)

            model.addMConstr(A, x, '=', b)

            model.optimize()

            if model.status == GRB.OPTIMAL: 
                ret_sum+=sum(x.X)
            elif model.status == GRB.INFEASIBLE : return float('inf')

        return ret_sum

    nfa_arr_len = len(nfa_s_array)

    total = 0
    while queue:
        total+=1
        current_distance, current_state = heapq.heappop(queue)

        # acceptance condition: check if all states are in acceptance etc
        if all([current_state[i] in nfa_s_array[i].acceptance_states for i in range(len(current_state))]): 
            print(total)
            return word[current_state]

        if current_distance > distance[current_state] : continue

        reachable = []
        for char in fa:
            new_states = []
            break_loop = False

            new_state = [None]*nfa_arr_len
            for i in range(nfa_arr_len):
                new_state[i] = nfa_s_array[i].transition_function(current_state[i], char)

                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))
            for new_state in l:
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""
                reachable.append((new_state, char))

        for (state, char) in reachable:
            # where state is an array of the current states of the nfas 
            # in nfa_array_single_state
            # associated by index
            
            # A* 
            # new_distance = distance[current_state] + 1 + heuristic(state)

            # best first search
            new_distance = heuristic(state)

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state))

    return None

def slow_dfs_buchi_intersection(buchi_array):

    # assert everything is a buchi automaton
    for buchi in buchi_array: assert isinstance(buchi, Buchi)

    # assert all buchi automata have the same alphabet
    fa = buchi_array[0].alphabet
    for buchi in buchi_array: assert fa == buchi.alphabet

    # create initial state for on the fly construction
    # states consist of (state info, states visited so far, acceptance states visited so far, trail)
    initial_state = tuple([tuple(buchi.initial_state for buchi in buchi_array), set(), set(), []])

    stack = [initial_state]
    discovered = set()

    buchi_arr_len = len(buchi_array)

    while stack:
        # print(stack)
        current_state, prev_states, acc_states, trail = stack.pop()
        discovered.add(current_state)

        # print(discovered)
        consider_break = acc_states != set()

        for char in fa:
            new_states = []
            break_loop = False

            new_state = [None]*buchi_arr_len
            for i in range(buchi_arr_len):
                new_state[i] = buchi_array[i].transition_function(current_state[i], char)

                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l: 
                # print(new_state)
                # takes advantage of fast termination
                if (consider_break and new_state in prev_states) or (new_state in acc_states and current_state == new_state):
                    trail.append((current_state, char))
                    trail.append((new_state, "*"))

                    final_state = new_state
                    in_cycle = False
                    word = ""
                    for state, char in trail:
                        if in_cycle:
                            if state == final_state : continue
                            word+=char
                        elif state == final_state: 
                            in_cycle = True
                            word+="("
                            word+=char
                        else:
                            word+=char

                    # indiciative of a cycle
                    word+=")*"
                    return word
                else:
                    if new_state in discovered : continue

                current_state_new = new_state
                prev_states_new = prev_states.copy()
                prev_states_new.add(current_state)
                acc_states_new = acc_states.copy()
                trail_new = trail.copy()
                trail_new.append((current_state, char))

                # hoping there's optimization here in the interpreter lmao
                if all([new_state[i] in buchi_array[i].acceptance_states for i in range(buchi_arr_len)]):
                    acc_states_new.add(new_state)

                stack_addition = (current_state_new, prev_states_new, acc_states_new, trail_new)

                stack.append(stack_addition)

    return None

def fast_dfs_buchi_intersection(buchi_array):

    fa = buchi_array[0].alphabet
    for buchi in buchi_array:
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

    initial_pairing = tuple([buchi.initial_state for buchi in buchi_array])
    cstr_state_mapping, rev_cstr_state_mapping = {}, {}
    cstr_state_mapping[initial_pairing] = 0 # initialize initial state as "0" (as in, state 0)
    rev_cstr_state_mapping[0] = initial_pairing

    word = {}
    initial_state = tuple([0, set(), set()])

    stack = [initial_state]
    acceptance_state_leading = {}
    word[0] = ""

    while stack:
        current_state, prev_states, acc_states = stack.pop() 

        # termination condition
        if (current_state in prev_states and
            (current_state in acc_states 
             or 
             any([ current_state in acceptance_state_leading[astate] for astate in acc_states])
             )
            ):

            word = word[current_state]

            print("counterexample: " + str(word))
            total = len(cstr_state_mapping)
            
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])
            print("percentage of states visited: ", round(total/total_possible*100, 3), "% -",str(total), "/", str(total_possible))

            return word

        elif current_state in prev_states : continue

        # skip condition : if our current distance is longer than established distance
        current_state_ite = rev_cstr_state_mapping[current_state]

        reachable = []
        for char in fa:
            new_states, break_loop = [], False

            new_state = [None] * len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state_ite[i], char)

                # if this is ever the case, we know the transition function led to nothing
                # thus this state is invalid
                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l:

                # if this is the case, we know we've already seen it before
                new_state_label = None
                if not new_state in cstr_state_mapping:
                    # print("new state!")
                    new_state_label = len(cstr_state_mapping)

                    cstr_state_mapping[new_state] = new_state_label
                    rev_cstr_state_mapping[new_state_label] = new_state

                    word[new_state_label] = ""

                else: 
                    new_state_label = cstr_state_mapping[new_state]

                # current_state_new = new_state
                prev_states_new = prev_states.copy()
                prev_states_new.add(current_state)
                acc_states_new = acc_states.copy()

                all_accepting = True
                for i in range(len(buchi_array)):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False

                if all_accepting:
                    acc_states_new.add(new_state_label)
                    acceptance_state_leading[new_state_label] = prev_states_new.copy()

                """
                if all([new_state[i] in buchi_array[i].acceptance_states for i in range(len(buchi_array))]):
                    acc_states_new.add(new_state)
                """

                stack_addition = (new_state_label, prev_states_new, acc_states_new)

                stack.append(stack_addition)
                word[new_state_label] = word[current_state] + char

    return None


def slow_djikstra_buchi_intersection(buchi_array):
    """
    We consider the following termination conditions:
    1, we run into a previously visited node (i.e., we find a cycle)
    2, we have previously visited an accepting state
    3, the node we ran into leads to an accepting state previously visited
    """

    # assert everything is a buchi automaton
    for buchi in buchi_array: assert isinstance(buchi, Buchi)

    # assert all buchi automata have the same alphabet
    fa = buchi_array[0].alphabet
    for buchi in buchi_array: assert fa == buchi.alphabet

    # create initial state for on the fly construction
    # states consist of (state info, states visited so far, acceptance states visited so far, trail)
    initial_state = tuple([tuple(buchi.initial_state for buchi in buchi_array), set(), set(), []])

    # initialize distance and word
    distance = {}
    word = {}

    # put first state onto the queue
    queue = []
    heapq.heappush(queue, (0, initial_state))
    distance[initial_state[0]] = 0
    word[initial_state[0]] = ""

    discovered = set()
    acceptance_state_leading = {}

    buchi_arr_len = len(buchi_array)

    total = 0
    while queue:
        total+=1

        # pop from heap
        current_distance, current_state_tuple = heapq.heappop(queue)

        # and get relevant things out of the tuple
        current_state, prev_states, acc_states, trail = current_state_tuple
        print(current_state, discovered)

        discovered.add(current_state)

        if current_distance > distance[current_state] : continue 

        # assemble all reachable states from the current
        reachable = []
        for char in fa:
            new_states = []
            break_loop = False

            new_state = [None]*buchi_arr_len
            for i in range(buchi_arr_len):
                new_state[i] = buchi_array[i].transition_function(current_state[i], char)

                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            # loop through all possible new states
            for new_state in l:

                # termination conditions (see 1,2,3 at the top of the function)
                if new_state in discovered and acc_states != set() and any([ new_state in acceptance_state_leading[astate] for astate in acc_states]):
                    trail.append((current_state, char))
                    trail.append((new_state, "*"))

                    final_state = new_state
                    in_cycle = False
                    word = ""

                    for state, char in trail:
                        if in_cycle:
                            if state == final_state : continue
                            word+=char
                        elif state == final_state: 
                            in_cycle = True
                            word+="("
                            word+=char
                        else:
                            word+=char

                    # indiciative of a cycle
                    print(total)
                    word+=")*"
                    return word
 
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""

                # create state tuple
                # that is, add to prev states, acc states, trail 
                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()
                # can be optimized... the for loop array is generated before
                if all([ new_state[i] in buchi_array[i].acceptance_states for i in range(buchi_arr_len) ]):
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state)
                    acceptance_state_leading[new_state] = new_prev_states.copy()

                trail_new = trail.copy()
                trail_new.append((current_state, char))

                t = (new_state, new_prev_states, new_acc_states, trail_new)

                reachable.append((t, char))

        for (state_t, char) in reachable:
            state = state_t[0]

            new_distance = distance[current_state] + 1

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state_t))

    return None

def djikstra_buchi_intersection(buchi_array):
    """
    idea: we use a slicker termination condition here to make this faster
    """

    fa = buchi_array[0].alphabet
    for buchi in buchi_array: 
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

    initial_state = tuple([tuple(buchi.initial_state for buchi in buchi_array), set(), set(), []])

    distance = {}
    word = {}

    queue = []
    heapq.heappush(queue, (0, initial_state))
    distance[initial_state[0]] = 0
    word[initial_state[0]] = ""

    discovered = set()
    acceptance_state_leading = {}

    buchi_arr_len = len(buchi_array)

    # keep track of total number of states visited
    total = 0
    while queue:
        total+=1

        current_distance, current_state_tuple = heapq.heappop(queue)

        current_state, prev_states, acc_states, trail = current_state_tuple 

        # the excellent termination condition
        if current_state in prev_states and (current_state in acc_states or any([ current_state in acceptance_state_leading[astate] for astate in acc_states])):
            trail.append((current_state, "*"))

            final_state = current_state
            in_cycle = False
            word = ""

            for state, char in trail:
                if in_cycle:
                    if state == final_state : continue
                    word+=char
                elif state == final_state: 
                    in_cycle = True
                    word+="("
                    word+=char
                else:
                    word+=char

            # indiciative of a cycle
            print(total)
            
            # print total possible states across all buchi for comparison
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])

            # print the percentage of states visited, as a percent, rounded off at 3 decimal places
            print("percentage of states visited: ", round(total/total_possible*100, 3), "%")

            word+=")*"
            return word

        discovered.add(current_state)

        if current_distance > distance[current_state] : continue

        reachable = []
        for char in fa:
            new_states = []
            break_loop = False

            new_state = [None]*buchi_arr_len
            for i in range(buchi_arr_len):
                new_state[i] = buchi_array[i].transition_function(current_state[i], char)
                # reminder that new_state[i] is an array of states

                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            buchi_arr_len = len(buchi_array)
            l = list(itertools.product(*new_state))

            # loop through all possible new states
            for new_state in l:
 
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""

                # create state tuple
                # that is, add to prev states, acc states, trail 
                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()
                # can be optimized... the for loop array is generated before

                all_accepting = True
                for i in range(buchi_arr_len):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False
                        break

                if all_accepting:
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state)
                    acceptance_state_leading[new_state] = new_prev_states.copy()

                trail_new = trail.copy()
                trail_new.append((current_state, char))

                t = (new_state, new_prev_states, new_acc_states, trail_new)

                reachable.append((t, char))

        for (state_t, char) in reachable:
            state = state_t[0]

            new_distance = distance[current_state] + 1

            distance[state] = new_distance
            word[state] = word[current_state] + char
            heapq.heappush(queue, (new_distance, state_t))

    return None

def fast_djikstra_buchi_intersection(buchi_array):
    """
    idea: we use,
        - slicker termination condition
        - state mappings 
    """

    fa = buchi_array[0].alphabet
    for buchi in buchi_array: 
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

    initial_pairing = tuple([buchi.initial_state for buchi in buchi_array])
    cstr_state_mapping, rev_cstr_state_mapping = {}, {}
    cstr_state_mapping[initial_pairing] = 0 # initialize initial state as "0" (as in, state 0)
    rev_cstr_state_mapping[0] = initial_pairing

    initial_state = tuple([0, set(), set(), []])

    distance, word, queue = {}, {}, []
    
    heapq.heappush(queue, (0, initial_state))
    # distance[0], word[0] = 0, ""
    distance[0] = 0

    acceptance_state_leading = {}

    while queue:
        current_distance, current_state_tuple = heapq.heappop(queue)
        current_state, prev_states, acc_states, trail = current_state_tuple 

        # termination condition
        if (current_state in prev_states and
            (current_state in acc_states 
             or 
             any([ current_state in acceptance_state_leading[astate] for astate in acc_states])
             )
            ):

            # print(word[current_state])
            trail.append((current_state, "*"))

            final_state = current_state
            in_cycle = False
            word = ""

            for state, char in trail:
                if in_cycle:
                    if state == final_state : continue
                    word+=char
                elif state == final_state: 
                    in_cycle = True
                    word+="("
                    word+=char
                else:
                    word+=char

            # indiciative of a cycle
            total = len(cstr_state_mapping)
            print(total)
            
            # print total possible states across all buchi for comparison
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])

            # print the percentage of states visited, as a percent, rounded off at 3 decimal places
            print("percentage of states visited: ", round(total/total_possible*100, 3), "%")

            word+=")*"
            return word
            """
            print("termination condition!")
            print(total)
            print(word[current_state])
            return 0
            """
        elif current_state in prev_states : continue

        # skip condition : if our current distance is longer than established distance
        if current_distance > distance[current_state] : continue
        
        current_state_ite = rev_cstr_state_mapping[current_state]

        # assemble new states
        reachable = []
        for char in fa:
            new_states, break_loop = [], False

            new_state = [None] * len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state_ite[i], char)

                # if this is ever the case, we know the transition function led to nothing
                # thus this state is invalid
                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l:

                # if this is the case, we know we've already seen it before
                new_state_label = None
                if not new_state in cstr_state_mapping:
                    # print("new state!")
                    new_state_label = len(cstr_state_mapping)

                    cstr_state_mapping[new_state] = new_state_label
                    rev_cstr_state_mapping[new_state_label] = new_state

                    distance[new_state_label] = float('inf')
                    word[new_state_label] = ""

                else: 
                    new_state_label = cstr_state_mapping[new_state]

                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()

                all_accepting = True
                for i in range(len(buchi_array)):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False
                        break

                if all_accepting:
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state_label)
                    acceptance_state_leading[new_state_label] = new_prev_states.copy()


                trail_new = trail.copy()
                trail_new.append((current_state, char))

                t = (new_state_label, new_prev_states, new_acc_states, trail_new)


                reachable.append((t, new_state, char))

        for (state_t, state, char) in reachable:
            state_label = state_t[0]

            new_distance = distance[current_state] + 1

            distance[state_label] = new_distance
            # word[state_label] = word[current_state] + char
            heapq.heappush(queue, (new_distance, state_t))

    return None

def fast_djikstra_buchi_intersection_notrail(buchi_array):
    """
    idea: we use,
        - slicker termination condition
        - state mappings 
    """

    fa = buchi_array[0].alphabet
    for buchi in buchi_array: 
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

    initial_pairing = tuple([buchi.initial_state for buchi in buchi_array])
    cstr_state_mapping, rev_cstr_state_mapping = {}, {}
    cstr_state_mapping[initial_pairing] = 0 # initialize initial state as "0" (as in, state 0)
    rev_cstr_state_mapping[0] = initial_pairing

    initial_state = tuple([0, set(), set()])

    distance, word, queue = {}, {}, []
    
    heapq.heappush(queue, (0, initial_state))
    distance[0], word[0] = 0, ""

    acceptance_state_leading = {}

    while queue:
        current_distance, current_state_tuple = heapq.heappop(queue)
        current_state, prev_states, acc_states = current_state_tuple 

        # termination condition
        if (current_state in prev_states and
            (current_state in acc_states 
             or 
             any([ current_state in acceptance_state_leading[astate] for astate in acc_states])
             )
            ):

            word = word[current_state]

            print("counterexample: " + str(word))
            total = len(cstr_state_mapping)
            
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])
            print("percentage of states visited: ", round(total/total_possible*100, 3), "% -",str(total), "/", str(total_possible))

            return word

        elif current_state in prev_states : continue

        # skip condition : if our current distance is longer than established distance
        if current_distance > distance[current_state] : continue
        
        current_state_ite = rev_cstr_state_mapping[current_state]

        # assemble new states
        reachable = []
        for char in fa:
            new_states, break_loop = [], False

            new_state = [None] * len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state_ite[i], char)

                # if this is ever the case, we know the transition function led to nothing
                # thus this state is invalid
                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l:

                # if this is the case, we know we've already seen it before
                new_state_label = None
                if not new_state in cstr_state_mapping:
                    # print("new state!")
                    new_state_label = len(cstr_state_mapping)

                    cstr_state_mapping[new_state] = new_state_label
                    rev_cstr_state_mapping[new_state_label] = new_state

                    distance[new_state_label] = float('inf')
                    word[new_state_label] = ""

                else: 
                    new_state_label = cstr_state_mapping[new_state]

                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()

                all_accepting = True
                for i in range(len(buchi_array)):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False
                        break

                if all_accepting:
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state_label)
                    acceptance_state_leading[new_state_label] = new_prev_states.copy()

                t = (new_state_label, new_prev_states, new_acc_states)

                reachable.append((t, new_state, char))

        for (state_t, state, char) in reachable:
            state_label = state_t[0]

            new_distance = distance[current_state] + 1

            distance[state_label] = new_distance
            word[state_label] = word[current_state] + char
            heapq.heappush(queue, (new_distance, state_t))

    return None


"""
the idea behind making a matrix ver is to think about 
integer disjunction. e.g., we have x >= 0, y >= 0,x+y=1
"""
def expr_z3_heuristic_buchi_intersection(buchi_array):
    """
    the heuristic here really has three states:
    1. before hitting the acceptance state
    2. when you hit the acceptance state
    3. after you hit the acceptance state

    based on these cases, the Ax = b matrix we're evaluating for the heuristic is different.

    however, they're all shared in the same goal: 
    we want to *underapproximate* the distance to the successful termination condition

    what is the termination condition, you may ask?
    1, we run into a previously visited node (i.e., we find a cycle)
    2, we have previously visited an accepting state
    3, the node we ran into leads to an accepting state previously visited
    """

    # assert eveything is a buchi automata
    for buchi in buchi_array: assert isinstance(buchi, Buchi)

    # assert all buchi automata have the same alphabet
    fa = buchi_array[0].alphabet
    for buchi in buchi_array: assert fa == buchi.alphabet

    initial_state = tuple([tuple(buchi.initial_state for buchi in buchi_array), set(), set(), []])

    distance = {}
    word = {}

    queue = []
    heapq.heappush(queue, (0, initial_state))
    distance[initial_state[0]] = 0
    word[initial_state[0]] = ""

    discovered = set()
    acceptance_state_leading = {}

    buchi_arr_len = len(buchi_array)

    # heurisitc required
    hitAcceptingStatesDict = {}

    # 1 : before hitting acceptance
    # 2 : in acceptance
    # 3 : after hitting acceptance
    heuristicBuchiState = {}

    transition_mapping_all = {}
    transition_mapping_reverse_all = {}

    incoming_transitions_all = {}
    outgoing_transitions_all = {}

    heuristic_wellconnected_states = {}

    processed_history = {}
    ret_sum_history = {}

    for buchi in buchi_array:
        """
        ok, here are the initialization things I need to do.
        1. get states whose states will be used to make the constraints via intersection
        2. create the incoming and outgoing transitions per state in each buchi
        3. create the transition mapping
        """

        hitAcceptingStatesDict[buchi] = set()

        """
        1. 
        """

        # yeah.. if the model is a promela model we can just do all states reachable from the initial state
        # otherwise we have to do intersection
        # can add a case if all the states are accepting (e.g., promela)
        """
        reachable_from_init = buchi.reachable_from(buchi.initial_state)
        to_accepting = set()
        for accepting in buchi.acceptance_states:
            to_accepting = to_accepting.union(buchi.reachable_to(accepting))
        usable_states = reachable_from_init.intersection(to_accepting)

        heuristic_wellconnected_states[buchi] = usable_states
        """
        usable_states = buchi.get_wellconnected_states()
        heuristic_wellconnected_states[buchi] = usable_states

        """
        2. 
        """

        incoming_transitions = {}
        outgoing_transitions = {}

        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        ts = []
        for (sA, i, sB) in buchi.transitions: 
            if sA in usable_states and sB in usable_states:
                ts.append((sA, i, sB))

        """
        3. 
        """

        transition_mapping = {}
        transition_mapping_reverse = {}
        transition_labels = ["t" + str(i) for i in range(len(ts))]

        count = 0
        for (sA, i, sB) in ts:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            transition_mapping_reverse[(sA, i, sB)] = transition_labels[count]

            outgoing_transitions[sA].add("t" + str(count))
            incoming_transitions[sB].add("t" + str(count))

            count = count + 1

        incoming_transitions_all[buchi] = incoming_transitions
        outgoing_transitions_all[buchi] = outgoing_transitions
        transition_mapping_all[buchi] = transition_mapping
        outgoing_transitions_all[buchi] = outgoing_transitions

        processed_history[buchi] = {}

        # do stuff... is there a linear time procedure to like, make it an accepting cycle thing

    def heuristic(state):
        ret_sum = 0
        buchi_current_states, prev_states, acc_states, trail = state
        for i in range(buchi_arr_len):

            buchi = buchi_array[i]  
            current_state = buchi_current_states[i]
            init_state = current_state

            incoming_transitions = incoming_transitions_all[buchi]
            outgoing_transitions = outgoing_transitions_all[buchi]

            # need to determine if we've visited any accepting states for the current "i" in buchi...
            # because we need to determine how we look for cycles
            hit_accepting_states = hitAcceptingStatesDict[buchi]

            hist = processed_history[buchi]
            if (init_state, tuple(hit_accepting_states)) in set(hist.keys()): 
                ret_sum += hist[(init_state, tuple(hit_accepting_states))]
                continue

            # no... we need to do things by individual acceptance state
            # it's clear we need to do like, a classification of the "state" for each individual acceptance state within the buchi.
            # especially when we're taking the minimum.
            """
            how taking the "or" over everything accepting is going to work.

            between the same initial input state:
            
            invariants: 
            - the regular flow constraints
            - the initial state flow constraints

            variants:
            - the flow constraints of the accepting

            ok, so we now hit this part of the function. what happens now?
            i guess we can use either gurobi or z3.

            this first try will be with z3 with smt formatting.

            steps:
            1. initialize the transition labels etc (can do that in pre-processing)
            2. load up z3 with init state and regular state flow constraints (the invariant)
            3. loop through all accepting states
                3a. determine, based on whether we've seen that accepting state,
                    if we should classify it as "1, 2, or 3"
                3b. based on the above, add the linear equation to an or constraint array
            4. solve, get the minimal, return the heuristic value
            """

            """
            2. 
            """
            
            # print(hit_accepting_states)

            transition_mapping = transition_mapping_all[buchi]
            usable_states = heuristic_wellconnected_states[buchi]


            s = Optimize()
            tosolve = {}
            for t in list(transition_mapping.keys()):
                tosolve[t] = Int(str(t))

            for t in tosolve.values():
                s.add(t >= 0)

            # minimize solution
            s.minimize(Sum(list(tosolve.values())))

            accepting_states_to_solve = set()

            for state in usable_states:
                if state in buchi.acceptance_states:
                    accepting_states_to_solve.add(state)
                elif state == init_state:
                    continue
                else: 
                    s.add(Sum([tosolve[t] for t in incoming_transitions[state]]) == Sum([tosolve[t] for t in outgoing_transitions[state]]))

            """
            ok, for each accepting we need:
            1. before hitting the acceptance state
            2. when you hit the acceptance state
            3. after you hit the acceptance state
            """

            # print("ESTABLISHED" in *hit_accepting_states)

            accept_or = []
            for accepting_state in accepting_states_to_solve:
                # here we determine the state...
                if accepting_state == init_state:
                    print("equal")
                    # print("case 2")
                    # case 2
                    # in this case, the accepting/init state should have one input, one output

                    accept_or.append(
                            And(
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]]) == 1),
                            (Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) == 1)
                            ))
                elif accepting_state in list({item for tup in hit_accepting_states for item in tup}):
                    print("in")
                    # print("case 3")
                    # case 3
                    # we've already seen the accepting state (we're in the cycle)

                    # init state and acceptance state

                    accept_or.append(
                            And(
                            (Sum([tosolve[t] for t in incoming_transitions[init_state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[init_state]])),
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]])  == Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) + 1)
                            ))
                else:
                    print("out")
                    # print("case 1")
                    # case 1
                    # we're not in the cycle.. we're doing the on-the-way path
                    # 3 constraints here: differential, out>1, int>2

                    # AND VER:
                    accept_or.append(
                            And(
                            (Sum([tosolve[t] for t in incoming_transitions[init_state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[init_state]])),
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]]) == Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) + 1),
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]]) == 2),
                            (Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) == 1)
                            ))

            # one of the big problems with this is that there will be just so many disjunctions in a Kripkie automata, who has all states set to accept ._.
            # print(*accept_or)
            s.add(Or(*accept_or)) 

            if s.check() == sat:
                m = s.model()
                solns = {var: m[var] for var in tosolve.values()}

                total_dist = 0
                for transition, value in solns.items():
                    # sA, i, sB = transition_mapping[str(transition)]
                    total_dist+=int(str(value)) # kinda slow, can be optimized 

                hist[(init_state, tuple(hit_accepting_states))] = total_dist


                ret_sum += total_dist

            else : return float('inf')

        return ret_sum

    # invoke queue
    total = 0
    while queue:
        total+=1

        current_distance, current_state_tuple = heapq.heappop(queue)

        current_state, prev_states, acc_states, trail = current_state_tuple
        discovered.add(current_state)

        # check if distance[current_state] has a value aleady.... if so, prev_states becomes the union

        if current_distance > distance[current_state] : continue 

        reachable = []
        for char in fa:
            new_states = []
            break_loop = False

            new_state = [None]*buchi_arr_len
            for i in range(buchi_arr_len):
                new_state[i] = buchi_array[i].transition_function(current_state[i], char)

                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            buchi_arr_len = len(buchi_array)
            l = list(itertools.product(*new_state))

            for new_state in l:

                if new_state in discovered and acc_states != set() and any([ new_state in acceptance_state_leading[astate] for astate in acc_states]):
                    trail.append((current_state, char))
                    trail.append((new_state, "*"))

                    final_state = new_state
                    in_cycle = False
                    word = ""

                    for state, char in trail:
                        if in_cycle:
                            if state == final_state : continue
                            word+=char
                        elif state == final_state: 
                            in_cycle = True
                            word+="("
                            word+=char
                        else:
                            word+=char

                    # indiciative of a cycle
                    word+=")*"
                    print(total)
                    return word
 
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""

                # that is, add to prev states, acc states, trail 
                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()
                # can be optimized... the for loop array is generated before
                # no, I *need* to optimize this
                addState = True
                for i in range(buchi_arr_len):
                    if new_state[i] in buchi_array[i].acceptance_states:
                        # do i add the new state?
                        hitAcceptingStatesDict[buchi_array[i]].add(new_state)

                    else:
                        addState = False

                if addState:
                    new_acc_states.add(new_state)
                    acceptance_state_leading[new_state] = new_prev_states.copy()

                """
                if all([ new_state[i] in buchi_array[i].acceptance_states for i in range(buchi_arr_len) ]):
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state)
                    acceptance_state_leading[new_state] = new_prev_states.copy()
                """

                trail_new = trail.copy()
                trail_new.append((current_state, char))

                t = (new_state, new_prev_states, new_acc_states, trail_new)

                reachable.append((t, char))

        for (state_t, char) in reachable:
            state = state_t[0]

            hval = heuristic(state_t)
            # print(hval, state)
            
            # print(hval)
            new_distance = distance[current_state] + hval + 1
            # new_distance = distance[current_state] + 1
            # new_distance = heuristic(state_t) 

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state_t))

    return None

def gurobi_heuristic_buchi_intersection(buchi_array):
    # assert eveything is a buchi automata
    for buchi in buchi_array: assert isinstance(buchi, Buchi)

    # assert all buchi automata have the same alphabet
    fa = buchi_array[0].alphabet
    for buchi in buchi_array: assert fa == buchi.alphabet

    initial_state = tuple([tuple(buchi.initial_state for buchi in buchi_array), set(), set(), []])

    distance = {}
    word = {}

    queue = []
    heapq.heappush(queue, (0, initial_state))
    distance[initial_state[0]] = 0
    word[initial_state[0]] = ""

    buchi_arr_len = len(buchi_array)

    discovered = set()

    # lots of dicts of dicts
    transition_mapping_all = {}
    transition_mapping_reverse_all = {}

    A_init_all = {}
    b_init_all = {}
    gurobi_model_all = {}
    x_all = {}

    # heurisitc required
    hitAcceptingStatesDict = {}

    # 1 : before hitting acceptance
    # 2 : in acceptance
    # 3 : after hitting acceptance
    heuristicBuchiState = {}

    acceptance_altered_matrices_A = {}
    acceptance_altered_matrices_b = {}
    acceptance_states_usable_all = {}

    incoming_transitions_all = {}
    outgoing_transitions_all = {}
    
    heuristic_wellconnected_states = {}
    state_mapping_all = {}

    acceptance_state_leading = {}

    processed_history = {}
    ret_sum_history = {}

    # initialization for each buchi 
    for buchi in buchi_array:
        """
        ok! for each buchi automata, what do we need to pre-process?
        1. the set of reachable states
        2. create incoming, outgoing, transitions, create their maps to the matrix values
        3. the matrices representing the flow constraints for the non-critical states
        4. the gurobi instances for each buchi automata
        5. doctor the matrices so, for each accepting state, we have a matrix with that acceptance state value removed
        """

        hitAcceptingStatesDict[buchi] = set()

        """
        1. 
        """

        """
        reachable_from_init = buchi.reachable_from(buchi.initial_state)
        to_accepting = set()
        for accepting in tqdm(buchi.acceptance_states):
            to_accepting = to_accepting.union(buchi.reachable_to(accepting))
        usable_states2 = reachable_from_init.intersection(to_accepting)
        """

        # heuristic_wellconnected_states[buchi] = usable_states
        usable_states = buchi.get_wellconnected_states()
        # print(usable_states2.difference(usable_states))
        heuristic_wellconnected_states[buchi] = usable_states

        """
        2. 
        """

        incoming_transitions = {}
        outgoing_transitions = {}

        # we only want to consider the acceptance states reachable from the init state
        acceptance_states_usable = set()

        state_mapping = {}
        state_mapping_reverse = {}

        sc = 0 
        for state in usable_states:
            if state in buchi.acceptance_states : acceptance_states_usable.add(state)
            # elif state != init_state : nonspecial_states.add(state)

            sc+=1

            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        ts = []
        # tsb = []
        for transition in buchi.transitions:
            # we know all items in tsb are in ts
            if transition[0] in usable_states and transition[2] in usable_states : ts.append(transition)
            # if (transition[0] in usable_states and transition[0] != init_state and transition[0] not in acceptance_states_usable
            #    and transition[2] in usable_states and transition[2] != init_state and transition[2] not in acceptance_states_usable):
            #    tsb.append(transition)

        transition_mapping = {}
        transition_mapping_reverse = {}
        transition_labels = ["t" + str(i) for i in range(len(ts))]


        count = 0
        for (sA, i, sB) in ts:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            transition_mapping_reverse[(sA, i, sB)] = transition_labels[count]

            outgoing_transitions[sA].add("t" + str(count))
            incoming_transitions[sB].add("t" + str(count))

            count = count + 1

        A_init = np.zeros( (len(usable_states), len(ts)) ) 
        b_init = np.zeros( (len(usable_states), 1) ) 

        transition_index_mapping = {}
        transition_index_mapping_reverse = {}

        tc = 0 # with the counter, we're going across the matrix from left to right and filling stuff in
        for (sA, i, sB) in ts:
            # we avoid setting flow constraints if a state is initial or accepting
            #if sA == init_state or sB == init_state: 
            #    continue
            
            #if sA in acceptance_states_usable or sB in acceptance_states_usable: 
            #    continue

            transition_index_mapping[tc] = (sA, i, sB)
            transition_index_mapping_reverse[(sA, i, sB)] = tc

            """
            # outgoing transitions
            A_init[state_mapping[sA]][tc]-=1

            # incoming transitions
            A_init[state_mapping[sB]][tc]+=1
            """

            tc+=1

        sc = 0
        for state in usable_states:

            state_mapping[state] = sc
            state_mapping_reverse[sc] = state
        
            # if state in acceptance_states_usable or state == buchi.initial_state : continue
            if state in acceptance_states_usable : continue



            inputs = incoming_transitions[state]
            output = outgoing_transitions[state]

            for i in inputs:
                a = transition_mapping[i]
                r = transition_index_mapping_reverse[a]
                A_init[sc][r]+=1

            for i in output:
                a = transition_mapping[i]
                r = transition_index_mapping_reverse[a]
                A_init[sc][r]-=1

            sc+=1

        model = Model('GAMING')
        num_constrs_m, num_vars = A_init.shape
        # print(A_init)

        model.setParam('OutputFlag', 0)

        x = model.addMVar(len(ts), vtype=GRB.INTEGER)


        model.setObjective(quicksum(x), GRB.MINIMIZE)

        transition_mapping_all[buchi] = transition_mapping
        transition_mapping_reverse_all[buchi] = transition_mapping_reverse

        acceptance_states_usable_all[buchi] = acceptance_states_usable
        incoming_transitions_all[buchi] = incoming_transitions
        outgoing_transitions_all[buchi] = outgoing_transitions
        state_mapping_all[buchi] = state_mapping

        A_init_all[buchi] = A_init
        b_init_all[buchi] = b_init
        gurobi_model_all[buchi] = model 
        x_all[buchi] = x

        processed_history[buchi] = {}

    def heuristic(state):
        ret_sum = 0
        buchi_current_states, prev_states, acc_states, trail = state

        if (tuple(buchi_current_states), tuple(acc_states)) in set(ret_sum_history.keys()):
            return ret_sum_history[(tuple(buchi_current_states), tuple(acc_states))]

        for i in range(buchi_arr_len):

            buchi = buchi_array[i]
            init_state = buchi_current_states[i]
            hit_accepting_states = hitAcceptingStatesDict[buchi]

            transition_mapping = transition_mapping_all[buchi]
            usable_states = heuristic_wellconnected_states[buchi]
            A_init = A_init_all[buchi]
            b_init = b_init_all[buchi]
            model = gurobi_model_all[buchi]
            incoming_transitions = incoming_transitions_all[buchi]
            outgoing_transitions = outgoing_transitions_all[buchi]
            acceptance_states_usable = acceptance_states_usable_all[buchi]
            x = x_all[buchi]
            hist = processed_history[buchi]
            # print(hist)
            if (init_state, tuple(hit_accepting_states)) in set(hist.keys()): 
                ret_sum += hist[(init_state, tuple(hit_accepting_states))]
                continue

            """
            altered_matrices_A = acceptance_altered_matrices_A[buchi]
            altered_matrices_b = acceptance_altered_matrices_b[buchi]

            acceptance_states_usable = acceptance_states_usable_all[buchi]
            """


            """
            OK, HERE ARE THE STEPS WE ARE GOING TO TAKE:
            1. initialize the matrix with spaces in the acceptance state bits and init state bits
                (can do this in pre-processing
            2. find the state we are in based on the init_state and the accepting state, for each accepting state
            3. deterine the states for each acceptance state
            4. add constraints 
                4.5. ADD CONSTRIANTS FOR OTHER ACCEPTING STATES 
            5. solve
            6. add total to ret_sum
            7. clear accepting states
            """

            """
            1. 
            """

            # first, assemble all indices for acceptance states
            indicators = model.addVars(len(acceptance_states_usable), vtype=GRB.BINARY)

            init_inputs = incoming_transitions[init_state]
            init_outputs = outgoing_transitions[init_state] 

            init_input_sum = quicksum(x[int(t[1:])] for t in init_inputs) #x1
            init_output_sum = quicksum(x[int(t[1:])] for t in init_outputs) #x2 

            indicator_proc = 0
            for acceptance_state in acceptance_states_usable:

                current_indicator = indicators[indicator_proc]
                
                acceptance_inputs = incoming_transitions[acceptance_state] 
                acceptance_outputs = outgoing_transitions[acceptance_state] 
                acceptance_input_sum = quicksum(x[int(t[1:])] for t in acceptance_inputs) #y1
                acceptance_output_sum = quicksum(x[int(t[1:])] for t in acceptance_outputs) #y2

                if acceptance_state == init_state:
                    # case 2: eq
                    print("init==accept")

                    snd_indicators = model.addVars(3 + len(acceptance_states_usable) - 1, vtype=GRB.BINARY, name="si")

                    M1 = len(init_inputs) * 3
                    M2 = len(init_outputs) * 3
                    M3 = (M1 + M2) * 2

                    model.addConstr(init_input_sum - 1 <= M1 * snd_indicators[0])
                    model.addConstr(1 - init_input_sum <= M1 * (1-snd_indicators[0]))
                    model.addConstr(init_output_sum - 1 <= M2 * snd_indicators[1])
                    model.addConstr(1 - init_output_sum <= M2 * (1-snd_indicators[1]))

                    model.addConstr(init_input_sum - init_output_sum <= M3 * (1-snd_indicators[2]))
                    model.addConstr(init_output_sum - init_input_sum <= M3 * (1-snd_indicators[2]))

                    model.addConstr(current_indicator >= snd_indicators[0] + snd_indicators[1] + snd_indicators[2] - 2)
                    model.addConstr(current_indicator <= snd_indicators[0])
                    model.addConstr(current_indicator <= snd_indicators[1])
                    model.addConstr(current_indicator <= snd_indicators[2])

                    diff = acceptance_states_usable.difference({acceptance_state})
                    curr_indicator = 3

                    for n_acc in diff:
                        a_input = quicksum(x[int(t[1:])] for t in incoming_transitions[n_acc])
                        a_output = quicksum(x[int(t[1:])] for t in outgoing_transitions[n_acc])
                        model.addConstr(a_input - a_output <= M3 * (1-snd_indicators[curr_indicator]))
                        model.addConstr(a_output - a_input <= M3 * (1-snd_indicators[curr_indicator]))

                        curr_indicator+=1

                elif acceptance_state in list({item for tup in hit_accepting_states for item in tup}):
                    print("in cycle")

                    # case 3: in cycle

                    snd_indicators = model.addVars(2 + len(acceptance_states_usable) - 1, vtype=GRB.BINARY, name="si")

                    M1 = len(init_inputs) * 3 
                    M2 = len(init_outputs) * 3
                    M3 = (M1 + M2) * 2

                    # 2nd order constraints
                    model.addConstr(init_output_sum - init_input_sum - 1 <= M1 * (1-snd_indicators[0]) )
                    model.addConstr(1-init_output_sum + init_input_sum <= M1 * (1 - snd_indicators[0]) )
                    model.addConstr(acceptance_input_sum - acceptance_output_sum - 1 <= M2 * (1 - snd_indicators[1]) ) 
                    model.addConstr(1-acceptance_input_sum + acceptance_output_sum <= M2 * (1 - snd_indicators[1]) ) 

                    # 1st order constraints
                    model.addConstr(current_indicator >= snd_indicators[0] + snd_indicators[1] - 1)
                    model.addConstr(current_indicator <= snd_indicators[0])
                    model.addConstr(current_indicator <= snd_indicators[1])

                    diff = acceptance_states_usable.difference({acceptance_state})
                    curr_indicator = 2

                    for n_acc in diff:
                        a_input = quicksum(x[int(t[1:])] for t in incoming_transitions[n_acc])
                        a_output = quicksum(x[int(t[1:])] for t in outgoing_transitions[n_acc])
                        model.addConstr(a_input - a_output <= M3 * (1-snd_indicators[curr_indicator]))
                        model.addConstr(a_output - a_input <= M3 * (1-snd_indicators[curr_indicator]))

                        curr_indicator+=1

                else: 
                    print("out of cycle")
                    # case 1: out of cycle

                    snd_indicators = model.addVars(4 + len(acceptance_states_usable) - 1, vtype=GRB.BINARY)
                    # snd_indicators = model.addVars(4, vtype=GRB.BINARY)

                    M1 = len(init_inputs) * 3 
                    M2 = len(init_outputs) * 3
                    M3 = len(acceptance_inputs) * 3
                    M4 = len(acceptance_outputs) * 3

                    model.addConstr(init_output_sum - init_input_sum - 1 <= M1 * (1-snd_indicators[0]) )
                    model.addConstr(1-init_output_sum + init_input_sum <= M1 * (1 - snd_indicators[0]) )
                    model.addConstr(acceptance_input_sum - acceptance_output_sum - 1 <= M2 * (1 - snd_indicators[1]) ) 
                    model.addConstr(1-acceptance_input_sum + acceptance_output_sum <= M2 * (1 - snd_indicators[1]) ) 

                    model.addConstr(acceptance_input_sum - 2 <= M3 * (snd_indicators[2]))
                    model.addConstr(2 - acceptance_input_sum <= M3 * (1 - snd_indicators[2]))
                    model.addConstr(acceptance_output_sum - 1 <= M4 * (snd_indicators[3]))
                    model.addConstr(1 - acceptance_output_sum <= M4 * (1 - snd_indicators[3]))

                    model.addConstr(current_indicator >= snd_indicators[0] + snd_indicators[1] + snd_indicators[2] + snd_indicators[3] - 3)
                    model.addConstr(current_indicator <= snd_indicators[0])
                    model.addConstr(current_indicator <= snd_indicators[1])
                    model.addConstr(current_indicator <= snd_indicators[2])
                    model.addConstr(current_indicator <= snd_indicators[3])

 
                    diff = acceptance_states_usable.difference({acceptance_state})
                    curr_indicator = 4

                    for n_acc in diff:
                        a_input = quicksum(x[int(t[1:])] for t in incoming_transitions[n_acc])
                        a_output = quicksum(x[int(t[1:])] for t in outgoing_transitions[n_acc])
                        model.addConstr(a_input - a_output <= M3 * (1-snd_indicators[curr_indicator]))
                        model.addConstr(a_output - a_input <= M3 * (1-snd_indicators[curr_indicator]))

                        curr_indicator+=1

                indicator_proc+=1

            model.addConstr(quicksum(indicators) == 1)
            
            """
            remove init state 
            """
            A_init_c = A_init.copy()
            A_init_c = np.delete(A_init_c, state_mapping[init_state], 0)
            b_init_c = b_init.copy()
            b_init_c = np.delete(b_init_c, state_mapping[init_state], 0)

            model.addMConstr(A_init_c, x, '=', b_init_c)

            model.update()
            model.optimize()

            if model.status == GRB.INFEASIBLE:
                model.remove(model.getConstrs())
                model.update()

                return float('inf')

            elif model.status == GRB.OPTIMAL:
                tv = 0
                for v in x:
                    ret_sum += v.x
                    tv+=v.x
            

                hist[(init_state, tuple(hit_accepting_states))] = tv    

                model.remove(model.getConstrs())
                model.update()

        ret_sum_history[(tuple(buchi_current_states), tuple(acc_states))] = ret_sum

        return ret_sum

    # invoke queue
    total = 0
    while queue:
        total+=1

        current_distance, current_state_tuple = heapq.heappop(queue)

        current_state, prev_states, acc_states, trail = current_state_tuple
        discovered.add(current_state)

        # check if distance[current_state] has a value aleady.... if so, prev_states becomes the union

        if current_distance > distance[current_state] : continue 

        reachable = []
        for char in fa:
            new_states = []
            break_loop = False

            new_state = [None]*buchi_arr_len
            for i in range(buchi_arr_len):
                new_state[i] = buchi_array[i].transition_function(current_state[i], char)

                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            buchi_arr_len = len(buchi_array)
            l = list(itertools.product(*new_state))

            for new_state in l:

                if new_state in discovered and acc_states != set() and any([ new_state in acceptance_state_leading[astate] for astate in acc_states]):
                    trail.append((current_state, char))
                    trail.append((new_state, "*"))

                    final_state = new_state
                    in_cycle = False
                    word = ""

                    for state, char in trail:
                        if in_cycle:
                            if state == final_state : continue
                            word+=char
                        elif state == final_state: 
                            in_cycle = True
                            word+="("
                            word+=char
                        else:
                            word+=char

                    # indiciative of a cycle
                    word+=")*"
                    print(total)
                    return word
 
                if not new_state in distance : distance[new_state] = float('inf')
                if not new_state in word : word[new_state] = ""

                # that is, add to prev states, acc states, trail 
                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()
                # can be optimized... the for loop array is generated before
                # no, I *need* to optimize this
                addState = True
                for i in range(buchi_arr_len):
                    if new_state[i] in buchi_array[i].acceptance_states:
                        # do i add the new state?
                        hitAcceptingStatesDict[buchi_array[i]].add(new_state)

                    else:
                        addState = False

                if addState:
                    new_acc_states.add(new_state)
                    acceptance_state_leading[new_state] = new_prev_states.copy()

                """
                if all([ new_state[i] in buchi_array[i].acceptance_states for i in range(buchi_arr_len) ]):
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state)
                    acceptance_state_leading[new_state] = new_prev_states.copy()
                """

                trail_new = trail.copy()
                trail_new.append((current_state, char))

                t = (new_state, new_prev_states, new_acc_states, trail_new)

                reachable.append((t, char))

        for (state_t, char) in reachable:
            state = state_t[0]
            # print(state)

            hval = heuristic(state_t)
            #print(hval, state)

            # new_distance = distance[current_state] + hval * 5
            new_distance = distance[current_state] + 1 + hval
            # new_distance = hval

            if new_distance < distance[state]:
                distance[state] = new_distance
                word[state] = word[current_state] + char
                heapq.heappush(queue, (new_distance, state_t))

    return None


def fast_z3_heuristic_buchi_intersection(buchi_array):

    # assert all buchi automata have the same alphabet
    fa = buchi_array[0].alphabet
    """
    for buchi in buchi_array: 
        assert fa == buchi.alphabet
        assert isinstance(buchi, Buchi)
    """
        # print(len(buchi.states))

    initial_pairing = tuple([buchi.initial_state for buchi in buchi_array])
    cstr_state_mapping, rev_cstr_state_mapping = {}, {}
    cstr_state_mapping[initial_pairing] = 0 # initialize initial state as "0" (as in, state 0)
    rev_cstr_state_mapping[0] = initial_pairing

    initial_state = tuple([0, set(), set()])

    distance, word, queue = {}, {}, []
    
    heapq.heappush(queue, (0, initial_state))
    distance[0], word[0] = 0, ""

    # heurisitc required
    hitAcceptingStatesDict = {}

    # 1 : before hitting acceptance
    # 2 : in acceptance
    # 3 : after hitting acceptance
    heuristicBuchiState = {}

    transition_mapping_all = {}
    transition_mapping_reverse_all = {}

    incoming_transitions_all = {}
    outgoing_transitions_all = {}

    heuristic_wellconnected_states = {}

    processed_history = {}
    ret_sum_history = {}

    acceptance_state_leading = {}

    for buchi in buchi_array:
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

        hitAcceptingStatesDict[buchi] = set()

        usable_states = buchi.get_wellconnected_states()
        heuristic_wellconnected_states[buchi] = usable_states

        incoming_transitions = {}
        outgoing_transitions = {}

        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        ts = []
        for (sA, i, sB) in buchi.transitions: 
            if sA in usable_states and sB in usable_states:
                ts.append((sA, i, sB))

        transition_mapping = {}
        transition_mapping_reverse = {}
        transition_labels = [i for i in range(len(ts))]

        count = 0
        for (sA, i, sB) in ts:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            transition_mapping_reverse[(sA, i, sB)] = transition_labels[count]

            outgoing_transitions[sA].add(count)
            incoming_transitions[sB].add(count)

            count = count + 1

        incoming_transitions_all[buchi] = incoming_transitions
        outgoing_transitions_all[buchi] = outgoing_transitions
        transition_mapping_all[buchi] = transition_mapping

        processed_history[buchi] = {}

    def heuristic(state):
        ret_sum = 0
        buchi_state_label, prev_states, acc_states = state
        buchi_current_states = rev_cstr_state_mapping[buchi_state_label]
        for i in range(len(buchi_array)):

            buchi = buchi_array[i]  
            current_state = buchi_current_states[i]
            init_state = current_state

            incoming_transitions = incoming_transitions_all[buchi]
            outgoing_transitions = outgoing_transitions_all[buchi]

            # need to determine if we've visited any accepting states for the current "i" in buchi...
            # because we need to determine how we look for cycles
            hit_accepting_states = hitAcceptingStatesDict[buchi]

            hist = processed_history[buchi]
            if (init_state, tuple(hit_accepting_states)) in set(hist.keys()): 
                ret_sum += hist[(init_state, tuple(hit_accepting_states))]
                continue

            transition_mapping = transition_mapping_all[buchi]
            usable_states = heuristic_wellconnected_states[buchi]

            s = Optimize()
            tosolve = {}
            for t in list(transition_mapping.keys()):
                tosolve[t] = Int(str(t))

            for t in tosolve.values():
                s.add(t >= 0)

            # minimize solution
            s.minimize(Sum(list(tosolve.values())))

            accepting_states_to_solve = set()

            for state in usable_states:
                if state in buchi.acceptance_states:
                    accepting_states_to_solve.add(state)
                elif state == init_state:
                    continue
                else: 
                    s.add(Sum([tosolve[t] for t in incoming_transitions[state]]) == Sum([tosolve[t] for t in outgoing_transitions[state]]))

            accept_or = []
            for accepting_state in accepting_states_to_solve:
                # here we determine the state...
                if accepting_state == init_state:
                    # print("equal")

                    accept_or.append(
                            And(
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]]) == 1),
                            (Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) == 1)
                            ))
                elif accepting_state in list({item for tup in hit_accepting_states for item in tup}):
                    # print("in")

                    accept_or.append(
                            And(
                            (Sum([tosolve[t] for t in incoming_transitions[init_state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[init_state]])),
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]])  == Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) + 1)
                            ))
                else:
                    # print("out")

                    # AND VER:
                    accept_or.append(
                            And(
                            (Sum([tosolve[t] for t in incoming_transitions[init_state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[init_state]])),
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]]) == Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) + 1),
                            (Sum([tosolve[t] for t in incoming_transitions[accepting_state]]) == 2),
                            (Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) == 1)
                            ))

            # print(*accept_or)
            s.add(Or(*accept_or)) 

            if s.check() == sat:
                m = s.model()
                solns = {var: m[var] for var in tosolve.values()}

                total_dist = 0
                for transition, value in solns.items():
                    # sA, i, sB = transition_mapping[str(transition)]
                    total_dist+=int(str(value)) # kinda slow, can be optimized 

                hist[(init_state, tuple(hit_accepting_states))] = total_dist

                ret_sum += total_dist

            else : return float('inf')
        
        return ret_sum

    while queue:
        current_distance, current_state_tuple = heapq.heappop(queue)
        current_state, prev_states, acc_states = current_state_tuple 

        # termination condition
        if (current_state in prev_states and
            (current_state in acc_states 
             or 
             any([ current_state in acceptance_state_leading[astate] for astate in acc_states])
             )
            ):

            word = word[current_state]

            print("counterexample: " + str(word))
            total = len(cstr_state_mapping)
            
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])
            print("percentage of states visited: ", round(total/total_possible*100, 3), "% -",str(total), "/", str(total_possible))

            return word

        elif current_state in prev_states : continue

        # skip condition : if our current distance is longer than established distance
        if current_distance > distance[current_state] : continue
        
        current_state_ite = rev_cstr_state_mapping[current_state]

        # assemble new states
        reachable = []
        for char in fa:
            new_states, break_loop = [], False

            new_state = [None] * len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state_ite[i], char)

                # if this is ever the case, we know the transition function led to nothing
                # thus this state is invalid
                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l:

                # if this is the case, we know we've already seen it before
                new_state_label = None
                if not new_state in cstr_state_mapping:
                    # print("new state!")
                    new_state_label = len(cstr_state_mapping)

                    cstr_state_mapping[new_state] = new_state_label
                    rev_cstr_state_mapping[new_state_label] = new_state

                    distance[new_state_label] = float('inf')
                    word[new_state_label] = ""

                else: 
                    new_state_label = cstr_state_mapping[new_state]

                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()

                all_accepting = True
                for i in range(len(buchi_array)):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False
                        break

                if all_accepting:
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state_label)
                    acceptance_state_leading[new_state_label] = new_prev_states.copy()

                t = (new_state_label, new_prev_states, new_acc_states)

                reachable.append((t, new_state, char))

        for (state_t, state, char) in reachable:
            state_label = state_t[0]


            hval = heuristic(state_t)

            # new_distance = distance[current_state] + 1 + hval
            new_distance = hval 

            distance[state_label] = new_distance
            word[state_label] = word[current_state] + char
            heapq.heappush(queue, (new_distance, state_t))

    return None

# it's actually nuts
def gurobi_zone_heuristic_buchi_intersection(buchi_array):
    """
    steps:
    1. initialize the standard things
    2. loop through each buchi automata
        2.1. get the heuristic wellconnected states for each buchi
        2.2. get the flow constraints for each buchi automata (linear time)
        2.3. with the flow constraints, for each acceptance state get the minimum (if any) acceptance cycle dist
            - we should store a list of acceptance states that'll actually matter
        2.4. throw the rest of the transitions into a matrix representation for gurobi
        2.5. simultanious BFS: establish distance between any state and the acceptance cycle
            2.5.1. start by assigning all acceptance states a distance value equal to the minimum acceptance cycle distance
            2.5.2. from each acceptance state, calculate the distance in waves, e.g., distance[new_state] = distance[state] + 1
            2.5.3. we only throw a new state on the stack if it is at least two less than the new state
            2.5.4. repeat until stack is empty (need to prove things about this algorithm...)
    3. define heuristic
        3.1. when we move through each automata we're intersecting, we keep track of acceptance states we've seen
        3.2. for all acceptance states we've visited, calculate the parikh image flow constraint-based distance using the matrix representation from the pre-processing
            3.2.note. we could use big M to turn the much lesser "or" here into a more efficient, single operation. but, given it's all over matrices it should be ok in practice.
            3.2.1. 
        3.3. to get the minimum distance, calculate minimum(all parikh image distances, distance[state])
    4. use Astar, move through the graph

    I think in the worst case, which is the all accepting case, this would be much much better. our disjunction is so much less nasty.
    this way we:
        - reduce having a three-cased system like in z3_heuristic_buchi_intersection
        - reduce our disjunction to something extremely managable. it's only over all acceptance states *visited*, so it's actually quite managable
    """

    fa = buchi_array[0].alphabet

    initial_pair = tuple([buchi.initial_state for buchi in buchi_array])
    cstr_state_mapping = { initial_pair : 0 }
    rev_cstr_state_mapping = { 0 : initial_pair }
    initial_state = tuple([0, set(), set()])

    distance = { 0 : 0 }
    word = { 0 : "" }
    queue = []
    heapq.heappush(queue, (0, initial_state))

    # all the dicts....
    all_wellconnected_states = {}
    all_incoming_transitions = {}
    all_outgoing_transitions = {}
    all_transition_mappings = {}
    all_rev_transition_mappings = {}
    all_transition_index_mappings = {}
    all_rev_transition_index_mappings = {}
    all_state_mappings = {}
    all_rev_state_mappings = {}
    all_state_distance_mappings = {}
    all_A = {}
    all_b = {}
    all_x = {}
    all_models = {}
    all_hit_acceptance_states = {}
    acceptance_state_leading = {}

    # and history
    processed_history = {}
    ret_sum_history = {}

    for buchi in buchi_array:
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

        """
        2.1. get the heuristic wellconnected states for each buchi
        """
        usable_states = buchi.get_wellconnected_states()
        all_wellconnected_states[buchi] = usable_states

        # also, get the usable acceptance states
        acceptance_states_usable = set()
        for astate in buchi.acceptance_states:
            if astate in usable_states : acceptance_states_usable.add(astate)

        if acceptance_states_usable == set() : return None

        """
        2.2. get the flow constraints for each buchi automata (linear time)
        2.3. with the flow constraints, for each acceptance state get the minimum (if any) acceptance cycle dist
            - we should store a list of acceptance states that'll actually matter
        2.4. throw the rest of the transitions into a matrix representation for gurobi
        """
        incoming_transitions = {}
        outgoing_transitions = {}

        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        ts = []
        for (sA, i, sB) in buchi.transitions:
            if sA in usable_states and sB in usable_states:
                ts.append((sA, i, sB))

        transition_mapping = {} # map transitions to integer values
        rev_transition_mapping = {} # map integer values to transitions
        transition_labels = [i for i in range(len(ts))]

        count = 0
        for (sA, i, sB) in ts:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            rev_transition_mapping[(sA, i, sB)] = transition_labels[count]

            outgoing_transitions[sA].add(count)
            incoming_transitions[sB].add(count)

            count = count + 1

        # note how the incoming and outgoing transition map to the labels

        # we want to create a matrix for all the flow constraints, 
        # and based on the current acceptance condition, add a value
        # constraining x

        A_init = np.zeros( (len(usable_states), len(ts)) ) 
        b_init = np.zeros( (len(usable_states), 1) ) 

        # need to create transition and state mappings to indices of matrix...
        transition_index_mapping = {} # map transition label to a column
        rev_transition_index_mapping = {} # map column to a transition label

        state_mapping = {} # map state to a row
        rev_state_mapping = {} # map row to a state

        tc = 0
        for (sA, i, sB) in ts:
            transition_index_mapping[tc] = rev_transition_mapping[(sA, i, sB)]
            rev_transition_index_mapping[rev_transition_mapping[(sA, i, sB)]] = tc
            tc+=1

        sc = 0
        for state in usable_states:
            state_mapping[state] = sc
            rev_state_mapping[sc] = state

            inputs = incoming_transitions[state]
            output = outgoing_transitions[state]

            for i in inputs:
                r = transition_index_mapping[i]
                A_init[sc][r]+=1

            for i in output:
                r = transition_index_mapping[i]
                A_init[sc][r]-=1

            sc+=1

        # create x for Ax=b
        model = Model('GAMING')
        model.setParam('OutputFlag', 0)
        x = model.addMVar(len(ts), vtype=GRB.INTEGER)
        model.setObjective(quicksum(x), GRB.MINIMIZE)

        """
        okay, now we need to:
        1. loop over all usable acceptance states
        2. alter A_init by deleting a row (after copying)
        3. constrain the right variables in x based on the inputs and outputs of the acceptance state
        4. throw things into the solver, and solve
        5. fill in some dict
        6. remove constraints, go again
        """
        # maps each state to the distance away from an accepting state, plus the accepting cycle
        state_distance_mapping = {} 
    
        for state in acceptance_states_usable:
            A_copy = np.delete(A_init, state_mapping[state], 0)
            b_copy = np.delete(b_init, state_mapping[state], 0)

            inputs = incoming_transitions[state]
            outputs = outgoing_transitions[state]
            # a reminder that these get *labels*

            c1 = model.addConstr(quicksum(x[transition_index_mapping[i]] for i in inputs) == 1,name="main1")
            c2 = model.addConstr(quicksum(x[transition_index_mapping[o]] for o in outputs) == 1,name="main2")

            c3 = model.addMConstr(A_copy, x, '=', b_copy, name="main3")
            model.update()
            model.optimize()
            
            if model.status == GRB.INFEASIBLE:
                pass
                # state_distance_mapping[state] = float('inf')
            elif model.status == GRB.OPTIMAL:
                state_distance_mapping[state] = sum(x.X)

            model.remove(model.getConstrByName("main1"))
            model.remove(model.getConstrByName("main2"))
            for constr in c3:
                model.remove(constr)

        model.update()

        # assert things are the same (I don't trust numpy)
        num_states_a, num_vars_A = A_init.shape
        num_states_b, num_vars_b = b_init.shape
        assert num_vars_A == len(ts)
        assert num_states_a == len(usable_states)
        assert num_states_b == len(usable_states)
        assert num_vars_b == 1

        """
        2.5. simultanious BFS: establish distance between any state and the acceptance cycle
            2.5.1. start by assigning all acceptance states a distance value equal to the minimum acceptance cycle distance
            2.5.2. from each acceptance state, calculate the distance in waves, e.g., distance[new_state] = distance[state] + 1
            2.5.3. we only throw a new state on the stack if it is at least two less than the new state
            2.5.4. repeat until stack is empty (need to prove things about this algorithm...)
        """

        qx = []
        for state in state_distance_mapping.keys() : qx.append(state)

        inf_states = usable_states.difference(set(state_distance_mapping.keys()))
        for state in inf_states : state_distance_mapping[state] = float('inf')

        while qx:
            curr_state = qx.pop(0)
            curr_dist = state_distance_mapping[curr_state]
            rev_reachable = buchi.reachability_object_reverse[curr_state]

            for (r_state, char) in rev_reachable:
                r_dist = state_distance_mapping[r_state]

                if curr_dist + 1 < r_dist:
                    state_distance_mapping[r_state] = curr_dist + 1
                    qx.append(r_state)
                    
        all_incoming_transitions[buchi] = incoming_transitions
        all_outgoing_transitions[buchi] = outgoing_transitions
        all_transition_mappings[buchi] = transition_mapping
        all_rev_transition_mappings[buchi] = rev_transition_mapping
        all_transition_index_mappings[buchi] = transition_index_mapping
        all_rev_transition_index_mappings[buchi] = rev_transition_index_mapping
        all_state_mappings[buchi] = state_mapping
        all_rev_state_mappings[buchi] = rev_state_mapping
        all_state_distance_mappings[buchi] = state_distance_mapping
        all_A[buchi] = A_init
        all_b[buchi] = b_init
        all_x[buchi] = x
        all_models[buchi] = model
        all_hit_acceptance_states[buchi] = set()

        processed_history[buchi] = {}

    def heuristic(state):
        ret_sum = 0
        buchi_state_label, prev_states, acc_states = state
        buchi_current_states = rev_cstr_state_mapping[buchi_state_label]

        # if (tuple(buchi_current_states), tuple(acc_states)) in set(ret_sum_history.keys()):
        #    return ret_sum_history[(tuple(buchi_current_states), tuple(acc_states))]

        for i in range(len(buchi_array)):

            buchi = buchi_array[i]
            init_state = buchi_current_states[i]

            # getting info from all the dicts...
            usable_states = all_wellconnected_states[buchi]
            incoming_transitions = all_incoming_transitions[buchi]
            outgoing_transitions = all_outgoing_transitions[buchi]
            transition_mapping = all_transition_mappings[buchi]
            rev_transition_mapping = all_rev_transition_mappings[buchi]
            transition_index_mapping = all_transition_index_mappings[buchi]
            rev_transition_index_mapping = all_rev_transition_index_mappings[buchi]
            state_mapping = all_state_mappings[buchi]
            rev_state_mapping = all_rev_state_mappings[buchi]
            state_distance_mapping = all_state_distance_mappings[buchi]
            A_init = all_A[buchi]
            b_init = all_b[buchi]
            x = all_x[buchi]
            model = all_models[buchi]

            hit_accepting_states = all_hit_acceptance_states[buchi]

            # checking history
            hist = processed_history[buchi]
            if (init_state, tuple(hit_accepting_states)) in set(hist.keys()): 
                ret_sum += hist[(init_state, tuple(hit_accepting_states))]
                continue

            """
            3. define heuristic
                3.1. when we move through each automata we're intersecting, we keep track of acceptance states we've seen
                3.2. for all acceptance states we've visited, calculate the parikh image flow constraint-based distance using the matrix representation from the pre-processing
                    3.2.note. we could use big M to turn the much lesser "or" here into a more efficient, single operation. but, given it's all over matrices it should be ok in practice.
                    3.2.1. 
                3.3. to get the minimum distance, calculate minimum(all parikh image distances, distance[state])
            """

            curr_distance = state_distance_mapping[init_state]

            if len(hit_accepting_states) == 0:
                ret_sum+=curr_distance
                hist[(init_state, tuple(hit_accepting_states))] = curr_distance
                continue

            to_del = []
            for astate in hit_accepting_states:
                if astate == init_state : continue
                to_del.append(state_mapping[astate])

            sorted(to_del, reverse=True)

            A_matrix = np.delete(A_init, to_del, axis=0)
            b_matrix = np.delete(b_init, to_del, axis=0)

            b_matrix[state_mapping[init_state]] = -1

            # create the "or" variables between the constraints of accepting states
            indicators = model.addVars(len(to_del), vtype=GRB.BINARY)
            
            counter = 0
            for astate in hit_accepting_states:
                if astate == init_state : continue

                indicator = indicators[counter]
                counter+=1

                inputs = incoming_transitions[astate]
                outputs = outgoing_transitions[astate]

                model.addConstr(quicksum(x[transition_index_mapping[i]] for i in inputs) - quicksum(x[transition_index_mapping[o]] for o in outputs) == indicator)

            model.addConstr(quicksum(indicators) == 1)
            model.addConstr(quicksum(x) + 1 <= curr_distance)
            model.addMConstr(A_matrix, x, '=', b_matrix)

            model.update()
            model.optimize()

            q = None
            if model.status == GRB.INFEASIBLE:
                ret_sum+=curr_distance
                q = curr_distance
            elif model.status == GRB.OPTIMAL:
                s = sum(x.X)
                ret_sum+=s
                q = s
        
            hist[(init_state, tuple(hit_accepting_states))] = q

            model.remove(model.getConstrs())
            # model.update()

        # ret_sum_history[(tuple(buchi_current_states), tuple(acc_states))] = ret_sum
        # print(ret_sum)

        return ret_sum

    print("running begin")

    while queue:
        current_distance, current_state_tuple = heapq.heappop(queue)
        current_state, prev_states, acc_states = current_state_tuple 

        # termination condition
        if (current_state in prev_states and
            (current_state in acc_states 
             or 
             any([ current_state in acceptance_state_leading[astate] for astate in acc_states])
             )
            ):

            word = word[current_state]

            print("counterexample: " + str(word))
            total = len(cstr_state_mapping)
            
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])
            print("percentage of states visited: ", round(total/total_possible*100, 3), "% -",str(total), "/", str(total_possible))

            return word

        elif current_state in prev_states : continue

        # skip condition : if our current distance is longer than established distance
        if current_distance > distance[current_state] : continue
        
        current_state_ite = rev_cstr_state_mapping[current_state]

        # assemble new states
        reachable = []
        for char in fa:
            new_states, break_loop = [], False

            new_state = [None] * len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state_ite[i], char)

                if new_state[i] in buchi_array[i].acceptance_states : all_hit_acceptance_states[buchi].add(new_state[i])

                # if this is ever the case, we know the transition function led to nothing
                # thus this state is invalid
                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l:

                # if this is the case, we know we've already seen it before
                new_state_label = None
                if not new_state in cstr_state_mapping:
                    # print("new state!")
                    new_state_label = len(cstr_state_mapping)

                    cstr_state_mapping[new_state] = new_state_label
                    rev_cstr_state_mapping[new_state_label] = new_state

                    distance[new_state_label] = float('inf')
                    word[new_state_label] = ""

                else: 
                    new_state_label = cstr_state_mapping[new_state]

                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()

                all_accepting = True
                for i in range(len(buchi_array)):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False
                        break

                if all_accepting:
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state_label)
                    acceptance_state_leading[new_state_label] = new_prev_states.copy()

                t = (new_state_label, new_prev_states, new_acc_states)

                reachable.append((t, new_state, char))

        for (state_t, state, char) in reachable:
            state_label = state_t[0]

            hval = heuristic(state_t)

            # new_distance = distance[current_state] + 1 + hval
            new_distance = hval 

            distance[state_label] = new_distance
            word[state_label] = word[current_state] + char
            heapq.heappush(queue, (new_distance, state_t))

    return None

# painfully unoptimized. todo, implement in C or rust or something
def z3_zone_heuristic_buchi_intersection(buchi_array):
    fa = buchi_array[0].alphabet

    initial_pair = tuple([buchi.initial_state for buchi in buchi_array])
    cstr_state_mapping = { initial_pair : 0 }
    rev_cstr_state_mapping = { 0 : initial_pair }
    initial_state = tuple([0, set(), set()])

    distance = { 0 : 0 }
    word = { 0 : "" }
    queue = []
    heapq.heappush(queue, (0, initial_state))

    # all the dicts....
    all_wellconnected_states = {}
    all_incoming_transitions = {}
    all_outgoing_transitions = {}
    all_transition_mappings = {}
    all_rev_transition_mappings = {}
    all_transition_index_mappings = {}
    all_rev_transition_index_mappings = {}
    all_state_mappings = {}
    all_rev_state_mappings = {}
    all_state_distance_mappings = {}
    all_A = {}
    all_b = {}
    all_x = {}
    all_models = {}
    all_hit_acceptance_states = {}
    acceptance_state_leading = {}

    # and history
    processed_history = {}
    ret_sum_history = {}

    for buchi in buchi_array:
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

        """
        2.1. get the heuristic wellconnected states for each buchi
        """
        usable_states = buchi.get_wellconnected_states()
        all_wellconnected_states[buchi] = usable_states

        # also, get the usable acceptance states
        acceptance_states_usable = set()
        for astate in buchi.acceptance_states:
            if astate in usable_states : acceptance_states_usable.add(astate)

        if acceptance_states_usable == set() : return None

        """
        2.2. get the flow constraints for each buchi automata (linear time)
        2.3. with the flow constraints, for each acceptance state get the minimum (if any) acceptance cycle dist
            - we should store a list of acceptance states that'll actually matter
        2.4. throw the rest of the transitions into a matrix representation for gurobi
        """
        incoming_transitions = {}
        outgoing_transitions = {}

        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        ts = []
        for (sA, i, sB) in buchi.transitions:
            if sA in usable_states and sB in usable_states:
                ts.append((sA, i, sB))

        transition_mapping = {} # map transitions to integer values
        rev_transition_mapping = {} # map integer values to transitions
        transition_labels = [i for i in range(len(ts))]

        count = 0
        for (sA, i, sB) in ts:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            rev_transition_mapping[(sA, i, sB)] = transition_labels[count]

            outgoing_transitions[sA].add(count)
            incoming_transitions[sB].add(count)

            count = count + 1

        # note how the incoming and outgoing transition map to the labels

        # we want to create a matrix for all the flow constraints, 
        # and based on the current acceptance condition, add a value
        # constraining x

        A_init = np.zeros( (len(usable_states), len(ts)) ) 
        b_init = np.zeros( (len(usable_states), 1) ) 

        # need to create transition and state mappings to indices of matrix...
        transition_index_mapping = {} # map transition label to a column
        rev_transition_index_mapping = {} # map column to a transition label

        state_mapping = {} # map state to a row
        rev_state_mapping = {} # map row to a state

        tc = 0
        for (sA, i, sB) in ts:
            transition_index_mapping[tc] = rev_transition_mapping[(sA, i, sB)]
            rev_transition_index_mapping[rev_transition_mapping[(sA, i, sB)]] = tc
            tc+=1

        sc = 0
        for state in usable_states:
            state_mapping[state] = sc
            rev_state_mapping[sc] = state

            inputs = incoming_transitions[state]
            output = outgoing_transitions[state]

            for i in inputs:
                r = transition_index_mapping[i]
                A_init[sc][r]+=1

            for i in output:
                r = transition_index_mapping[i]
                A_init[sc][r]-=1

            sc+=1

        # create x for Ax=b
        model = Model('GAMING')
        model.setParam('OutputFlag', 0)
        x = model.addMVar(len(ts), vtype=GRB.INTEGER)
        model.setObjective(quicksum(x), GRB.MINIMIZE)

        """
        okay, now we need to:
        1. loop over all usable acceptance states
        2. alter A_init by deleting a row (after copying)
        3. constrain the right variables in x based on the inputs and outputs of the acceptance state
        4. throw things into the solver, and solve
        5. fill in some dict
        6. remove constraints, go again
        """
        # maps each state to the distance away from an accepting state, plus the accepting cycle
        state_distance_mapping = {} 
    
        for state in acceptance_states_usable:
            A_copy = np.delete(A_init, state_mapping[state], 0)
            b_copy = np.delete(b_init, state_mapping[state], 0)

            inputs = incoming_transitions[state]
            outputs = outgoing_transitions[state]
            # a reminder that these get *labels*

            c1 = model.addConstr(quicksum(x[transition_index_mapping[i]] for i in inputs) == 1,name="main1")
            c2 = model.addConstr(quicksum(x[transition_index_mapping[o]] for o in outputs) == 1,name="main2")

            c3 = model.addMConstr(A_copy, x, '=', b_copy, name="main3")
            model.update()
            model.optimize()
            
            if model.status == GRB.INFEASIBLE:
                pass
                # state_distance_mapping[state] = float('inf')
            elif model.status == GRB.OPTIMAL:
                state_distance_mapping[state] = sum(x.X)

            model.remove(model.getConstrByName("main1"))
            model.remove(model.getConstrByName("main2"))
            for constr in c3:
                model.remove(constr)

        model.update()

        # assert things are the same (I don't trust numpy)
        num_states_a, num_vars_A = A_init.shape
        num_states_b, num_vars_b = b_init.shape
        assert num_vars_A == len(ts)
        assert num_states_a == len(usable_states)
        assert num_states_b == len(usable_states)
        assert num_vars_b == 1

        """
        2.5. simultanious BFS: establish distance between any state and the acceptance cycle
            2.5.1. start by assigning all acceptance states a distance value equal to the minimum acceptance cycle distance
            2.5.2. from each acceptance state, calculate the distance in waves, e.g., distance[new_state] = distance[state] + 1
            2.5.3. we only throw a new state on the stack if it is at least two less than the new state
            2.5.4. repeat until stack is empty (need to prove things about this algorithm...)
        """

        qx = []
        for state in state_distance_mapping.keys() : qx.append(state)

        inf_states = usable_states.difference(set(state_distance_mapping.keys()))
        for state in inf_states : state_distance_mapping[state] = float('inf')

        while qx:
            curr_state = qx.pop(0)
            curr_dist = state_distance_mapping[curr_state]
            rev_reachable = buchi.reachability_object_reverse[curr_state]

            for (r_state, char) in rev_reachable:
                r_dist = state_distance_mapping[r_state]

                if curr_dist + 1 < r_dist:
                    state_distance_mapping[r_state] = curr_dist + 1
                    qx.append(r_state)
                    
        all_incoming_transitions[buchi] = incoming_transitions
        all_outgoing_transitions[buchi] = outgoing_transitions
        all_transition_mappings[buchi] = transition_mapping
        all_rev_transition_mappings[buchi] = rev_transition_mapping
        all_transition_index_mappings[buchi] = transition_index_mapping
        all_rev_transition_index_mappings[buchi] = rev_transition_index_mapping
        all_state_mappings[buchi] = state_mapping
        all_rev_state_mappings[buchi] = rev_state_mapping
        all_state_distance_mappings[buchi] = state_distance_mapping
        all_A[buchi] = A_init
        all_b[buchi] = b_init
        all_x[buchi] = x
        all_models[buchi] = model
        all_hit_acceptance_states[buchi] = set()

        processed_history[buchi] = {}

    def heuristic(state):
        ret_sum = 0
        buchi_state_label, prev_states, acc_states = state
        buchi_current_states = rev_cstr_state_mapping[buchi_state_label]

        #if (tuple(buchi_current_states), tuple(acc_states)) in set(ret_sum_history.keys()):
        #    return ret_sum_history[(tuple(buchi_current_states), tuple(acc_states))]

        for i in range(len(buchi_array)):

            buchi = buchi_array[i]
            init_state = buchi_current_states[i]

            # getting info from all the dicts...
            usable_states = all_wellconnected_states[buchi]
            incoming_transitions = all_incoming_transitions[buchi]
            outgoing_transitions = all_outgoing_transitions[buchi]
            transition_mapping = all_transition_mappings[buchi]
            rev_transition_mapping = all_rev_transition_mappings[buchi]
            transition_index_mapping = all_transition_index_mappings[buchi]
            rev_transition_index_mapping = all_rev_transition_index_mappings[buchi]
            state_mapping = all_state_mappings[buchi]
            rev_state_mapping = all_rev_state_mappings[buchi]
            state_distance_mapping = all_state_distance_mappings[buchi]
            A_init = all_A[buchi]
            b_init = all_b[buchi]
            x = all_x[buchi]
            model = all_models[buchi]

            hit_accepting_states = all_hit_acceptance_states[buchi]

            # checking history
            hist = processed_history[buchi]
            if (init_state, tuple(hit_accepting_states)) in set(hist.keys()): 
                ret_sum += hist[(init_state, tuple(hit_accepting_states))]
                continue

            """
            3. define heuristic
                3.1. when we move through each automata we're intersecting, we keep track of acceptance states we've seen
                3.2. for all acceptance states we've visited, calculate the parikh image flow constraint-based distance using the matrix representation from the pre-processing
                    3.2.note. we could use big M to turn the much lesser "or" here into a more efficient, single operation. but, given it's all over matrices it should be ok in practice.
                    3.2.1. 
                3.3. to get the minimum distance, calculate minimum(all parikh image distances, distance[state])
            """

            curr_distance = state_distance_mapping[init_state]

            if len(hit_accepting_states) == 0:
                ret_sum+=curr_distance
                hist[(init_state, tuple(hit_accepting_states))] = curr_distance
                continue

            s = Optimize()
            for t in list(transition_mapping.keys()):
                tosolve[t] = Int(str(t))

            for t in tosolve.values():
                s.add(t >= 0)

            s.minimize(Sum(list(tosolve.values())))

            accepting_states_to_solve = set()

            for state in usable_states:
                if state in hit_accepting_states:
                    accepting_states_to_solve.add(state)
                elif state == init_state:
                    continue
                else: 
                    s.add(Sum([tosolve[t] for t in incoming_transitions[state]]) == Sum([tosolve[t] for t in outgoing_transitions[state]]))

            accept_or = []
            for accepting_state in accepting_states_to_solve:
                accept_or.append(
                        And(
                        (Sum([tosolve[t] for t in incoming_transitions[init_state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[init_state]])),
                        (Sum([tosolve[t] for t in incoming_transitions[accepting_state]])  == Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) + 1)
                        ))

            s.add(Or(*accept_or))

            if s.check() == sat:
                print("hit!")
                m = s.model()
                solns = {var: m[var] for var in tosolve.values()}

                total_dist = 0
                for transition, value in solns.items():
                    # sA, i, sB = transition_mapping[str(transition)]
                    total_dist+=int(str(value)) # kinda slow, can be optimized 

                hist[(init_state, tuple(hit_accepting_states))] = total_dist

                ret_sum += total_dist

            else : return float('inf')

        # ret_sum_history[(tuple(buchi_current_states), tuple(acc_states))] = ret_sum
        # print(ret_sum)

        return ret_sum

    print("running begin")

    while queue:
        current_distance, current_state_tuple = heapq.heappop(queue)
        current_state, prev_states, acc_states = current_state_tuple 

        # termination condition
        if (current_state in prev_states and
            (current_state in acc_states 
             or 
             any([ current_state in acceptance_state_leading[astate] for astate in acc_states])
             )
            ):

            word = word[current_state]

            print("counterexample: " + str(word))
            total = len(cstr_state_mapping)
            
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])
            print("percentage of states visited: ", round(total/total_possible*100, 3), "% -",str(total), "/", str(total_possible))

            return word

        elif current_state in prev_states : continue

        # skip condition : if our current distance is longer than established distance
        if current_distance > distance[current_state] : continue
        
        current_state_ite = rev_cstr_state_mapping[current_state]

        # assemble new states
        reachable = []
        for char in fa:
            new_states, break_loop = [], False

            new_state = [None] * len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state_ite[i], char)

                if new_state[i] in buchi_array[i].acceptance_states : all_hit_acceptance_states[buchi].add(new_state[i])

                # if this is ever the case, we know the transition function led to nothing
                # thus this state is invalid
                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l:

                # if this is the case, we know we've already seen it before
                new_state_label = None
                if not new_state in cstr_state_mapping:
                    # print("new state!")
                    new_state_label = len(cstr_state_mapping)

                    cstr_state_mapping[new_state] = new_state_label
                    rev_cstr_state_mapping[new_state_label] = new_state

                    distance[new_state_label] = float('inf')
                    word[new_state_label] = ""

                else: 
                    new_state_label = cstr_state_mapping[new_state]

                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()

                all_accepting = True
                for i in range(len(buchi_array)):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False
                        break

                if all_accepting:
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state_label)
                    acceptance_state_leading[new_state_label] = new_prev_states.copy()

                t = (new_state_label, new_prev_states, new_acc_states)

                reachable.append((t, new_state, char))

        for (state_t, state, char) in reachable:
            state_label = state_t[0]

            hval = heuristic(state_t)

            # new_distance = distance[current_state] + 1 + hval
            # because our heuristic is very accurate, best first search works 
            # well here
            new_distance = hval 

            distance[state_label] = new_distance
            word[state_label] = word[current_state] + char
            heapq.heappush(queue, (new_distance, state_t))

    return None

def z3_zone_heuristic_buchi_intersection_threading(buchi_array):
    fa = buchi_array[0].alphabet

    initial_pair = tuple([buchi.initial_state for buchi in buchi_array])
    cstr_state_mapping = { initial_pair : 0 }
    rev_cstr_state_mapping = { 0 : initial_pair }
    initial_state = tuple([0, set(), set()])

    distance = { 0 : 0 }
    word = { 0 : "" }
    queue = []
    heapq.heappush(queue, (0, initial_state))

    # all the dicts....
    all_wellconnected_states = {}
    all_incoming_transitions = {}
    all_outgoing_transitions = {}
    all_transition_mappings = {}
    all_rev_transition_mappings = {}
    all_transition_index_mappings = {}
    all_rev_transition_index_mappings = {}
    all_state_mappings = {}
    all_rev_state_mappings = {}
    all_state_distance_mappings = {}
    all_A = {}
    all_b = {}
    all_x = {}
    all_models = {}
    all_hit_acceptance_states = {}
    acceptance_state_leading = {}

    # and history
    processed_history = {}
    ret_sum_history = {}

    for buchi in buchi_array:
        assert isinstance(buchi, Buchi)
        assert fa == buchi.alphabet

        """
        2.1. get the heuristic wellconnected states for each buchi
        """
        usable_states = buchi.fast_get_wellconnected_states()
        all_wellconnected_states[buchi] = usable_states

        # also, get the usable acceptance states
        acceptance_states_usable = set()
        for astate in buchi.acceptance_states:
            if astate in usable_states : acceptance_states_usable.add(astate)

        if acceptance_states_usable == set() : return None

        """
        2.2. get the flow constraints for each buchi automata (linear time)
        2.3. with the flow constraints, for each acceptance state get the minimum (if any) acceptance cycle dist
            - we should store a list of acceptance states that'll actually matter
        2.4. throw the rest of the transitions into a matrix representation for gurobi
        """
        incoming_transitions = {}
        outgoing_transitions = {}

        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        ts = []
        for (sA, i, sB) in buchi.transitions:
            if sA in usable_states and sB in usable_states:
                ts.append((sA, i, sB))

        transition_mapping = {} # map transitions to integer values
        rev_transition_mapping = {} # map integer values to transitions
        transition_labels = [i for i in range(len(ts))]

        count = 0
        for (sA, i, sB) in ts:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            rev_transition_mapping[(sA, i, sB)] = transition_labels[count]

            outgoing_transitions[sA].add(count)
            incoming_transitions[sB].add(count)

            count = count + 1

        # note how the incoming and outgoing transition map to the labels

        # we want to create a matrix for all the flow constraints, 
        # and based on the current acceptance condition, add a value
        # constraining x

        A_init = np.zeros( (len(usable_states), len(ts)) ) 
        b_init = np.zeros( (len(usable_states), 1) ) 

        # need to create transition and state mappings to indices of matrix...
        transition_index_mapping = {} # map transition label to a column
        rev_transition_index_mapping = {} # map column to a transition label

        state_mapping = {} # map state to a row
        rev_state_mapping = {} # map row to a state

        tc = 0
        for (sA, i, sB) in ts:
            transition_index_mapping[tc] = rev_transition_mapping[(sA, i, sB)]
            rev_transition_index_mapping[rev_transition_mapping[(sA, i, sB)]] = tc
            tc+=1

        sc = 0
        for state in usable_states:
            state_mapping[state] = sc
            rev_state_mapping[sc] = state

            inputs = incoming_transitions[state]
            output = outgoing_transitions[state]

            for i in inputs:
                r = transition_index_mapping[i]
                A_init[sc][r]+=1

            for i in output:
                r = transition_index_mapping[i]
                A_init[sc][r]-=1

            sc+=1

        # create x for Ax=b
        model = Model('GAMING')
        model.setParam('OutputFlag', 0)
        x = model.addMVar(len(ts), vtype=GRB.INTEGER)
        model.setObjective(quicksum(x), GRB.MINIMIZE)

        """
        okay, now we need to:
        1. loop over all usable acceptance states
        2. alter A_init by deleting a row (after copying)
        3. constrain the right variables in x based on the inputs and outputs of the acceptance state
        4. throw things into the solver, and solve
        5. fill in some dict
        6. remove constraints, go again
        """
        # maps each state to the distance away from an accepting state, plus the accepting cycle
        state_distance_mapping = {} 
    
        for state in acceptance_states_usable:
            A_copy = np.delete(A_init, state_mapping[state], 0)
            b_copy = np.delete(b_init, state_mapping[state], 0)

            inputs = incoming_transitions[state]
            outputs = outgoing_transitions[state]
            # a reminder that these get *labels*

            c1 = model.addConstr(quicksum(x[transition_index_mapping[i]] for i in inputs) == 1,name="main1")
            c2 = model.addConstr(quicksum(x[transition_index_mapping[o]] for o in outputs) == 1,name="main2")

            c3 = model.addMConstr(A_copy, x, '=', b_copy, name="main3")
            model.update()
            model.optimize()
            
            if model.status == GRB.INFEASIBLE:
                pass
                # state_distance_mapping[state] = float('inf')
            elif model.status == GRB.OPTIMAL:
                state_distance_mapping[state] = sum(x.X)

            model.remove(model.getConstrByName("main1"))
            model.remove(model.getConstrByName("main2"))
            for constr in c3:
                model.remove(constr)

        model.update()

        # assert things are the same (I don't trust numpy)
        num_states_a, num_vars_A = A_init.shape
        num_states_b, num_vars_b = b_init.shape
        assert num_vars_A == len(ts)
        assert num_states_a == len(usable_states)
        assert num_states_b == len(usable_states)
        assert num_vars_b == 1

        """
        2.5. simultanious BFS: establish distance between any state and the acceptance cycle
            2.5.1. start by assigning all acceptance states a distance value equal to the minimum acceptance cycle distance
            2.5.2. from each acceptance state, calculate the distance in waves, e.g., distance[new_state] = distance[state] + 1
            2.5.3. we only throw a new state on the stack if it is at least two less than the new state
            2.5.4. repeat until stack is empty (need to prove things about this algorithm...)
        """

        qx = []
        for state in state_distance_mapping.keys() : qx.append(state)

        inf_states = usable_states.difference(set(state_distance_mapping.keys()))
        for state in inf_states : state_distance_mapping[state] = float('inf')

        while qx:
            curr_state = qx.pop(0)
            curr_dist = state_distance_mapping[curr_state]
            rev_reachable = buchi.reachability_object_reverse[curr_state]

            for (r_state, char) in rev_reachable:
                r_dist = state_distance_mapping[r_state]

                if curr_dist + 1 < r_dist:
                    state_distance_mapping[r_state] = curr_dist + 1
                    qx.append(r_state)
                    
        all_incoming_transitions[buchi] = incoming_transitions
        all_outgoing_transitions[buchi] = outgoing_transitions
        all_transition_mappings[buchi] = transition_mapping
        all_rev_transition_mappings[buchi] = rev_transition_mapping
        all_transition_index_mappings[buchi] = transition_index_mapping
        all_rev_transition_index_mappings[buchi] = rev_transition_index_mapping
        all_state_mappings[buchi] = state_mapping
        all_rev_state_mappings[buchi] = rev_state_mapping
        all_state_distance_mappings[buchi] = state_distance_mapping
        all_A[buchi] = A_init
        all_b[buchi] = b_init
        all_x[buchi] = x
        all_models[buchi] = model
        all_hit_acceptance_states[buchi] = set()

        processed_history[buchi] = {}

    def heuristic(state):
        ret_sum = 0
        buchi_state_label, prev_states, acc_states = state
        buchi_current_states = rev_cstr_state_mapping[buchi_state_label]

        def automaton_calculations(i, buchi_current_states, ret_sum):
            buchi = buchi_array[i]
            init_state = buchi_current_states[i]

            # getting info from all the dicts...
            usable_states = all_wellconnected_states[buchi]
            incoming_transitions = all_incoming_transitions[buchi]
            outgoing_transitions = all_outgoing_transitions[buchi]
            transition_mapping = all_transition_mappings[buchi]
            rev_transition_mapping = all_rev_transition_mappings[buchi]
            transition_index_mapping = all_transition_index_mappings[buchi]
            rev_transition_index_mapping = all_rev_transition_index_mappings[buchi]
            state_mapping = all_state_mappings[buchi]
            rev_state_mapping = all_rev_state_mappings[buchi]
            state_distance_mapping = all_state_distance_mappings[buchi]
            A_init = all_A[buchi]
            b_init = all_b[buchi]
            x = all_x[buchi]
            model = all_models[buchi]

            hit_accepting_states = all_hit_acceptance_states[buchi]

            # checking history
            hist = processed_history[buchi]
            if (init_state, tuple(hit_accepting_states)) in set(hist.keys()): 
                ret_sum += hist[(init_state, tuple(hit_accepting_states))]
                return ret_sum

            """
            3. define heuristic
                3.1. when we move through each automata we're intersecting, we keep track of acceptance states we've seen
                3.2. for all acceptance states we've visited, calculate the parikh image flow constraint-based distance using the matrix representation from the pre-processing
                    3.2.note. we could use big M to turn the much lesser "or" here into a more efficient, single operation. but, given it's all over matrices it should be ok in practice.
                    3.2.1. 
                3.3. to get the minimum distance, calculate minimum(all parikh image distances, distance[state])
            """

            curr_distance = state_distance_mapping[init_state]

            if len(hit_accepting_states) == 0:
                ret_sum+=curr_distance
                hist[(init_state, tuple(hit_accepting_states))] = curr_distance
                return ret_sum

            s = Optimize()
            for t in list(transition_mapping.keys()):
                tosolve[t] = Int(str(t))

            for t in tosolve.values():
                s.add(t >= 0)

            s.minimize(Sum(list(tosolve.values())))

            accepting_states_to_solve = set()

            for state in usable_states:
                if state in hit_accepting_states:
                    accepting_states_to_solve.add(state)
                elif state == init_state:
                    continue
                else: 
                    s.add(Sum([tosolve[t] for t in incoming_transitions[state]]) == Sum([tosolve[t] for t in outgoing_transitions[state]]))

            accept_or = []
            for accepting_state in accepting_states_to_solve:
                accept_or.append(
                        And(
                        (Sum([tosolve[t] for t in incoming_transitions[init_state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[init_state]])),
                        (Sum([tosolve[t] for t in incoming_transitions[accepting_state]])  == Sum([tosolve[t] for t in outgoing_transitions[accepting_state]]) + 1)
                        ))

            s.add(Or(*accept_or))

            if s.check() == sat:
                print("hit!")
                m = s.model()
                solns = {var: m[var] for var in tosolve.values()}

                total_dist = 0
                for transition, value in solns.items():
                    # sA, i, sB = transition_mapping[str(transition)]
                    total_dist+=int(str(value)) # kinda slow, can be optimized 

                hist[(init_state, tuple(hit_accepting_states))] = total_dist

                ret_sum += total_dist

            else : return float('inf')

            return ret_sum

        threads = []
        results = [0] * len(buchi_array)  # To store results from threads

        for i in range(len(buchi_array)):
            # Create a thread for the calculation on the i-th automaton
            t = Thread(target=automaton_calculations, args=(i, buchi_current_states, results))
            threads.append(t)
            t.start()

        # Wait for all threads to finish
        for t in threads:
            t.join()

        return sum(results)

    print("running begin")

    while queue:
        current_distance, current_state_tuple = heapq.heappop(queue)
        current_state, prev_states, acc_states = current_state_tuple 

        # termination condition
        if (current_state in prev_states and
            (current_state in acc_states 
             or 
             any([ current_state in acceptance_state_leading[astate] for astate in acc_states])
             )
            ):

            word = word[current_state]

            print("counterexample: " + str(word))
            total = len(cstr_state_mapping)
            
            total_possible = prod([ len(buchi.states) for buchi in buchi_array ])
            print("percentage of states visited: ", round(total/total_possible*100, 3), "% -",str(total), "/", str(total_possible))

            return word

        elif current_state in prev_states : continue

        # skip condition : if our current distance is longer than established distance
        if current_distance > distance[current_state] : continue
        
        current_state_ite = rev_cstr_state_mapping[current_state]

        # assemble new states
        reachable = []
        for char in fa:
            new_states, break_loop = [], False

            new_state = [None] * len(buchi_array)
            for i in range(len(buchi_array)):
                new_state[i] = buchi_array[i].transition_function(current_state_ite[i], char)

                if new_state[i] in buchi_array[i].acceptance_states : all_hit_acceptance_states[buchi].add(new_state[i])

                # if this is ever the case, we know the transition function led to nothing
                # thus this state is invalid
                if new_state[i] == []:
                    break_loop = True
                    break

            if break_loop : continue

            new_state = tuple(new_state)
            l = list(itertools.product(*new_state))

            for new_state in l:

                # if this is the case, we know we've already seen it before
                new_state_label = None
                if not new_state in cstr_state_mapping:
                    # print("new state!")
                    new_state_label = len(cstr_state_mapping)

                    cstr_state_mapping[new_state] = new_state_label
                    rev_cstr_state_mapping[new_state_label] = new_state

                    distance[new_state_label] = float('inf')
                    word[new_state_label] = ""

                else: 
                    new_state_label = cstr_state_mapping[new_state]

                new_prev_states = prev_states.copy()
                new_prev_states.add(current_state)

                new_acc_states = acc_states.copy()

                all_accepting = True
                for i in range(len(buchi_array)):
                    if not new_state[i] in buchi_array[i].acceptance_states:
                        all_accepting = False
                        break

                if all_accepting:
                    # we have an accepting state (no way bro)
                    new_acc_states.add(new_state_label)
                    acceptance_state_leading[new_state_label] = new_prev_states.copy()

                t = (new_state_label, new_prev_states, new_acc_states)

                reachable.append((t, new_state, char))

        for (state_t, state, char) in reachable:
            state_label = state_t[0]

            hval = heuristic(state_t)

            # new_distance = distance[current_state] + 1 + hval
            new_distance = hval 

            distance[state_label] = new_distance
            word[state_label] = word[current_state] + char
            heapq.heappush(queue, (new_distance, state_t))

    return None
