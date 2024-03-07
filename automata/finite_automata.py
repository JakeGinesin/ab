#!/usr/bin/env python3

from pydot import Dot, Edge, Node
from collections import deque
from typing import Set, Tuple, List 
import heapq

"""
PARENT CLASS FOR NFA, DFA, BUCHI

follows from the classical definition of a finite automata
which is a 5-tuple, (Q, Sigma, T, q0, F)
where T : Q x Sigma -> P(Q) (where P(Q) is the powerset of Q)
"""
class Finite_Automata(object):
    def __init__(
        self,
        states : set,
        initial_state : str, # this is really just a state label, not necessarily a state object
        acceptance_states : set,
        alphabet : set,
        transitions : List[Tuple[str, str, str]] = []):

        self.states = sorted(list(set(states)))
        self.initial_state = initial_state
        self.alphabet = sorted(list(set(alphabet)))
        self.acceptance_states = sorted(list(set(acceptance_states)))
        self.transitions = transitions

        assert self.initial_state in self.states, "initial state not in states"

        for state in acceptance_states: 
            assert(state in self.states)

    def construct_transition_function_hashmap(self) -> None:
        """ constructs a hashmap with all the transition relations, e.g. T : State -> Alphabet x State """
        self.transition_function_hashmap = {}

        for state in self.states:
            self.transition_function_hashmap[state] = {}
            for char in self.alphabet:
                self.transition_function_hashmap[state][char] = set()

        for (sA, i, sB) in self.transitions:
            assert(sA in self.states)
            assert(sB in self.states)
            assert(i in self.alphabet)
            self.transition_function_hashmap[sA][i].add(sB)

    def construct_reachability_hashmap(self) -> None:
        """ constructs a hashmap with all reachability relations, 
        e.g. what states are immediately reachable from the current state """
        self.reachability_hashmap = {}
        self.reachability_hashmap_reverse = {}

        for state in self.states:
            self.reachability_hashmap[state] = set()
            self.reachability_hashmap_reverse[state] = set()

        for (sA, i, sB) in self.transitions:
            assert(sA in self.states)
            assert(sB in self.states)
            assert(i in self.alphabet)
            self.reachability_hashmap[sA].add((sB, i))
            self.reachability_hashmap_reverse[sB].add((sA, i))

    def transition_function(self, state : str, char : str) -> set[str]:
        """ finite automata transition relation, e.g. state, char -> set of states"""
        if not hasattr(self, "transition_fnuction_hashmap") : self.construct_transition_function_hashmap()
        assert state in self.states, "transition_function: state called not in NFA"
        assert char in self.alphabet, "transition_function: char not in alphabet"
        return list(self.transition_function_hashmap[state][char])

    def show_diagram(self, path : str = "automata.png") -> None:
        """ creates a diagram of this finite Automata """
        graph = Dot(graph_type='digraph', rankdir='LR')
        nodes = {}

        for state in self.states:
            if state in self.acceptance_states:
                if state == self.initial_state:
                    nodes[state] = Node(state, shape='doublecircle', color='green')
                else:
                    nodes[state] = Node(state, shape='doublecircle')
            else:
                if state == self.initial_state:
                    nodes[state] = Node(state, shape='circle', color='green')
                else:
                    nodes[state] = Node(state, shape='circle')

            graph.add_node(nodes[state])

        for (sA, i, sB) in self.transitions:
            if i == '':
                graph.add_edge(Edge(nodes[sA], nodes[sB], label='ε'))
            else:
                graph.add_edge(Edge(nodes[sA], nodes[sB], label=i))

        graph.write_png(path)

    def reachable_from(self, state : str) -> set[str]:
        """ returns a set of all reachable states from the selected state """
        assert state in self.states, "all_reachable_from: state not in automata"
        if not hasattr(self, "reachability_hashmap") : self.construct_reachability_hashmap()
        queue = []
        queue.append(state)
        visited = set()
        while queue:
            current_state = queue.pop(0)
            reachable = self.reachability_hashmap[current_state]
            visited.add(current_state)
            for (state, char) in reachable:
                if state not in visited : queue.append(state)

        return visited

    def reachable_to(self, state : str) -> set[str]:
        """ returns a set of all states that eventually reach the current state """
        assert state in self.states, "all_reachable_from: state not in NFA"
        if not hasattr(self, "reachability_hashmap_reverse") : self.construct_reachability_hashmap()
        queue = []
        queue.append(state)
        visited = set()
        while queue:
            current_state = queue.pop(0)
            reachable = self.reachability_hashmap_reverse[current_state]
            visited.add(current_state)
            for (state, char) in reachable:
                if state not in visited : queue.append(state)

        return visited

    def get_wellconnected_states(self) -> set[str]:
        """ 
        returns a list of states that are reachable from the initial state intersected with
        the states that lead to an accepting state in linear time.

        in other words, returns the states that are relevant to the language of the automata
        """
        if not hasattr(self, "reachability_hashmap") : self.construct_reachability_hashmap()
        visited = set()
        prev = {state: set() for state in self.states}
        ret = set()

        queue = deque([self.initial_state])

        while queue:
            current_state = queue.popleft()
            reachable = self.reachability_hashmap[current_state]

            visited.add(current_state)

            for (state, char) in reachable:
                prev[state] |= prev[current_state]
                prev[state].add(current_state)

                if state not in visited and state not in queue:
                    queue.append(state)

                if state in self.acceptance_states or state in ret:
                    ret |= prev[state]
                    ret.add(state)

        return ret

    def dijkstra_shortest_distance(self, init_state : str, final_state : str) -> (int, str):
        """
        djikstra search, returning the shortest distance from the first state to the second state, 
        along with the associated word
        """
        assert init_state in self.states, "dijkstra: init_state not in NFA"
        assert final_state in self.states, "dijkstra: final_state not in NFA"

        distance = {}
        word = {}
        for state in self.states:
            distance[state] = float('inf')
            word[state] = ""

        queue = []
        heapq.heappush(queue, (0, init_state))
        distance[init_state] = 0

        while queue:
            current_distance, current_state = heapq.heappop(queue)

            if current_distance > distance[current_state]:
                continue

            reachable = self.reachability_hashmap[current_state]
            for (state, char) in reachable:
                new_distance = distance[current_state] + 1  

                # if the new distance is less than the stored distance, update
                if new_distance < distance[state]:
                    distance[state] = new_distance
                    word[state] = word[current_state] + char
                    heapq.heappush(queue, (new_distance, state))

        return distance[final_state], word[final_state]
