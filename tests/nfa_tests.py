from automata.nfa import NFA
from automata.dfa import DFA
from algorithms.intersection_emptiness import *

nfa1_transitions = [
    ('s1', 'a', 's2'),
    ('s2', 'b', 's4'),
    ('s4', 'c', 's3'),
    ('s2', 'c', 's3')
]

nfa1 = NFA(["s1", "s2", "s3", "s4"], "s1", ["s3"], ["a","b","c"], nfa1_transitions) 

assert nfa1.test_word(["a","c"],"s1") == {'s3'}
assert nfa1.test_word(["a"], "s1") == set()
assert nfa1.test_word(["a","b"], "s1") == set()
assert nfa1.test_word(["a","b","c"], "s1") == {'s3'}
assert nfa1.test_word(["c"], "s2") == {'s3'}

assert nfa1.reachable_from("s1") == {'s2', 's1', 's4', 's3'}
assert nfa1.reachable_from("s2") == {'s2', 's3', 's4'}
assert nfa1.reachable_from("s3") == {'s3'}
assert nfa1.reachable_from("s4") == {'s3', 's4'}

assert nfa1.dijkstra_shortest_distance("s1","s3") == (2, 'ac')

assert dijkstra_nfa_intersection_emptiness([nfa1, nfa1]) != None

nfa2_transitions = [
    ('s1','a','s2'),
    ('s1','a','s3'),
    ('s3','b','s2'),
    ('s2','c','s4'),
    ('s4','a','s3')
]

nfa2 = NFA(["s1", "s2", "s3", "s4"], "s1", ["s4"], ["a","b","c"], nfa2_transitions) 
assert dijkstra_nfa_intersection_emptiness([nfa1, nfa2]) != None
assert z3_dijkstra_nfa_intersection_emptiness([nfa1, nfa2]) != None

nfa3 = NFA(["s1", "s2", "s3", "s4"], "s1", [], ["a","b","c"], nfa2_transitions) 
assert dijkstra_nfa_intersection_emptiness([nfa1, nfa3]) == None

nfa4 = NFA(["s1", "s2", "s3", "s4"], "s1", ["s2","s3"], ["a","b","c"], nfa2_transitions) 
assert dijkstra_nfa_intersection_emptiness([nfa4.to_single_acceptance_state(), nfa4]) != None
assert z3_dijkstra_nfa_intersection_emptiness([nfa4.to_single_acceptance_state(), nfa4]) != None

nfa5_transitions = [
    ('s1', 'a', 's2'),
    ('s2', 'b', 's3'),
    ('s1', 'c', 's4'), 
    ('s4', 'd', 's5'),
    ('s5', 'e', 's3')
]

nfa5 = NFA(["s1", "s2", "s3", "s4", "s5"], "s1", ["s3"], ["a","b","c","d","e"], nfa5_transitions)
assert nfa5.minimal_parikh_image() == {'a': 1, 'b': 1, 'c': 0, 'd': 0, 'e': 0, '': 0}

nfa6_transitions = [
    ('s1', 'c', 's4'),
    ('s4', 'd', 's2'),
    ('s2', 'e', 's4'),
    ('s1', 'a', 's2'),
    ('s2', 'b', 's3'),
    ('s4', 'f', 's3')
]

nfa6 = NFA(["s1", "s2", "s3", "s4"], "s1", ["s3"], ["a","b","c","d","e","f"], nfa6_transitions)
nfa6_t = nfa6.minimal_parikh_image()
assert nfa6_t == {'a': 1, 'b': 1, 'c': 0, 'd': 0, 'e': 0, 'f': 0, '': 0} or nfa6_t == {'a': 0, 'b': 0, 'c': 1, 'd': 0, 'e': 0, 'f': 1, '': 0}

nfa8_transitions = [
    ('s1', 'c', 's4'),
    ('s4', 'd', 's2'),
    ('s2', 'e', 's4'),
    ('s1', 'a', 's2'),
    ('s2', 'b', 's3'),
    ('s4', 'f', 's3'),
    ('s7', 'd', 's8'),
    ('s8', 'd', 's7'),
    ('s1', 'e', 's7'),
    ('s2', 'e', 's7')
]

nfa8 = NFA(["s1", "s2", "s3", "s4", "s7", "s8"], "s1", ["s3"], ["a","b","c","d","e","f"], nfa8_transitions)
nfa8_t = nfa8.minimal_parikh_image()
assert nfa8_t == {'a': 1, 'b': 1, 'c': 0, 'd': 0, 'e': 0, 'f': 0, '': 0} or nfa8_t == {'a': 0, 'b': 0, 'c': 1, 'd': 0, 'e': 0, 'f': 1, '': 0}
