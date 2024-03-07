from algorithms.intersection_emptiness import *
from automata.buchi import Buchi

buchi1_transitions = [
    ('s1', 'a', 's2'),
    ('s2', 'b', 's4'),
    ('s4', 'c', 's3'),
    ('s3', 'c', 's2')
]

buchi1 = Buchi(["s1", "s2", "s3", "s4"], "s1", ["s3"], ["a","b","c"], buchi1_transitions)

buchi2_transitions = [
    ('s1','a','s2'),
    ('s1','a','s3'),
    ('s3','b','s2'),
    ('s2','c','s4'),
    ('s4','a','s3')
]

buchi2 = Buchi(["s1", "s2", "s3", "s4"], "s1", ["s4"], ["a","b","c"], buchi2_transitions)

# assert dijkstra_buchi_intersection_emptiness([buchi1,buchi2]) == None
assert dijkstra_buchi_intersection_emptiness([buchi2, buchi2]) != None

buchi3_transitions = [
    ('s1','a','s2'),
    ('s2','b','s2')
]

buchi3 = Buchi(["s1","s2"],"s1",["s2"], ["a","b"], buchi3_transitions)
assert dijkstra_buchi_intersection_emptiness([buchi3, buchi3]) != None

buchi_4_transitions = [
    ("s1", "a", "s2"),
    ("s2", "b", "s3"), 
    ("s1", "a", "s4")
]

buchi_4_states = ["s1", "s2", "s3", "s4"]
buchi_4_alphabet = ["a", "b"]
buchi_4_initial_state = "s1"
buchi_4_accepting_states = ["s3"]
buchi_4 = Buchi(buchi_4_states, buchi_4_initial_state, buchi_4_accepting_states, buchi_4_alphabet, buchi_4_transitions)

assert buchi_4.get_wellconnected_states() == {"s1", "s2", "s3"}


buchi_5_transitions = [
    ("s1", "a", "s2"),
    ("s2", "b", "s3"), 
    ("s1", "a", "s4"), 
    ("s3", "b", "s5")
]

buchi_5_states = ["s1", "s2", "s3", "s4", "s5"]
buchi_5_alphabet = ["a", "b"]
buchi_5_initial_state = "s1"
buchi_5_accepting_states = ["s3"]
buchi_5 = Buchi(buchi_5_states, buchi_5_initial_state, buchi_5_accepting_states, buchi_5_alphabet, buchi_5_transitions)

assert buchi_5.get_wellconnected_states() == {"s1", "s2", "s3"}

buchi_6_transitions = [
    ("s1", "a", "s2"),
    ("s2", "b", "s3"), 
    ("s1", "a", "s4"), 
    ("s3", "b", "s5")
]

buchi_6_states = ["s1", "s2", "s3", "s4", "s5"]
buchi_6_alphabet = ["a", "b"]
buchi_6_initial_state = "s4"
buchi_6_accepting_states = ["s3"]
buchi_6 = Buchi(buchi_6_states, buchi_6_initial_state, buchi_6_accepting_states, buchi_6_alphabet, buchi_6_transitions)

assert buchi_6.get_wellconnected_states() == set()

buchi_7_transitions = [
    ("s1", "a", "s2"),
    ("s2", "b", "s3"), 
    ("s3", "a", "s1"), 
]

buchi_7_states = ["s1", "s2", "s3"]
buchi_7_alphabet = ["a", "b"]
buchi_7_initial_state = "s1"
buchi_7_accepting_states = ["s2"]
buchi_7 = Buchi(buchi_7_states, buchi_7_initial_state, buchi_7_accepting_states, buchi_7_alphabet, buchi_7_transitions)

assert buchi_7.get_wellconnected_states() == {'s3', 's1', 's2'}

buchi_8_states = ["s1", "s2", "s3", "s4"]

buchi_8_transitions = [
    ("s1", "a", "s2"),
    ("s2", "b", "s3"),
    ("s3", "c", "s4"),
    ("s4", "d", "s3")
]

buchi_8_alphabet = ["a", "b", "c", "d"]
buchi_8_initial_state = "s1"
buchi_8_accepting_states = ["s4"]

buchi_8 = Buchi(buchi_8_states, buchi_8_initial_state, buchi_8_accepting_states, buchi_8_alphabet, buchi_8_transitions)
assert buchi_8.minimal_parikh_image() == {'a': 1, 'b': 1, 'c': 2, 'd': 1, '': 0}
assert buchi_8.parikh_image_inclusion(buchi_8.minimal_parikh_image())
assert not buchi_8.parikh_image_inclusion({'a': 1, 'b': 1, 'c': 1, 'd': 0, '': 0})

buchi_9_states = ["s1", "s2", "s3", "s4"]

buchi_9_transitions = [
    ("s1", "a", "s2"),
    ("s2", "b", "s3"),
    ("s3", "c", "s4"),
    ("s4", "d", "s3")
]

buchi_9_alphabet = ["a", "b", "c", "d"]
buchi_9_initial_state = "s1"
buchi_9_accepting_states = ["s3"]

buchi_9 = Buchi(buchi_9_states, buchi_9_initial_state, buchi_9_accepting_states, buchi_9_alphabet, buchi_9_transitions)
assert buchi_9.minimal_parikh_image() == {'a': 1, 'b': 1, 'c': 1, 'd': 1, '': 0}

