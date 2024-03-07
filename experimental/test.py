from utils import *
from algorithms.intersection_emptiness import *
b = generate_random_fa(5,1,3,10)
assert dijkstra_buchi_intersection_emptiness([b.complement(), b]) == None
