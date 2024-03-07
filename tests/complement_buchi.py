from utils import *
from algorithms.intersection_emptiness import *
b = generate_random_fa(10,1,3,20)
print(dijkstra_buchi_intersection_emptiness([b.complement(), b]))
assert dijkstra_buchi_intersection_emptiness([b.complement(), b]) == None
