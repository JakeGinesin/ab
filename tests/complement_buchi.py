from utils import *
from algorithms.intersection_emptiness import *
b = generate_random_fa(5,1,3,10)
# print(dijkstra_buchi_intersection_emptiness([b.complement(), b]))
bcomp = b.complement()
if not dijkstra_buchi_intersection_emptiness([bcomp, b]) == None:
    print("failed")
    #b.show_diagram("b.png")
    #bcomp.show_diagram("bcomp.png")
    print(dijkstra_buchi_intersection_emptiness([bcomp, b]))
    #b.to_ba_file("b.ba")
    #bcomp.to_ba_file("bcomp.ba")
    

else : print("good")

