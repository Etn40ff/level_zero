# At the moment this file is just a recording of an experimental session to compute generalized minors

# First of all attach test_conjecture.py

# Find some g-vectors you want to test, for example this should give a fail with two g-vectors: (-1, -1, 1, -1, 2) and (-1, -1, 0, 1, 0)
# test_conjecture_on_type(['B',4,1])

# Configuration
C = CartanMatrix(['B',4,1])
c = (0, 1, 2, 3, 4)
g_vect = (-1, -1, 1, -1, 2)
#g_vect = (-1, -1, 0, 1, 0)

# Setup
B = b_matrix(C, c)
Bdp = block_matrix([[B],[identity_matrix(B.ncols())],[identity_matrix(B.ncols())]])
A = ClusterAlgebra(Bdp)
A.find_g_vector(g_vect)
R = RootSystem(C.cartan_type())
La = R.weight_space(extended = True).basis()
V = crystals.LSPaths(sum([x*y for x,y in zip(g_vect,La)]))


# truncate to a sufficiently large depth
W = V.subcrystal(max_depth=10)

# plot
# The chosen sidest weight is highlited
# arrows are colored with colors from rainbow(len( W.index_set())) in order
# x_i correspond to the incoming arrow of colour rainbow[i]
# x_{\overline{i}} is outgoing with the same convention
G = W.digraph()
G.show(method='js',vertex_labels=False, edge_partition=[ [e for e in G.edges() if e[2] == i] for i in W.index_set()],charge=-500, vertex_partition=[W.module_generators])

