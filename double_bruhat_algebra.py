from __future__ import print_function
from tropical_cluster_algebra import TropicalClusterAlgebra
#from cluster_algebra import ClusterAlgebra
from root_system import RootSystem as RootSystemDylan

class DoubleBruhatAlgebra(SageObject):

    def __init__(self, data):

        Q = ClusterQuiver(data)
        self._n = Q.n()
        self._B = Q.b_matrix()

        self._c = coxeter_element(self._B)
        self._A = CartanMatrix(2-self._B.apply_map(abs))

        self._RS = RootSystem(self._A.cartan_type())

        self._La = self._RS.weight_space().basis()

        # generic element
        self._g = [ ('-', i) for i in self._c ]
        self._g.append(('h',None))
        self._g += [ ('+', i) for i in reversed(self._c) ]

        # ambient ring
        self._R = PolynomialRing(QQ, sorted(flatten([('h%s'%i,'r%s'%i,'t%s'%i) for i in xrange(self._n)])))

        self._crystal_weight_dict = dict()
        self._admissible_paths_dict = dict()
        self._step_dict = dict()


    def find_steps(self, crystal):
        crystal_weight = crystal.module_generators[0].weight()
        if crystal_weight not in self._crystal_weight_dict:
            self.find_repeated_weights(crystal)

        if crystal_weight not in self._step_dict:
            self._step_dict[crystal_weight] = dict()
        for v in crystal:
            #compute natural steps inside the crystal
            for i in range(self._n):
                if (v, '+', i) not in self._step_dict[crystal_weight]:
                    self._step_dict[crystal_weight][(v, '+', i)] = dict()
                for r in range(v.epsilon(i) + 1):
                    next_node = v.e_string([i]*r)
                    step_factor = binomial(v.phi(i)+r,r) * self._R.gens()[2*self._n+i]**r
                    if (i, r) not in self._step_dict[crystal_weight][(v, '+', i)]:
                        self._step_dict[crystal_weight][(v, '+', i)][(i,r)] = [(next_node, step_factor)]
                    else:
                        self._step_dict[crystal_weight][(v, '+', i)][(i,r)].append((next_node, step_factor))
                if (v, '-', i) not in self._step_dict[crystal_weight]:
                    self._step_dict[crystal_weight][(v, '-', i)] = dict()
                for s in range(v.phi(i) + 1):
                    next_node = v.f_string([i]*s)
                    step_factor = binomial(v.epsilon(i)+s,s) * self._R.gens()[self._n+i]**s
                    if (i, s) not in self._step_dict[crystal_weight][(v, '-', i)]:
                        self._step_dict[crystal_weight][(v, '-', i)][(i,s)] = [(next_node, step_factor)]
                    else:
                        self._step_dict[crystal_weight][(v, '-', i)][(i,s)].append((next_node, step_factor))
            wt = v.weight()
            if wt in self._crystal_weight_dict[crystal_weight]:
                for w in self._crystal_weight_dict[crystal_weight][wt]:
                    if w != v:
                        #compute jump steps inside the crystal`
                        for i in range(self._n):
                            for j in range(self._n):
                                if v.epsilon(i) != 0 and v.phi(i) != 0 and w.epsilon(j) != 0 and w.phi(j) != 0:
                                    post_v = v.e_string([i]*v.epsilon(i))
                                    pre_v = v.f_string([i]*v.phi(i))
                                    post_w = w.e_string([j]*w.epsilon(j))
                                    pre_w = w.f_string([j]*w.phi(j))
                                    has_roof = False
                                    if  pre_w.phi(i) != 0 and pre_v.phi(j) != 0:
                                        coeff = pre_v.phi(j) * self._R.gens()[2*self._n+i] * self._R.gens()[2*self._n+j]
                                        for r in range(pre_w.phi(i)):
                                            for s in range(pre_v.phi(j)):
                                                if pre_w.f_string([i]*(r+1)) == pre_v.f_string([j]*(s+1)):
                                                    has_roof = True
                                                    if (pre_v, '+', i) not in self._step_dict[crystal_weight]:
                                                        self._step_dict[crystal_weight][(pre_v, '+', i)] = dict()
                                                    key = (j, w.epsilon(j))
                                                    if key not in self._step_dict[crystal_weight][(pre_v, '+', i)]:
                                                        self._step_dict[crystal_weight][(pre_v, '+', i)][key] = [(post_w, coeff)]
                                                    else:
                                                        self._step_dict[crystal_weight][(pre_v, '+', i)][key].append((post_w, coeff))
                                    if  post_w.epsilon(i) != 0 and post_v.epsilon(j) != 0:
                                        coeff = post_v.epsilon(j) * self._R.gens()[self._n+i] * self._R.gens()[self._n+j]
                                        for r in range(post_w.epsilon(i)):
                                            for s in range(post_v.epsilon(j)):
                                                if post_w.e_string([i]*(r+1)) == post_v.e_string([j]*(s+1)):
                                                    has_roof = True
                                                    if (post_v, '-', i) not in self._step_dict[crystal_weight]:
                                                        self._step_dict[crystal_weight][(post_v, '-', i)] = dict()
                                                    key = (j, w.phi(j))
                                                    if key not in self._step_dict[crystal_weight][(post_v, '-', i)]:
                                                        self._step_dict[crystal_weight][(post_v, '-', i)][key] = [(pre_w, coeff)]
                                                    else:
                                                        self._step_dict[crystal_weight][(post_v, '-', i)][key].append((pre_w, coeff))
                                    if has_roof:
                                        #parallel steps may be possible
                                        for k in range(self._n):
                                            if k != i and k != j:
                                                if v.phi(k) != 0 or w.epsilon(k) != 0:
                                                    for r in range(v.phi(k) + 1):
                                                        for s in range(w.epsilon(k) + 1):
                                                            if r + s != 0:
                                                                prev_node = v.f_string([k]*r)
                                                                next_node = w.e_string([k]*s)
                                                                coeff = binomial(prev_node.phi(k)+r+s, r+s) * self._R.gens()[2*self._n+k]**(r+s)
                                                                if (prev_node, '+', k) not in self._step_dict[crystal_weight]:
                                                                    self._step_dict[crystal_weight][(prev_node, '+', k)] = dict()
                                                                if (k, r+s) not in self._step_dict[crystal_weight][(prev_node, '+', k)]:
                                                                    self._step_dict[crystal_weight][(prev_node, '+', k)][(k, r+s)] = [(next_node, coeff)]
                                                                else:
                                                                    self._step_dict[crystal_weight][(prev_node, '+', k)][(k, r+s)].append((next_node, coeff))
                                                                coeff = binomial(next_node.epsilon(k)+r+s, r+s) * self._R.gens()[self._n+k]**(r+s)
                                                                if (next_node, '-', k) not in self._step_dict[crystal_weight]:
                                                                    self._step_dict[crystal_weight][(next_node, '-', k)] = dict()
                                                                if (k, r+s) not in self._step_dict[crystal_weight][(next_node, '-', k)]:
                                                                    self._step_dict[crystal_weight][(next_node, '-', k)][(k, r+s)] = [(prev_node, coeff)]
                                                                else:
                                                                    self._step_dict[crystal_weight][(next_node, '-', k)][(k, r+s)].append((prev_node, coeff))


    def find_repeated_weights(self, crystal):
        crystal_weight = crystal.module_generators[0].weight()
        self._crystal_weight_dict[crystal_weight] = dict()
        tmp_dict = dict()
        for v in crystal:
            wt = v.weight()
            if wt in tmp_dict:
                tmp_dict[wt].append(v)
            else:
                tmp_dict[wt] = [v]
        for wt in tmp_dict:
            if len(tmp_dict[wt]) > 1:
                self._crystal_weight_dict[crystal_weight][wt] = tmp_dict[wt]


    def _recursive_path_graph(self, g_pos, ls_path, crystal, G = None, path_so_far = []):

        if not G:
            G = DiGraph(weighted = True, loops = True)

        if g_pos == -1:
            if ls_path == crystal.module_generators[0]:
                for e in path_so_far:
                    G.add_edge(e)
                self._admissible_paths_dict[crystal.module_generators[0].weight()].append(path_so_far)
            return G

        crystal_weight = crystal.module_generators[0].weight()
        if crystal_weight not in self._step_dict:
            self.find_steps(crystal)

        cpt , i = self._g[g_pos]
        wt = ls_path.weight()

        if cpt == 'h':
            monomial = prod([ h**a for (h,a) in zip(self._R.gens()[:self._n],vector(wt))])
            new_edge = [(ls_path,ls_path,('','','','h',monomial))]
            return self._recursive_path_graph(g_pos-1, ls_path, crystal, G=G, path_so_far=path_so_far+new_edge)

        if (ls_path, cpt, i) in self._step_dict[crystal_weight]:
            for j, s in self._step_dict[crystal_weight][(ls_path, cpt, i)]:
                step_list = self._step_dict[crystal_weight][(ls_path, cpt, i)][(j, s)]
#               if cpt == '-':
#                   correction_factor = len(step_list)
#               else:
#                   correction_factor = 1
                correction_factor = 1 #len(step_list)
                for p in range(self._n):
                    if g_pos-p > -1 and (cpt, j) == self._g[g_pos-p]:
                        for new_ls_path, coeff in step_list:
                            new_edge = []
                            if j == i and s > 0:
                                new_edge = [(ls_path, new_ls_path, (cpt, i, s, 'step', coeff/correction_factor))]
                            elif s > 0:
                                new_edge = [(ls_path, new_ls_path, (cpt, i, j, 'jump', coeff/correction_factor))]
                            G = self._recursive_path_graph(g_pos-1-p, new_ls_path, crystal, G=G, path_so_far=path_so_far+new_edge)
        return G


    def path_graph(self, weight, to_show=True):
        # TODO: improve choice of max_depth
        if type(weight) == tuple:
            weight = sum([x*y for x,y in zip(weight,self._La)])
        V = crystals.LSPaths(weight).subcrystal(max_depth=self._n**4)
        if weight not in self._admissible_paths_dict:
            self._admissible_paths_dict[weight] = []
            G = self._recursive_path_graph(len(self._g)-1, V.module_generators[0], V)
        else:
            G = DiGraph(weighted = True, loops = True)
            for path in self._admissible_paths_dict[weight]:
                for e in path:
                    G.add_edge(e)

        if to_show:
            G.show(method='js', link_distance=300, vertex_labels=False, edge_labels=True, charge=-1000, vertex_partition=[V.module_generators])

        return G


    def generalized_minor(self, weight):
        if type(weight) == tuple:
            weight = sum([x*y for x,y in zip(weight,self._La)])
        if not weight in self._admissible_paths_dict:
            self.path_graph(weight, to_show=False)
        output = 0
        for path in self._admissible_paths_dict[weight]:
            term = 1
            for e in path:
                #print("edge=",e)
                term *= e[2][4]
            output += term
        return output


def test_conjecture_on_type(cartan_type):
    r"""
    Test our conjecture on all the acyclic matrices whose cartan companion is of cartan_type
    """
    A = CartanType(cartan_type).cartan_matrix()
    n = A.ncols()

    # HACK: type D is not well understood by QuiverMutationType so we pass it along if we can
    try:
        qm = QuiverMutationType(cartan_type)
    except:
        qm = None

    # list all matrices to test
    to_check = []
    for c in Permutations(range(n)):
        B = b_matrix(A,c)
        if B not in [ BB for (BB,cc) in to_check]:
            to_check.append((B,tuple(c)))

    # run test
    if all(test_conjecture_on_matrix(B, coxeter=c, mutation_type=qm, cartan_type=CartanType(cartan_type)) for (B,c) in to_check):
        return True

    return False

def test_conjecture_on_matrix(b_matrix, mutation_type=None, coxeter=None, cartan_type=None):
    r"""
    Test our conjecture on a single b_matrix
    """
    print(b_matrix)
    if not coxeter:
        coxeter = coxeter_element(b_matrix)

    # blackbox to produce the vectors we need
    T = TropicalClusterAlgebra(b_matrix,mutation_type=mutation_type)
    if cartan_type:
        T.cartan_type.set_cache(cartan_type)
    #if not (T.is_affine() and T.is_acyclic()) :
    #    raise NotImplementedError("This code works only with acyclic affine exchange matrices")

    print("################################################################################\n")
    print("Testing Dynkin diagram\n"+T.cartan_type().ascii_art()+"\ncoxeter = " + str(coxeter))
    print

    # blackbox to produce cluster variables
    n = b_matrix.ncols()
    Bdp = block_matrix([[b_matrix],[identity_matrix(n)],[identity_matrix(n)]])
    A = ClusterAlgebra(Bdp)

    # blackbox to produce minors
    R = DoubleBruhatAlgebra(b_matrix)

    # find regular_g_vectors
    if R._A.is_affine():
        tubes = flatten(T.affine_tubes())
        regular_g_vectors = map(lambda x: tuple(vector(T.to_weight(x))),tubes)
    elif R._A.is_finite():
        A.explore_to_depth(10)
        regular_g_vectors = A.g_vectors_so_far()

    # change of coordinates
    substitution = dict()
    for i in range(n):
        substitution[A.ambient().gens()[i]] = R._R.gen(i)
        substitution[A.ambient().gens()[i+n]] = R._R.gen(i+n)*R._R.gen(i)*prod([R._R.gen(j)**(-max(b_matrix[j,i],0)) for j in range(n)])
        substitution[A.ambient().gens()[i+2*n]] = R._R.gen(i+2*n)*R._R.gen(i)*prod([R._R.gen(j)**(-max(b_matrix[j,i],0)) for j in range(n)])

    got_problems = False
    for gvect in regular_g_vectors:
        print(str(gvect))
        #print(str(gvect) + "\tin the orbit of the finite type dominant weight " + str((0,)+tuple(R._level_zero_dominant_conjugate(R._g_to_weight(gvect)))[:b_matrix.ncols()-1]) + "\t: ", end="")
        A.find_g_vector(gvect)
        variable = A.cluster_variable(gvect).lift().subs(substitution)
        minor = R.generalized_minor(gvect)
        if variable == minor:
            print("Pass")
        else:
            print("Fail")
            print([m/variable.denominator() for m in variable.numerator().monomials() if m not in minor.numerator().monomials()])
            print([m/minor.denominator() for m in minor.numerator().monomials() if m not in variable.numerator().monomials()])
            print(minor-variable)
            print(minor)
            got_problems = True
    print
    return not got_problems

def b_matrix(A,c):
    r"""
    Produce the exchange matrix associated to the cartan matrix A and the coxeter element c
    """
    n = A.ncols()
    B = matrix(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if c.index(i) < c.index(j):
                B[i,j]=-A[i,j]
            else:
                B[i,j]=A[i,j]
    return B

def coxeter_element(b_matrix):
    r"""
    Returns a list expressing the coxeter element corresponding to b_matrix
    (twisted) reflections are applied from top of the list, for example
    [2, 1, 0] correspond to s_2s_1s_0
    Sources == non positive columns == leftmost letters
    """
    zero_vector = vector([0 for x in range(b_matrix.ncols())])
    coxeter = []
    columns = b_matrix.columns()
    source = None
    for j in range(b_matrix.ncols()):
        for i in range(b_matrix.ncols()):
            if all(x <=0 for x in columns[i]) and columns[i] != zero_vector:
                source = i
                break
        if source == None:
            if b_matrix != matrix(b_matrix.ncols()):
                raise ValueError("Unable to find a Coxeter element representing self._b_matrix")
            coxeter += [ x for x in range(b_matrix.ncols()) if x not in coxeter]
            break
        coxeter.append(source)
        columns[source] = zero_vector
        b_matrix = matrix(columns).transpose()
        b_matrix[source] = zero_vector
        columns = b_matrix.columns()
        source = None
    return tuple(coxeter)

