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
        self._octagon_dict = dict()


    def find_octagons(self, crystal):
        crystal_weight = crystal.module_generators[0].weight()
        if crystal_weight not in self._crystal_weight_dict:
            self.find_repeated_weights(crystal)

        self._octagon_dict[crystal_weight] = dict()
        for v in crystal:
            wt = v.weight()
            if wt in self._crystal_weight_dict[crystal_weight]:
                for w in self._crystal_weight_dict[crystal_weight][wt]:
                    if w != v:
                        for i in range(self._n):
                            for j in range(self._n):
                                if v.epsilon(i) != 0 and v.phi(i) != 0 and w.epsilon(j) != 0 and w.phi(j) != 0:
                                    post_v = v.e_string([i]*v.epsilon(i))
                                    pre_v = v.f_string([i]*v.phi(i))
                                    #post_w = w.e_string([j]*w.epsilon(j))
                                    pre_w = w.f_string([j]*w.phi(j))
                                    if pre_w.phi(i) != 0 and pre_w.f(i) == pre_v.f(j):
                                        #after jumping, it is not necessary to follow the steps which are used show existence of an octagon
                                        for k in range(self._n):
                                            if post_v.epsilon(k) != 0 and w.epsilon(k) != 0:
                                                coeff = post_v.epsilon(k) * self._R.gens()[2*self._n+i] * self._R.gens()[2*self._n+k]
                                                post_w = w.e_string([k]*w.epsilon(k))
                                                if (pre_v, '+', i) in self._octagon_dict[crystal_weight]:
                                                    self._octagon_dict[crystal_weight][(pre_v, '+', i)].append((k, post_w, coeff))
                                                else:
                                                    self._octagon_dict[crystal_weight][(pre_v, '+', i)] = [(k, post_w, coeff)]
                                                coeff = self._R.gens()[self._n+i] * self._R.gens()[self._n+k] #probably needs a scalar coefficient in general
                                                if (post_w, 'i', k) in self._octagon_dict[crystal_weight]:
                                                    self._octagon_dict[crystal_weight][(post_w, '-', k)].append((i, pre_v, coeff))
                                                else:
                                                    self._octagon_dict[crystal_weight][(post_w, '-', k)] = [(i, pre_v, coeff)]
                                    #if post_w.epsilon(i) != 0 and post_w.e(i) == post_v.e(j):
                                    #    coeff = pre_v.phi(j) * self._R.gens()[self._n+i] * self._R.gens()[self._n+j]
                                    #    if (post_v, '-', i) in self._octagon_dict[crystal_weight]:
                                    #        self._octagon_dict[crystal_weight][(post_v, '-', i)].append((j, pre_w, coeff))
                                    #    else:
                                    #        self._octagon_dict[crystal_weight][(post_v, '-', i)] = [(j, pre_w, coeff)]


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



    def _recursive_path_graph(self, g_pos, ls_path, crystal, G = None, jump_history = dict(), path_so_far = []):

        if not G:
            G = DiGraph(weighted = True, loops = True)

        if g_pos == -1:
            if ls_path == crystal.module_generators[0]:
                for e in path_so_far:
#                    if e[2][3] == 'jump' and (e[2][0] == '+' and (e[1],e[0],('-',e[2][2],e[2][1],'jump',e[2][4])) not in path_so_far or e[2][0] == '-' and (e[1],e[0],('+',e[2][2],e[2][1],'jump',e[2][4])) not in path_so_far):
#                        break
                    G.add_edge(e)
                self._admissible_paths_dict[crystal.module_generators[0].weight()].append(path_so_far)
            return G

        crystal_weight = crystal.module_generators[0].weight()
        if crystal_weight not in self._crystal_weight_dict:
            self.find_repeated_weights(crystal)

        if crystal_weight not in self._octagon_dict:
            self.find_octagons(crystal)

        cpt , i = self._g[g_pos]
        wt = ls_path.weight()

        if cpt == 'h':
            monomial = prod([ h**a for (h,a) in zip(self._R.gens()[:self._n],vector(wt))])
            return self._recursive_path_graph(g_pos-1, ls_path, crystal, G=G, jump_history=jump_history, path_so_far=path_so_far+[(ls_path,ls_path,('','','','h',monomial))])
        elif cpt == '-':
            step = lambda v,i: v.f(i)
            step_string = lambda v,l: v.f_string(l)
            k = lambda v,i: v.phi(i)
            l = lambda v,i: v.epsilon(i)
            gen = self._R.gens()[self._n+i]
        else:
            step = lambda v,i: v.e(i)
            step_string = lambda v,l: v.e_string(l)
            k = lambda v,i: v.epsilon(i)
            l = lambda v,i: v.phi(i)
            gen = self._R.gens()[2*self._n+i]

        for r in range(k(ls_path, i) + 1):
            new_ls_path = step_string(ls_path, [i]*r)
            if r == 0:
                new_edge = []
            else:
                monomial = binomial( l(ls_path, i) + r, r) * gen**r
                new_edge = [(ls_path, new_ls_path, (cpt, i, r, '', monomial))]
            G = self._recursive_path_graph(g_pos-1, new_ls_path, crystal, G=G, jump_history=jump_history, path_so_far=path_so_far+new_edge)

        if (ls_path, cpt, i) in self._octagon_dict[crystal_weight]:
            for j, new_ls_path, coeff in self._octagon_dict[crystal_weight][(ls_path, cpt, i)]:
                for p in range(self._n):
                    if g_pos-1-p > 0 and (cpt, j) == self._g[g_pos-1-p]:
                        new_edge = [(ls_path, new_ls_path, (cpt, i, j, 'jump', coeff))]
                        G = self._recursive_path_graph(g_pos-2-p, new_ls_path, crystal, G=G, jump_history=jump_history, path_so_far=path_so_far+new_edge)

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

    def _recursive_evaluation(self, g, ls_path, crystal, jump_history=dict()):

        if not g:
            if ls_path == crystal.module_generators[0]:
                return 1
            else:
                return 0

        crystal_weight = crystal.module_generators[0].weight()
        if not self._crystal_weight_dict[crystal_weight]:
            self._crystal_weight_dict[crystal_weight] = dict()
            for v in crystal:
                wt = v.weight()
                if wt in self._crystal_weight_dict[crystal_weight]:
                    self._crystal_weight_dict[crystal_weight][wt].append(v)
                else:
                    self._crystal_weight_dict[crystal_weight][wt] = [v]

        new_g = copy(g)
        cpt , i = new_g.pop()
        wt = ls_path.weight()

        if cpt == 'h':
            return self._recursive_evaluation(new_g, copy(ls_path), crystal, jump_history=jump_history) * prod([ h**a for (h,a) in zip(self._R.gens()[:self._n],vector(wt))])
        elif cpt == '-':
            step = lambda v,i: v.f(i)
            step_op = lambda v,i: v.e(i)
            step_string = lambda v,l: v.f_string(l)
            step_op_string = lambda v,l: v.e_string(l)
            k = lambda v,i: v.phi(i)
            l = lambda v,i: v.epsilon(i)
            gen = self._R.gens()[self._n+i]
        else:
            step = lambda v,i: v.e(i)
            step_op = lambda v,i: v.f(i)
            step_string = lambda v,l: v.e_string(l)
            step_op_string = lambda v,l: v.f_string(l)
            k = lambda v,i: v.epsilon(i)
            l = lambda v,i: v.phi(i)
            gen = self._R.gens()[2*self._n+i]

        jump_list = [ls_path]
        steps_list1 = [ j for j in range(self._n) if k(ls_path,j) != 0 ]
        op_steps_list = [ j for j in range(self._n) if l(ls_path,j) != 0 ]
        for v in self._crystal_weight_dict[crystal_weight][wt]:
            if v != ls_path:
                steps_list2 = [ j for j in range(self._n) if k(v,j) != 0 ]
                for r in steps_list1:
                    for s in steps_list2:
                        w1 = step_string(v,[s]*k(v,s))
                        w2 = step_string(ls_path,[r]*k(ls_path,r))
                        if k(w2,s) != 0 and step(w2,s) == step(w1,r) and v not in jump_list:
                            jump_list.append(v)
                for j in op_steps_list:
                    w1 = step_op_string(v,[i]*l(v,i))
                    w2 = step_op_string(ls_path,[j]*l(ls_path,j))
                    if l(w2,i) != 0 and step_op(w2,i) == step_op(w1,j) and v not in jump_list:
                        jump_list.append(v)
        if ls_path in jump_history:
            jump_list = [ls_path]
            for v in jump_history[ls_path]:
                if v not in jump_list:
                    jump_list.append(v)

        #jump_list = jump_dict[wt]

        output = 0

        for v in jump_list:
            new_ls_path = copy(v)
            current_k = k(v,i)
            current_l = l(v,i)
            for j in range(current_k+1):
                if v == ls_path:
                    output += self._recursive_evaluation(new_g, new_ls_path, crystal, jump_history=jump_history) * gen**j * binomial(current_l+j,current_l)
                elif j > 0:
                    new_jump_history = copy(jump_history)
                    if v not in new_jump_history:
                        new_jump_history[v] = [ls_path]
                    else:
                        new_jump_history[v].append(ls_path)
                    old_output = output
                    output += self._recursive_evaluation(new_g, new_ls_path, crystal, jump_history=new_jump_history) * gen**j * binomial(current_l+j,current_l)

                if j < current_k:
                    new_ls_path = step(new_ls_path, i)

        return output


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
    if all(test_conjecture_on_matrix(B, coxeter=c, mutation_type=qm) for (B,c) in to_check):
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

    # find regular_g_vectors
    #tubes = flatten(T.affine_tubes())
    #regular_g_vectors = map(lambda x: tuple(vector(T.to_weight(x))),tubes)
    A.explore_to_depth(10)
    regular_g_vectors = A.g_vectors_so_far()

    # blackbox to produce minors
    R = DoubleBruhatAlgebra(b_matrix)

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

