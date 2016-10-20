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

        self._La = self._RS.weight_space(extended = True).basis()

        # generic element
        self._g = [ ('-', i) for i in self._c ]
        self._g.append(('h',None))
        self._g += [ ('+', i) for i in reversed(self._c) ]

        # ambient ring
        self._R = PolynomialRing(QQ, sorted(flatten([('h%s'%i,'r%s'%i,'t%s'%i) for i in xrange(self._n)])))

        self._jump_locations = []


    def _recursive_path_graph(self, g, ls_path, crystal, jump_dict=None, G = None):

        if not G:
            G = DiGraph(weighted=True, loops=True)

        if not g:
            return G

        if not jump_dict:
            jump_dict = dict()
            for v in crystal:
                wt = v.weight()
                if wt in jump_dict:
                    jump_dict[wt].append(v)
                else:
                    jump_dict[wt] = [v]

        new_g = copy(g)
        cpt , i = new_g.pop()
        wt = ls_path.weight()

        if cpt == 'h':
#            G.add_edge((ls_path,ls_path,'h'))
            return self._recursive_path_graph(new_g, ls_path, crystal, jump_dict=jump_dict, G=G)
        elif cpt == '-':
            step = lambda v,i: v.f(i)
            k = lambda v,i: v.phi(i)
        else:
            step = lambda v,i: v.e(i)
            k = lambda v,i: v.epsilon(i)

        jump_list = jump_dict[wt]

        for v in jump_list:
#            if v != ls_path:
#                G.add_edge((ls_path,v,'jump'))
            new_ls_path = copy(v)
            current_k = k(v,i)
            for j in range(current_k+1):
                if v == ls_path or j > 0:
                    G = self._recursive_path_graph(new_g, new_ls_path, crystal, jump_dict=jump_dict, G=G)
                if j < current_k:
                    new_ls_path = step(new_ls_path, i)
#                    if cpt == '+':
                    G.add_edge((v,new_ls_path,(cpt,i,j+1)))
        return G

    def path_graph(self, weight):
        La = self._RS.weight_space(extended = True).basis()
        # TODO: improve choice of max_depth
        V = crystals.LSPaths(sum([x*y for x,y in zip(weight,La)])).subcrystal(max_depth=self._n**2)
        v = V.module_generators[0]
        G = self._recursive_path_graph(self._g, v, V)
        sinks = G.sinks()
        while sinks:
            for w in sinks:
                G.delete_vertex(w)
            sinks = G.sinks()
        return (G, v)

    def _recursive_evaluation(self, g, ls_path, crystal, jump_dict=None, jump_history=dict()):
        #if len(jump_history) > 1:
        #    print("lots of jumps " + str(len(jump_history)))

        if not g:
            if ls_path == crystal.module_generators[0]:
                return 1
            else:
                return 0

        if not jump_dict:
            jump_dict = dict()
            for v in crystal:
                wt = v.weight()
                if wt in jump_dict:
                    jump_dict[wt].append(v)
                else:
                    jump_dict[wt] = [v]

        new_g = copy(g)
        cpt , i = new_g.pop()
        wt = ls_path.weight()

        if cpt == 'h':
            return self._recursive_evaluation(new_g, copy(ls_path), crystal, jump_dict=jump_dict, jump_history=jump_history) * prod([ h**a for (h,a) in zip(self._R.gens()[:self._n],vector(wt))])
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
        for v in jump_dict[wt]:
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
            jump_list.append(jump_history[ls_path])

        #jump_list = jump_dict[wt]

        output = 0

        for v in jump_list:
#            if v != ls_path:
#                G.add_edge((ls_path,v,'jump'))
            new_ls_path = copy(v)
            current_k = k(v,i)
            current_l = l(v,i)
            old_output = output
            for j in range(current_k+1):
                if v == ls_path:
                    shift = 0
                    #if new_ls_path in jump_history:
                    #    shift = k(new_ls_path, i)
                    output += self._recursive_evaluation(new_g, new_ls_path, crystal, jump_dict=jump_dict, jump_history=jump_history) * gen**j * binomial(current_l-shift+j,j)
                elif j > 0:
                    new_jump_history = copy(jump_history)
                    new_jump_history[v] = ls_path
                    output += self._recursive_evaluation(new_g, new_ls_path, crystal, jump_dict=jump_dict, jump_history=new_jump_history) * gen**j * binomial(current_l+j,j)
                    if output != old_output and ls_path not in self._jump_locations:
                        self._jump_locations.append(ls_path)

                if j < current_k:
                    new_ls_path = step(new_ls_path, i)
#                    if cpt == '+':

        return output


    def _recursive_evaluation2(self, g, ls_path, crystal, jump_dict=None, just_jumped = False):

        if not jump_dict:
            jump_dict = dict()
            for v in crystal:
                wt = v.weight()
                if wt in jump_dict:
                    jump_dict[wt].append(v)
                else:
                    jump_dict[wt] = [v]

        if not g:
            if ls_path == crystal.module_generators[0]:
                return 1
            else:
                return 0

        new_g = copy(g)
        cpt , i = new_g.pop()
        wt = ls_path.weight()

        if cpt == 'h':
            return self._recursive_evaluation(new_g, copy(ls_path), crystal, jump_dict=jump_dict) * prod([ h**a for (h,a) in zip(self._R.gens()[:self._n],vector(wt))])
        elif cpt == '-':
            step = lambda v,i: v.f(i)
            k = lambda v,i: v.phi(i)
            l = lambda v,i: v.epsilon(i)
            gen = self._R.gens()[self._n+i]
        else:
            step = lambda v,i: v.e(i)
            k = lambda v,i: v.epsilon(i)
            l = lambda v,i: v.phi(i)
            gen = self._R.gens()[2*self._n+i]

        if just_jumped:
            jump_list = [ls_path]
        else:
            jump_list = jump_dict[wt]

        output = 0

        for v in jump_list:
            new_ls_path = copy(v)
            current_k = k(v,i)
            current_l = l(v,i)
            for j in range(current_k+1):
                output += self._recursive_evaluation(new_g, new_ls_path, crystal, jump_dict=jump_dict, just_jumped = v != ls_path) * gen**j * binomial(current_l+j,current_l)
                #if j < current_k:
                new_ls_path = step(new_ls_path, i)
                #if new_ls_path == None:
                #    break
        return output

    def generalized_minor(self, weight):
        La = self._RS.weight_space(extended = True).basis()
        # TODO: improve choice of max_depth
        V = crystals.LSPaths(sum([x*y for x,y in zip(weight,La)])).subcrystal(max_depth=self._n**2)
        v = V.module_generators[0]
        return self._recursive_evaluation(self._g, v, V)


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
    if not (T.is_affine() and T.is_acyclic()) :
        raise NotImplementedError("This code works only with acyclic affine exchange matrices")

    print("################################################################################\n")
    print("Testing Dynkin diagram\n"+T.cartan_type().ascii_art()+"\ncoxeter = " + str(coxeter))
    print

    # find regular_g_vectors
    tubes = flatten(T.affine_tubes())
    regular_g_vectors = map(lambda x: vector(T.to_weight(x)),tubes)

    # blackbox to produce cluster variables
    n = b_matrix.ncols()
    Bdp = block_matrix([[b_matrix],[identity_matrix(n)],[identity_matrix(n)]])
    A = ClusterAlgebra(Bdp)

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
            #print(minor-variable)
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

