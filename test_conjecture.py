from __future__ import print_function
from tropical_cluster_algebra import TropicalClusterAlgebra
from cluster_algebra import ClusterAlgebra
from root_system import RootSystem 

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

    problems = False
    for gvect in regular_g_vectors:
        print(str(gvect) + "\tin the orbit of the finite type dominant weight " + str((0,)+tuple(R._level_zero_dominant_conjugate(R._g_to_weight(gvect)))[:b_matrix.ncols()-1]) + "\t: ", end="")
        A.find_cluster_variable(gvect)
        variable = A.cluster_variable(gvect).lift().subs(substitution)
        minor = R.level_zero_minor(gvect)
        if variable == minor:
            print("Pass")
        else:
            print("Fail")
            problems = True
    print
    return not problems

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

class DoubleBruhatAlgebra(SageObject):
    r"""
    This code assumes that the affine node is the one labelled by 0
    """

    def __init__(self, data):

        # process input 
        Q = ClusterQuiver(data)
        # This is commented because QuiverMutationType is not able to identify type D
        #if not Q.mutation_type().is_affine():
        #    raise ValueError("Input must be of affine type")
        n = Q.n()
        B = Q.b_matrix()

        # extended Cartan matrix       
        self._A = block_diagonal_matrix(2-matrix(n,map(abs,B.list())),matrix(1))
        self._A[0,n] = 1
        self._A[n,0] = 1

        # extended symmetrizer
        self._D = diagonal_matrix(self._A.is_symmetrizable(return_diag=True))
        self._D *= self._D.denominator()
        
        # extended root system
        self._rs = RootSystem(n+1,self._A,self._D)

        # finite Cartan matrix
        self._A_f = self._A[1:n,1:n]

        # finite symmetrizer
        self._D_f = self._D[1:n,1:n]

        # finite root system
        self._rs_f = RootSystem(n-1,self._A_f,self._D_f)
        
        # generic element factorization 
        # elements are of the form (s,i) with s one of '-', 'h', '+' and i in range(n)
        # the second entry is irrelevant if the first entry is 'h'
        c = coxeter_element(B)
        self._g = [ ('-', i) for i in c ]
        self._g.append(('h',None))
        self._g += [ ('+', i) for i in reversed(c) ]

        # ambient ring
        self._R = PolynomialRing(QQ, sorted(flatten([('h%s'%i,'r%s'%i,'t%s'%i) for i in xrange(n)])))

        # 
        self._n = n
        self._B = B

    def _recursive_evaluation(self, g,  wt1, wt2, highest_wt):
        if not g:
            if wt1 == wt2:
                return 1
            else:
                return 0

        new_g = copy(g)
        cpt , i = new_g.pop()
        new_wt1 = copy(wt1)
        
        if cpt == 'h':
            return self._recursive_evaluation(new_g, new_wt1, wt2, highest_wt) * prod([ h**a for (h,a) in zip(self._R.gens()[:self._n],self._rs.weightify(wt1)[:self._n])])

        eps = 1 if cpt == '+' else -1
        alpha = eps * self._rs._simple_roots[i] 
        pairing = self._rs.pairing(alpha, wt1)

        output = 0
        j = 0    
        while self._level_zero_weight_multiplicity(highest_wt, new_wt1) != 0:
            k,n = self._alpha_string(wt1,highest_wt,alpha)
            if cpt == '+':
                output += self._recursive_evaluation(new_g, new_wt1, wt2, highest_wt) * self._R.gens()[self._n+i]**j * binomial(n-k+j,n-k)
            else:
                output += self._recursive_evaluation(new_g, new_wt1, wt2, highest_wt) * self._R.gens()[2*self._n+i]**j * binomial(n-k+j,n-k)
            j += 1
            new_wt1 += alpha
        return output

    def _truncate_weight(self,wt):
        return sum([self._rs.weightify(wt)[i]*self._rs_f.fundamental_weight(i-1) for i in xrange(1,self._n)])

    def _level_zero_weight_multiplicity(self, highest_wt, wt):
        # return multiplicity of wt in level zero representation indexed by dominant finite-type highest_wt
        return self._rs_f.weight_multiplicity(highest_wt,self._truncate_weight(wt))

    def _alpha_string(self, wt, highest_wt, alpha):
        #determines the length of the alpha string containing wt
        current_wt = copy(wt)
        current_wt_mult = self._level_zero_weight_multiplicity(highest_wt, current_wt)
        initial_wt_mult = current_wt_mult
        num_steps = 0
        while current_wt_mult != 0:
            current_wt += alpha
            num_steps += 1
            current_wt_mult = self._level_zero_weight_multiplicity(highest_wt, current_wt)
        string_length = self._rs.pairing(alpha,current_wt)
        #print "sl_2 weight=",string_length
        return (num_steps-1,string_length-2)

    def _g_to_weight(self,gvect):
        return sum([gvect[i]*self._rs.fundamental_weight(i) for i in xrange(self._n)])

    def _level_zero_dominant_conjugate(self, wt):
        # wt is an element of the finite-type weight subspace of the affine weight space
        wt = self._rs.weightify(wt)
        trunc_wt = self._truncate_weight(wt)
        for w in self._rs_f.Weyl_group():
            Weyl_wt = w*trunc_wt           
            if self._rs_f.is_dominant(Weyl_wt):
                return self._rs_f.weightify(Weyl_wt)
        return self._rs_f._zero()                     

    def level_zero_minor(self, g_vect):
        # This assumes that g_vect is a level zero weight
        wt = self._g_to_weight(g_vect)
        highest_wt = self._level_zero_dominant_conjugate(wt)
        return self._recursive_evaluation(self._g, wt, wt, highest_wt)
