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
    A = ClusterAlgebra(b_matrix, principal_coefficients=True)

    # blackbox to produce minors
    R = DoubleBruhatAlgebra(b_matrix)

    initial_cluster = R._initial_cluster

    for gvect in regular_g_vectors:
        print(str(gvect) + "\tin the orbit of the finite type dominant weight " + str((0,)+tuple(R.level_zero_dominant_conjugate(R.g_to_weight(gvect)))[:b_matrix.ncols()-1]) + "\t: ", end="")
        A.find_cluster_variable(gvect)
        variable = A.cluster_variable(gvect).lift().subs(initial_cluster)
        minor = R.generic_evaluation(R._double_coxeter,R.g_to_weight(gvect))
        if variable == minor:
            print("Pass")
        else:
            print("Fail")
            #return False
    print
    return True

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

    def __init__(self, b_matrix):
        """
        We assume that b_matrix is an orientation of an affine Dynkin diagram,
        and that the first column is the affine root.
        """

        if not b_matrix.is_skew_symmetrizable():
            raise ValueError("The input must be a skew symmetrizable integer matrix")
        self._B = copy(b_matrix)
        self._rank = self._B.ncols()

        # create extended affine root system
        self._ext_Cartan_mat = block_diagonal_matrix(2-matrix(self._rank,map(abs,self._B.list())),matrix(1))
        self._ext_Cartan_mat[0,self._rank] = 1
        self._ext_Cartan_mat[self._rank,0] = 1
        self._ext_symm_mat = diagonal_matrix(self._ext_Cartan_mat.is_symmetrizable(return_diag=True))
        self._ext_symm_mat *= self._ext_symm_mat.denominator()
        self._RootSystem = RootSystem(self._rank+1,self._ext_Cartan_mat,self._ext_symm_mat)

        # create finite sub root system
        self._sub_Cartan_mat = self._ext_Cartan_mat[1:self._rank,1:self._rank]
        self._sub_symm_mat = self._ext_symm_mat[1:self._rank,1:self._rank]
        self._sub_RootSystem = RootSystem(self._rank-1,self._sub_Cartan_mat,self._sub_symm_mat)

        # create Coxeter words and element
        self._coxeter_word = self.coxeter()
        self._coxeter_element = prod([self._RootSystem._simple_reflections[i] for i in self._coxeter_word])
        self._double_coxeter = [(i,-1) for i in self._coxeter_word]
        cv = list(self._coxeter_word)
        cv.reverse()
        self._double_coxeter += [(i,1) for i in cv]

        # create ambient polynomial ring
        self._parameter_polynomial_ring = PolynomialRing(QQ,['t%s'%i for i in xrange(self._rank)]+['u%s'%i for i in xrange(self._rank)])
        self._polygens = self._parameter_polynomial_ring.gens()

        # create cluster algebra with principal coefficients
        self._cluster_algebra = ClusterAlgebra(block_matrix([[self._B],[identity_matrix(self._rank)]]))

        # compute generic evaluations of principal coefficients
        self._coefficients = []
        temp_coeff = []
        for i in xrange(self._rank):
            # next two lines depend on the implementation of weights by RootSystem
            root = self._RootSystem._fundamental_weights[i]-self._coxeter_element*self._RootSystem._fundamental_weights[i]
            coeff = self._polygens[self._rank+i]**(-1)*prod([self._polygens[j]**root[self._rank+1+j] for j in xrange(self._rank)]) 
            temp_coeff.append(coeff)
        for j in xrange(self._rank):
            coeff = temp_coeff[j]
            for i in self._coxeter_word:
                if i == j:
                    break
                else:
                    coeff *= temp_coeff[i]**self._ext_Cartan_mat[i,j]
            self._coefficients.append(coeff)

        # specify cluster variables to generalized minors
        clgens = self._cluster_algebra.ambient().gens()
        self._initial_cluster = dict([(clgens[i],self._polygens[self._rank+i]**(-1)) for i in xrange(self._rank)]+[(clgens[self._rank+i],self._coefficients[i]) for i in xrange(self._rank)])

    def coxeter(self):
        r"""
        Returns a list expressing the coxeter element corresponding to self._B
        (twisted) reflections are applied from top of the list, for example
        [2, 1, 0] correspond to s_2s_1s_0

        Sources == non positive columns == leftmost letters
        """
        zero_vector = vector([0 for x in range(self._rank)])
        coxeter = []
        B = copy(self._B)
        columns = B.columns()
        source = None
        for j in range(self._rank):
            for i in range(self._rank):
                if all(x <=0 for x in columns[i]) and columns[i] != zero_vector:
                    source = i
                    break
            if source == None:
                if B != matrix(self._rank):
                    raise ValueError("Unable to find a Coxeter element representing self._B")
                coxeter += [ x for x in range(self._rank) if x not in coxeter]
                break
            coxeter.append(source)
            columns[source] = zero_vector
            B = matrix(columns).transpose()
            B[source] = zero_vector
            columns = B.columns()
            source = None
        return tuple(coxeter)

    def g_to_weight(self,gvect):
        return sum([gvect[i]*self._RootSystem.fundamental_weight(i) for i in xrange(self._rank)])

    def truncate_weight(self,wt):
        return sum([self._RootSystem.weightify(wt)[i]*self._sub_RootSystem.fundamental_weight(i-1) for i in xrange(1,self._rank)])

    def level_zero_weight_multiplicity(self, highest_wt, wt):
        # return multiplicity of wt in level zero representation indexed by dominant finite-type highest_wt
        return self._sub_RootSystem.weight_multiplicity(highest_wt,self.truncate_weight(wt))

    def validate_weight(self, xlist, wt1, wt2, highest_wt, alpha):
        # check whether there is an ambiguity in the next step of generic_evaluation
        current_wt = copy(wt1)
        current_wt_mult = self.level_zero_weight_multiplicity(highest_wt, current_wt)
        initial_wt_mult = current_wt_mult
        while current_wt_mult != 0:
            if current_wt_mult < initial_wt_mult:
                print("There was an ambiguity.")
                print("initial_wt_mult = ", initial_wt_mult)
                print("current_wt_mult = ", current_wt_mult)
                print("current_wt = ", current_wt)
                print("alpha = ", alpha)
                print("xlist = ", xlist)
                print("wt1 = ", wt1)
                print("wt2 = ", wt2)
            current_wt += alpha
            current_wt_mult = self.level_zero_weight_multiplicity(highest_wt, current_wt)

    def alpha_string(self, wt, highest_wt, alpha):
        #determines the length of the alpha string containing wt
        #might fail if validate_weight gives a complaint
        current_wt = copy(wt)
        current_wt_mult = self.level_zero_weight_multiplicity(highest_wt, current_wt)
        initial_wt_mult = current_wt_mult
        num_steps = 0
        while current_wt_mult != 0:
            current_wt += alpha
            num_steps += 1
            current_wt_mult = self.level_zero_weight_multiplicity(highest_wt, current_wt)
        string_length = self._RootSystem.pairing(alpha,current_wt)
        #print "sl_2 weight=",string_length
        return (num_steps-1,string_length-2)

    def level_zero_dominant_conjugate(self, wt):
        # wt is an element of the finite-type weight subspace of the affine weight space
        wt = self._RootSystem.weightify(wt)
        trunc_wt = self.truncate_weight(wt)
        for w in self._sub_RootSystem.Weyl_group():
            Weyl_wt = w*trunc_wt
            if self._sub_RootSystem.is_dominant(Weyl_wt):
                return self._sub_RootSystem.weightify(Weyl_wt)
        return self._sub_RootSystem._zero()

    def generic_evaluation(self, xlist, wt1, wt2 = None, highest_wt = None):
        if wt2 == None:
            wt2 = copy(wt1)
        if highest_wt == None:
            highest_wt = self.level_zero_dominant_conjugate(wt2)
        if xlist == []:
            if wt1 == wt2:
                return 1
            else:
                return 0
        new_xlist = copy(xlist)
        i, eps = new_xlist.pop()
        alpha = eps * self._RootSystem._simple_roots[i]
        pairing = self._RootSystem.pairing(alpha, wt1)
        self.validate_weight(xlist, wt1, wt2, highest_wt, sign(pairing)*alpha if pairing != 0 else alpha)
        output = 0
        j = 0
        new_wt1 = copy(wt1)
        while self.level_zero_weight_multiplicity(highest_wt, new_wt1) != 0:
            k,n = self.alpha_string(wt1,highest_wt,alpha)
            if eps > 0:
                #if self._polygens[i]**j * binomial(n-k+j,n-k) != 1:
                #    print self._RootSystem.weightify(wt1), self._polygens[i]**j * binomial(n-k+j,n-k), n-k+j, n-k, k, n
                # this records the action of the matrix [[1,t],[0,1]]
                output += self.generic_evaluation(new_xlist, new_wt1, wt2, highest_wt) * self._polygens[i]**j * binomial(n-k+j,n-k)
            else:
                #if self._polygens[self._rank + i]**(pairing + j) * binomial(n-k+j,n-k) != 1:
                #    print self._RootSystem.weightify(wt1), self._polygens[self._rank + i]**(pairing + j) * binomial(n-k+j,n-k), n-k+j,n-k, k, n
                # this records the action of the matrix [[u^{-1},0],[1,u]] = [[1,0],[u,1]]*[[u^{-1},0],[0,u]]
                output += self.generic_evaluation(new_xlist, new_wt1, wt2, highest_wt) * self._polygens[self._rank + i]**(pairing + j) * binomial(n-k+j,n-k)
            j += 1
            new_wt1 += alpha
        return output

    def compare_constructions(self,glist):
        """
        Input: A list of g-vectors
        Output: A comparison of the cluster variables with these g-vectors (evaluated in the parameter ring) and the corresponding
                sidest weight minors evaluated at a generic point of the reduced double Bruhat cell
        """
        for gvect in glist:
            self._cluster_algebra.find_cluster_variable(gvect)
            cl_minor = self._cluster_algebra.cluster_variable(gvect).lift().subs(self._initial_cluster)
            gen_minor = self.generic_evaluation(self._double_coxeter,self.g_to_weight(gvect))
            if cl_minor == gen_minor:
                print(str(gvect)+": True")
            else:
                print(str(gvect)+": False")
                #print("  Cluster minor=",cl_minor)
                #print("  Generalized minor=",gen_minor)
                print("  Diff=",expand(factor(cl_minor-gen_minor)))

