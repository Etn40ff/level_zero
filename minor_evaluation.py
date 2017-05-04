class MinorEvaluation(SageObject):
    
    def __init__(self, X, element=None):
        self._C = CartanType(X)
        self._P = self._C.root_system().weight_space()
        self._W = self._P.weyl_group()
        if element == None:
            self.set_generic_element()
        elif set(element) == set(self._C.index_set()):
            self.set_coxeter_element(element)
        else:
            self.set_element(element)

    def set_element(self, element):
        self._g = tuple(element)

    def set_generic_element(self):
        w0 = self._W.long_element().reduced_word()
        element = [ ('f', i, var('t%s'%k)) for (k,i) in enumerate(w0) ]
        element += [ ('h', k, var('h%s'%k)) for k in set(self._C.index_set()) ]
        element += reversed([ ('e', i, var('s%s'%k)) for (k,i) in enumerate(w0) ])
        self.set_element(element)

    def set_coxeter_element(self, coxeter):
        element = [ ('f', i, var('t%s'%i)) for i in coxeter ]
        element += [ ('h', k, var('h%s'%k)) for k in set(self._C.index_set()) ]
        element += [ ('e', i, var('s%s'%i)) for i in reversed(coxeter) ]
        self.set_element(element)

    def minor(self, la, mu=None, element=None, crystal_generator=None):
        # WARNING: this function gives wrong answers precisely to the extent that crystals are not the same thing as representations

        if not element:
            element = self._g

        # compute principal coefficients by default
        if not mu:
            mu = la

        # sanitize input
        la = self._P.from_vector(vector(la))
        mu = self._P.from_vector(vector(mu))

        if not crystal_generator:
            crystal_generator = la.to_dominant_chamber()
        # NOTE TO SELF: here crystals.LSPaths builds the wrong crystal if crystal_generator is an element of self._P. Is this a bug in sage?
        # WORKAROUND: use vectors instead
        # EXAMPLE:
        # sage: C = CartanType(['A',2])
        # sage: P = C.root_system().weight_space()
        # sage: W = P.weyl_group()
        # sage: la = P.from_vector(vector((1,0)))
        # sage: la
        # Lambda[1]
        # sage: crystals.LSPaths(C,la)
        # The crystal of LS paths of type ['A', 1] and weight Lambda[2]
        crystal = crystals.LSPaths(self._C, vector(crystal_generator))

        # The following line assumes that la is a 1-dimensional weight space in the representation
        initial_path = [v for v in crystal if v.weight() == la][0]
        paths = [(initial_path,1)]
        for (l,i,x) in reversed(self._g):
            new_paths = []
            for (wt,mon) in paths:
                if l in ['e','f']:
                    fwd = getattr(wt, 'epsilon' if l == 'e' else 'phi')(i)
                    bwd = getattr(wt, 'epsilon' if l == 'f' else 'phi')(i)
                    for k in range(fwd+1):
                        new_paths.append((wt,binomial(bwd+k,k)*mon*x**k))
                        wt = getattr(wt, l)(i)
                else: # l='h'
                    new_paths.append((wt,mon*x**wt.weight()[i]))
            paths = new_paths
        return sum([ mon for (wt,mon) in paths if wt.weight() == mu ])
                    
