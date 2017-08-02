from itertools import ifilter

class MinorEvaluation(SageObject):
    
    def __init__(self, X, element=None):
        self._C = CartanType(X)
        self._P = self._C.root_system().weight_space(extended = self._C.is_affine())
        if element == None:
            self.set_generic_element()
        elif set(element) == set(self._C.index_set()):
            self.set_coxeter_element(element)
        else:
            self.set_element(element)

    def set_element(self, element):
        self._g = tuple(element)

    def set_generic_element(self):
        w0 = self._P.weyl_group().long_element().reduced_word()
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

        # compute principal minors by default
        if not mu:
            mu = la

        # sanitize input
        la = self._P.from_vector(vector(la))
        mu = self._P.from_vector(vector(mu))
        
        if crystal_generator:
            # this assumes that the representation is highest weight (because of how sage.combinat.crystals.crystals.CrystalBacktracker is implemented)
            crystal_generator = self._P.from_vector(vector(crystal_generator)).to_dominant_chamber()
            crystal = crystals.LSPaths(crystal_generator)
            # this picks the first vertex of the crystal having weight la. The assumption is that there is only one such.
            initial_path = ifilter(lambda v: v.weight() == la, crystal).next()
        else: # by default compute minors in the cyclic irreducible representations having the required weight as extremal 
            crystal = crystals.LSPaths(la)
            initial_path = crystal((la,))
        
        paths = [(initial_path,1)]
        for (l,i,x) in reversed(self._g):
            new_paths = []
            for (wt,mon) in paths:
                if l == 'h':
                    new_paths.append((wt,mon*x**wt.weight()[i]))
                else: # l is either 'e' or 'f'
                    fwd = getattr(wt, 'epsilon' if l == 'e' else 'phi')(i)
                    bwd = getattr(wt, 'epsilon' if l == 'f' else 'phi')(i)
                    for k in range(fwd+1):
                        new_paths.append((wt,binomial(bwd+k,k)*mon*x**k))
                        wt = getattr(wt, l)(i)
            paths = new_paths
        return sum([ mon for (wt,mon) in paths if wt.weight() == mu ])

    def minor_a2sidest(self, la, mu=None, element=None):

        if not element:
            element = self._g

        # compute principal minors by default
        if not mu:
            mu = copy(la)

        # sanitize input
        la = tuple(la)
        mu = tuple(mu)
        
        crystal = A2LoopCrystal(la)
        initial_path = crystal(la)
        
        paths = [(initial_path,1)]
        for (l,i,x) in reversed(self._g):
            new_paths = []
            for (wt,mon) in paths:
                if l == 'h':
                    new_paths.append((wt,mon*x**wt.weight()[i]))
                else: # l is either 'e' or 'f'
                    fwd = getattr(wt, 'epsilon' if l == 'e' else 'phi')(i)
                    bwd = getattr(wt, 'epsilon' if l == 'f' else 'phi')(i)
                    for k in range(fwd+1):
                        new_paths.append((wt,binomial(bwd+k,k)*mon*x**k))
                        wt = getattr(wt, l)(i)
            paths = new_paths
        return sum([ mon for (wt,mon) in paths if wt.weight() == mu ])

from sage.structure.element import Element
class A2LoopCrystalVertex(Element):

    def __init__(self, parent, pt):
        self._pt = tuple(pt)
        Element.__init__(self, parent)

    def e(self, i):
        if i == 0 and self.epsilon(0) > 0:
            return self.parent()((self._pt[0]+2,self._pt[1]-2))
        if i == 1 and self.epsilon(1) > 0: 
            return self.parent()((self._pt[0],self._pt[1]+2))

    def epsilon(self, i):
        if i == 0:
            return self.parent()._m-max(0, self.parent()._m - self._pt[1])/2
        if i == 1:
            return self.parent()._m-max(0, self.parent()._m + self._pt[1])/2

    def f(self, i):
        if i == 0 and self.phi(0) > 0:
            return self.parent()((self._pt[0]-2,self._pt[1]+2))
        if i == 1 and self.phi(1) > 0: 
            return self.parent()((self._pt[0],self._pt[1]-2))

    def phi(self, i):
        if i == 0:
            return max(0, self.parent()._m - self._pt[1])/2
        if i == 1:
            return max(0, self.parent()._m + self._pt[1])/2

    def _repr_(self):
        return str(self._pt)

    def weight(self):
        return copy(self._pt)

class A2LoopCrystal(Parent):
    r"""
    Construct a crystal for the sidest weight representation associated to wt
    """
    
    #Element = A2LoopCrystalVertex

    def __init__(self, wt):
        Parent.__init__(self)
        self._m = abs(wt[1])
        self._wt = tuple(wt)

    def _element_constructor_(self, pt):
        return A2LoopCrystalVertex(self, pt)
        
    def _repr_(self):
        return "The sidest weight representation generated by %s"%str(self._wt)

