#class AffineA1LevelZeroRepresentation(CombinatorialFreeModule):
#
#    def __init__(self, base, factors):
#        basis_keys = [ cartesian_product([ZZ,tuple(range(-d,d+1,2))]) for (d,a) in factors ]
#        basis_keys = cartesian_product(basis_keys)
#        CombinatorialFreeModule.__init__(self, base, basis_keys)
#        self._twists = tuple(a for (d,a) in factors)
#        self._factors = tuple(d for (d,a) in factors)
#
#    def e(self, i, k, t, v):
#        if v.length() > 1:
#            return sum( coeff * self.e(i, k, t, self.basis()[key]) for (key,coeff) in v)
#        u = self(0)
#        key = v.support()[0]
#        if i == 0:
#            for comp in IntegerVectors(k,len(key)):
#                tmp_key = tuple( (a+c,b-2*c) for ((a,b),c) in zip(key,comp) )
#                if tmp_key in self.basis().keys():
#                    u += prod((self._twists[j]*t)**comp[j]*binomial((self._factors[j]+key[j][1])/2,comp[j]) for j in range(len(key)))*self.basis()[tmp_key]
#        if i == 1:
#            for comp in IntegerVectors(k,len(key)):
#                tmp_key = tuple( (a,b+2*c) for ((a,b),c) in zip(key,comp) )
#                if tmp_key in self.basis().keys():
#                    u += prod(t**comp[j]*binomial((self._factors[j]-key[j][1])/2,comp[j]) for j in range(len(key)))*self.basis()[tmp_key]
#        return u*v.coefficients()[0]
#
#    def f(self, i, k, t, v):
#        if v.length() > 1:
#            return sum( coeff * self.f(i, k, t, self.basis()[key]) for (key,coeff) in v)
#        u = self(0)
#        key = v.support()[0]
#        if i == 0:
#            for comp in IntegerVectors(k,len(key)):
#                tmp_key = tuple( (a-c,b+2*c) for ((a,b),c) in zip(key,comp) )
#                if tmp_key in self.basis().keys():
#                    u += prod((self._twists[j]**(-1)*t)**comp[j]*binomial((self._factors[j]-key[j][1])/2,comp[j]) for j in range(len(key)))*self.basis()[tmp_key]
#        if i == 1:
#            for comp in IntegerVectors(k,len(key)):
#                tmp_key = tuple( (a,b-2*c) for ((a,b),c) in zip(key,comp) )
#                if tmp_key in self.basis().keys():
#                    u += prod(t**comp[j]*binomial((self._factors[j]+key[j][1])/2,comp[j]) for j in range(len(key)))*self.basis()[tmp_key]
#        return u*v.coefficients()[0]
#
#    def h(self, i, t, v):
#        if v.length() > 1:
#            return sum( coeff * self.h(i, t, self.basis()[key]) for (key,coeff) in v)
#        return t**sum( (l if i==1 else -l) for (k,l) in v.support()[0]) * v
#
#    def E(self, i, t, v):
#        k = 1
#        w = self(0)
#        u = self.e(i, k, t, v)
#        while u != self(0):
#            w += u
#            k += 1
#            u = self.e(i, k, t, v)
#        return v+w
#
#    def F(self, i, t, v):
#        k = 1
#        w = self(0)
#        u = self.f(i, k, t, v)
#        while u != self(0):
#            w += u
#            k += 1
#            u = self.f(i, k, t, v)
#        return v+w

class AffineA1LevelZeroRepresentation(CombinatorialFreeModule):

    def __init__(self, base, factors):
        basis_keys = [ tuple(range(-d,d+1,2)) for (d,a) in factors ]
        basis_keys = cartesian_product([ZZ]+basis_keys)
        CombinatorialFreeModule.__init__(self, base, basis_keys)
        self._twists = tuple(a for (d,a) in factors)
        self._factors = tuple(d for (d,a) in factors)
        self._nfactors = len(factors)

    def e(self, i, k, t, v):
        if v.length() > 1:
            return sum( coeff * self.e(i, k, t, self.basis()[key]) for (key,coeff) in v)
        u = self(0)
        key = v.support_of_term()
        if i == 0:
            for comp in IntegerVectors(k,self._nfactors):
                tmp_key = (key[0]+k,) + tuple( a-2*c for (a,c) in zip(key[1:],comp))
                if tmp_key in self.basis().keys():
                    u += t**k * prod(self._twists[j]**comp[j] * binomial((self._factors[j]+key[j+1])/2,comp[j]) for j in range(self._nfactors)) * self.basis()[tmp_key]
        if i == 1:
            for comp in IntegerVectors(k,self._nfactors):
                tmp_key = (key[0],) + tuple( a+2*c for (a,c) in zip(key[1:],comp))
                if tmp_key in self.basis().keys():
                    u += t**k * prod(binomial((self._factors[j]-key[j+1])/2,comp[j]) for j in range(self._nfactors)) * self.basis()[tmp_key]
        return u*v.coefficients()[0]

    def f(self, i, k, t, v):
        if v.length() > 1:
            return sum( coeff * self.f(i, k, t, self.basis()[key]) for (key,coeff) in v)
        u = self(0)
        key = v.support_of_term()
        if i == 0:
            for comp in IntegerVectors(k,self._nfactors):
                tmp_key = (key[0]-k,) + tuple( a+2*c for (a,c) in zip(key[1:],comp))
                if tmp_key in self.basis().keys():
                    u += t**k * prod(self._twists[j]**(-comp[j]) * binomial((self._factors[j]-key[j+1])/2,comp[j]) for j in range(self._nfactors)) * self.basis()[tmp_key]
        if i == 1:
            for comp in IntegerVectors(k,self._nfactors):
                tmp_key = (key[0],) + tuple( a-2*c for (a,c) in zip(key[1:],comp))
                if tmp_key in self.basis().keys():
                    u += t**k * prod(binomial((self._factors[j]+key[j+1])/2,comp[j]) for j in range(self._nfactors)) * self.basis()[tmp_key]
        return u*v.coefficients()[0]

    def h(self, i, t, v):
        if v.length() > 1:
            return sum( coeff * self.h(i, t, self.basis()[key]) for (key,coeff) in v)
        return t**( (1 if i == 1 else -1) * sum( v.support_of_term()[1:] )) * v

    def E(self, i, t, v):
        k = 1
        w = self(0)
        u = self.e(i, k, t, v)
        while u != self(0):
            w += u
            k += 1
            u = self.e(i, k, t, v)
        return v+w

    def F(self, i, t, v):
        k = 1
        w = self(0)
        u = self.f(i, k, t, v)
        while u != self(0):
            w += u
            k += 1
            u = self.f(i, k, t, v)
        return v+w

