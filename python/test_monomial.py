'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-18
Copyright © Department of Physics, Tsinghua University.  All rights reserved
'''

import _polymanip as cpolymanip

class PolymanipIndexVec:
    vec: cpolymanip.PolymanipIndexVec
    size: int
    def __init__(self, vec:list):
        self.vec = cpolymanip.PolymanipIndexVec(vec)
        self.size = len(vec)
    def get(self):
        return self.vec

class PolymanipScalarVec:
    vec: cpolymanip.PolymanipScalarVec
    size: int
    def __init__(self, vec:list):
        self.vec = cpolymanip.PolymanipScalarVec(vec)
        self.size = len(vec)
    def get(self):
        return self.vec

class Monomial:
    coeff: complex
    dim:int
    _soul: cpolymanip._Monomial
    def __init__(self, **kargs):
        if('cmonomial' in kargs.keys()):
            self._soul = kargs['cmonomial']
            self.coeff = self._soul._coeff
        else:
            coeff = kargs['coeff']
            order_var = kargs['order_var']
            self.coeff = coeff
            self._soul = cpolymanip._Monomial(coeff, order_var.get())
        self.dim = self._soul.dim

    def __LT__(self, another):
        return (self._soul > another._soul)
    
    def __LE__(self, another):
        return (self._soul >= another._soul)
    
    def __mul__(self, another):
        new_soul = self._soul * another._soul
        return Monomial(cmonomial=new_soul)
    
    def __pow__(self, n):
        new_soul = self._soul ** n
        return Monomial(cmonomial=new_soul)
    
    def eval(self, vec:PolymanipScalarVec):
        return self._soul.eval(vec.get())
    
    def get_order(self, j):
        return self._soul.get_order(j)
    
    def get_str(self,name_list:list):
        name_str = str(self.coeff)
        for j in range(self.dim):
            if(self.get_order(j)):
                name_str += " %s^{%d} "%(name_list[j], self.get_order(j))
        return name_str
        
if __name__ == '__main__':
    coeff = 1 + 1j
    order_var = PolymanipIndexVec([2,1,3])
    f = Monomial(coeff = coeff, order_var = order_var)
    g = f**2
    var_str = ['x', 'y', 'z']
    print(f.get_str(var_str))
    print(g.get_str(var_str))
    print(f.eval(PolymanipScalarVec([3.0,2.0,2.0])))
    