import sympy as sym

def test_comp():
    x, y, z = sym.symbols('x,y,z')
    u,v,w = sym.symbols('u,v,w')

    f = (3+2j) * x**2 * y *z - 2 * x*z

    g = [u + v, u - v, u + 2*v**2]

    f_subs = f.subs([(x,g[0]), (y, g[1]), (z, g[2])])
    f_subs = f_subs.expand()
    print(f_subs)

if __name__ == '__main__':
    test_comp()