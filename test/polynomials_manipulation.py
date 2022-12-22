import sympy as sym

def test_comp():
    x, y, z = sym.symbols('x,y,z')
    u,v,w = sym.symbols('u,v,w')

    f = (3+2j) * x**2 * y *z - 2 * x*z

    g = [u + v, u - v, u + 2*v**2]

    f_subs = f.subs([(x,g[0]), (y, g[1]), (z, g[2])])
    f_subs = f_subs.expand()
    print(f_subs)

def test_der():
    x, y, z = sym.symbols('x,y,z')
    u,v,w = sym.symbols('u,v,w')

    f = (3+2j) * x**2 * y *z - 2.0j * x*z + 3.0*(y**2) *(z**2)
    # f = x

    g = [x + y, x - z, y + 2*z**2]

    h = sym.diff(f, x) * g[0] + sym.diff(f, y) * g[1] + sym.diff(f, z) * g[2]
    print(h.expand())
 
if __name__ == '__main__':
    test_comp()
    # test_der()