all: obj/polys.o obj/poly_utils.o obj/invariant_manifold_solver.o

obj/polys.o:
	c++ -O2 -Wall -shared -fPIC -c src/polys.cpp -o obj/polys.o 

obj/poly_utils.o:
	c++ -O2 -Wall -shared -fPIC -c src/poly_utils.cpp -o obj/poly_utils.o

obj/invariant_manifold_solver.o:
	c++ -O2 -Wall -shared -fPIC -c src/invariant_manifold_solver.cpp -o obj/invariant_manifold_solver.o