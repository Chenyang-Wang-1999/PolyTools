all: cc_all rc_all rr_all
cc_all: obj/polys_cc.o obj/invariant_manifold_solver_cc.o obj/poly_utils.o
rc_all: obj/polys_rc.o obj/invariant_manifold_solver_rc.o obj/poly_utils.o
rr_all: obj/polys_rr.o obj/invariant_manifold_solver_rr.o obj/poly_utils.o

DEFINATIONS=
CC=g++

obj/polys_cc.o:
	$(CC) -O3 ${DEFINATIONS} -DSCALAR_MODE=0 -Wall -shared -fPIC -c src/polys.cpp -o obj/polys_cc.o 

obj/polys_rc.o:
	$(CC) -O3 ${DEFINATIONS} -DSCALAR_MODE=1 -Wall -shared -fPIC -c src/polys.cpp -o obj/polys_rc.o 

obj/polys_rr.o:
	$(CC) -O3 ${DEFINATIONS} -DSCALAR_MODE=2 -Wall -shared -fPIC -c src/polys.cpp -o obj/polys_rr.o 

obj/poly_utils.o:
	$(CC) -O3 ${DEFINATIONS} -Wall -shared -fPIC -c src/poly_utils.cpp -o obj/poly_utils.o

obj/invariant_manifold_solver_cc.o:
	$(CC) -O3 ${DEFINATIONS} -DSCALAR_MODE=0 -Wall -shared -fPIC -c src/invariant_manifold_solver.cpp -o obj/invariant_manifold_solver_cc.o

obj/invariant_manifold_solver_rc.o:
	$(CC) -O3 ${DEFINATIONS} -DSCALAR_MODE=1 -Wall -shared -fPIC -c src/invariant_manifold_solver.cpp -o obj/invariant_manifold_solver_rc.o

obj/invariant_manifold_solver_rr.o:
	$(CC) -O3 ${DEFINATIONS} -DSCALAR_MODE=2 -Wall -shared -fPIC -c src/invariant_manifold_solver.cpp -o obj/invariant_manifold_solver_rr.o
