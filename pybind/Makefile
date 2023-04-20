PYTHONARCH=python3
# PYTHONARCH=/bin/python3
INCLUDE+=$(shell $(PYTHONARCH)-config --includes)
INCLUDE+=-I/usr/include/eigen3
SUFFIX=$(shell $(PYTHONARCH)-config --extension-suffix)


all: _invariant_manifold_cc$(SUFFIX)  _invariant_manifold_rc$(SUFFIX) _invariant_manifold_rr$(SUFFIX)

CC=g++

_invariant_manifold_cc$(SUFFIX): ../obj/*cc.o ../obj/poly_utils.o
	$(CC) -O3 -Wall -shared  -fPIC -DSCALAR_MODE=0 $(INCLUDE) python_invariant_manifold.cpp ../obj/*cc.o ../obj/poly_utils.o -o _invariant_manifold_cc$(SUFFIX)

_invariant_manifold_rc$(SUFFIX): ../obj/*rc.o ../obj/poly_utils.o
	$(CC) -O3 -Wall -shared  -fPIC -DSCALAR_MODE=1 $(INCLUDE) python_invariant_manifold.cpp ../obj/*rc.o ../obj/poly_utils.o -o _invariant_manifold_rc$(SUFFIX)

_invariant_manifold_rr$(SUFFIX): ../obj/*rr.o ../obj/poly_utils.o
	$(CC) -O3 -Wall -shared  -fPIC -DSCALAR_MODE=2 $(INCLUDE) python_invariant_manifold.cpp ../obj/*rr.o ../obj/poly_utils.o -o _invariant_manifold_rr$(SUFFIX)

_jet_transport_rr$(SUFFIX): ../obj/*rr.o ../obj/poly_utils.o
	$(CC) -O3 -Wall -shared  -fPIC -DSCALAR_MODE=2 $(INCLUDE) python_jet_transport.cpp ../obj/*rr.o ../obj/poly_utils.o -o _jet_transport_rr$(SUFFIX)
	mv _jet_transport_rr$(SUFFIX) ../python/py_jet_transport
