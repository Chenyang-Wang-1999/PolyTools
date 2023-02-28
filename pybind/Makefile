PYTHONARCH=python3
# PYTHONARCH=/bin/python3
INCLUDE+=$(shell $(PYTHONARCH)-config --includes)
INCLUDE+=-I/usr/include/eigen3
SUFFIX=$(shell $(PYTHONARCH)-config --extension-suffix)


all: _invariant_manifold_cc$(SUFFIX)  _invariant_manifold_rc$(SUFFIX) _invariant_manifold_rr$(SUFFIX)

_invariant_manifold_cc$(SUFFIX): ../obj/*cc.o ../obj/poly_utils.o
	c++ -O2 -Wall -shared  -fPIC -DSCALAR_MODE=0 $(INCLUDE) python_invariant_manifold.cpp ../obj/*cc.o ../obj/poly_utils.o -o _invariant_manifold_cc$(SUFFIX)

_invariant_manifold_rc$(SUFFIX): ../obj/*rc.o ../obj/poly_utils.o
	c++ -O2 -Wall -shared  -fPIC -DSCALAR_MODE=1 $(INCLUDE) python_invariant_manifold.cpp ../obj/*rc.o ../obj/poly_utils.o -o _invariant_manifold_rc$(SUFFIX)

_invariant_manifold_rr$(SUFFIX): ../obj/*rr.o ../obj/poly_utils.o
	c++ -O2 -Wall -shared  -fPIC -DSCALAR_MODE=2 $(INCLUDE) python_invariant_manifold.cpp ../obj/*rr.o ../obj/poly_utils.o -o _invariant_manifold_rr$(SUFFIX)
