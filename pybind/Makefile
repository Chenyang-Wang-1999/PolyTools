PYTHONARCH=python3
# PYTHONARCH=/bin/python3
INCLUDE+=$(shell $(PYTHONARCH)-config --includes)
INCLUDE+=-I/usr/include/eigen3
SUFFIX=$(shell $(PYTHONARCH)-config --extension-suffix)

PROJNAME=invariant_manifold

all: _$(PROJNAME)$(SUFFIX)

_$(PROJNAME)$(SUFFIX): ../obj/*.o
	c++ -O2 -Wall -shared  -fPIC $(INCLUDE) python_$(PROJNAME).cpp ../obj/*.o -o _$(PROJNAME)$(SUFFIX)
