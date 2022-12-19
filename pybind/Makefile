PYTHONARCH=python3
# PYTHONARCH=/bin/python3
INCLUDE+=$(shell $(PYTHONARCH)-config --includes)
# INCLUDE+=-I/usr/include/eigen3
SUFFIX=$(shell $(PYTHONARCH)-config --extension-suffix)

PROJNAME=polymanip

_$(PROJNAME)$(SUFFIX): ../obj/polys.o
	c++ -O3 -Wall -shared  -fPIC $(INCLUDE) python_$(PROJNAME).cpp ../obj/polys.o -o _$(PROJNAME)$(SUFFIX)
	cp _$(PROJNAME)$(SUFFIX) ../python
