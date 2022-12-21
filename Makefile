all: obj/poly_utils.o obj/polys.o

obj/polys.o:
	c++ -Wall -c src/polys.cpp -o obj/polys.o 

obj/poly_utils.o:
	c++ -Wall -c src/poly_utils.cpp -o obj/poly_utils.o