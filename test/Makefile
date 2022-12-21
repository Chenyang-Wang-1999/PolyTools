
test_homogen_basics.exe:  ../obj/polys.o
	c++ -Wall test_homogen_basics.cpp ../obj/polys.o -o test_homogen_basics.exe 


test_poly_mul_sweeper.exe:  ../obj/polys.o
	c++ -Wall test_poly_mul_sweeper.cpp ../obj/polys.o -o test_poly_mul_sweeper.exe 


test_remove_zeros.exe:  ../obj/polys.o
	c++ -Wall test_remove_zeros.cpp ../obj/polys.o -o test_remove_zeros.exe 

test_eval_and_mul.exe:  ../obj/polys.o
	c++ -Wall -DDEBUG test_eval_and_mul.cpp ../obj/polys.o -o test_eval_and_mul.exe 

test_mul_seq.exe:  ../obj/polys.o
	c++ -Wall -DDEBUG test_mul_seq.cpp ../obj/polys.o -o test_mul_seq.exe 

test_term_comp.exe:  ../obj/polys.o ../obj/poly_utils.o
	c++ -Wall -DDEBUG test_term_comp.cpp ../obj/polys.o ../obj/poly_utils.o -o test_term_comp.exe 

test_order_iteration.exe:  ../obj/poly_utils.o
	c++ -Wall -DDEBUG test_order_iteration.cpp ../obj/poly_utils.o -o test_order_iteration.exe 

test_matrix_mul.exe:  ../obj/polys.o ../obj/poly_utils.o
	c++ -Wall -DDEBUG test_matrix_mul.cpp ../obj/polys.o ../obj/poly_utils.o -o test_matrix_mul.exe 
