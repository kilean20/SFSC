CC = g++
OPTS = -O3

run_pass_sc_quad: statvec.cpp statvec.h pass_sc_equad.cpp pass_sc_equad.h test_scQuad.cpp 
	$(CC) $(OPTS) -o $@

clean:
	rm -f *.o 

.PHONY: clean
