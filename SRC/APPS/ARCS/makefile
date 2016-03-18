SRC=../../../../ORBITER/SRC
DISCRETA_CPP = /usr/bin/g++
CPPFLAGS = -Wall -I$(SRC)/LIB -I$(SRC)/LIB/DISCRETA
LIB = \
	$(SRC)/LIB/TOP_LEVEL/libtoplevel.a \
	$(SRC)/LIB/DISCRETA/libdiscreta.a \
	$(SRC)/LIB/SNAKES_AND_LADDERS/libsnakesandladders.a  \
	$(SRC)/LIB/ACTION/libaction.a \
	$(SRC)/LIB/INCIDENCE/libincidence.a \
	$(SRC)/LIB/GALOIS/libgalois.a 
LFLAGS = -lm


all: arcs_main.out \
	test_hyperoval.out 


arcs_main.out:  arcs_main.o arc_generator.o ;
	$(DISCRETA_CPP)  arcs_main.o arc_generator.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

test_hyperoval.out:  test_hyperoval.o;
	$(DISCRETA_CPP)  test_hyperoval.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@




arcs_main.o: arcs_main.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c arcs_main.C -o $@
arc_generator.o: arc_generator.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c arc_generator.C -o $@
test_hyperoval.o: test_hyperoval.C; $(DISCRETA_CPP) $(CPPFLAGS) -c test_hyperoval.C -o $@


clean:
	- rm *.o *.out *~

