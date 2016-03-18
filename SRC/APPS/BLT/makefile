SRC=../..
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



all: \
	blt_main.out 




blt_main.out:  blt_main.o blt_set.o blt_set2.o ;
	$(DISCRETA_CPP)  blt_main.o blt_set.o blt_set2.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@


clean:
	- rm *.o *.a *.out *~

blt_main.o: blt_main.C blt.h; $(DISCRETA_CPP) $(CPPFLAGS) -c blt_main.C -o $@
blt_set.o: blt_set.C blt.h ; $(DISCRETA_CPP) $(CPPFLAGS) -c blt_set.C -o $@
blt_set2.o: blt_set2.C blt.h ; $(DISCRETA_CPP) $(CPPFLAGS) -c blt_set2.C -o $@




