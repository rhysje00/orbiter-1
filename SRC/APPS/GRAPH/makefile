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



all: graph.out

graph.out:  graph.o graph_generator.o ;
	$(DISCRETA_CPP)  graph.o graph_generator.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

clean:
	- rm *.o *.a *.out *~

graph.o: graph.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c graph.C -o $@
graph_generator.o: graph_generator.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c graph_generator.C -o $@




