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



all: all_cliques.out \
	all_cycles.out \
	all_rainbow_cliques.out \
	draw_colored_graph.out \
	layered_graph_main.out \
	treedraw.out \
	plot_decomposition_matrix.out  \
	dlx.out \
	solve_diophant.out \
	scheduler.out 




all_cliques.out:  all_cliques.o ;
	$(DISCRETA_CPP)  all_cliques.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

all_cycles.out:  all_cycles.o ;
	$(DISCRETA_CPP)  all_cycles.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

all_rainbow_cliques.out:  all_rainbow_cliques.o ;
	$(DISCRETA_CPP)  all_rainbow_cliques.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

draw_colored_graph.out:  draw_colored_graph.o ;
	$(DISCRETA_CPP)  draw_colored_graph.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

layered_graph_main.out:  layered_graph_main.o ;
	$(DISCRETA_CPP)  layered_graph_main.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

treedraw.out:  treedraw.o ;
	$(DISCRETA_CPP)  treedraw.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

plot_decomposition_matrix.out:  plot_decomposition_matrix.o;
	$(DISCRETA_CPP)  plot_decomposition_matrix.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

dlx.out:  dlx.o;
	$(DISCRETA_CPP)  dlx.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

solve_diophant.out:  solve_diophant.o;
	$(DISCRETA_CPP)  solve_diophant.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

scheduler.out:  scheduler.o;
	$(DISCRETA_CPP)  scheduler.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@




clean:
	- rm *.o *.a *.out *~

all_cliques.o: all_cliques.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c all_cliques.C -o $@
all_cycles.o: all_cycles.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c all_cycles.C -o $@
all_rainbow_cliques.o: all_rainbow_cliques.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c all_rainbow_cliques.C -o $@
draw_colored_graph.o: draw_colored_graph.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c draw_colored_graph.C -o $@
layered_graph_main.o: layered_graph_main.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c layered_graph_main.C -o $@
treedraw.o: treedraw.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c treedraw.C -o $@
plot_decomposition_matrix.o: plot_decomposition_matrix.C; $(DISCRETA_CPP) $(CPPFLAGS) -c plot_decomposition_matrix.C -o $@
dlx.o: dlx.C; $(DISCRETA_CPP) $(CPPFLAGS) -c dlx.C -o $@
solve_diophant.o: solve_diophant.C; $(DISCRETA_CPP) $(CPPFLAGS) -c solve_diophant.C -o $@
scheduler.o: scheduler.C; $(DISCRETA_CPP) $(CPPFLAGS) -c scheduler.C -o $@




