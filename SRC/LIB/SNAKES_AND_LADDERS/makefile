DISCRETA_CPP = /usr/bin/g++
CPPFLAGS = -Wall -I.. -I../GALOIS
LIB = 
LFLAGS = -lm
CPP_PREFIX=c_
UNIX_TO_DOS = ~/bin/unix2dos -n 

MODULES = \
	compute_stabilizer.o \
	extension.o \
	generator.o \
	generator_init.o \
	generator_classify.o \
	generator_combinatorics.o \
	generator_draw.o \
	generator_io.o \
	generator_trace.o \
	oracle.o \
	oracle_downstep.o \
	oracle_downstep_subspace_action.o \
	oracle_io.o \
	oracle_upstep.o \
	oracle_upstep_subspace_action.o \
	recognize.o \
	set_stabilizer_compute.o \
	snakes_and_ladders_global.o \
	upstep_work.o \
	upstep_work_subspace_action.o \
	upstep_work_trace.o 


		

LIB = libsnakesandladders.a

libsnakesandladders.a: $(MODULES) ;
	ar ruvs $(LIB) $(MODULES)

clean:
	- rm *.o *.a *~

# BASIC_FILES:
compute_stabilizer.o: compute_stabilizer.C; $(DISCRETA_CPP) $(CPPFLAGS) -c compute_stabilizer.C -o $@
extension.o: extension.C; $(DISCRETA_CPP) $(CPPFLAGS) -c extension.C -o $@
generator.o: generator.C; $(DISCRETA_CPP) $(CPPFLAGS) -c generator.C -o $@
generator_init.o: generator_init.C; $(DISCRETA_CPP) $(CPPFLAGS) -c generator_init.C -o $@
generator_classify.o: generator_classify.C; $(DISCRETA_CPP) $(CPPFLAGS) -c generator_classify.C -o $@
generator_combinatorics.o: generator_combinatorics.C; $(DISCRETA_CPP) $(CPPFLAGS) -c generator_combinatorics.C -o $@
generator_draw.o: generator_draw.C; $(DISCRETA_CPP) $(CPPFLAGS) -c generator_draw.C -o $@
generator_io.o: generator_io.C; $(DISCRETA_CPP) $(CPPFLAGS) -c generator_io.C -o $@
generator_trace.o: generator_trace.C; $(DISCRETA_CPP) $(CPPFLAGS) -c generator_trace.C -o $@
oracle.o: oracle.C; $(DISCRETA_CPP) $(CPPFLAGS) -c oracle.C -o $@
oracle_downstep.o: oracle_downstep.C; $(DISCRETA_CPP) $(CPPFLAGS) -c oracle_downstep.C -o $@
oracle_downstep_subspace_action.o: oracle_downstep_subspace_action.C; $(DISCRETA_CPP) $(CPPFLAGS) -c oracle_downstep_subspace_action.C -o $@
oracle_io.o: oracle_io.C; $(DISCRETA_CPP) $(CPPFLAGS) -c oracle_io.C -o $@
oracle_upstep.o: oracle_upstep.C; $(DISCRETA_CPP) $(CPPFLAGS) -c oracle_upstep.C -o $@
oracle_upstep_subspace_action.o: oracle_upstep_subspace_action.C; $(DISCRETA_CPP) $(CPPFLAGS) -c oracle_upstep_subspace_action.C -o $@
recognize.o: recognize.C; $(DISCRETA_CPP) $(CPPFLAGS) -c recognize.C -o $@
set_stabilizer_compute.o: set_stabilizer_compute.C; $(DISCRETA_CPP) $(CPPFLAGS) -c set_stabilizer_compute.C -o $@
snakes_and_ladders_global.o: snakes_and_ladders_global.C; $(DISCRETA_CPP) $(CPPFLAGS) -c snakes_and_ladders_global.C -o $@
upstep_work.o: upstep_work.C; $(DISCRETA_CPP) $(CPPFLAGS) -c upstep_work.C -o $@
upstep_work_subspace_action.o: upstep_work_subspace_action.C; $(DISCRETA_CPP) $(CPPFLAGS) -c upstep_work_subspace_action.C -o $@
upstep_work_trace.o: upstep_work_trace.C; $(DISCRETA_CPP) $(CPPFLAGS) -c upstep_work_trace.C -o $@



cpp:
	$(UNIX_TO_DOS) compute_stabilizer.C $(CPP_PREFIX)compute_stabilizer.cpp
	$(UNIX_TO_DOS) extension.C $(CPP_PREFIX)extension.cpp
	$(UNIX_TO_DOS) generator.C $(CPP_PREFIX)generator.cpp
	$(UNIX_TO_DOS) generator_init.C $(CPP_PREFIX)generator_init.cpp
	$(UNIX_TO_DOS) generator_classify.C $(CPP_PREFIX)generator_classify.cpp
	$(UNIX_TO_DOS) generator_combinatorics.C $(CPP_PREFIX)generator_combinatorics.cpp
	$(UNIX_TO_DOS) generator_draw.C $(CPP_PREFIX)generator_draw.cpp
	$(UNIX_TO_DOS) generator_io.C $(CPP_PREFIX)generator_io.cpp
	$(UNIX_TO_DOS) generator_trace.C $(CPP_PREFIX)generator_trace.cpp
	$(UNIX_TO_DOS) oracle.C $(CPP_PREFIX)oracle.cpp
	$(UNIX_TO_DOS) oracle_downstep.C $(CPP_PREFIX)oracle_downstep.cpp
	$(UNIX_TO_DOS) oracle_downstep_subspace_action.C $(CPP_PREFIX)oracle_downstep_subspace_action.cpp
	$(UNIX_TO_DOS) oracle_io.C $(CPP_PREFIX)oracle_io.cpp
	$(UNIX_TO_DOS) oracle_upstep.C $(CPP_PREFIX)oracle_upstep.cpp
	$(UNIX_TO_DOS) oracle_upstep_subspace_action.C $(CPP_PREFIX)oracle_upstep_subspace_action.cpp
	$(UNIX_TO_DOS) recognize.C $(CPP_PREFIX)recognize.cpp
	$(UNIX_TO_DOS) set_stabilizer_compute.C $(CPP_PREFIX)set_stabilizer_compute.cpp
	$(UNIX_TO_DOS) snakes_and_ladders_global.C $(CPP_PREFIX)snakes_and_ladders_global.cpp
	$(UNIX_TO_DOS) upstep_work.C $(CPP_PREFIX)upstep_work.cpp
	$(UNIX_TO_DOS) upstep_work_subspace_action.C $(CPP_PREFIX)upstep_work_orbits_on_subspaces.cpp
	$(UNIX_TO_DOS) upstep_work_trace.C $(CPP_PREFIX)upstep_work_trace.cpp
	$(UNIX_TO_DOS) snakesandladders.h ../WINDOWS/snakesandladders.h
	mv $(CPP_PREFIX)* ../WINDOWS
	
	
