DISCRETA_CPP = /usr/bin/g++
DISCRETA_CC = /usr/bin/gcc
CPPFLAGS = -Wall -I.. -I../GALOIS
CCFLAGS = -Wall -I.. -I../GALOIS
LIB = 
LFLAGS = -lm
CPP_PREFIX=d_
UNIX_TO_DOS = ~/bin/unix2dos -n 

MODULES = \
	dynamic_memory.o \
	geo_parameter.o \
	geometric_tests.o \
	inc_gen_global.o \
	incidence_global.o \
	mckay.o \
	packing.o \
	point_line.o \
	refine_rows.o \
	refine_columns.o \
	refine_3design.o \
	tdo_scheme.o \
	tdo_data.o 





LIB = libincidence.a

libincidence.a: $(MODULES) ;
	ar ruvs $(LIB) $(MODULES)

clean:
	- rm *.o *.a *~

# BASIC_FILES:
dynamic_memory.o: dynamic_memory.C   ; $(DISCRETA_CPP) $(CPPFLAGS) -c dynamic_memory.C -o $@
geo_parameter.o: geo_parameter.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c geo_parameter.C -o $@
geometric_tests.o: geometric_tests.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c geometric_tests.C -o $@
inc_gen_global.o: inc_gen_global.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c inc_gen_global.C -o $@
incidence_global.o: incidence_global.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c incidence_global.C -o $@
mckay.o: mckay.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c mckay.C -o $@
packing.o: packing.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c packing.C -o $@
point_line.o: point_line.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c point_line.C -o $@
refine_rows.o: refine_rows.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c refine_rows.C -o $@
refine_columns.o: refine_columns.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c refine_columns.C -o $@
refine_3design.o: refine_3design.C   ; $(DISCRETA_CPP) $(CPPFLAGS) -c refine_3design.C -o $@
tdo_scheme.o: tdo_scheme.C   ; $(DISCRETA_CPP) $(CPPFLAGS) -c tdo_scheme.C -o $@
tdo_data.o: tdo_data.C   ; $(DISCRETA_CPP) $(CPPFLAGS) -c tdo_data.C -o $@






cpp:
	$(UNIX_TO_DOS) dynamic_memory.C $(CPP_PREFIX)dynamic_memory.cpp
	$(UNIX_TO_DOS) geo_parameter.C $(CPP_PREFIX)geo_parameter.cpp
	$(UNIX_TO_DOS) geometric_tests.C $(CPP_PREFIX)geometric_tests.cpp
	$(UNIX_TO_DOS) inc_gen_global.C $(CPP_PREFIX)inc_gen_global.cpp
	$(UNIX_TO_DOS) incidence_global.C $(CPP_PREFIX)incidence_global.cpp
	$(UNIX_TO_DOS) mckay.C $(CPP_PREFIX)mckay.cpp
	$(UNIX_TO_DOS) nauty_interface.C $(CPP_PREFIX)nauty_interface.cpp
	$(UNIX_TO_DOS) packing.C $(CPP_PREFIX)packing.cpp
	$(UNIX_TO_DOS) point_line.C $(CPP_PREFIX)plane.cpp
	$(UNIX_TO_DOS) refine_3design.C $(CPP_PREFIX)refine_3design.cpp
	$(UNIX_TO_DOS) refine_columns.C $(CPP_PREFIX)refine_columns.cpp
	$(UNIX_TO_DOS) refine_rows.C $(CPP_PREFIX)refine_rows.cpp
	$(UNIX_TO_DOS) tdo_data.C $(CPP_PREFIX)tdo_data.cpp
	$(UNIX_TO_DOS) tdo_scheme.C $(CPP_PREFIX)tdo_scheme.cpp
	$(UNIX_TO_DOS) incidence.h ../WINDOWS/incidence.h
	mv $(CPP_PREFIX)* ../WINDOWS
	

