DISCRETA_CPP = /usr/bin/g++
CPPFLAGS = -Wall -I.. -I../GALOIS
LIB = 
LFLAGS = -lm
CPP_PREFIX=f_
UNIX_TO_DOS = ~/bin/unix2dos -n 

MODULES = \
	analyze_group.o \
	choose_points_or_lines.o \
	decomposition.o \
	elliptic_curve.o \
	exact_cover.o \
	exact_cover_arguments.o \
	extra.o \
	factor_group.o \
	incidence_structure.o \
	isomorph.o \
	isomorph_arguments.o \
	isomorph_database.o \
	isomorph_files.o \
	isomorph_global.o \
	isomorph_testing.o \
	isomorph_trace.o \
	knarr.o \
	orbit_of_sets.o \
	orbit_of_subspaces.o \
	orbit_rep.o \
	polar.o \
	projective_space.o \
	recoordinatize.o \
	representatives.o \
	search_blocking_set.o \
	singer_cycle.o \
	subspace_orbits.o \
	translation_plane.o \
	translation_plane2.o \
	translation_plane_via_andre_model.o \
	w3q.o \
	young.o
		

LIB = libtoplevel.a

libtoplevel.a: $(MODULES) ;
	ar ruvs $(LIB) $(MODULES)

clean:
	- rm *.o *.a *~

# BASIC_FILES:
analyze_group.o: analyze_group.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c analyze_group.C -o $@
choose_points_or_lines.o: choose_points_or_lines.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c choose_points_or_lines.C -o $@
decomposition.o: decomposition.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c decomposition.C -o $@
elliptic_curve.o: elliptic_curve.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c elliptic_curve.C -o $@
exact_cover.o: exact_cover.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c exact_cover.C -o $@
exact_cover_arguments.o: exact_cover_arguments.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c exact_cover_arguments.C -o $@
extra.o: extra.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c extra.C -o $@
factor_group.o: factor_group.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c factor_group.C -o $@
incidence_structure.o: incidence_structure.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c incidence_structure.C -o $@
isomorph.o: isomorph.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c isomorph.C -o $@
isomorph_arguments.o: isomorph_arguments.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c isomorph_arguments.C -o $@
isomorph_database.o: isomorph_database.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c isomorph_database.C -o $@
isomorph_files.o: isomorph_files.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c isomorph_files.C -o $@
isomorph_global.o: isomorph_global.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c isomorph_global.C -o $@
isomorph_testing.o: isomorph_testing.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c isomorph_testing.C -o $@
isomorph_trace.o: isomorph_trace.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c isomorph_trace.C -o $@
knarr.o: knarr.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c knarr.C -o $@
orbit_of_sets.o: orbit_of_sets.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c orbit_of_sets.C -o $@
orbit_of_subspaces.o: orbit_of_subspaces.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c orbit_of_subspaces.C -o $@
orbit_rep.o: orbit_rep.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c orbit_rep.C -o $@
polar.o: polar.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c polar.C -o $@
projective_space.o: projective_space.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c projective_space.C -o $@
recoordinatize.o: recoordinatize.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c recoordinatize.C -o $@
representatives.o: representatives.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c representatives.C -o $@
search_blocking_set.o: search_blocking_set.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c search_blocking_set.C -o $@
singer_cycle.o: singer_cycle.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c singer_cycle.C -o $@
subspace_orbits.o: subspace_orbits.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c subspace_orbits.C -o $@
translation_plane.o: translation_plane.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c translation_plane.C -o $@
translation_plane2.o: translation_plane2.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c translation_plane2.C -o $@
translation_plane_via_andre_model.o: translation_plane_via_andre_model.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c translation_plane_via_andre_model.C -o $@
w3q.o: w3q.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c w3q.C -o $@
young.o: young.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c young.C -o $@



cpp:
	$(UNIX_TO_DOS) analyze_group.C $(CPP_PREFIX)analyze_group.cpp
	$(UNIX_TO_DOS) buekenhout_metz.C $(CPP_PREFIX)buekenhout_metz.cpp
	$(UNIX_TO_DOS) choose_points_or_lines.C $(CPP_PREFIX)choose_points_or_lines.cpp
	$(UNIX_TO_DOS) decomposition.C $(CPP_PREFIX)decomposition.cpp
	$(UNIX_TO_DOS) elliptic_curve.C $(CPP_PREFIX)elliptic_curve.cpp
	$(UNIX_TO_DOS) exact_cover.C $(CPP_PREFIX)exact_cover.cpp
	$(UNIX_TO_DOS) exact_cover_arguments.C $(CPP_PREFIX)exact_cover_arguments.cpp
	$(UNIX_TO_DOS) extra.C $(CPP_PREFIX)extra.cpp
	$(UNIX_TO_DOS) factor_group.C $(CPP_PREFIX)factor_group.cpp
	$(UNIX_TO_DOS) incidence_structure.C $(CPP_PREFIX)incidence_structure.cpp
	$(UNIX_TO_DOS) isomorph.C $(CPP_PREFIX)isomorph.cpp
	$(UNIX_TO_DOS) isomorph_arguments.C $(CPP_PREFIX)isomorph_arguments.cpp
	$(UNIX_TO_DOS) isomorph_database.C $(CPP_PREFIX)isomorph_database.cpp
	$(UNIX_TO_DOS) isomorph_files.C $(CPP_PREFIX)isomorph_files.cpp
	$(UNIX_TO_DOS) isomorph_global.C $(CPP_PREFIX)isomorph_global.cpp
	$(UNIX_TO_DOS) isomorph_testing.C $(CPP_PREFIX)isomorph_testing.cpp
	$(UNIX_TO_DOS) isomorph_trace.C $(CPP_PREFIX)isomorph_trace.cpp
	$(UNIX_TO_DOS) knarr.C $(CPP_PREFIX)knarr.cpp
	$(UNIX_TO_DOS) orbit_of_sets.C $(CPP_PREFIX)orbit_of_sets.cpp
	$(UNIX_TO_DOS) orbit_of_subspaces.C $(CPP_PREFIX)orbit_of_subspaces.cpp
	$(UNIX_TO_DOS) orbit_rep.C $(CPP_PREFIX)orbit_rep.cpp
	$(UNIX_TO_DOS) polar.C $(CPP_PREFIX)polar.cpp
	$(UNIX_TO_DOS) projective_space.C $(CPP_PREFIX)projective_space.cpp
	$(UNIX_TO_DOS) recoordinatize.C $(CPP_PREFIX)recoordinatize.cpp
	$(UNIX_TO_DOS) representatives.C $(CPP_PREFIX)representatives.cpp
	$(UNIX_TO_DOS) search_blocking_set.C $(CPP_PREFIX)search_blocking_set.cpp
	$(UNIX_TO_DOS) singer_cycle.C $(CPP_PREFIX)singer_cycle.cpp
	$(UNIX_TO_DOS) subspace_orbits.C $(CPP_PREFIX)subspace_orbits.cpp
	$(UNIX_TO_DOS) translation_plane.C $(CPP_PREFIX)translation_plane.cpp
	$(UNIX_TO_DOS) translation_plane2.C $(CPP_PREFIX)translation_plane2.cpp
	$(UNIX_TO_DOS) translation_plane_via_andre_model.C $(CPP_PREFIX)translation_plane_via_andre_model.cpp
	$(UNIX_TO_DOS) w3q.C $(CPP_PREFIX)w3q.cpp
	$(UNIX_TO_DOS) young.C $(CPP_PREFIX)young.cpp
	$(UNIX_TO_DOS) top_level.h ../WINDOWS/top_level.h
	mv $(CPP_PREFIX)* ../WINDOWS


