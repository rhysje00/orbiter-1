DISCRETA_CPP = /usr/bin/g++
DISCRETA_CC = /usr/bin/gcc
CPPFLAGS = -Wall
LIB = 
LFLAGS = -lm
CPP_PREFIX=a_
UNIX_TO_DOS = ~/bin/unix2dos -n 

MODULES = \
	INT_matrix.o \
	INT_vector.o \
	a_domain.o \
	andre_construction.o \
	andre_construction_line_element.o \
	andre_construction_point_element.o \
	brick_domain.o \
	buekenhout_metz.o \
	classify.o \
	clique_finder.o \
	colored_graph.o \
	combinatorics.o \
	data.o \
	data_file.o \
	decomposition.o \
	diophant.o \
	dlx.o \
	draw.o \
	fancy_set.o \
	file_output.o \
	finite_ring.o \
	finite_field.o \
	finite_field_tables.o \
	finite_field_linear_algebra.o \
	finite_field_representations.o \
	galois_global.o \
	geometric_object.o \
	geometric_operations.o \
	gl_classes.o \
	graph_layer.o \
	graph_node.o \
	grassmann.o \
	grassmann_embedded.o \
	hermitian.o \
	hjelmslev.o \
	incidence_structure.o \
	klein_correspondence.o \
	layered_graph.o \
	layered_graph_draw_options.o \
	longinteger_object.o \
	longinteger_domain.o \
	mp_graphics.o \
	memory.o \
	memory_object.o \
	mindist.o \
	nauty_interface.o \
	nauty.o \
	nautil.o \
	naugraph.o \
	norm_tables.o \
	null_polarity_generator.o \
	number_theory.o \
	orthogonal.o \
	orthogonal_points.o \
	partitionstack.o \
	plot.o \
	projective.o \
	projective_space.o \
	rainbow_cliques.o \
	rank_checker.o \
	set_of_sets.o \
	sorting.o \
	spreadsheet.o \
	subfield_structure.o \
	super_fast_hash.o \
	tensor.o \
	tree.o \
	tree_node.o \
	unipoly.o \
	unusual.o \
	util.o \
	vector_hashing.o


LIB = libgalois.a

libgalois.a: $(MODULES) ;
	ar ruvs $(LIB) $(MODULES)

clean:
	- rm *.o *.a *~



INT_matrix.o: INT_matrix.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c INT_matrix.C -o $@
INT_vector.o: INT_vector.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c INT_vector.C -o $@
a_domain.o: a_domain.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c a_domain.C -o $@
andre_construction.o: andre_construction.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c andre_construction.C -o $@
andre_construction_line_element.o: andre_construction_line_element.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c andre_construction_line_element.C -o $@
andre_construction_point_element.o: andre_construction_point_element.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c andre_construction_point_element.C -o $@
brick_domain.o: brick_domain.C; $(DISCRETA_CPP) $(CPPFLAGS) -c brick_domain.C -o $@
buekenhout_metz.o: buekenhout_metz.C; $(DISCRETA_CPP) $(CPPFLAGS) -c buekenhout_metz.C -o $@
classify.o: classify.C; $(DISCRETA_CPP) $(CPPFLAGS) -c classify.C -o $@
clique_finder.o: clique_finder.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c clique_finder.C -o $@
colored_graph.o: colored_graph.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c colored_graph.C -o $@
combinatorics.o: combinatorics.C; $(DISCRETA_CPP) $(CPPFLAGS) -c combinatorics.C -o $@
data.o: data.C data_BLT.C data_TP.C data_DH.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c data.C -o $@
data_file.o: data_file.C; $(DISCRETA_CPP) $(CPPFLAGS) -c data_file.C -o $@
decomposition.o: decomposition.C; $(DISCRETA_CPP) $(CPPFLAGS) -c decomposition.C -o $@
diophant.o: diophant.C; $(DISCRETA_CPP) $(CPPFLAGS) -O3 -c diophant.C -o $@
dlx.o: dlx.C; $(DISCRETA_CPP) $(CPPFLAGS) -O3 -c dlx.C -o $@
draw.o: draw.C; $(DISCRETA_CPP) $(CPPFLAGS) -c draw.C -o $@
fancy_set.o: fancy_set.C; $(DISCRETA_CPP) $(CPPFLAGS) -c fancy_set.C -o $@
file_output.o: file_output.C; $(DISCRETA_CPP) $(CPPFLAGS) -c file_output.C -o $@
finite_ring.o: finite_ring.C; $(DISCRETA_CPP) $(CPPFLAGS) -c finite_ring.C -o $@
finite_field.o: finite_field.C; $(DISCRETA_CPP) $(CPPFLAGS) -c finite_field.C -o $@
finite_field_tables.o: finite_field_tables.C; $(DISCRETA_CPP) $(CPPFLAGS) -c finite_field_tables.C -o $@
finite_field_linear_algebra.o: finite_field_linear_algebra.C; $(DISCRETA_CPP) $(CPPFLAGS) -c finite_field_linear_algebra.C -o $@
finite_field_representations.o: finite_field_representations.C; $(DISCRETA_CPP) $(CPPFLAGS) -c finite_field_representations.C -o $@
galois_global.o: galois_global.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c galois_global.C -o $@
geometric_object.o: geometric_object.C; $(DISCRETA_CPP) $(CPPFLAGS) -c geometric_object.C -o $@
geometric_operations.o: geometric_operations.C; $(DISCRETA_CPP) $(CPPFLAGS) -c geometric_operations.C -o $@
gl_classes.o: gl_classes.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c gl_classes.C -o $@
graph_layer.o: graph_layer.C; $(DISCRETA_CPP) $(CPPFLAGS) -c graph_layer.C -o $@
graph_node.o: graph_node.C; $(DISCRETA_CPP) $(CPPFLAGS) -c graph_node.C -o $@
grassmann.o: grassmann.C; $(DISCRETA_CPP) $(CPPFLAGS) -c grassmann.C -o $@
grassmann_embedded.o: grassmann_embedded.C; $(DISCRETA_CPP) $(CPPFLAGS) -c grassmann_embedded.C -o $@
hermitian.o: hermitian.C; $(DISCRETA_CPP) $(CPPFLAGS) -c hermitian.C -o $@
hjelmslev.o: hjelmslev.C; $(DISCRETA_CPP) $(CPPFLAGS) -c hjelmslev.C -o $@
incidence_structure.o: incidence_structure.C; $(DISCRETA_CPP) $(CPPFLAGS) -c incidence_structure.C -o $@
klein_correspondence.o: klein_correspondence.C; $(DISCRETA_CPP) $(CPPFLAGS) -c klein_correspondence.C -o $@
layered_graph.o: layered_graph.C; $(DISCRETA_CPP) $(CPPFLAGS) -c layered_graph.C -o $@
layered_graph_draw_options.o: layered_graph_draw_options.C; $(DISCRETA_CPP) $(CPPFLAGS) -c layered_graph_draw_options.C -o $@
longinteger_object.o: longinteger_object.C; $(DISCRETA_CPP) $(CPPFLAGS) -c longinteger_object.C -o $@
longinteger_domain.o: longinteger_domain.C; $(DISCRETA_CPP) $(CPPFLAGS) -c longinteger_domain.C -o $@
mp_graphics.o: mp_graphics.C; $(DISCRETA_CPP) $(CPPFLAGS) -c mp_graphics.C -o $@
memory.o: memory.C; $(DISCRETA_CPP) $(CPPFLAGS) -c memory.C -o $@
memory_object.o: memory_object.C; $(DISCRETA_CPP) $(CPPFLAGS) -c memory_object.C -o $@
mindist.o: mindist.C; $(DISCRETA_CPP) $(CPPFLAGS) -c mindist.C -o $@
nauty_interface.o: nauty_interface.C nauty.h ; $(DISCRETA_CPP) $(CPPFLAGS) -c nauty_interface.C -o $@
nauty.o: nauty.c nauty.h; $(DISCRETA_CC) $(CCFLAGS) -c nauty.c -o $@
nautil.o: nautil.c nauty.h; $(DISCRETA_CC) $(CCFLAGS) -c nautil.c -o $@
naugraph.o: naugraph.c nauty.h; $(DISCRETA_CC) $(CCFLAGS) -c naugraph.c -o $@
norm_tables.o: norm_tables.C; $(DISCRETA_CPP) $(CPPFLAGS) -c norm_tables.C -o $@
null_polarity_generator.o: null_polarity_generator.C; $(DISCRETA_CPP) $(CPPFLAGS) -c null_polarity_generator.C -o $@
number_theory.o: number_theory.C; $(DISCRETA_CPP) $(CPPFLAGS) -c number_theory.C -o $@
orthogonal.o: orthogonal.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c orthogonal.C -o $@
orthogonal_points.o: orthogonal_points.C; $(DISCRETA_CPP) $(CPPFLAGS) -c orthogonal_points.C -o $@
partitionstack.o: partitionstack.C; $(DISCRETA_CPP) $(CPPFLAGS) -c partitionstack.C -o $@
plot.o: plot.C; $(DISCRETA_CPP) $(CPPFLAGS) -c plot.C -o $@
projective.o: projective.C; $(DISCRETA_CPP) $(CPPFLAGS) -c projective.C -o $@
projective_space.o: projective_space.C; $(DISCRETA_CPP) $(CPPFLAGS) -c projective_space.C -o $@
rainbow_cliques.o: rainbow_cliques.C; $(DISCRETA_CPP) $(CPPFLAGS) -c rainbow_cliques.C -o $@
rank_checker.o: rank_checker.C; $(DISCRETA_CPP) $(CPPFLAGS) -c rank_checker.C -o $@
set_of_sets.o: set_of_sets.C; $(DISCRETA_CPP) $(CPPFLAGS) -c set_of_sets.C -o $@
sorting.o: sorting.C; $(DISCRETA_CPP) $(CPPFLAGS) -c sorting.C -o $@
spreadsheet.o: spreadsheet.C; $(DISCRETA_CPP) $(CPPFLAGS) -c spreadsheet.C -o $@
subfield_structure.o: subfield_structure.C; $(DISCRETA_CPP) $(CPPFLAGS) -c subfield_structure.C -o $@
super_fast_hash.o: super_fast_hash.C  ; $(DISCRETA_CPP) $(CPPFLAGS) -c super_fast_hash.C -o $@
tensor.o: tensor.C; $(DISCRETA_CPP) $(CPPFLAGS) -c tensor.C -o $@
tree.o: tree.C; $(DISCRETA_CPP) $(CPPFLAGS) -c tree.C -o $@
tree_node.o: tree_node.C; $(DISCRETA_CPP) $(CPPFLAGS) -c tree_node.C -o $@
unipoly.o: unipoly.C; $(DISCRETA_CPP) $(CPPFLAGS) -c unipoly.C -o $@
unusual.o: unusual.C; $(DISCRETA_CPP) $(CPPFLAGS) -c unusual.C -o $@
util.o: util.C; $(DISCRETA_CPP) $(CPPFLAGS) -c util.C -o $@
vector_hashing.o: vector_hashing.C; $(DISCRETA_CPP) $(CPPFLAGS) -c vector_hashing.C -o $@



cpp:
	$(UNIX_TO_DOS) INT_matrix.C $(CPP_PREFIX)INT_matrix.cpp
	$(UNIX_TO_DOS) INT_vector.C $(CPP_PREFIX)INT_vector.cpp
	$(UNIX_TO_DOS) a_domain.C $(CPP_PREFIX)a_domain.cpp
	$(UNIX_TO_DOS) andre_construction.C $(CPP_PREFIX)andre_construction.cpp
	$(UNIX_TO_DOS) andre_construction_line_element.C $(CPP_PREFIX)andre_construction_line_element.cpp
	$(UNIX_TO_DOS) andre_construction_point_element.C $(CPP_PREFIX)andre_construction_point_element.cpp
	$(UNIX_TO_DOS) brick_domain.C $(CPP_PREFIX)brick_domain.cpp
	$(UNIX_TO_DOS) buekenhout_metz.C $(CPP_PREFIX)buekenhout_metz.cpp
	$(UNIX_TO_DOS) classify.C $(CPP_PREFIX)classify.cpp
	$(UNIX_TO_DOS) clique_finder.C $(CPP_PREFIX)clique_finder.cpp
	$(UNIX_TO_DOS) colored_graph.C $(CPP_PREFIX)colored_graph.cpp
	$(UNIX_TO_DOS) combinatorics.C $(CPP_PREFIX)combinatorics.cpp
	$(UNIX_TO_DOS) data.C $(CPP_PREFIX)data.cpp
	$(UNIX_TO_DOS) data_BLT.C $(CPP_PREFIX)data_BLT.cpp
	$(UNIX_TO_DOS) data_DH.C $(CPP_PREFIX)data_DH.cpp
	$(UNIX_TO_DOS) data_TP.C $(CPP_PREFIX)data_TP.cpp
	$(UNIX_TO_DOS) data_file.C $(CPP_PREFIX)data_file.cpp
	$(UNIX_TO_DOS) data_hyperovals.C $(CPP_PREFIX)data_hyperovals.cpp
	$(UNIX_TO_DOS) data_packings_PG_3_3.C $(CPP_PREFIX)data_packings_PG_3_3.cpp
	$(UNIX_TO_DOS) decomposition.C $(CPP_PREFIX)decomposition.cpp
	$(UNIX_TO_DOS) diophant.C $(CPP_PREFIX)diophant.cpp
	$(UNIX_TO_DOS) dlx.C $(CPP_PREFIX)dlx.cpp
	$(UNIX_TO_DOS) draw.C $(CPP_PREFIX)draw.cpp
	$(UNIX_TO_DOS) fancy_set.C $(CPP_PREFIX)fancy_set.cpp
	$(UNIX_TO_DOS) file_output.C $(CPP_PREFIX)file_output.cpp
	$(UNIX_TO_DOS) finite_field.C $(CPP_PREFIX)finite_field.cpp
	$(UNIX_TO_DOS) finite_field_linear_algebra.C $(CPP_PREFIX)finite_field_linear_algebra.cpp
	$(UNIX_TO_DOS) finite_field_representations.C $(CPP_PREFIX)finite_field_representations.cpp
	$(UNIX_TO_DOS) finite_field_tables.C $(CPP_PREFIX)finite_field_tables.cpp
	$(UNIX_TO_DOS) finite_ring.C $(CPP_PREFIX)finite_ring.cpp
	$(UNIX_TO_DOS) galois_global.C $(CPP_PREFIX)galois_global.cpp
	$(UNIX_TO_DOS) geometric_object.C $(CPP_PREFIX)geometric_object.cpp
	$(UNIX_TO_DOS) geometric_operations.C $(CPP_PREFIX)geometric_operations.cpp
	$(UNIX_TO_DOS) gl_classes.C $(CPP_PREFIX)gl_classes.cpp
	$(UNIX_TO_DOS) graph_layer.C $(CPP_PREFIX)graph_layer.cpp
	$(UNIX_TO_DOS) graph_node.C $(CPP_PREFIX)graph_node.cpp
	$(UNIX_TO_DOS) grassmann.C $(CPP_PREFIX)grassmann.cpp
	$(UNIX_TO_DOS) grassmann_embedded.C $(CPP_PREFIX)grassmann_embedded.cpp
	$(UNIX_TO_DOS) hermitian.C $(CPP_PREFIX)hermitian.cpp
	$(UNIX_TO_DOS) hjelmslev.C $(CPP_PREFIX)hjelmslev.cpp
	$(UNIX_TO_DOS) incidence_structure.C $(CPP_PREFIX)incidence_structure.cpp
	$(UNIX_TO_DOS) klein_correspondence.C $(CPP_PREFIX)klein_correspondence.cpp
	$(UNIX_TO_DOS) layered_graph.C $(CPP_PREFIX)layered_graph.cpp
	$(UNIX_TO_DOS) layered_graph_draw_options.C $(CPP_PREFIX)layered_graph_draw_options.cpp
	$(UNIX_TO_DOS) longinteger_domain.C $(CPP_PREFIX)longinteger_domain.cpp
	$(UNIX_TO_DOS) longinteger_object.C $(CPP_PREFIX)longinteger_object.cpp
	$(UNIX_TO_DOS) memory.C $(CPP_PREFIX)memory.cpp
	$(UNIX_TO_DOS) memory_object.C $(CPP_PREFIX)memory_object.cpp
	$(UNIX_TO_DOS) mp_graphics.C $(CPP_PREFIX)mp_graphics.cpp
	$(UNIX_TO_DOS) nauty.h ../WINDOWS/nauty.h
	$(UNIX_TO_DOS) nauty.c ../WINDOWS/$(CPP_PREFIX)nauty.c
	$(UNIX_TO_DOS) nautil.c ../WINDOWS/$(CPP_PREFIX)nautil.c
	$(UNIX_TO_DOS) naugraph.c ../WINDOWS/$(CPP_PREFIX)naugraph.c
	$(UNIX_TO_DOS) nauty_interface.C $(CPP_PREFIX)nauty_interface.cpp
	$(UNIX_TO_DOS) norm_tables.C $(CPP_PREFIX)norm_tables.cpp
	$(UNIX_TO_DOS) null_polarity_generator.C $(CPP_PREFIX)null_polarity_generator.cpp
	$(UNIX_TO_DOS) number_theory.C $(CPP_PREFIX)number_theory.cpp
	$(UNIX_TO_DOS) orthogonal.C $(CPP_PREFIX)orthogonal.cpp
	$(UNIX_TO_DOS) orthogonal_points.C $(CPP_PREFIX)orthogonal_points.cpp
	$(UNIX_TO_DOS) partitionstack.C $(CPP_PREFIX)partitionstack.cpp
	$(UNIX_TO_DOS) plot.C $(CPP_PREFIX)plot.cpp
	$(UNIX_TO_DOS) projective.C $(CPP_PREFIX)projective.cpp
	$(UNIX_TO_DOS) projective_space.C $(CPP_PREFIX)projective_space.cpp
	$(UNIX_TO_DOS) rainbow_cliques.C $(CPP_PREFIX)rainbow_cliques.cpp
	$(UNIX_TO_DOS) rank_checker.C $(CPP_PREFIX)rank_checker.cpp
	$(UNIX_TO_DOS) set_of_sets.C $(CPP_PREFIX)set_of_sets.cpp
	$(UNIX_TO_DOS) sorting.C $(CPP_PREFIX)sorting.cpp
	$(UNIX_TO_DOS) spreadsheet.C $(CPP_PREFIX)spreadsheet.cpp
	$(UNIX_TO_DOS) subfield_structure.C $(CPP_PREFIX)subfield_structure.cpp
	$(UNIX_TO_DOS) super_fast_hash.C $(CPP_PREFIX)super_fast_hash.cpp
	$(UNIX_TO_DOS) tensor.C $(CPP_PREFIX)tensor.cpp
	$(UNIX_TO_DOS) tree.C $(CPP_PREFIX)tree.cpp
	$(UNIX_TO_DOS) tree_node.C $(CPP_PREFIX)tree_node.cpp
	$(UNIX_TO_DOS) unipoly.C $(CPP_PREFIX)unipoly.cpp
	$(UNIX_TO_DOS) unusual.C $(CPP_PREFIX)unusual.cpp
	$(UNIX_TO_DOS) util.C $(CPP_PREFIX)util.cpp
	$(UNIX_TO_DOS) vector_hashing.C $(CPP_PREFIX)vector_hashing.cpp
	$(UNIX_TO_DOS) galois.h ../WINDOWS/galois.h
	mv $(CPP_PREFIX)* ../WINDOWS


