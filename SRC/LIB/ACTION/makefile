DISCRETA_CPP = /usr/bin/g++
CPPFLAGS = -Wall -I.. -I../GALOIS
LIB = 
LFLAGS = -lm
CPP_PREFIX=b_
UNIX_TO_DOS = ~/bin/unix2dos -n 

MODULES = \
	action.o \
	action_by_conjugation.o \
	action_by_representation.o \
	action_by_restriction.o \
	action_by_right_multiplication.o \
	action_by_subfield_structure.o \
	action_cb.o \
	action_global.o \
	action_induce.o \
	action_init.o \
	action_on_andre.o \
	action_on_bricks.o \
	action_on_determinant.o \
	action_on_cosets.o \
	action_on_factor_space.o \
	action_on_grassmannian.o \
	action_on_k_subsets.o \
	action_on_orthogonal.o \
	action_on_sets.o \
	action_on_spread_set.o \
	action_on_wedge_product.o \
	backtrack.o \
	desarguesian_spread.o \
	group.o \
	interface.o \
	interface_matrix_group.o \
	interface_perm_group.o \
	linear_group.o \
	linear_group_description.o \
	matrix_group.o \
	page_storage.o \
	perm_group.o \
	product_action.o \
	schreier.o \
	schreier_sims.o \
	schreier_vector.o \
	sims.o \
	sims2.o \
	sims_global.o \
	strong_generators.o \
	strong_generators_groups.o \
	union_find.o \
	union_find_on_k_subsets.o \
	vector_ge.o 

LIB = libaction.a

libaction.a: $(MODULES) ;
	ar ruvs $(LIB) $(MODULES)

clean:
	- rm *.o *.a *~



#files for liborbiter:
action.o: action.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action.C -o $@
action_induce.o: action_induce.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_induce.C -o $@
action_init.o: action_init.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_init.C -o $@
action_on_andre.o: action_on_andre.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_andre.C -o $@
action_on_bricks.o: action_on_bricks.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_bricks.C -o $@
action_cb.o: action_cb.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_cb.C -o $@
action_global.o: action_global.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_global.C -o $@
action_on_k_subsets.o: action_on_k_subsets.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_k_subsets.C -o $@
action_by_restriction.o: action_by_restriction.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_by_restriction.C -o $@
action_by_right_multiplication.o: action_by_right_multiplication.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_by_right_multiplication.C -o $@
action_by_conjugation.o: action_by_conjugation.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_by_conjugation.C -o $@
action_by_representation.o: action_by_representation.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_by_representation.C -o $@
action_by_subfield_structure.o: action_by_subfield_structure.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_by_subfield_structure.C -o $@
action_on_determinant.o: action_on_determinant.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_determinant.C -o $@
action_on_cosets.o: action_on_cosets.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_cosets.C -o $@
action_on_factor_space.o: action_on_factor_space.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_factor_space.C -o $@
action_on_grassmannian.o: action_on_grassmannian.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_grassmannian.C -o $@
action_on_orthogonal.o: action_on_orthogonal.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_orthogonal.C -o $@
action_on_sets.o: action_on_sets.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_sets.C -o $@
action_on_spread_set.o: action_on_spread_set.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_spread_set.C -o $@
action_on_wedge_product.o: action_on_wedge_product.C; $(DISCRETA_CPP) $(CPPFLAGS) -c action_on_wedge_product.C -o $@
backtrack.o: backtrack.C; $(DISCRETA_CPP) $(CPPFLAGS) -c backtrack.C -o $@
desarguesian_spread.o: desarguesian_spread.C; $(DISCRETA_CPP) $(CPPFLAGS) -c desarguesian_spread.C -o $@
group.o: group.C; $(DISCRETA_CPP) $(CPPFLAGS) -c group.C -o $@
interface.o: interface.C; $(DISCRETA_CPP) $(CPPFLAGS) -c interface.C -o $@
interface_matrix_group.o: interface_matrix_group.C; $(DISCRETA_CPP) $(CPPFLAGS) -c interface_matrix_group.C -o $@
interface_perm_group.o: interface_perm_group.C; $(DISCRETA_CPP) $(CPPFLAGS) -c interface_perm_group.C -o $@
linear_group.o: linear_group.C; $(DISCRETA_CPP) $(CPPFLAGS) -c linear_group.C -o $@
linear_group_description.o: linear_group_description.C; $(DISCRETA_CPP) $(CPPFLAGS) -c linear_group_description.C -o $@
matrix_group.o: matrix_group.C; $(DISCRETA_CPP) $(CPPFLAGS) -c matrix_group.C -o $@
page_storage.o: page_storage.C; $(DISCRETA_CPP) $(CPPFLAGS) -c page_storage.C -o $@
perm_group.o: perm_group.C; $(DISCRETA_CPP) $(CPPFLAGS) -c perm_group.C -o $@
product_action.o: product_action.C; $(DISCRETA_CPP) $(CPPFLAGS) -c product_action.C -o $@
schreier.o: schreier.C; $(DISCRETA_CPP) $(CPPFLAGS) -c schreier.C -o $@
schreier_sims.o: schreier_sims.C; $(DISCRETA_CPP) $(CPPFLAGS) -c schreier_sims.C -o $@
schreier_vector.o: schreier_vector.C; $(DISCRETA_CPP) $(CPPFLAGS) -c schreier_vector.C -o $@
sims.o: sims.C; $(DISCRETA_CPP) $(CPPFLAGS) -c sims.C -o $@
sims2.o: sims2.C; $(DISCRETA_CPP) $(CPPFLAGS) -c sims2.C -o $@
sims_global.o: sims_global.C; $(DISCRETA_CPP) $(CPPFLAGS) -c sims_global.C -o $@
strong_generators.o: strong_generators.C; $(DISCRETA_CPP) $(CPPFLAGS) -c strong_generators.C -o $@
strong_generators_groups.o: strong_generators_groups.C; $(DISCRETA_CPP) $(CPPFLAGS) -c strong_generators_groups.C -o $@
union_find.o: union_find.C; $(DISCRETA_CPP) $(CPPFLAGS) -c union_find.C -o $@
union_find_on_k_subsets.o: union_find_on_k_subsets.C; $(DISCRETA_CPP) $(CPPFLAGS) -c union_find_on_k_subsets.C -o $@
vector_ge.o: vector_ge.C; $(DISCRETA_CPP) $(CPPFLAGS) -c vector_ge.C -o $@


cpp:
	$(UNIX_TO_DOS) action.C $(CPP_PREFIX)action.cpp
	$(UNIX_TO_DOS) action_by_conjugation.C $(CPP_PREFIX)action_by_conjugation.cpp
	$(UNIX_TO_DOS) action_by_representation.C $(CPP_PREFIX)action_by_representation.cpp
	$(UNIX_TO_DOS) action_by_restriction.C $(CPP_PREFIX)action_by_restriction.cpp
	$(UNIX_TO_DOS) action_by_right_multiplication.C $(CPP_PREFIX)action_by_right_multiplication.cpp
	$(UNIX_TO_DOS) action_by_subfield_structure.C $(CPP_PREFIX)action_by_subfield_structure.cpp
	$(UNIX_TO_DOS) action_on_andre.C $(CPP_PREFIX)action_on_andre.cpp
	$(UNIX_TO_DOS) action_on_bricks.C $(CPP_PREFIX)action_on_bricks.cpp
	$(UNIX_TO_DOS) action_cb.C $(CPP_PREFIX)action_cb.cpp
	$(UNIX_TO_DOS) action_induce.C $(CPP_PREFIX)action_induce.cpp
	$(UNIX_TO_DOS) action_init.C $(CPP_PREFIX)action_init.cpp
	$(UNIX_TO_DOS) action_on_determinant.C $(CPP_PREFIX)action_on_determinant.cpp
	$(UNIX_TO_DOS) action_on_cosets.C $(CPP_PREFIX)action_on_cosets.cpp
	$(UNIX_TO_DOS) action_on_factor_space.C $(CPP_PREFIX)action_on_factor_space.cpp
	$(UNIX_TO_DOS) action_on_grassmannian.C $(CPP_PREFIX)action_on_grassmannian.cpp
	$(UNIX_TO_DOS) action_on_orthogonal.C $(CPP_PREFIX)action_on_orthogonal.cpp
	$(UNIX_TO_DOS) action_on_sets.C $(CPP_PREFIX)action_on_sets.cpp
	$(UNIX_TO_DOS) action_on_k_subsets.C $(CPP_PREFIX)action_on_k_subsets.cpp
	$(UNIX_TO_DOS) action_on_wedge_product.C $(CPP_PREFIX)action_on_wedge_product.cpp
	$(UNIX_TO_DOS) backtrack.C $(CPP_PREFIX)backtrack.cpp
	$(UNIX_TO_DOS) desarguesian_spread.C $(CPP_PREFIX)desarguesian_spread.cpp
	$(UNIX_TO_DOS) group.C $(CPP_PREFIX)group.cpp
	$(UNIX_TO_DOS) interface.C $(CPP_PREFIX)interface.cpp
	$(UNIX_TO_DOS) interface_matrix_group.C $(CPP_PREFIX)interface_matrix_group.cpp
	$(UNIX_TO_DOS) interface_perm_group.C $(CPP_PREFIX)interface_perm_group.cpp
	$(UNIX_TO_DOS) linear_group.C $(CPP_PREFIX)linear_group.cpp
	$(UNIX_TO_DOS) linear_group_description.C $(CPP_PREFIX)linear_group_description.cpp
	$(UNIX_TO_DOS) matrix_group.C $(CPP_PREFIX)matrix_group.cpp
	$(UNIX_TO_DOS) page_storage.C $(CPP_PREFIX)page_storage.cpp
	$(UNIX_TO_DOS) perm_group.C $(CPP_PREFIX)perm_group.cpp
	$(UNIX_TO_DOS) product_action.C $(CPP_PREFIX)product_action.cpp
	$(UNIX_TO_DOS) schreier.C $(CPP_PREFIX)schreier.cpp
	$(UNIX_TO_DOS) schreier_sims.C $(CPP_PREFIX)schreier_sims.cpp
	$(UNIX_TO_DOS) schreier_vector.C $(CPP_PREFIX)schreier_vector.cpp
	$(UNIX_TO_DOS) sims.C $(CPP_PREFIX)sims.cpp
	$(UNIX_TO_DOS) sims2.C $(CPP_PREFIX)sims2.cpp
	$(UNIX_TO_DOS) sims_global.C $(CPP_PREFIX)sims_global.cpp
	$(UNIX_TO_DOS) strong_generators.C $(CPP_PREFIX)strong_generators.cpp
	$(UNIX_TO_DOS) strong_generators_groups.C $(CPP_PREFIX)strong_generators_groups.cpp
	$(UNIX_TO_DOS) union_find.C $(CPP_PREFIX)union_find.cpp
	$(UNIX_TO_DOS) union_find_on_k_subsets.C $(CPP_PREFIX)union_find_on_k_subsets.cpp
	$(UNIX_TO_DOS) vector_ge.C $(CPP_PREFIX)vector_ge.cpp
	$(UNIX_TO_DOS) action.h ../WINDOWS/action.h
	mv $(CPP_PREFIX)* ../WINDOWS

	

