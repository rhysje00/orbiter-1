// action_init.C
//
// Anton Betten
// 1/1/2009

#include "galois.h"
#include "action.h"

void action::init_orthogonal_group(INT epsilon, 
	INT n, finite_field *F, 
	INT f_on_points, INT f_on_lines, INT f_on_points_and_lines, 
	INT f_semilinear, 
	INT f_basis, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	action *A;
	orthogonal *O;
	action_on_orthogonal *AO;
	INT q = F->q;

	if (f_v) {
		cout << "action::init_orthogonal_group" << endl;
		}
	A = new action;
	if (f_vv) {
		cout << "action::init_orthogonal_group before A->init_projective_group" << endl;
		}
	A->init_projective_group(n, F, f_semilinear, TRUE /* f_basis */, verbose_level - 2);

	O = new orthogonal;
	if (f_vv) {
		cout << "action::init_orthogonal_group before O->init" << endl;
		}
	O->init(epsilon, n, F, 0 /*verbose_level*/);

	AO = new action_on_orthogonal;
	if (f_vv) {
		cout << "action::init_orthogonal_group before AO->init" << endl;
		}
	AO->init(A, O, f_on_points, f_on_lines, f_on_points_and_lines, verbose_level - 2);
	
	type_G = action_on_orthogonal_t;
	G.AO = AO;

	f_has_subaction = TRUE;
	subaction = A;
	degree = AO->degree;
	low_level_point_size = A->low_level_point_size;
	elt_size_in_INT = A->elt_size_in_INT;
	coded_elt_size_in_char = A->coded_elt_size_in_char;
	allocate_element_data();
	init_function_pointers_induced_action();


	sprintf(group_prefix, "O%ld_%ld_%ld", epsilon, n, q);
	sprintf(label, "O^%s(%ld,%ld)", plus_minus_string(epsilon), n, q);
	sprintf(label_tex, "O^{%s}(%ld,%ld)", plus_minus_string(epsilon), n, q);



	if (f_basis) {
		longinteger_object target_go;
	
		order_PO_epsilon(f_semilinear, epsilon, n - 1, F->q, target_go, verbose_level);

		if (f_vv) {
			cout << "action::init_orthogonal_group before create_orthogonal_group" << endl;
			}
		create_orthogonal_group(A /*subaction*/, 
			TRUE /* f_has_target_go */, target_go, 
			callback_choose_random_generator_orthogonal, 
			verbose_level - 2);
		if (f_vv) {
			cout << "action::init_orthogonal_group create_orthogonal_group done" << endl;
			}
		}

#if 0
	if (f_vv) {
		cout << "action::init_orthogonal_group before induced_action_on_orthogonal" << endl;
		}
	induced_action_on_orthogonal(A, 
		AO, 
		FALSE /* f_induce_action */, NULL, 
		verbose_level - 2);
#endif


	if (f_v) {
		cout << "action::init_orthogonal_group done" << endl;
		}
}

void action::init_BLT(finite_field *F, INT f_basis, INT f_init_hash_table, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT p, hh, epsilon, n;
	INT f_semilinear = FALSE;
	
	if (f_v) {
		cout << "action::init_BLT q=" << F->q << " f_init_hash_table=" << f_init_hash_table << endl;
		cout << "f_basis=" << f_basis << endl;
		cout << "verbose_level=" << verbose_level << endl;
		}
	is_prime_power(F->q, p, hh);
	if (hh > 1)
		f_semilinear = TRUE;
	else
		f_semilinear = FALSE;
	epsilon = 0;
	n = 5;


	if (f_v) {
		cout << "action::init_BLT before init_orthogonal_group" << endl;
		}
	init_orthogonal_group(epsilon, n, F, 
		TRUE /* f_on_points */, FALSE /* f_on_lines */, FALSE /* f_on_points_and_lines */, 
		f_semilinear, 
		f_basis, 
		verbose_level - 2);
	if (f_v) {
		cout << "action::init_BLT after init_orthogonal_group" << endl;
		}



	if (f_v) {
		cout << "action::init_BLT computing lex least base" << endl;
		}
	lex_least_base_in_place(verbose_level - 2);
	if (f_v) {
		cout << "action::init_BLT computing lex least base done" << endl;
		cout << "base: ";
		INT_vec_print(cout, base, base_len);
		cout << endl;
		}

	if (f_v) {
		print_base();
		}


	if (f_has_strong_generators) {
		if (f_v) {
			cout << "action::init_BLT strong generators have been computed" << endl;
			}
		if (f_vv) {
			Strong_gens->print_generators();
			}
		}
	else {
		cout << "action::init_BLT we don't have strong generators" << endl;
		exit(1);
		}

#if 0
	if (f_init_hash_table) {
		matrix_group *M;
		orthogonal *O;

		M = subaction->G.matrix_grp;
		O = M->O;
	
		if (f_v) {
			cout << "calling init_hash_table_parabolic" << endl;
			}
		init_hash_table_parabolic(*O->F, 4, 0 /* verbose_level */);
		}
#endif

	if (f_v) {
		print_info();
		}
	if (f_v) {
		cout << "action::init_BLT done" << endl;
		}
}

void action::init_group_from_strong_generators(vector_ge *gens, sims *K, 
	INT given_base_length, INT *given_base,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	sims *G;
	longinteger_object G_order;
	INT i;
	
	
	if (f_v) {
		cout << "action::init_group_from_strong_generators" << endl;
		}
	strcpy(label, "from sgs");
	strcpy(label_tex, "from sgs");
	if (f_vv) {
		cout << "generators are" << endl;
		gens->print(cout);
		cout << endl;
		}
	

	if (f_vv) {
		cout << "calling allocate_base_data, initial base:";
		INT_vec_print(cout, given_base, given_base_length);
		cout << " of length " << given_base_length << endl;
		}
	allocate_base_data(given_base_length);
	base_len = given_base_length;
	
	// init base:
	for (i = 0; i < base_len; i++) {
		base[i] = given_base[i];
		}



	if (f_vv) {
		cout << "action::init_group_from_strong_generators, now trying to set up the group from the given generators" << endl;
		}
	
	G = NEW_OBJECT(sims);
	
	G->init(this);
	G->init_trivial_group(verbose_level - 1);
	G->group_order(G_order);
	
	G->build_up_group_from_generators(K, gens, 
		FALSE, NULL, /* target_go */
		FALSE /* f_override_choose_next_base_point */,
		NULL, 
		verbose_level);
	
	G->group_order(G_order);


	if (f_vvv) {
		//G.print(TRUE);
		//G.print_generator_depth_and_perm();
		}

	if (f_v) {
		cout << "init_group_from_strong_generators: found a group of order " << G_order << endl;
		if (f_vv) {
			cout << "transversal lengths:" << endl;
			INT_vec_print(cout, G->orbit_len, base_len);
			cout << endl;
			}
		}
	
	if (f_vv) {
		cout << "init_sims" << endl;
		}
	init_sims(G, 0/*verbose_level - 1*/);
	if (f_vv) {
		cout << "after init_sims" << endl;
		}
	compute_strong_generators_from_sims(0/*verbose_level - 2*/);

	if (f_v) {
		print_info();
		}
}

void action::init_projective_special_group(INT n, finite_field *F, INT f_semilinear, INT f_basis, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "action::init_projective_special_group" << endl;
		cout << "n=" << n << " q=" << F->q << endl;
		cout << "f_semilinear=" << f_semilinear << endl;
		cout << "f_basis=" << f_basis << endl;
		}

	init_projective_group(n, F, f_semilinear, f_basis, verbose_level);

	{
		action A_on_det;
		longinteger_object go;
		strong_generators *gens;
		sims *Sims2;
		
		gens = new strong_generators;
		
		if (f_v) {
			cout << "action::init_projective_special_group computing intersection with special linear group" << endl;
			}
		A_on_det.induced_action_on_determinant(Sims, verbose_level);
		if (f_v) {
			cout << "action::init_projective_special_group  induced_action_on_determinant finished" << endl;
			A_on_det.Kernel->group_order(go);
			cout << "action::init_projective_special_group  intersection has order " << go << endl;
			}
		gens->init_from_sims(A_on_det.Kernel, verbose_level - 1);
		
		
		Sims2 = gens->create_sims(verbose_level - 1);

		delete gens;
		init_sims(Sims2, verbose_level);
	}

	if (f_v) {
		cout << "action::init_projective_special_group done" << endl;
		}
}

void action::init_projective_group(INT n, finite_field *F, INT f_semilinear, INT f_basis, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);
	matrix_group *M;
	
	if (f_v) {
		cout << "action::init_projective_group" << endl;
		cout << "n=" << n << " q=" << F->q << endl;
		cout << "f_semilinear=" << f_semilinear << endl;
		cout << "f_basis=" << f_basis << endl;
		}

	M = NEW_OBJECT(matrix_group);



	type_G = matrix_group_t;
	G.matrix_grp = M;
	f_allocated = TRUE;

	f_is_linear = TRUE;
	dimension = n;
		
	if (f_v) {
		cout << "action::init_projective_group before M->init_projective_group" << endl;
		}
	M->init_projective_group(n, F, f_semilinear, this, verbose_level - 1);


	low_level_point_size = M->low_level_point_size;
	if (f_v) {
		cout << "action::init_projective_group low_level_point_size=" << low_level_point_size<< endl;
		}
	strcpy(label, M->label);
	strcpy(label_tex, M->label_tex);
	if (f_v) {
		cout << "action::init_projective_group label=" << label << endl;
		}

	degree = M->degree;
	make_element_size = M->elt_size_INT_half;

	init_function_pointers_matrix_group();

	elt_size_in_INT = M->elt_size_INT;
	coded_elt_size_in_char = M->char_per_elt;
	allocate_element_data();

	strcpy(group_prefix, label);

	if (f_basis) {
		setup_linear_group_from_strong_generators(M, verbose_level);
		}
	if (f_v) {
		cout << "action::init_projective_group, finished setting up " << group_prefix;
		cout << ", a permutation group of degree " << degree << " ";
		cout << "and of order ";
		print_group_order(cout);
		cout << endl;
		//cout << "make_element_size=" << make_element_size << endl;
		//cout << "base_len=" << base_len << endl;
		//cout << "f_semilinear=" << f_semilinear << endl;
		}	
}

void action::init_affine_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT f_basis, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);
	matrix_group *M;
	
	if (f_v) {
		cout << "action::init_affine_group" << endl;
		cout << "n=" << n << " q=" << F->q << endl;
		cout << "f_semilinear=" << f_semilinear << endl;
		cout << "f_basis=" << f_basis << endl;
		}

	M = NEW_OBJECT(matrix_group);



	type_G = matrix_group_t;
	G.matrix_grp = M;
	f_allocated = TRUE;

	f_is_linear = TRUE;
	dimension = n;
		
	M->init_affine_group(n, F, f_semilinear, /*f_no_translations,*/ this, verbose_level - 1);


	low_level_point_size = M->low_level_point_size;
	if (f_v) {
		cout << "action::init_affine_group low_level_point_size=" << low_level_point_size<< endl;
		}
	strcpy(label, M->label);
	strcpy(label_tex, M->label_tex);
	if (f_v) {
		cout << "action::init_affine_group label=" << label << endl;
		}

	degree = M->degree;
	make_element_size = M->elt_size_INT_half;

	init_function_pointers_matrix_group();

	elt_size_in_INT = M->elt_size_INT;
	coded_elt_size_in_char = M->char_per_elt;
	allocate_element_data();

	strcpy(group_prefix, label);

	if (f_basis) {
		setup_linear_group_from_strong_generators(M, verbose_level);
		}
	if (f_v) {
		cout << "action::init_affine_group, finished setting up " << group_prefix;
		cout << ", a permutation group of degree " << degree << " ";
		cout << "and of order ";
		print_group_order(cout);
		cout << endl;
		//cout << "make_element_size=" << make_element_size << endl;
		//cout << "base_len=" << base_len << endl;
		//cout << "f_semilinear=" << f_semilinear << endl;
		}	
}

void action::init_general_linear_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT f_basis, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);
	matrix_group *M;
	
	if (f_v) {
		cout << "action::init_general_linear_group" << endl;
		cout << "n=" << n << " q=" << F->q << endl;
		cout << "f_semilinear=" << f_semilinear << endl;
		cout << "f_basis=" << f_basis << endl;
		}

	M = NEW_OBJECT(matrix_group);



	type_G = matrix_group_t;
	G.matrix_grp = M;
	f_allocated = TRUE;

	f_is_linear = TRUE;
	dimension = n;
		
	M->init_general_linear_group(n, F, f_semilinear, this, verbose_level - 1);


	low_level_point_size = M->low_level_point_size;
	if (f_v) {
		cout << "action::init_general_linear_group low_level_point_size=" << low_level_point_size<< endl;
		}
	strcpy(label, M->label);
	strcpy(label_tex, M->label_tex);
	if (f_v) {
		cout << "action::init_general_linear_group label=" << label << endl;
		}

	degree = M->degree;
	make_element_size = M->elt_size_INT_half;

	init_function_pointers_matrix_group();

	elt_size_in_INT = M->elt_size_INT;
	coded_elt_size_in_char = M->char_per_elt;
	allocate_element_data();

	strcpy(group_prefix, label);

	if (f_basis) {
		setup_linear_group_from_strong_generators(M, verbose_level);
		}
	if (f_v) {
		cout << "action::init_general_linear_group, finished setting up " << group_prefix;
		cout << ", a permutation group of degree " << degree << " ";
		cout << "and of order ";
		print_group_order(cout);
		cout << endl;
		//cout << "make_element_size=" << make_element_size << endl;
		//cout << "base_len=" << base_len << endl;
		//cout << "f_semilinear=" << f_semilinear << endl;
		}	
}

void action::setup_linear_group_from_strong_generators(matrix_group *M, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);


	if (f_v) {
		cout << "action::setup_linear_group_from_strong_generators setting up a basis" << endl;
		cout << "action::setup_linear_group_from_strong_generators before init_matrix_group_strong_generators_builtin" << endl;
		}
	init_matrix_group_strong_generators_builtin(M, verbose_level - 2);
		// see below
	if (f_v) {
		cout << "action::setup_linear_group_from_strong_generators after init_matrix_group_strong_generators_builtin" << endl;
		}

	if (f_v) {
		cout << "action::setup_linear_group_from_strong_generators before S->compute_base_orbits_known_length" << endl;
		}
	sims *S;
	//longinteger_object go1;

	S = NEW_OBJECT(sims);

	S->init(this);
	S->init_generators(*Strong_gens->gens, 0/*verbose_level*/);
	S->compute_base_orbits_known_length(transversal_length, verbose_level - 1);
	//S->group_order(go1);


	if (f_v) {
		cout << "action::setup_linear_group_from_strong_generators before init_sims" << endl;
		}

	init_sims(S, verbose_level);

	if (f_v) {
		cout << "action::setup_linear_group_from_strong_generators after init_sims" << endl;
		}
}

void action::init_matrix_group_strong_generators_builtin(matrix_group *M, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT n, q;
	finite_field *F;
	INT *data;
	INT size, nb_gens;
	
	if (f_v) {
		cout << "action::init_matrix_group_strong_generators_builtin" << endl;
		}
	F = M->GFq;
	n = M->n;
	q = F->q;
	if (f_v) {
		cout << "action::init_matrix_group_strong_generators_builtin computing strong generators builtin group" << endl;
		cout << "n=" << n << endl;
		cout << "q=" << q << endl;
		cout << "p=" << F->p << endl;
		cout << "e=" << F->e << endl;
		cout << "f_semilinear=" << M->f_semilinear << endl;
		}

	if (M->f_projective) {
		strong_generators_for_projective_linear_group(n, F, 
			M->f_semilinear, 
			data, size, nb_gens, 
			verbose_level - 1);
			// in GALOIS/projective.C
		}
	else if (M->f_affine) {
		strong_generators_for_affine_linear_group(n, F, 
			M->f_semilinear, 
			data, size, nb_gens, 
			verbose_level - 1);
			// in GALOIS/projective.C
		}
	else if (M->f_general_linear) {
		strong_generators_for_general_linear_group(n, F, 
			M->f_semilinear, 
			data, size, nb_gens, 
			verbose_level - 1);
			// in GALOIS/projective.C
		}
	else {
		cout << "action::init_matrix_group_strong_generators_builtin unknown group type" << endl;
		exit(1);
		}



	Strong_gens = NEW_OBJECT(strong_generators);
	Strong_gens->init_from_data(this, data, nb_gens, size, transversal_length, verbose_level - 1);

	FREE_INT(data);
	
	if (f_v) {
		cout << "action::init_matrix_group_strong_generators_builtin computing strong generators builtin group finished" << endl;
		}
}

void action::init_permutation_group(INT degree, INT verbose_level)
{
	INT page_length_log = PAGE_LENGTH_LOG;
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	perm_group *P;
	
	if (f_v) {
		cout << "action::init_permutation_group, degree=" << degree << endl;
		}
	sprintf(group_prefix, "Perm%ld", degree);
	sprintf(label, "Perm%ld", degree);
	sprintf(label_tex, "Perm%ld", degree);
	P = NEW_OBJECT(perm_group);
	type_G = perm_group_t;
	G.perm_grp = P;
	f_allocated = TRUE;
	
	if (f_v) {
		cout << "action::init_permutation_group before P->init" << endl;
		}

	P->init(degree, page_length_log, verbose_level);

	if (f_v) {
		cout << "action::init_permutation_group after P->init" << endl;
		}
	
	init_function_pointers_permutation_group();

	
	elt_size_in_INT = P->elt_size_INT;
	coded_elt_size_in_char = P->char_per_elt;
	if (f_vv) {
		cout << "elt_size_in_INT = " << elt_size_in_INT << endl;
		cout << "coded_elt_size_in_char = " << coded_elt_size_in_char << endl;
		}
	allocate_element_data();
	action::degree = degree;
	make_element_size = degree;
	
	if (f_v) {
		cout << "action::init_permutation_group, finished" << endl;
		cout << "a permutation group of degree " << action::degree << endl;
		cout << "and of order ";
		print_group_order(cout);
		cout << endl;
		print_info();
		}
	
}

void action::init_permutation_group_from_generators(INT degree, 
	INT nb_gens, INT *gens, 
	INT given_base_length, INT *given_base,
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 1);
	INT f_vvv = (verbose_level >= 3);
	//INT f_vvvv = (verbose_level >= 4);
	sims *G;
	longinteger_object go, G_order, quo, rem;
	INT drop_out_level, image, base_point;
	INT nb_generators_found;
	INT *idx_generators_found;
	INT *generators_found;
	INT level, i, j, k, f_added;
	
	
	if (f_v) {
		cout << "action::init_permutation_group_from_generators degree=" << degree << " nb_gens=" << nb_gens << " given_base_length=" << given_base_length << endl;
		}
	sprintf(group_prefix, "Perm%ld", degree);
	sprintf(label, "Perm%ld", degree);
	sprintf(label_tex, "Perm%ld", degree);
	if (f_vv) {
		cout << "action::init_permutation_group_from_generators the " << nb_gens << " generators are" << endl;
		for (i = 0; i < nb_gens; i++) {
			cout << i << " : ";
			if (degree < 20) {
				perm_print(cout, gens + i * degree, degree);
				}
			else {
				cout << "too large to print";
				}
			cout << endl;
			}
		}
	
	if (f_vv) {
		cout << "action::init_permutation_group_from_generators calling init_permutation_group" << endl;
		}
	init_permutation_group(degree, verbose_level /* - 2 */);
	if (f_vv) {
		cout << "action::init_permutation_group_from_generators after init_permutation_group" << endl;
		}

	if (f_vv) {
		cout << "action::init_permutation_group_from_generators calling allocate_base_data" << endl;
		cout << "initial base:";
		INT_vec_print(cout, given_base, given_base_length);
		cout << " of length " << given_base_length << endl;
		}
	allocate_base_data(given_base_length);
	base_len = given_base_length;
	
	// init base:
	for (i = 0; i < base_len; i++) {
		base[i] = given_base[i];
		}



	if (f_vv) {
		cout << "action::init_permutation_group_from_generators, now trying to set up the group from the given generators" << endl;
		}
	
	G = NEW_OBJECT(sims);
	
	G->init(this);
	G->init_trivial_group(verbose_level - 1);
	G->group_order(G_order);
	
	idx_generators_found = NEW_INT(nb_gens);
	generators_found = NEW_INT(nb_gens * degree);
	
	for (level = given_base_length - 1; level >= 0; level--) {
		base_point = given_base[level];
		if (f_vv) {
			cout << "action::init_permutation_group_from_generators level " << level << " base point " << base_point << endl;
			}
		perm_group_find_strong_generators_at_level(level, degree, 
			given_base_length, given_base,
			nb_gens, gens, nb_generators_found, idx_generators_found);
		if (f_vv) {
			cout << "action::init_permutation_group_from_generators found " << nb_generators_found << " generators:";
			INT_vec_print(cout, idx_generators_found, nb_generators_found);
			cout << endl;
			}
		for (i = 0; i < nb_generators_found; i++) {
			k = idx_generators_found[i];
			for (j = 0; j < degree; j++) {
				generators_found[i * degree + j] = gens[k * degree + j];
				}
			}
		if (f_vvv) {
			cout << "action::init_permutation_group_from_generators they are" << endl;
			for (i = 0; i < nb_generators_found; i++) {
				cout << i << " : ";
				perm_print(cout, generators_found + i * degree, degree);
				cout << endl;
				}
			}
		{
		schreier O;
	
		//cout << "calling O.init()" << endl;
	
		O.init(this);

		//cout << "calling O.init_generators()" << endl;
		O.init_generators(nb_generators_found, generators_found);
		
		if (FALSE) {
			cout << "action::init_permutation_group_from_generators generators in schreier" << endl;
			O.print_generators();
			}

		if (f_vv) {
			cout << "action::init_permutation_group_from_generators computing orbit of point " << base_point << endl;
			}
		O.compute_point_orbit(base_point, 0);
		if (f_vv) {
			cout << "action::init_permutation_group_from_generators point " << base_point << " lies in an orbit of size " << O.orbit_len[0] << endl;
			if (FALSE) {
				O.print(cout);
				O.print_tables(cout, FALSE);
				}
			}
		for (j = 0; j < O.orbit_len[0]; j++) {
			if (FALSE) {
				cout << "action::init_permutation_group_from_generators level " << level << " coset rep " << j << endl;
				}
			O.coset_rep(j);
			if (FALSE) {
				perm_print(cout, O.cosetrep, degree);
				//A.element_print(O.cosetrep, cout);
				cout << endl;
				}
			if (G->strip(O.cosetrep, Elt2 /* residue */, drop_out_level, image, 0 /*verbose_level - 1*/)) {
				if (FALSE) {
					cout << "action::init_permutation_group_from_generators element strips through." << endl;
					cout << "if it is the identity element, that's OK," << endl;
					cout << "otherwise please add another base point," << endl;
					cout << "a point which is moved by the residue" << endl;
					element_print(Elt2, cout);
					}
				f_added = FALSE;
				}
			else {
				f_added = TRUE;
				G->add_generator_at_level(Elt2, drop_out_level, 0 /*verbose_level - 1*/);
				}

			G->group_order(G_order);


			if (f_v && f_added) {
				cout << "action::init_permutation_group_from_generators level " << level << " coset " << j << " group of order increased to " << G_order << endl;
				}
			if (f_vv) {
				cout << "action::init_permutation_group_from_generators level " << level << " coset " << j << " found a group of order " << G_order << endl;
				}
			}
		} // end of schreier
		
		} // next level
	if (f_v) {
		cout << "action::init_permutation_group_from_generators found a group of order " << G_order << endl;
		}

	if (f_v) {
		cout << "action::init_permutation_group_from_generators: found a group of order " << G_order << endl;
		}
		
	FREE_INT(idx_generators_found);
	FREE_INT(generators_found);

	
#if 0
	strong_generators = NEW_OBJECT(vector_ge);
	tl = NEW_INT(base_len);
	f_has_strong_generators = TRUE;
	
	if (f_v) {
		cout << "action::init_permutation_group_from_generators: extracting strong generators" << endl;
		}
	G->extract_strong_generators_in_order(*strong_generators, tl, 0/*verbose_level - 2*/);
	if (f_v) {
		cout << "action::init_permutation_group_from_generators: extracting strong generators done" << endl;
		}
	if (FALSE) {
		cout << "strong generators are:" << endl;
		strong_generators->print(cout);
		cout << endl;
		}
	if (f_vv) {
		cout << "action::init_permutation_group_from_generators transversal lengths:" << endl;
		for (i = 0; i < base_len; i++) {
			cout << tl[i] << " ";
			}
		cout << endl;
		}


	if (f_vvv) {
		//G.print(TRUE);
		//G.print_generator_depth_and_perm();
		}
#endif
	
	if (f_vv) {
		cout << "action::init_permutation_group_from_generators init_sims" << endl;
		}
	init_sims(G, 0/*verbose_level - 1*/);

	if (f_vv) {
		cout << "action::init_permutation_group_from_generators after init_sims" << endl;
		}



	if (f_vv) {
		cout << "action::init_permutation_group_from_generators compute_strong_generators_from_sims" << endl;
		}
	compute_strong_generators_from_sims( 0 /*verbose_level - 2*/);
	if (f_vv) {
		cout << "action::init_permutation_group_from_generators after_strong_generators_from_sims" << endl;
		}
	if (f_v) {
		print_info();
		}
	if (f_v) {
		cout << "action::init_permutation_group_from_generators done" << endl;
		}
}

void action::init_affine_group(INT n, INT q, INT f_translations, 
	INT f_semilinear, INT frobenius_power, 
	INT f_multiplication, INT multiplication_order, INT verbose_level)
{
	INT nb_gens, degree;
	INT *gens;
	INT given_base_length;
	INT *given_base;
	finite_field F;
	
	sprintf(label, "AGL_%ld_%ld", n, q);
	sprintf(label_tex, "AGL(%ld,%ld)", n, q);

	F.init(q, verbose_level - 1);
	
	
	F.affine_generators(n, f_translations, 
		f_semilinear, frobenius_power, 
		f_multiplication, multiplication_order, 
		nb_gens, degree, gens, 
		given_base_length, given_base);

	init_permutation_group_from_generators(degree, 
		nb_gens, gens, 
		given_base_length, given_base,
		verbose_level);

	FREE_INT(gens);
	FREE_INT(given_base);
}

void action::init_affine_grid_group(INT q1, INT q2, 
	INT f_translations1, INT f_translations2, 
	INT f_semilinear1, INT frobenius_power1, 
	INT f_semilinear2, INT frobenius_power2, 
	INT f_multiplication1, INT multiplication_order1, 
	INT f_multiplication2, INT multiplication_order2, 
	INT f_diagonal, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT n = 1;
	INT nb_gens1, degree1;
	INT nb_gens2, degree2;
	INT nb_gens3, degree3;
	INT *gens1;
	INT *gens2;
	INT *gens3;
	INT given_base_length1;
	INT given_base_length2;
	INT given_base_length3;
	INT *given_base1;
	INT *given_base2;
	INT *given_base3;
	finite_field F1, F2;
	
	sprintf(label, "AGL_1_%ld_x_AGL_1_%ld", q1, q2);
	sprintf(label_tex, "AGL(1,%ld) x AGL(1,%ld)", q1, q2);

	F1.init(q1, verbose_level - 1);
	F2.init(q2, verbose_level - 1);
	
	if (f_v) {
		cout << "computing generators first group" << endl;
		}
	F1.affine_generators(n, f_translations1, 
		f_semilinear1, frobenius_power1, 
		f_multiplication1, multiplication_order1, 
		nb_gens1, degree1, gens1, 
		given_base_length1, given_base1);
	if (f_v) {
		cout << "computing generators second group" << endl;
		}
	F2.affine_generators(n, f_translations2, 
		f_semilinear2, frobenius_power2, 
		f_multiplication2, multiplication_order2, 
		nb_gens2, degree2, gens2, 
		given_base_length2, given_base2);

	if (f_v) {
		cout << "computing generators for the direct product" << endl;
		}
	if (f_diagonal) {
		perm_group_generators_direct_product(1, degree1, degree2, degree3, 
			nb_gens1, nb_gens2, nb_gens3, 
			gens1, gens2, gens3, 
			given_base_length1, given_base_length2, given_base_length3, 
			given_base1, given_base2, given_base3);
		}
	else {
		perm_group_generators_direct_product(0, degree1, degree2, degree3, 
			nb_gens1, nb_gens2, nb_gens3, 
			gens1, gens2, gens3, 
			given_base_length1, given_base_length2, given_base_length3, 
			given_base1, given_base2, given_base3);
		}
	
	init_permutation_group_from_generators(degree3, 
		nb_gens3, gens3, 
		given_base_length3, given_base3,
		verbose_level);
	if (f_v) {
		cout << "elt_size_in_INT=" << elt_size_in_INT << endl;
		}

	FREE_INT(gens1);
	FREE_INT(gens2);
	FREE_INT(gens3);
	FREE_INT(given_base1);
	FREE_INT(given_base2);
	FREE_INT(given_base3);
}

void action::init_symmetric_group(INT degree, INT verbose_level)
{
	INT nb_gens, *gens;
	INT given_base_length, *given_base;
	INT i, j;
	
	sprintf(label, "Sym_%ld", degree);
	sprintf(label_tex, "Sym(%ld)", degree);

	make_element_size = degree;
	nb_gens = degree - 1;
	given_base_length = degree - 1;
	gens = NEW_INT(nb_gens * degree);
	given_base = NEW_INT(given_base_length);
	for (i = 0; i < nb_gens; i++) {
		for (j = 0; j < degree; j++) {
			gens[i * degree + j] = j;
			}
		gens[i * degree + i] = i + 1;
		gens[i * degree + i + 1] = i;
		}
	for (i = 0; i < given_base_length; i++) {
		given_base[i] = i;
		}
	init_permutation_group_from_generators(degree, 
		nb_gens, gens, 
		given_base_length, given_base,
		verbose_level);
	FREE_INT(gens);
	FREE_INT(given_base);
}

void action::null_function_pointers()
{
	ptr_element_image_of = NULL;
	ptr_element_image_of_low_level = NULL;
	ptr_element_linear_entry_ij = NULL;
	ptr_element_linear_entry_frobenius = NULL;
	ptr_element_one = NULL;
	ptr_element_is_one = NULL;
	ptr_element_unpack = NULL;
	ptr_element_pack = NULL;
	ptr_element_retrieve = NULL;
	ptr_element_store = NULL;
	ptr_element_mult = NULL;
	ptr_element_invert = NULL;
	ptr_element_move = NULL;
	ptr_element_dispose = NULL;
	ptr_element_print = NULL;
	ptr_element_print_quick = NULL;
	ptr_element_print_latex = NULL;
	ptr_element_print_verbose = NULL;
	ptr_element_print_for_make_element = NULL;
	ptr_print_point = NULL;
}

void action::init_function_pointers_matrix_group()
{
	ptr_element_image_of = matrix_group_element_image_of;
	ptr_element_image_of_low_level = matrix_group_element_image_of_low_level;
	ptr_element_linear_entry_ij = matrix_group_element_linear_entry_ij;
	ptr_element_linear_entry_frobenius = matrix_group_element_linear_entry_frobenius;
	ptr_element_one = matrix_group_element_one;
	ptr_element_is_one = matrix_group_element_is_one;
	ptr_element_unpack = matrix_group_element_unpack;
	ptr_element_pack = matrix_group_element_pack;
	ptr_element_retrieve = matrix_group_element_retrieve;
	ptr_element_store = matrix_group_element_store;
	ptr_element_mult = matrix_group_element_mult;
	ptr_element_invert = matrix_group_element_invert;
	ptr_element_move = matrix_group_element_move;
	ptr_element_dispose = matrix_group_element_dispose;
	ptr_element_print = matrix_group_element_print;
	ptr_element_print_quick = matrix_group_element_print_quick;
	ptr_element_print_latex = matrix_group_element_print_latex;
	ptr_element_print_verbose = matrix_group_element_print_verbose;
	ptr_element_print_for_make_element = matrix_group_element_print_for_make_element;
	ptr_print_point = matrix_group_print_point;
}

void action::init_function_pointers_permutation_group()
{
	ptr_element_image_of = perm_group_element_image_of;
	ptr_element_image_of_low_level = NULL;
	ptr_element_linear_entry_ij = NULL;
	ptr_element_linear_entry_frobenius = NULL;
	ptr_element_one = perm_group_element_one;
	ptr_element_is_one = perm_group_element_is_one;
	ptr_element_unpack = perm_group_element_unpack;
	ptr_element_pack = perm_group_element_pack;
	ptr_element_retrieve = perm_group_element_retrieve;
	ptr_element_store = perm_group_element_store;
	ptr_element_mult = perm_group_element_mult;
	ptr_element_invert = perm_group_element_invert;
	ptr_element_move = perm_group_element_move;
	ptr_element_dispose = perm_group_element_dispose;
	ptr_element_print = perm_group_element_print;
	ptr_element_print_quick = perm_group_element_print; // no quick version here!
	ptr_element_print_latex = perm_group_element_print_latex;
	ptr_element_print_verbose = perm_group_element_print_verbose;
	ptr_element_print_for_make_element = perm_group_element_print_for_make_element;
	ptr_print_point = perm_group_print_point;
}

void action::init_function_pointers_induced_action()
{
	//ptr_get_transversal_rep = induced_action_get_transversal_rep;
	ptr_element_image_of = induced_action_element_image_of;
	ptr_element_image_of_low_level = induced_action_element_image_of_low_level;
	ptr_element_linear_entry_ij = NULL;
	ptr_element_linear_entry_frobenius = NULL;
	ptr_element_one = induced_action_element_one;
	ptr_element_is_one = induced_action_element_is_one;
	ptr_element_unpack = induced_action_element_unpack;
	ptr_element_pack = induced_action_element_pack;
	ptr_element_retrieve = induced_action_element_retrieve;
	ptr_element_store = induced_action_element_store;
	ptr_element_mult = induced_action_element_mult;
	ptr_element_invert = induced_action_element_invert;
	ptr_element_move = induced_action_element_move;
	ptr_element_dispose = induced_action_element_dispose;
	ptr_element_print = induced_action_element_print;
	ptr_element_print_quick = induced_action_element_print_quick;
	ptr_element_print_latex = induced_action_element_print_latex;
	ptr_element_print_verbose = induced_action_element_print_verbose;
	ptr_element_print_for_make_element = induced_action_element_print_for_make_element;
	ptr_print_point = induced_action_print_point;
}

void action::create_sims(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	sims *S;
	
	if (f_v) {
		cout << "action::create_sims" << endl;
		}
	if (!f_has_strong_generators) {
		cout << "action::create_sims we need strong generators" << endl;
		exit(1);
		}

	S = Strong_gens->create_sims(verbose_level - 1);
	
	init_sims(S, verbose_level);
	if (f_v) {
		cout << "action::create_sims done" << endl;
		}
}



void action::create_orthogonal_group(action *subaction, 
	INT f_has_target_group_order, longinteger_object &target_go, 
	void (* callback_choose_random_generator)(INT iteration, 
		INT *Elt, void *data, INT verbose_level), 
	INT verbose_level)
{
	//verbose_level = 10;
	
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);

	if (f_v) {
		cout << "action::create_orthogonal_group we will use a schreier_sims object" << endl;
		}
	
	{
	schreier_sims ss;

	ss.init(this, verbose_level - 1);

	ss.interested_in_kernel(subaction, verbose_level - 1);
	
	if (f_has_target_group_order) {
		ss.init_target_group_order(target_go, verbose_level - 1);
		}
	
	ss.init_random_process(
		callback_choose_random_generator,   // see action_global.C
		&ss, 
		verbose_level - 1);
	
	ss.create_group(verbose_level);
	
	init_sims(ss.G, verbose_level);
	f_has_kernel = TRUE;
	Kernel = ss.K;

	ss.K = NULL;
	ss.G = NULL;
	
	//init_transversal_reps_from_stabilizer_chain(G, verbose_level - 2);
	if (f_v) {
		cout << "action::create_orthogonal_group after init_sims, calling compute_strong_generators_from_sims" << endl;
		}
	compute_strong_generators_from_sims(verbose_level - 2);
	if (f_v) {
		cout << "action::create_orthogonal_group done" << endl;
		}
	}
}


