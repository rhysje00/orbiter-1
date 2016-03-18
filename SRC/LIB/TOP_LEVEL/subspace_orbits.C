// subspace_orbits.C
// 
// Anton Betten
//
// started:    January 25, 2010
// moved here: March 29, 2012
// 
//
//

#include "orbiter.h"


subspace_orbits::subspace_orbits()
{
	LG = NULL;
	tmp_M = NULL;
	tmp_M2 = NULL;
	tmp_M3 = NULL;
	base_cols = NULL;
	v = NULL;
	w = NULL;
	Gen = NULL;
	
	f_print_generators = FALSE;
	f_has_extra_test_func = FALSE;
}

subspace_orbits::~subspace_orbits()
{
	if (tmp_M) {
		FREE_INT(tmp_M);
		}
	if (tmp_M2) {
		FREE_INT(tmp_M2);
		}
	if (tmp_M3) {
		FREE_INT(tmp_M3);
		}
	if (base_cols) {
		FREE_INT(base_cols);
		}
	if (v) {
		FREE_INT(v);
		}
	if (w) {
		FREE_INT(w);
		}
	if (Gen) {
		delete Gen;
		}
}

void subspace_orbits::init(int argc, const char **argv, 
	linear_group *LG, INT depth, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "subspace_orbits::init" << endl;
		}

	subspace_orbits::LG = LG;
	subspace_orbits::depth = depth;
	n = LG->vector_space_dimension;
	F = LG->F;
	q = F->q;

	if (f_v) {
		cout << "subspace_orbits::init n=" << n << " q=" << q << endl;
		}


	
	tmp_M = NEW_INT(n * n);
	tmp_M2 = NEW_INT(n * n);
	tmp_M3 = NEW_INT(n * n);
	base_cols = NEW_INT(n);
	v = NEW_INT(n);
	w = NEW_INT(n);
	Gen = new generator;

	if (f_v) {
		cout << "subspace_orbits::init before Gen->read_arguments" << endl;
		}

	Gen->read_arguments(argc, argv, 0);

	if (f_v) {
		cout << "subspace_orbits::init after Gen->read_arguments" << endl;
		}

	if (f_v) {
		cout << "subspace_orbits::init LG->prefix=" << LG->prefix << endl;
		}

	sprintf(Gen->fname_base, "%s", LG->prefix);
	
	
	Gen->depth = depth;

	if (f_v) {
		cout << "subspace_orbits::init before init_group" << endl;
		}

	init_group(verbose_level);

	if (f_v) {
		cout << "subspace_orbits::init after init_group" << endl;
		}


	if (f_v) {
		cout << "subspace_orbits::init done" << endl;
		}
}


void subspace_orbits::init_group(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);



	if (f_print_generators) {
		INT f_print_as_permutation = FALSE;
		INT f_offset = TRUE;
		INT offset = 1;
		INT f_do_it_anyway_even_for_big_degree = TRUE;
		INT f_print_cycles_of_length_one = TRUE;
		
		cout << "subspace_orbits->init_group printing generators for the group:" << endl;
		LG->Strong_gens->gens->print(cout, f_print_as_permutation, 
			f_offset, offset, 
			f_do_it_anyway_even_for_big_degree, 
			f_print_cycles_of_length_one);
		}

	if (f_v) {
		cout << "subspace_orbits->init_group before Gen->init" << endl;
		}
	Gen->init(LG->A_linear, LG->A2, LG->Strong_gens, Gen->depth, verbose_level);

#if 0
	Gen->init_check_func(
		subspace_orbits_test_func, 
		this /* candidate_check_data */);
#endif
	Gen->init_early_test_func(
		subspace_orbits_early_test_func, 
		this /*void *data */,  
		verbose_level);

	//Gen->init_incremental_check_func(
		//check_mindist_incremental, 
		//this /* candidate_check_data */);

	Gen->init_vector_space_action(n, 
		F, 
		subspace_orbits_rank_point_func, 
		subspace_orbits_unrank_point_func, 
		this, 
		verbose_level);
#if 0
	Gen->f_print_function = TRUE;
	Gen->print_function = print_set;
	Gen->print_function_data = this;
#endif	

	INT nb_oracle_nodes = 1000;
	
	if (f_v) {
		cout << "subspace_orbits->init_group before Gen->init_oracle" << endl;
		}
	Gen->init_oracle(nb_oracle_nodes, verbose_level - 1);
	if (f_v) {
		cout << "subspace_orbits->init_group calling Gen->init_root_node" << endl;
		}
	Gen->root[0].init_root_node(Gen, verbose_level - 1);
	
	schreier_depth = Gen->depth;
	f_use_invariant_subset_if_available = FALSE;
	f_implicit_fusion = FALSE;
	f_debug = FALSE;
	if (f_v) {
		cout << "subspace_orbits->init_group done" << endl;
		}
}


#if 0
void subspace_orbits::read_data_file(INT depth_completed, const BYTE *fname_data_file, 
	INT f_exportmagma, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_recompute_schreier = TRUE;
	
	if (f_v) {
		cout << "subspace_orbits::read_data_file" << endl;
		}
	Gen->read_data_file(depth_completed, 
		fname_data_file, 
		verbose_level - 1);
	
	// ignore the last level: Schreier vectors have not yet been computed
	depth_completed--;
	
	
	if (f_v) {
		cout << "read_data_file: after reading file " << fname_data_file << endl;
		cout << "depth_completed = " << depth_completed << endl;
		}

	if (f_recompute_schreier) {
		if (f_v) {
			cout << "recomputing Schreier vectors" << endl;
			}
		Gen->recreate_schreier_vectors_up_to_level(
			depth_completed, 
			TRUE /* f_compact */, 
			MINIMUM(verbose_level, 1));
		}
	if (f_v) {
		cout << "read_data_file: recreated Schreier vectors" << endl;
		}


	if (f_exportmagma) {
		INT level;

		level = depth_completed + 1;

		wedge_product_export_magma(Gen, 
			n, q, vector_space_dimension, 
			level, verbose_level);
			// SNAKES_AND_LADDERS/snakes_and_ladders_global.C

		

		} // exportmagma
}
#endif


void subspace_orbits::compute_orbits(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT t0 = os_ticks();
	
	if (f_v) {
		cout << "subspace_orbits::compute_orbits calling generator_main" << endl;
		cout << "A=";
		Gen->A->print_info();
		cout << "A2=";
		Gen->A2->print_info();
		}
	Gen->main(t0, 
		schreier_depth, 
		f_use_invariant_subset_if_available, 
		f_implicit_fusion, 
		f_debug, 
		verbose_level - 1);
	
	INT nb_orbits;
	
	if (f_v) {
		cout << "subspace_orbits::compute_orbits done with generator_main" << endl;
		}
	nb_orbits = Gen->nb_orbits_at_level(depth);
	if (f_v) {
		cout << "subspace_orbits::compute_orbits we found " << nb_orbits << " orbits at depth " << depth << endl;
		}
}

void subspace_orbits::unrank_set_to_M(INT len, INT *S)
{
	INT i;
	
	for (i = 0; i < len; i++) {
		PG_element_unrank_modified(*F, tmp_M + i * n, 1, n, S[i]);
		}
}



void subspace_orbits::Kramer_Mesner_matrix(INT t, INT k, INT f_print_matrix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "subspace_orbits::Kramer_Mesner_matrix t=" << t << " k=" << k << ":" << endl;
		}


	matrix Mtk;

	compute_Kramer_Mesner_matrix(Gen, 
		t, k, Mtk, TRUE /* f_subspaces */, q, verbose_level - 2);
		// in DISCRETA/discreta_global.C

	if (f_v) {
		cout << "The Kramer Mesner matrix has size " << Mtk.s_m() << " x " << Mtk.s_n() << endl;
		//cout << Mtk << endl;
		}

	if (f_print_matrix) {
		cout << "The Kramer Mesner matrix has size " << Mtk.s_m() << " x " << Mtk.s_n() << endl;
		cout << Mtk << endl;
		}

	if (f_v) {
		cout << "creating diophant:" << endl;
		}

	diophant *D;

	matrix_to_diophant(Mtk, D, verbose_level);

	if (f_v) {
		cout << "diophant has been created" << endl;
		}

	BYTE fname[1000];

	sprintf(fname, "%s_KM_%ld_%ld.system", Gen->fname_base, t, k);
	cout << "saving diophant under the name " << fname << endl;
	D->save_in_general_format(fname, verbose_level);

#if 0

	cout << "eqn 6:" << endl;
		for (j = 0; j < nb_cols; j++) {
			if (Mtk.s_iji(6, j)) {
				cout << Mtk.s_iji(6, j) << " in col " << j << " : ";
				}
			}
	cout << endl;
#endif


	cout << "closing diophant:" << endl;
	delete D;
	
	
#if 0
	Mtk_sup_to_inf(Gen, t, k, Mtk, Mtk_inf, verbose_level - 2);
	cout << "M_{" << t << "," << k << "} inf has been computed" << endl;
	//Mtk_inf.print(cout);
	
	//cout << endl;
	//cout << endl;
	
	INT nb_t_orbits;
	INT nb_k_orbits;
	INT first_t, first_k;
	INT len, rep, size;
	INT set1[1000];
	//INT set2[1000];
	
	first_t = Gen->first_oracle_node_at_level[t];
	first_k = Gen->first_oracle_node_at_level[k];
	nb_t_orbits = Mtk_inf.s_m();
	nb_k_orbits = Mtk_inf.s_n();
	for (j = 0; j < nb_k_orbits; j++) {
		cout << "   ";
		}
	cout << "| ";
	cout << " t-orbit orbit_lenth" << endl;
	for (i = 0; i < nb_t_orbits; i++) {
		len = Gen->orbit_length_as_INT(i, t);
		//cout << "i=" << i << " len=" << len << endl;
		Gen->get_set(first_t + i, set1, size);
		if (size != t) {
			cout << "size != t" << endl;
			exit(1);
			}
		for (j = 0; j < nb_k_orbits; j++) {
			a = Mtk_inf.s_iji(i, j);
			cout << setw(2) << a << " ";
			}
		cout << "| ";
		cout << setw(3) << i << " " << setw(3) << len << " ";
		if (t == 1) {
			rep = set1[0];
			schreier Schreier;

			Schreier.init(Gen->A2);
			Schreier.init_generators_by_hdl(Gen->root[0].nb_strong_generators, Gen->root[0].hdl_strong_generators, verbose_level - 1);
			Schreier.compute_point_orbit(rep, 0 /* verbose_level */);
			if (Schreier.orbit_len[0] != len) {
				cout << "Schreier.orbit_len[0] != len" << endl;
				exit(1);
				}
			INT *pts;
			INT len1;

			pts = NEW_INT(len);
			Schreier.get_orbit(0 /* orbit_idx */, pts, len1, 0 /* verbose_level */);
			
			//cout << "{";
			INT_vec_print(cout, pts, len);
			//cout << "}";
			FREE_INT(pts);
			}
		cout << endl;
		}
	cout << "t-orbits, t=" << t << " :" << endl;
	cout << "i : orbit_length of i-th orbit" << endl;
	for (i = 0; i < nb_t_orbits; i++) {
		len = Gen->orbit_length_as_INT(i, t);
		cout << i << " : " << len << endl;
		}
	cout << "k-orbits, k=" << k << " :" << endl;
	cout << "i : orbit_length of i-th orbit" << endl;
	for (i = 0; i < nb_k_orbits; i++) {
		len = Gen->orbit_length_as_INT(i, k);
		cout << i << " : " << len << endl;
		}
#endif
}

INT subspace_orbits::test_dim_C_cap_Cperp_property(INT len, INT *S, INT d)
{
	//INT i;
	
#if 0
	cout << "subspace_orbits::test_dim_C_Cperp_property" << endl;
	cout << "Set ";
	INT_vec_print(cout, S, len);
	cout << endl;
#endif

	unrank_set_to_M(len, S);

#if 0
	cout << "coordinate matrix:" << endl;
	print_integer_matrix_width(cout, tmp_M, len, n, n, F->log10_of_q);
#endif

	INT k = len;
	INT k3;


	F->perp_standard_with_temporary_data(n, k, tmp_M, 
		tmp_M2, tmp_M3, base_cols, 
		0 /*verbose_level*/);

	//cout << "C perp:" << endl;
	//print_integer_matrix_width(cout, tmp_M + k * n, n - k, n, n, P->F->log10_of_q);


	F->intersect_subspaces(n, k, tmp_M, n - k, tmp_M + k * n, 
		k3, tmp_M2, 0 /*verbose_level*/);

	//cout << "\\dim (C \\cap C^\\bot) = " << k3 << endl;

	//cout << "basis for C \\cap C^\\bot:" << endl;
	//print_integer_matrix_width(cout, tmp_M2, k3, n, n, P->F->log10_of_q);

	if (k3 == d) {
		return TRUE;
		}
	else {
		return FALSE;
		}
}

INT subspace_orbits::compute_minimum_distance(INT len, INT *S)
{
	INT i, d = 0;
	
#if 0
	cout << "subspace_orbits::compute_minimum_distance" << endl;
	cout << "Set ";
	INT_vec_print(cout, S, len);
	cout << endl;
#endif

	unrank_set_to_M(len, S);

#if 0
	cout << "coordinate matrix:" << endl;
	print_integer_matrix_width(cout, tmp_M, len, n, n, P->F->log10_of_q);
#endif


	INT *weights;

	weights = NEW_INT(n + 1);
	F->code_projective_weight_enumerator(n, len, 
		tmp_M, // [k * n]
		weights, // will be allocated [N]
		0 /*verbose_level*/);


	//cout << "projective weights: " << endl;
	for (i = 1; i <= n; i++) {
		if (weights[i]) {
			d = i;
			break;
			}
		}
	FREE_INT(weights);
	return d;
}

void subspace_orbits::print_set(INT len, INT *S)
{
	INT i;
	
	cout << "subspace_orbits::print_set" << endl;
	cout << "Set ";
	INT_vec_print(cout, S, len);
	cout << endl;

	unrank_set_to_M(len, S);


	cout << "coordinate matrix:" << endl;
	print_integer_matrix_width(cout, tmp_M, len, n, n, F->log10_of_q);


	INT *weights;

	weights = NEW_INT(n + 1);
	F->code_projective_weight_enumerator(n, len, 
		tmp_M, // [k * n]
		weights, // will be allocated [N]
		0 /*verbose_level*/);


	cout << "projective weights: " << endl;
	for (i = 0; i <= n; i++) {
		if (weights[i] == 0) {
			continue;
			}
		cout << i << " : " << weights[i] << endl;
		}
	FREE_INT(weights);

	INT k = len;
	INT k3;


	F->perp_standard_with_temporary_data(n, k, tmp_M, 
		tmp_M2, tmp_M3, base_cols, 
		0 /*verbose_level*/);

	cout << "C perp:" << endl;
	print_integer_matrix_width(cout, tmp_M + k * n, n - k, n, n, F->log10_of_q);


	F->intersect_subspaces(n, k, tmp_M, n - k, tmp_M + k * n, 
		k3, tmp_M2, 0 /*verbose_level*/);

	cout << "\\dim (C \\cap C^\\bot) = " << k3 << endl;

	cout << "basis for C \\cap C^\\bot:" << endl;
	print_integer_matrix_width(cout, tmp_M2, k3, n, n, F->log10_of_q);
	
	
}

INT subspace_orbits::test_set(INT len, INT *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT ret = TRUE;
	INT rk;
	
	if (f_v) {
		cout << "subspace_orbits::test_set" << endl;
		cout << "Testing set ";
		INT_vec_print(cout, S, len);
		cout << endl;
		}
	unrank_set_to_M(len, S);
	if (f_vv) {
		cout << "coordinate matrix:" << endl;
		print_integer_matrix_width(cout, tmp_M, len, n, n, F->log10_of_q);
		}
	rk = F->Gauss_simple(tmp_M, len, n, base_cols, 0 /*verbose_level - 2*/);
	if (f_v) {
		cout << "the matrix has rank " << rk << endl;
		}
	if (rk < len) {
		ret = FALSE;
		}
	if (ret) {
		if (f_has_extra_test_func) {
			ret = (*extra_test_func)(this, len, S, extra_test_func_data, verbose_level);
			}
		}

	if (ret) {
		if (f_v) {
			cout << "OK" << endl;
			}
		}
	else {
		if (f_v) {
			cout << "not OK" << endl;
			}
		}
	return ret;
}

INT subspace_orbits::test_minimum_distance(INT len, INT *S, INT mindist, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT ret = TRUE;
	INT i, h, wt, N, k;
	INT *msg;
	INT *word;
	INT *M;
	
	if (f_v) {
		cout << "subspace_orbits::test_minimum_distance" << endl;
		cout << "Testing set ";
		INT_vec_print(cout, S, len);
		cout << endl;
		}
	k = len;
	M = tmp_M;
	unrank_set_to_M(len, S);
	if (f_vv) {
		cout << "coordinate matrix:" << endl;
		print_integer_matrix_width(cout, M, len, n, n, F->log10_of_q);
		}
	N = nb_PG_elements(k - 1, q);
	msg = v;
	word = w;
	for (h = 0; h < N; h++) {
		PG_element_unrank_modified(*F, msg, 1, k, h);
		//AG_element_unrank(q, msg, 1, k, h);
		F->mult_vector_from_the_left(msg, M, word, k, n);
		wt = 0;
		for (i = 0; i < n; i++) {
			if (word[i]) {
				wt++;
				}
			}
		if (wt < mindist) {
			ret = FALSE;
			break;
			}
		}
	if (f_v) {
		if (ret) {
			cout << "is OK" << endl;
			}
		else {
			cout << "is not OK" << endl;
			}
		}
	return ret;
}

INT subspace_orbits::test_if_self_orthogonal(INT len, INT *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT ret;
	INT i, j, a;
	INT *M;
	
	if (f_v) {
		cout << "subspace_orbits::test_if_self_orthogonal" << endl;
		cout << "Testing set ";
		INT_vec_print(cout, S, len);
		cout << endl;
		}
	M = tmp_M;

	unrank_set_to_M(len, S);
	if (f_vv) {
		cout << "coordinate matrix:" << endl;
		print_integer_matrix_width(cout, M, len, n, n, F->log10_of_q);
		}

	ret = TRUE;
	for (i = 0; i < len; i++) {
		for (j = i; j < len; j++) {
			a = F->dot_product(n, M + i * n, M + j * n);
			if (a) {
				ret = FALSE;
				break;
				}
			}
		if (j < len) {
			break;
			}
		}
	if (f_v) {
		if (ret) {
			cout << "is OK" << endl;
			}
		else {
			cout << "is not OK" << endl;
			}
		}
	return ret;
}




// ####################################################################################
// global functions:
// ####################################################################################


INT subspace_orbits_rank_point_func(INT *v, void *data)
{
	subspace_orbits *G;
	generator *gen;
	INT rk;
	
	G = (subspace_orbits *) data;
	gen = G->Gen;
	PG_element_rank_modified(*gen->F, v, 1, gen->vector_space_dimension, rk);
	return rk;
}

void subspace_orbits_unrank_point_func(INT *v, INT rk, void *data)
{
	subspace_orbits *G;
	generator *gen;
	
	G = (subspace_orbits *) data;
	gen = G->Gen;
	PG_element_unrank_modified(*gen->F, v, 1, gen->vector_space_dimension, rk);
}

void subspace_orbits_early_test_func(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	//verbose_level = 1;

	subspace_orbits *SubOrb;
	INT f_v = (verbose_level >= 1);
	INT i;

	SubOrb = (subspace_orbits *) data;

	if (f_v) {
		cout << "subspace_orbits_early_test_func" << endl;
		cout << "testing " << nb_candidates << " candidates" << endl;
		}
	nb_good_candidates = 0;
	for (i = 0; i < nb_candidates; i++) {
		S[len] = candidates[i];
		if (SubOrb->test_set(len + 1, S, verbose_level - 1)) {
			good_candidates[nb_good_candidates++] = candidates[i];
			}
		}
	if (f_v) {
		cout << "subspace_orbits_early_test_func" << endl;
		cout << "Out of " << nb_candidates << " candidates, " << nb_good_candidates << " survive" << endl;
		}
}




