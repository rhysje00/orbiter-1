// code_generator.C
//
// Anton Betten
//
// moved here from codes.C: May 18, 2009
//
// December 30, 2003

#include "orbiter.h"
#include "discreta.h"
#include "codes.h"

// ##################################################################################################
// start of class code_generator
// ##################################################################################################


void code_generator::read_arguments(int argc, const char **argv)
{
	INT i;
	INT f_n = FALSE;
	INT f_k = FALSE;
	INT f_q = FALSE;
	INT f_d = FALSE;
	
	
	if (argc <= 4) {
		print_usage();
		exit(1);
		}

	gen->read_arguments(argc, argv, 0);

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-lex") == 0) {
			f_lex = TRUE;
			cout << "-lex " << endl;
			}
		else if (strcmp(argv[i], "-debug") == 0) {
			f_debug = TRUE;
			cout << "-debug " << endl;
			}
		else if (strcmp(argv[i], "-schreier_depth") == 0) {
			schreier_depth = atoi(argv[++i]);
			cout << "-schreier_depth " << schreier_depth << endl;
			}
		else if (strcmp(argv[i], "-n") == 0) {
			f_n = TRUE;
			n = atoi(argv[++i]);
			cout << "-n " << n << endl;
			}
		else if (strcmp(argv[i], "-k") == 0) {
			f_k = TRUE;
			k = atoi(argv[++i]);
			cout << "-k " << k << endl;
			}
		else if (strcmp(argv[i], "-q") == 0) {
			f_q = TRUE;
			q = atoi(argv[++i]);
			cout << "-q " << q << endl;
			}
		else if (strcmp(argv[i], "-d") == 0) {
			f_d = TRUE;
			d = atoi(argv[++i]);
			cout << "-d " << d << endl;
			}
		else if (strcmp(argv[i], "-draw_poset") == 0) {
			f_draw_poset = TRUE;
			cout << "-draw_poset " << endl;
			}
		else if (strcmp(argv[i], "-print_data_structure") == 0) {
			f_print_data_structure = TRUE;
			cout << "-print_data_structure " << endl;
			}
		else if (strcmp(argv[i], "-list") == 0) {
			f_list = TRUE;
			cout << "-list" << endl;
			}
		else if (strcmp(argv[i], "-table_of_nodes") == 0) {
			f_table_of_nodes = TRUE;
			cout << "-table_of_nodes" << endl;
			}
		}
	
	if (!f_n) {
		cout << "Please use option -n <n> to specify n" << endl;
		exit(1);
		}
	if (!f_k) {
		cout << "Please use option -k <k> to specify k" << endl;
		exit(1);
		}
	if (!f_q) {
		cout << "Please use option -q <q> to specify q" << endl;
		exit(1);
		}
	if (!f_d) {
		cout << "Please use option -d <d> to specify d" << endl;
		exit(1);
		}
	cout << "n=" << n << endl;
	cout << "k=" << k << endl;
	cout << "q=" << q << endl;
	cout << "d=" << d << endl;
		
	f_irreducibility_test = TRUE;
	
	INT p, h;
	
	is_prime_power(q, p, h);
	if (h > 1) {
		f_semilinear = TRUE;
		}
	else {
		f_semilinear = FALSE;
		}

}

code_generator::code_generator()
{
	null();
}

code_generator::~code_generator()
{
	freeself();
}

void code_generator::null()
{
	gen = NULL;
	F = NULL;
	A = NULL;
	schreier_depth = 1000;
	f_list = FALSE;
	f_table_of_nodes = FALSE;
	f_use_invariant_subset_if_available = TRUE;
	f_debug = FALSE;
	f_lex = FALSE;
	f_draw_poset = FALSE;
	f_print_data_structure = FALSE;
}

void code_generator::freeself()
{
	if (A) {
		delete A;
		}
	if (F) {
		delete F;
		}
	if (gen) {
		delete gen;
		}
	null();
}

void code_generator::init(int argc, const char **argv)
{
	F = new finite_field;
	A = new action;
	gen = new generator;
	INT f_basis = TRUE;
	

	read_arguments(argc, argv);
	
	INT verbose_level = gen->verbose_level;
	
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "code_generator::init" << endl;
		}
	nmk = n - k;

	sprintf(gen->fname_base, "codes_%ld_%ld_%ld_%ld", n, k, q, d);
	
	if (f_v) {
		cout << "code_generator::init initializing finite field of order " << q << endl;
		}
	F->init(q, 0);
	if (f_v) {
		cout << "code_generator::init initializing finite field of order " << q << " done" << endl;
		}

	if (f_v) {
		cout << "code_generator::init calling init_projective_group, dimension = " << nmk << endl;
		}
	
	A->init_projective_group(nmk, F, 
		f_semilinear, 
		f_basis, 
		verbose_level - 2);
	
	if (f_v) {
		cout << "code_generator::init finished with init_projective_group" << endl;
		}
	if (f_v) {
		cout << "code_generator::init degree = " << A->degree << endl;
		}
	

	
	gen->depth = n;
	
	if (f_v) {
		cout << "code_generator::init group set up, calling gen->init" << endl;
		cout << "A->f_has_strong_generators=" << A->f_has_strong_generators << endl;
		}
	
	gen->init(A, A, A->Strong_gens, gen->depth /* sz */, verbose_level);

#if 0
	if (f_v) {
		cout << "code_generator::init group set up, calling gen->init_check_func" << endl;
		}

	gen->init_check_func(check_mindist, this /* candidate_check_data */);

	if (f_v) {
		cout << "code_generator::init group set up, calling gen->init_incremental_check_func" << endl;
		}

	gen->init_incremental_check_func(check_mindist_incremental, this /* candidate_check_data */);
#endif
	if (f_v) {
		cout << "code_generator::init group set up, calling gen->init_early_test_func" << endl;
		}
	gen->init_early_test_func(
		check_mindist_early_test_func, 
		this,  
		verbose_level);
	//gen->f_its_OK_to_not_have_an_early_test_func = TRUE;

	
	rc.init(F, nmk, n, d);

	if (FALSE && gen->A->degree < 1000) {
		cout << "the elements of PG(" << n - k - 1 << "," << F->q << ") are:" << endl;
		display_all_PG_elements(n - k - 1, *F);
		}

	gen->f_print_function = TRUE;
	gen->print_function = print_code;
	gen->print_function_data = this;
	

	INT nb_oracle_nodes = ONE_MILLION;
	
	if (f_v) {
		cout << "code_generator::init group set up, calling gen->init_oracle" << endl;
		}

	gen->init_oracle(nb_oracle_nodes, verbose_level - 1);

	if (f_v) {
		cout << "code_generator::init group set up, calling gen->root[0].init_root_node" << endl;
		}

	gen->root[0].init_root_node(gen, gen->verbose_level - 2);
}


void code_generator::print(INT len, INT *S)
{
	INT i, j;
	
	if (len == 0) {
		return;
		}
	cout << "generator matrix:" << endl;
	for (j = 0; j < len; j++) {
		PG_element_unrank_modified(*F, rc.M1 + j, len /* stride */, nmk /* len */, S[j]);
		}
	print_integer_matrix(cout, rc.M1, nmk, len);

	INT_matrix_print_tex(cout, rc.M1, nmk, len);

	INT *weights;
	INT *M = rc.M1;
	INT n = len;
	INT k = nmk;

	weights = NEW_INT(n + 1);
	F->code_projective_weight_enumerator(n, k, 
		M, // [k * n]
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

	
}

void code_generator::main()
{
	INT depth;
	INT f_embedded = TRUE;
	INT f_sideways = TRUE;
	INT verbose_level = 0;

	depth = gen->main(t0, 
		schreier_depth, 
		f_use_invariant_subset_if_available, 
		f_lex, 
		f_debug, 
		gen->verbose_level);

	if (f_table_of_nodes) {
		INT *Table;
		INT nb_rows, nb_cols;

		gen->get_table_of_nodes(Table, nb_rows, nb_cols, verbose_level);
	
		INT_matrix_write_csv("data.csv", Table, nb_rows, nb_cols);


		FREE_INT(Table);
		}

	if (f_list) {
		INT f_show_stab = TRUE, f_show_whole_orbit = FALSE;
		
		gen->list_all_orbits_at_level(depth, 
			TRUE, 
			print_code, 
			this, 
			f_show_stab, f_show_whole_orbit);

		INT d;
		for (d = 0; d < 3; d++) {
			gen->print_schreier_vectors_at_depth(d, verbose_level);
			}
		}

	if (f_draw_poset) {
		gen->draw_poset(gen->fname_base, depth, 0 /* data1 */, f_embedded, f_sideways, gen->verbose_level);
		}
	if (f_print_data_structure) {
		gen->print_data_structure_tex(depth, gen->verbose_level);
		}
}




// ##################################################################################################
// callback functions
// ##################################################################################################

void check_mindist_early_test_func(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level)
{
	code_generator *cg = (code_generator *) data;
	INT i, j, node, f, l, pt, nb_good_orbits;
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "check_mindist_early_test_func S=";
		INT_vec_print(cout, S, len);
		cout << " testing " << nb_candidates << " candidates" << endl;
		//INT_vec_print(cout, candidates, nb_candidates);
		//cout << endl;
		}
	node = cg->gen->find_oracle_node_for_set(len, S, FALSE /* f_tolerant */, 0);
	oracle *O;
	O = cg->gen->root + node;

	if (f_v) {
		cout << "check_mindist_early_test_func for ";
		O->print_set(cg->gen);
		cout << endl;
		}

	schreier Schreier;

	Schreier.init(cg->gen->A2);
	Schreier.init_generators_by_hdl(O->nb_strong_generators, O->hdl_strong_generators, 0);
	Schreier.orbits_on_invariant_subset_fast(nb_candidates, candidates, 0/*verbose_level*/);

	if (f_v) {
		cout << "after Schreier.compute_all_orbits_on_invariant_subset, we found " << Schreier.nb_orbits << " orbits" << endl;
		}
	nb_good_candidates = 0;
	nb_good_orbits = 0;
	for (i = 0; i < Schreier.nb_orbits; i++) {
		f = Schreier.orbit_first[i];
		l = Schreier.orbit_len[i];
		pt = Schreier.orbit[f];
		S[len] = pt;
		if (cg->rc.check_rank_last_two_are_fixed(len + 1, 
			S, verbose_level - 1)) {
			for (j = 0; j < l; j++) {
				pt = Schreier.orbit[f + j];
				good_candidates[nb_good_candidates++] = pt;
				}	
			nb_good_orbits++;
			}
		}

	INT_vec_heapsort(good_candidates, nb_good_candidates);
	if (f_v) {
		cout << "after Schreier.compute_all_orbits_on_invariant_subset, we found " << nb_good_candidates << " good candidates in " << nb_good_orbits << " good orbits" << endl;
		}
	
#if 0
	nb_good_candidates = 0;
	for (i = 0; i < nb_candidates; i++) {
		S[len] = candidates[i];
		if (cg->rc.check_rank_last_two_are_fixed(len + 1, 
			S, verbose_level - 1)) {		
			good_candidates[nb_good_candidates++] = candidates[i];
			}
		if ((i % 1000) == 0 && i) {
			INT t1, dt;
	
			t1 = os_ticks();
			dt = t1 - cg->gen->t0;
			cout << "Time ";
			time_check_delta(cout, dt);
			cout << " : " << i << endl;
			}
		}
#endif
	
	if (f_v) {
		cout << "check_mindist_early_test_func nb_good_candidates=" << nb_good_candidates << endl;
		}
}

INT check_mindist(INT len, INT *S, void *data, INT verbose_level)
{
	code_generator *cg = (code_generator *) data;
	INT f_OK = TRUE;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		}
	if (!cg->rc.check_rank(len, S, verbose_level - 1)) {
		return FALSE;
		}

	
	if (f_OK) {
		if (f_v) {
			cout << "OK" << endl;
			}
		return TRUE;
		}
	else {
		return FALSE;
		}
}

INT check_mindist_incremental(INT len, INT *S, void *data, INT verbose_level)
{
	code_generator *cg = (code_generator *) data;
	INT f_OK = TRUE;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "checking set ";
		print_set(cout, len, S);
		cout << " (incrementally)";
		}
	if (!cg->rc.check_rank_last_two_are_fixed(len, S, verbose_level - 1)) {
		return FALSE;
		}

	
	if (f_OK) {
		if (f_v) {
			cout << "OK" << endl;
			}
		return TRUE;
		}
	else {
		return FALSE;
		}
}

void print_code(INT len, INT *S, void *data)
{
	code_generator *cg = (code_generator *) data;
	
	cg->print(len, S);
}



