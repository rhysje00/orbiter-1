// projective_space.C
//
// Anton Betten
// Jan 17, 2010

#include "galois.h"

#define MAX_NUMBER_OF_LINES 100000



INT projective_space::cntr_new = 0;
INT projective_space::cntr_objects = 0;
INT projective_space::f_debug_memory = FALSE;

void *projective_space::operator new(size_t bytes)
{
	cntr_new++;
	cntr_objects++;
	if (f_debug_memory) {
		cout << "projective_space::operator new bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	return malloc(bytes);
}

void *projective_space::operator new[](size_t bytes)
{
	INT n;
	
	n = bytes / sizeof(projective_space);
	cntr_new++;
	cntr_objects += n;
	if (f_debug_memory) {
		cout << "projective_space::operator new[] n=" << n 
			<< " bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	return malloc(bytes);
}

void projective_space::operator delete(void *ptr, size_t bytes)
{
	if (f_debug_memory) {
		cout << "projective_space::operator delete bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects--;
	return free(ptr);
}

void projective_space::operator delete[](void *ptr, size_t bytes)
{
	INT n;
	
	n = bytes / sizeof(projective_space);
	if (f_debug_memory) {
		cout << "projective_space::operator delete[] n=" << n 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects -= n;
	return free(ptr);
}

projective_space::projective_space()
{
	null();
};

projective_space::~projective_space()
{
	freeself();
}

void projective_space::null()
{
#if 0
	A = NULL;
	A2 = NULL;
	A_lines = NULL;
#endif
	//SG = NULL;
	Grass_lines = NULL;
	F = NULL;
	Go = NULL;
	//A_gens = NULL; // not allocated, only a reference!
	//A2_gens = NULL; // not allocated, only a reference!
	incidence_bitvec = NULL;
	Line_through_two_points = NULL;
	Line_intersection = NULL;
	Lines = NULL;
	Lines_on_point = NULL;
	Polarity_point_to_hyperplane = NULL;
	Polarity_hyperplane_to_point = NULL;
	v = NULL;
	w = NULL;
}

void projective_space::freeself()
{
	INT f_v = FALSE;

	if (f_v) {
		cout << "projective_space::freeself" << endl;
		}
#if 0
	if (A) {
		if (f_v) {
			cout << "projective_space::freeself deleting A" << endl;
			}
		FREE_OBJECT(A);
		}
	if (A2) {
		if (f_v) {
			cout << "projective_space::freeself deleting A2" << endl;
			}
		FREE_OBJECT(A2);
		}
#endif
	if (Go) {
		if (f_v) {
			cout << "projective_space::freeself deleting Go" << endl;
			}
		FREE_OBJECT(Go);
		}
	if (v) {
		if (f_v) {
			cout << "projective_space::freeself deleting v" << endl;
			}
		FREE_INT(v);
		}
	if (w) {
		FREE_INT(w);
		}

#if 0
	// A_lines should not be freed since it is part of A2
	if (A_lines) {
		if (f_v) {
			cout << "projective_space::freeself deleting A_lines" << endl;
			}
		delete A_lines;
		}
#endif
	if (Grass_lines) {
		if (f_v) {
			cout << "projective_space::freeself deleting Grass_lines" << endl;
			}
		FREE_OBJECT(Grass_lines);
		}
	if (incidence_bitvec) {
		FREE_UBYTE(incidence_bitvec);
		}
	if (Line_through_two_points) {
		FREE_INT(Line_through_two_points);
		}
	if (Line_intersection) {
		FREE_INT(Line_intersection);
		}
	if (Lines) {
		FREE_INT(Lines);
		}
	if (Lines_on_point) {
		if (f_v) {
			cout << "projective_space::freeself deleting Lines_on_point" << endl;
			}
		FREE_INT(Lines_on_point);
		}
	if (Polarity_point_to_hyperplane) {
		FREE_INT(Polarity_point_to_hyperplane);
		}
	if (Polarity_hyperplane_to_point) {
		FREE_INT(Polarity_hyperplane_to_point);
		}
	null();
	if (f_v) {
		cout << "projective_space::freeself done" << endl;
		}
}

void projective_space::init(INT n, finite_field *F, 
	INT f_init_incidence_structure, 
	INT verbose_level)
// n is projective dimension
{
	INT f_v = (verbose_level >= 1);

	projective_space::n = n;
	projective_space::F = F;
	projective_space::q = F->q;
	
	if (f_v) {
		cout << "projective_space::init PG(" << n << "," << q << ")" << endl;
		cout << "f_init_incidence_structure=" << f_init_incidence_structure << endl;
		}

	v = NEW_INT(n + 1);
	w = NEW_INT(n + 1);

	Grass_lines = NEW_OBJECT(grassmann);
	Grass_lines->init(n + 1, 2, F, verbose_level - 2);


	N_points = generalized_binomial(n + 1, 1, q);
	if (f_v) {
		cout << "projective_space::init N_points=" << N_points << endl;
		}
	N_lines = generalized_binomial(n + 1, 2, q);
	if (f_v) {
		cout << "projective_space::init N_lines=" << N_lines << endl;
		}
	r = generalized_binomial(n, 1, q);
	if (f_v) {
		cout << "projective_space::init r=" << r << endl;
		}
	k = q + 1;
	if (f_v) {
		cout << "projective_space::init k=" << k << endl;
		}



	if (f_init_incidence_structure) {
		if (f_v) {
			cout << "projective_space::init calling init_incidence_structure" << endl;
			}
		init_incidence_structure(verbose_level);
		if (f_v) {
			cout << "projective_space::init init_incidence_structure done" << endl;
			}
		}
	else {
		if (f_v) {
			cout << "projective_space::init we don't initialize the incidence structure data" << endl;
			}
		}
	if (f_v) {
		
		cout << "projective_space::init n=" << n << " q=" << q << " done" << endl;
		}
}

#if 0
void projective_space::init_group(INT f_with_line_action, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		cout << "projective_space::init_group" << endl;
		cout << "f_with_line_action=" << f_with_line_action << endl;
		}


	if (f_v) {
		if (F->e > 1) {
			cout << "projective_space::init_group the minimum polynomial is:" << endl;
			F->print_minimum_polynomial(F->p, F->polynomial);
			cout << endl;
			}
		}
	if (f_vv) {
		cout << "projective_space::init_group the finite field tables are:" << endl;
		F->print_tables();
		F->print_add_mult_tables();
		}


	A = NEW_OBJECT(action);
	
	Go = NEW_OBJECT(longinteger_object);
	if (f_v) {
		cout << "projective_space::init_group before init_projective_group" << endl;
		}

	A->init_projective_group(n + 1, 
		F, 
		f_semilinear, 
		f_basis, 
		verbose_level - 5);

	if (f_v) {
		cout << "projective_space::init_group after init_projective_group" << endl;
		}


	
	A->print_base();
	A->group_order(*Go);
	if (f_v) {
		cout << "projective_space::init_group group order = " << *Go << endl;
		}
	
	if (!A->f_has_strong_generators) {
		cout << "projective_space::init_group induced action does not have strong generators" << endl;
		}

	if (f_with_line_action) {
		init_line_action(verbose_level);
		}

	if (f_v) {
		cout << "projective_space::init_group computing strong generators done" << endl;
		}
}

void projective_space::init_line_action(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "projective_space::init_line_action" << endl;
		}
	A2 = NEW_OBJECT(action);

	A_lines = NEW_OBJECT(action_on_grassmannian);

	A_lines->init(*A, Grass_lines, verbose_level - 5);
	
	
	if (f_v) {
		cout << "projective_space::init_line_action action on grassmannian established" << endl;
		}

	if (f_v) {
		cout << "projective_space::init_line_action initializing A2" << endl;
		}
	INT f_induce_action = TRUE;
	sims S;
	longinteger_object go1;

	S.init(A);
	S.init_generators(*A->Strong_gens->gens, 0/*verbose_level*/);
	S.compute_base_orbits_known_length(A->transversal_length, 0/*verbose_level - 1*/);
	S.group_order(go1);
	if (f_v) {
		cout << "projective_space::init_line_action group order " << go1 << endl;
		}
	
	if (f_v) {
		cout << "projective_space::init_line_action initializing action on grassmannian" << endl;
		}
	A2->induced_action_on_grassmannian(A, A_lines, 
		f_induce_action, &S, verbose_level);
	if (f_v) {
		cout << "projective_space::init_line_action initializing A2 done" << endl;
		A2->print_info();
		}

	if (f_v) {
		cout << "projective_space::init_line_action computing strong generators" << endl;
		}
	if (!A2->f_has_strong_generators) {
		cout << "projective_space::init_line_action induced action does not have strong generators" << endl;
		}
	if (f_v) {
		cout << "projective_space::init_line_action done" << endl;
		}
}

void projective_space::init_without_group(INT n, INT f_init_incidence_structure, INT verbose_level)
// n is projective dimension
{
	//INT f_v = (verbose_level >= 1);
	
	init(n, F, 
		FALSE /* f_init_group*/, FALSE /*f_init_line_action*/,
		f_init_incidence_structure, 
		FALSE, FALSE, 
		verbose_level);
}
#endif

void projective_space::init_incidence_structure(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = FALSE; //(verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);

	INT i, j, a, b, i1, i2, j1, j2;
	
	if (f_v) {
		cout << "projective_space::init_incidence_structure" << endl;
		}


	if (N_lines < MAX_NUMBER_OF_LINES) {
		if (f_v) {
			cout << "projective_space::init_incidence_structure allocating Incidence (bitvector)" << endl;
			}
		INT len = (N_points * N_lines + 7) >> 3;
		if (f_v) {
			cout << "projective_space::init_incidence_structure allocating Incidence (bitvector) len = " << len << endl;
			}
		//Incidence = NEW_INT(N_points * N_lines);
		//INT_vec_zero(Incidence, N_points * N_points);
		incidence_bitvec = NEW_UBYTE(len);
		for (i = 0; i < len; i++) {
			incidence_bitvec[i] = 0;
			}

		}
	else {
		cout << "projective_space::init_incidence_structure: N_lines too big, we do not initialize the incidence matrix" << endl;
		//return;
		incidence_bitvec = NULL;
		}




	if (f_v) {
		cout << "projective_space::init_incidence_structure allocating Lines" << endl;
		}
	Lines = NEW_INT(N_lines * k);

	if (f_v) {
		cout << "projective_space::init_incidence_structure allocating Lines_on_point" << endl;
		cout << "projective_space::init_incidence_structure allocating r=" << r << endl;
		}
	Lines_on_point = NEW_INT(N_points * r);



	if (N_points * N_points < ONE_MILLION) {

		if (f_v) {
			cout << "projective_space::init_incidence_structure allocating Line_through_two_points" << endl;
			}
		Line_through_two_points = NEW_INT(N_points * N_points);
		}
	else {
		Line_through_two_points = NULL;
		}

	if (N_lines * N_lines < ONE_MILLION) {
		if (f_v) {
			cout << "projective_space::init_incidence_structure allocating Line_intersection" << endl;
			}
		Line_intersection = NEW_INT(N_lines * N_lines);
		INT_vec_zero(Line_through_two_points, N_points * N_points);
		for (i = 0; i < N_lines * N_lines; i++) {
			Line_intersection[i] = -1;
			}
		}
	else {
		Line_intersection = NULL;
		}	

	
	if (f_v) {
		cout << "number of points = " << N_points << endl;
		}
	if (f_vv) {
		for (i = 0; i < N_points; i++) {
			PG_element_unrank_modified(*F, v, 1, n + 1, i);
			cout << "point " << i << " : ";
			INT_vec_print(cout, v, n + 1);
			cout << " = ";
			F->INT_vec_print(cout, v, n + 1);

			PG_element_normalize_from_front(*F, v, 1, n + 1);
			cout << " = ";
			INT_vec_print(cout, v, n + 1);

		
			cout << " = ";
			F->INT_vec_print(cout, v, n + 1);

			
			cout << endl;
			}
		}
	if (f_v) {
		cout << "number of lines = " << N_lines << endl;
		}
	if (f_v) {
		cout << "computing lines..." << endl;
		}

	INT *R;

	R = NEW_INT(N_points);
	INT_vec_zero(R, N_points);
	
	for (i = 0; i < N_lines; i++) {
		if ((i % 1000000) == 0) {
			cout << "Line " << i << " / " << N_lines << ":" << endl;
			}
		Grass_lines->unrank_INT(i, 0/*verbose_level - 4*/);
		if (FALSE) {
			print_integer_matrix_width(cout, Grass_lines->M, 2, n + 1, n + 1, F->log10_of_q + 1);
			}
		j = Grass_lines->rank_INT(0/*verbose_level - 4*/);
		if (j != i) {
			cout << "rank yields " << j << " != " << i << endl;
			exit(1);
			}
		for (a = 0; a < k; a++) {
			PG_element_unrank_modified(*F, v, 1, 2, a);
			F->mult_matrix(v, Grass_lines->M, w, 1, 2, n + 1);
			PG_element_rank_modified(*F, w, 1, n + 1, b);
			if (incidence_bitvec) {
				incidence_m_ii(b, i, 1);
				}
			Lines[i * k + a] = b;
			Lines_on_point[b * r + R[b]] = i;
			R[b]++;
			}
		if (f_vv) {
			cout << "line " << i << ":" << endl;
			print_integer_matrix_width(cout, Grass_lines->M, 2, n + 1, n + 1, F->log10_of_q + 1);
			cout << "points on line " << i << " : ";
			INT_vec_print(cout, Lines + i * k, k);
			cout << endl;
			}
		
		}
	for (i = 0; i < N_points; i++) {
		if (R[i] != r) {
			cout << "R[i] != r" << endl;
			exit(1);
			}
		}

	FREE_INT(R);

	if (n == 2) {
		if (f_v) {
			cout << "computing polarity information..." << endl;
			}
		Polarity_point_to_hyperplane = NEW_INT(N_points);
		Polarity_hyperplane_to_point = NEW_INT(N_points);
		for (i = 0; i < N_lines; i++) {
			INT *A, a;
			A = NEW_INT((n + 1) * (n + 1));
			Grass_lines->unrank_INT(i, 0 /*verbose_level - 4*/);
			for (j = 0; j < 2 * (n + 1); j++) {
				A[j] = Grass_lines->M[j];
				}
			if (FALSE) {
				print_integer_matrix_width(cout, A, 2, n + 1, n + 1, F->log10_of_q + 1);
				}
			F->perp_standard(n + 1, 2, A, 0);
			if (FALSE) {
				print_integer_matrix_width(cout, A, n + 1, n + 1, n + 1, F->log10_of_q + 1);
				}
			PG_element_rank_modified(*F, A + 2 * (n + 1), 1, n + 1, a);
			if (f_vv) {
				cout << "line " << i << " is ";
				INT_vec_print(cout, A + 2 * (n + 1), n + 1);
				cout << "^\\perp = " << a << "^\\perp" << endl;
				}
			FREE_INT(A);
			Polarity_point_to_hyperplane[a] = i;
			Polarity_hyperplane_to_point[i] = a;
			}
		if (FALSE /* f_vv */) {
			cout << "i : pt_to_hyperplane[i] : hyperplane_to_pt[i]" << endl;
			for (i = 0; i < N_lines; i++) {
				cout << setw(4) << i << " " 
					<< setw(4) << Polarity_point_to_hyperplane[i] << " " 
					<< setw(4) << Polarity_hyperplane_to_point[i] << endl;
				}
			}
		}
	


#if 0
	if (f_v) {
		cout << "computing Lines_on_point..." << endl;
		}
	for (i = 0; i < N_points; i++) {
		if ((i % 1000) == 0) {
			cout << "point " << i << " / " << N_points << ":" << endl;
			}
		a = 0;
		for (j = 0; j < N_lines; j++) {	
			if (is_incident(i, j)) {
				Lines_on_point[i * r + a] = j;
				a++;
				}
			}
		if (FALSE /*f_vv */) {
			cout << "lines on point " << i << " = ";
			PG_element_unrank_modified(*F, w, 1, n + 1, i);
			INT_vec_print(cout, w, n + 1);
			cout << " : ";
			INT_vec_print(cout, Lines_on_point + i * r, r);
			cout << endl;
			}
		}
	if (f_v) {
		cout << "computing Lines_on_point done" << endl;
		}
#endif

	if (FALSE) {
		//cout << "incidence matrix:" << endl;
		//print_integer_matrix_width(cout, Incidence, N_points, N_lines, N_lines, 1);
		cout << "Lines:" << endl;
		print_integer_matrix_width(cout, Lines, N_lines, k, k, 3);
		cout << "Lines_on_point:" << endl;
		print_integer_matrix_width(cout, Lines_on_point, N_points, r, r, 3);
		}
	

	if (Line_through_two_points) {
		if (f_v) {
			cout << "computing Line_through_two_points..." << endl;
			}
		for (i1 = 0; i1 < N_points; i1++) {
			for (a = 0; a < r; a++) {
				j = Lines_on_point[i1 * r + a];
				for (b = 0; b < k; b++) {
					i2 = Lines[j * k + b];
					if (i2 == i1)
						continue;
					Line_through_two_points[i1 * N_points + i2] = j;
					Line_through_two_points[i2 * N_points + i1] = j;
					}
				}
			}
		if (FALSE) {
			cout << "line through points i,j is" << endl;
			for (i = 0; i < N_points; i++) {
				for (j = i + 1; j < N_points; j++) {
					cout << i << " , " << j << " : " << Line_through_two_points[i * N_points + j] << endl;
					}
				}
			//cout << "Line_through_two_points:" << endl;
			//print_integer_matrix_width(cout, Line_through_two_points, N_points, N_points, N_points, 2);
			}
		}

	if (Line_intersection) {
		if (f_v) {
			cout << "computing Line_intersection..." << endl;
			}
		for (j1 = 0; j1 < N_lines; j1++) {
			for (a = 0; a < k; a++) {
				i = Lines[j1 * k + a];
				for (b = 0; b < r; b++) {
					j2 = Lines_on_point[i * r + b];
					if (j2 == j1)
						continue;
					Line_intersection[j1 * N_lines + j2] = i;
					Line_intersection[j2 * N_lines + j1] = i;
					}
				}
			}
		if (FALSE) {
			cout << "point of intersection of lines i,j is" << endl;
			for (i = 0; i < N_lines; i++) {
				for (j = i + 1; j < N_lines; j++) {
					cout << i << " , " << j << " : " << Line_intersection[i * N_lines + j] << endl;
					}
				}
			//cout << "Line_intersection:" << endl;
			//print_integer_matrix_width(cout, Line_intersection, N_lines, N_lines, N_lines, 2);
			}
		}
	if (f_v) {
		cout << "projective_space::init_incidence_structure done" << endl;
		}
}

void projective_space::make_incidence_matrix(INT &m, INT &n, INT *&Inc, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, h, j;

	if (f_v) {
		cout << "projective_space::make_incidence_matrix" << endl;
		}
	m = N_points;
	n = N_lines;
	Inc = NEW_INT(m * n);
	INT_vec_zero(Inc, m * n);
	for (i = 0; i < N_points; i++) {
		for (h = 0; h < r; h++) {
			j = Lines_on_point[i * r + h];
			Inc[i * n + j] = 1;
			}
		}
	if (f_v) {
		cout << "projective_space::make_incidence_matrix done" << endl;
		}
}

INT projective_space::is_incident(INT pt, INT line)
{
	INT a;

	if (incidence_bitvec == 0) {
		cout << "projective_space::is_incident incidence_bitvec == 0" << endl;
		exit(1);
		}
	a = pt * N_lines + line;
	return bitvector_s_i(incidence_bitvec, a);
}

void projective_space::incidence_m_ii(INT pt, INT line, INT a)
{
	INT b;

	if (incidence_bitvec == 0) {
		cout << "projective_space::incidence_m_ii incidence_bitvec == 0" << endl;
		exit(1);
		}
	b = pt * N_lines + line;
	bitvector_m_ii(incidence_bitvec, b, a);
}

void projective_space::make_incidence_structure_and_partition(incidence_structure *&Inc, 
	partitionstack *&Stack, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *M;
	INT i, j, h;

	if (f_v) {
		cout << "projective_space::make_incidence_structure_and_partition" << endl;
		cout << "N_points=" << N_points << endl;
		cout << "N_lines=" << N_lines << endl;
		}
	Inc = new incidence_structure;

	
	if (f_v) {
		cout << "projective_space::make_incidence_structure_and_partition allocating M of size " << N_points * N_lines << endl;
		}
	M = NEW_INT(N_points * N_lines);
	if (f_v) {
		cout << "projective_space::make_incidence_structure_and_partition after allocating M of size " << N_points * N_lines << endl;
		}
	INT_vec_zero(M, N_points * N_lines);

	if (Lines_on_point == NULL) {
		cout << "projective_space::make_incidence_structure_and_partition Lines_on_point == NULL" << endl;
		exit(1);
		}
	for (i = 0; i < N_points; i++) {
		for (h = 0; h < r; h++) {
			j = Lines_on_point[i * r + h];
			M[i * N_lines + j] = 1;
			}
		}
	if (f_v) {
		cout << "projective_space::make_incidence_structure_and_partition before Inc->init_by_matrix" << endl;
		}
	Inc->init_by_matrix(N_points, N_lines, M, verbose_level - 1);
	if (f_v) {
		cout << "projective_space::make_incidence_structure_and_partition after Inc->init_by_matrix" << endl;
		}
	FREE_INT(M);


	Stack = new partitionstack;
	Stack->allocate(N_points + N_lines, 0 /* verbose_level */);
	Stack->subset_continguous(N_points, N_lines);
	Stack->split_cell(0 /* verbose_level */);
	Stack->sort_cells();
	
	if (f_v) {
		cout << "projective_space::make_incidence_structure_and_partition done" << endl;
		}
}

INT projective_space::nb_rk_k_subspaces_as_INT(INT k)
{
	longinteger_domain D;
	longinteger_object aa;
	INT N;
	INT d = n + 1;

	D.q_binomial(aa, d, k, q, 0/*verbose_level*/);
	N = aa.as_INT();
	return N;
}

void projective_space::print_all_points()
{
	INT *v;
	INT i;

	v = NEW_INT(n + 1);
	cout << "All points in PG(" << n << "," << q << "):" << endl;
	for (i = 0; i < N_points; i++) {
		unrank_point(v, i);
		cout << setw(3) << i << " : ";
		INT_vec_print(cout, v, n + 1);
		cout << endl;
		}
}

INT projective_space::rank_point(INT *v)
{
	INT b;
	
	PG_element_rank_modified(*F, v, 1, n + 1, b);
	return b;
}

void projective_space::unrank_point(INT *v, INT rk)
{	
	PG_element_unrank_modified(*F, v, 1, n + 1, rk);
}

INT projective_space::rank_line(INT *basis)
{
	INT b, i;
	
	for (i = 0; i < 2 * (n + 1); i++) {
		Grass_lines->M[i] = basis[i];
		}
	b = Grass_lines->rank_INT(0/*verbose_level - 4*/);
	return b;
}

void projective_space::unrank_line(INT *basis, INT rk)
{	
	INT i;
	
	Grass_lines->unrank_INT(rk, 0/*verbose_level - 4*/);
	for (i = 0; i < 2 * (n + 1); i++) {
		basis[i] = Grass_lines->M[i];
		}
}

INT projective_space::line_through_two_points(INT p1, INT p2)
{
	INT b;
	
	unrank_point(Grass_lines->M, p1);
	unrank_point(Grass_lines->M + n + 1, p2);
	b = Grass_lines->rank_INT(0/*verbose_level - 4*/);
	return b;
}

INT projective_space::test_if_lines_are_disjoint(INT l1, INT l2)
{
	if (Lines) {
		return test_if_sets_are_disjoint(Lines + l1 * k, Lines + l2 * k, k, k);
		}
	else {
		return test_if_lines_are_disjoint_from_scratch(l1, l2);
		}
}

INT projective_space::test_if_lines_are_disjoint_from_scratch(INT l1, INT l2)
{
	INT *Mtx;
	INT m, rk;

	m = n + 1;
	Mtx = NEW_INT(4 * m);
	Grass_lines->unrank_INT_here(Mtx, l1, 0/*verbose_level - 4*/);
	Grass_lines->unrank_INT_here(Mtx + 2 * m, l2, 0/*verbose_level - 4*/);
	rk = F->Gauss_easy(Mtx, 4, m);
	FREE_INT(Mtx);
	if (rk == 4) {
		return TRUE;
		}
	else {
		return FALSE;
		}
}

INT projective_space::line_intersection(INT l1, INT l2)
// works only for projective planes, i.e., n = 2
{
	INT *Mtx1;
	INT *Mtx3;
	INT i, b;
	
	if (n != 2) {
		cout << "projective_space::line_intersection n != 2" << endl;
		exit(1);
		}
	Mtx1 = NEW_INT(3 * 3);
	Mtx3 = NEW_INT(3 * 3);
	
	Grass_lines->unrank_INT(l1, 0/*verbose_level - 4*/);
	for (i = 0; i < 2 * 3; i++) {
		Mtx1[i] = Grass_lines->M[i];
		}
	F->perp_standard(3, 2, Mtx1, 0);
	for (i = 0; i < 3; i++) {
		Mtx3[i] = Mtx1[6 + i];
		}
	
	Grass_lines->unrank_INT(l2, 0/*verbose_level - 4*/);
	for (i = 0; i < 2 * 3; i++) {
		Mtx1[i] = Grass_lines->M[i];
		}
	F->perp_standard(3, 2, Mtx1, 0);
	for (i = 0; i < 3; i++) {
		Mtx3[3 + i] = Mtx1[6 + i];
		}
	F->perp_standard(3, 2, Mtx3, 0);
	b = rank_point(Mtx3 + 6);

	FREE_INT(Mtx1);
	FREE_INT(Mtx3);
	
	return b;
}

INT projective_space::arc_test(INT *input_pts, INT nb_pts, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *Pts;
	INT *Mtx;
	INT set[3];
	INT ret = TRUE;
	INT h, i, N;

	if (f_v) {
		cout << "projective_space::arc_test" << endl;
		}
	if (n != 2) {
		cout << "projective_space::arc_test n != 2" << endl;
		exit(1);
		}
	Pts = NEW_INT(nb_pts * 3);
	Mtx = NEW_INT(3 * 3);
	for (i = 0; i < nb_pts; i++) {
		unrank_point(Pts + i * 3, input_pts[i]);
		}
	N = INT_n_choose_k(nb_pts, 3);
	for (h = 0; h < N; h++) {
		unrank_k_subset(h, set, nb_pts, 3);
		INT_vec_copy(Pts + set[0] * 3, Mtx, 3);
		INT_vec_copy(Pts + set[1] * 3, Mtx + 3, 3);
		INT_vec_copy(Pts + set[2] * 3, Mtx + 6, 3);
		if (F->rank_of_matrix(Mtx, 3, 0 /* verbose_level */) < 3) {
			if (f_v) {
				cout << "Points P_" << set[0] << ", P_" << set[1] << " and P_" << set[2] << " are collinear" << endl;
				}
			ret = FALSE;
			}
		}

	FREE_INT(Pts);
	FREE_INT(Mtx);
	if (f_v) {
		cout << "projective_space::arc_test done" << endl;
		}
	return ret;
}

INT projective_space::determine_conic_in_plane(INT *input_pts, INT nb_pts, INT *six_coeffs, INT verbose_level)
// returns FALSE is the rank of the coefficient matrix is not 5. TRUE otherwise.
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *coords; // [nb_pts * 3];
	INT *system; // [nb_pts * 6];
	INT kernel[6 * 6];
	INT base_cols[6];
	INT i, x, y, z, rk;
	INT kernel_m, kernel_n;

	if (f_v) {
		cout << "projective_space::determine_conic_in_plane" << endl;
		}
	if (n != 2) {
		cout << "projective_space::determine_conic_in_plane n != 2" << endl;
		exit(1);
		}
	if (nb_pts < 5) {
		cout << "projective_space::determine_conic_in_plane need at least 5 points" << endl;
		exit(1);
		}

	if (!arc_test(input_pts, nb_pts, verbose_level)) {
		cout << "projective_space::determine_conic_in_plane some 3 of the points are collinear" << endl;
		exit(1);
		}


	coords = NEW_INT(nb_pts * 3);
	system = NEW_INT(nb_pts * 6);
	for (i = 0; i < nb_pts; i++) {
		unrank_point(coords + i * 3, input_pts[i]);
		}
	if (f_vv) {
		cout << "projective_space::determine_conic_in_plane points:" << endl;
		print_integer_matrix_width(cout, coords, nb_pts, 3, 3, F->log10_of_q);
		}
	for (i = 0; i < nb_pts; i++) {
		x = coords[i * 3 + 0];
		y = coords[i * 3 + 1];
		z = coords[i * 3 + 2];
		system[i * 6 + 0] = F->mult(x, x);
		system[i * 6 + 1] = F->mult(y, y);
		system[i * 6 + 2] = F->mult(z, z);
		system[i * 6 + 3] = F->mult(x, y);
		system[i * 6 + 4] = F->mult(x, z);
		system[i * 6 + 5] = F->mult(y, z);
		}
	if (f_v) {
		cout << "projective_space::determine_conic_in_plane system:" << endl;
		print_integer_matrix_width(cout, system, nb_pts, 6, 6, F->log10_of_q);
		}



	rk = F->Gauss_simple(system, nb_pts, 6, base_cols, verbose_level - 2);
	if (rk != 5) {
		if (f_v) {
			cout << "projective_space::determine_conic_in_plane system underdetermined" << endl;
			}
		return FALSE;
		}
	F->matrix_get_kernel(system, 5, 6, base_cols, rk, 
		kernel_m, kernel_n, kernel);
	if (f_v) {
		cout << "projective_space::determine_conic_in_plane conic:" << endl;
		print_integer_matrix_width(cout, kernel, 1, 6, 6, F->log10_of_q);
		}
	for (i = 0; i < 6; i++) {
		six_coeffs[i] = kernel[i];
		}
	FREE_INT(coords);
	FREE_INT(system);
	return TRUE;
}

void projective_space::determine_quadric_in_solid(INT *nine_pts_or_more, INT nb_pts, INT *ten_coeffs, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *coords; // [nb_pts * 4]
	INT *system; // [nb_pts * 10]
	INT kernel[10 * 10];
	INT base_cols[10];
	INT i, x, y, z, w, rk;
	INT kernel_m, kernel_n;

	if (f_v) {
		cout << "projective_space::determine_quadric_in_solid" << endl;
		}
	if (n != 3) {
		cout << "projective_space::determine_quadric_in_solid n != 3" << endl;
		exit(1);
		}
	if (nb_pts < 9) {
		cout << "projective_space::determine_quadric_in_solid you need to give at least 9 points" << endl;
		exit(1);
		}
	coords = NEW_INT(nb_pts * 4);
	system = NEW_INT(nb_pts * 10);
	for (i = 0; i < nb_pts; i++) {
		unrank_point(coords + i * 4, nine_pts_or_more[i]);
		}
	if (f_vv) {
		cout << "projective_space::determine_quadric_in_solid points:" << endl;
		print_integer_matrix_width(cout, coords, nb_pts, 4, 4, F->log10_of_q);
		}
	for (i = 0; i < nb_pts; i++) {
		x = coords[i * 4 + 0];
		y = coords[i * 4 + 1];
		z = coords[i * 4 + 2];
		w = coords[i * 4 + 3];
		system[i * 10 + 0] = F->mult(x, x);
		system[i * 10 + 1] = F->mult(y, y);
		system[i * 10 + 2] = F->mult(z, z);
		system[i * 10 + 3] = F->mult(w, w);
		system[i * 10 + 4] = F->mult(x, y);
		system[i * 10 + 5] = F->mult(x, z);
		system[i * 10 + 6] = F->mult(x, w);
		system[i * 10 + 7] = F->mult(y, z);
		system[i * 10 + 8] = F->mult(y, w);
		system[i * 10 + 9] = F->mult(z, w);
		}
	if (f_v) {
		cout << "projective_space::determine_quadric_in_solid system:" << endl;
		print_integer_matrix_width(cout, system, nb_pts, 10, 10, F->log10_of_q);
		}



	rk = F->Gauss_simple(system, nb_pts, 10, base_cols, verbose_level - 2);
	if (rk != 9) {
		cout << "projective_space::determine_quadric_in_solid system underdetermined" << endl;
		cout << "rk=" << rk << endl;
		exit(1);
		}
	F->matrix_get_kernel(system, 9, 10, base_cols, rk, 
		kernel_m, kernel_n, kernel);
	if (f_v) {
		cout << "projective_space::determine_quadric_in_solid conic:" << endl;
		print_integer_matrix_width(cout, kernel, 1, 10, 10, F->log10_of_q);
		}
	for (i = 0; i < 10; i++) {
		ten_coeffs[i] = kernel[i];
		}
}

void projective_space::conic_points_brute_force(INT *six_coeffs, 
	INT *points, INT &nb_points, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT v[3];
	INT i, a;

	if (f_v) {
		cout << "projective_space::conic_points_brute_force" << endl;
		}
	nb_points = 0;
	for (i = 0; i < N_points; i++) {
		unrank_point(v, i);
		a = F->evaluate_conic_form(six_coeffs, v);
		if (f_vv) {
			cout << "point " << i << " = ";
			INT_vec_print(cout, v, 3);
			cout << " gives a value of " << a << endl;
			}
		if (a == 0) {
			if (f_vv) {
				cout << "point " << i << " = ";
				INT_vec_print(cout, v, 3);
				cout << " lies on the conic" << endl;
				}
			points[nb_points++] = i;
			}
		}
	if (f_v) {
		cout << "projective_space::conic_points_brute_force done, we found " << nb_points << " points" << endl;
		}
	if (f_vv) {
		cout << "They are : ";
		INT_vec_print(cout, points, nb_points);
		cout << endl;
		}
}

void projective_space::quadric_points_brute_force(INT *ten_coeffs, 
	INT *points, INT &nb_points, INT verbose_level)
// quadric in PG(3,q)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT v[3];
	INT i, a;

	if (f_v) {
		cout << "projective_space::quadric_points_brute_force" << endl;
		}
	nb_points = 0;
	for (i = 0; i < N_points; i++) {
		unrank_point(v, i);
		a = F->evaluate_quadric_form_in_PG_three(ten_coeffs, v);
		if (f_vv) {
			cout << "point " << i << " = ";
			INT_vec_print(cout, v, 3);
			cout << " gives a value of " << a << endl;
			}
		if (a == 0) {
			if (f_vv) {
				cout << "point " << i << " = ";
				INT_vec_print(cout, v, 4);
				cout << " lies on the quadric" << endl;
				}
			points[nb_points++] = i;
			}
		}
	if (f_v) {
		cout << "projective_space::quadric_points_brute_force done, we found " << nb_points << " points" << endl;
		}
	if (f_vv) {
		cout << "They are : ";
		INT_vec_print(cout, points, nb_points);
		cout << endl;
		}
}

void projective_space::conic_points(INT *five_pts, INT *six_coeffs, 
	INT *points, INT &nb_points, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT Gram_matrix[9];
	INT Basis[9];
	INT Basis2[9];
	INT v[3], w[3];
	INT i, j, l, a, av, ma, b, bv, t;

	
	if (f_v) {
		cout << "projective_space::conic_points" << endl;
		}
	if (n != 2) {
		cout << "projective_space::conic_points n != 2" << endl;
		exit(1);
		}
	Gram_matrix[0 * 3 + 0] = F->add(six_coeffs[0], six_coeffs[0]);
	Gram_matrix[1 * 3 + 1] = F->add(six_coeffs[1], six_coeffs[1]);
	Gram_matrix[2 * 3 + 2] = F->add(six_coeffs[2], six_coeffs[2]);
	Gram_matrix[0 * 3 + 1] = Gram_matrix[1 * 3 + 0] = six_coeffs[3];
	Gram_matrix[0 * 3 + 2] = Gram_matrix[2 * 3 + 0] = six_coeffs[4];
	Gram_matrix[1 * 3 + 2] = Gram_matrix[2 * 3 + 1] = six_coeffs[5];
	if (f_vv) {
		cout << "projective_space::conic_points Gram matrix:" << endl;
		print_integer_matrix_width(cout, Gram_matrix, 3, 3, 3, F->log10_of_q);
		}
	
	unrank_point(Basis, five_pts[0]);
	for (i = 1; i < 5; i++) {
		unrank_point(Basis + 3, five_pts[i]);
		a = F->evaluate_bilinear_form(3, Basis, Basis + 3, Gram_matrix);
		if (a) {
			break;
			}
		}
	if (i == 5) {
		cout << "projective_space::conic_points did not find non-orthogonal vector" << endl;
		exit(1);
		}
	if (a != 1) {
		av = F->inverse(a);
		for (i = 0; i < 3; i++) {
			Basis[3 + i] = F->mult(av, Basis[3 + i]);
			}
		}
	if (f_v) {	
		cout << "projective_space::conic_points Hyperbolic pair:" << endl;
		print_integer_matrix_width(cout, Basis, 2, 3, 3, F->log10_of_q);
		}
	F->perp(3, 2, Basis, Gram_matrix);
	if (f_v) {	
		cout << "projective_space::conic_points perp:" << endl;
		print_integer_matrix_width(cout, Basis, 3, 3, 3, F->log10_of_q);
		}
	a = F->evaluate_conic_form(six_coeffs, Basis + 6);
	if (f_v) {	
		cout << "projective_space::conic_points form value = " << a << endl;
		}
	if (a == 0) {
		cout << "projective_space::conic_points the form is degenerate or we are in characteristic zero" << endl;
		exit(1);
		}
	l = F->log_alpha(a);
	if ((l % 2) == 0) {
		j = l / 2;
		b = F->alpha_power(j);
		bv = F->inverse(b);
		for (i = 0; i < 3; i++) {
			Basis[6 + i] = F->mult(bv, Basis[6 + i]);
			}
		a = F->evaluate_conic_form(six_coeffs, Basis + 6);
		if (f_v) {	
			cout << "form value = " << a << endl;
			}
		}
	for (i = 0; i < 3; i++) {
		Basis2[3 + i] = Basis[6 + i];
		}
	for (i = 0; i < 3; i++) {
		Basis2[0 + i] = Basis[0 + i];
		}
	for (i = 0; i < 3; i++) {
		Basis2[6 + i] = Basis[3 + i];
		}
	if (f_v) {	
		cout << "Basis2:" << endl;
		print_integer_matrix_width(cout, Basis2, 3, 3, 3, F->log10_of_q);
		}
	// Now the form is a^{-1}y_1^2 = y_0y_2 
	// (or, equivalently, a^{-1}y_1^2 - y_0y_2 = 0)
	// and  the quadratic form on (0,1,0) in y-coordinates is a.
	// 
	// In the y-coordinates, the points on this conic are 
	// (1,0,0) and (t^2,t,-a) for t \in GF(q).
	// In the x-coordinates, the points are 
	// (1,0,0) * Basis2 and (t^2,t,-a) * Basis2 for t \in GF(q).

	v[0] = 1;
	v[1] = 0;
	v[2] = 0;

	F->mult_vector_from_the_left(v, Basis2, w, 3, 3);
	if (f_v) {	
		cout << "vector corresponding to 100:" << endl;
		INT_vec_print(cout, w, 3);
		cout << endl;
		}
	b = rank_point(w);
	points[0] = b;
	nb_points = 1;

	ma = F->negate(a);
	
	for (t = 0; t < F->q; t++) {
		v[0] = F->mult(t, t);
		v[1] = t;
		v[2] = ma;
		F->mult_vector_from_the_left(v, Basis2, w, 3, 3);
		if (f_v) {	
			cout << "vector corresponding to t=" << t << ":" << endl;
			INT_vec_print(cout, w, 3);
			cout << endl;
			}
		b = rank_point(w);
		points[nb_points++] = b;
		}
	if (f_vv) {	
		cout << "projective_space::conic_points conic points:" << endl;
		INT_vec_print(cout, points, nb_points);
		cout << endl;
		}
	
}

void projective_space::find_tangent_lines_to_conic(INT *six_coeffs, 
	INT *points, INT nb_points, 
	INT *tangents, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	//INT v[3];
	INT Basis[9];
	INT Gram_matrix[9];
	INT i;
	
	if (f_v) {
		cout << "projective_space::find_tangent_lines_to_conic" << endl;
		}
	if (n != 2) {
		cout << "projective_space::find_tangent_lines_to_conic n != 2" << endl;
		exit(1);
		}
	Gram_matrix[0 * 3 + 0] = F->add(six_coeffs[0], six_coeffs[0]);
	Gram_matrix[1 * 3 + 1] = F->add(six_coeffs[1], six_coeffs[1]);
	Gram_matrix[2 * 3 + 2] = F->add(six_coeffs[2], six_coeffs[2]);
	Gram_matrix[0 * 3 + 1] = Gram_matrix[1 * 3 + 0] = six_coeffs[3];
	Gram_matrix[0 * 3 + 2] = Gram_matrix[2 * 3 + 0] = six_coeffs[4];
	Gram_matrix[1 * 3 + 2] = Gram_matrix[2 * 3 + 1] = six_coeffs[5];
	
	for (i = 0; i < nb_points; i++) {
		unrank_point(Basis, points[i]);
		F->perp(3, 1, Basis, Gram_matrix);
		if (f_vv) {	
			cout << "perp:" << endl;
			print_integer_matrix_width(cout, Basis, 3, 3, 3, F->log10_of_q);
			}
		tangents[i] = rank_line(Basis + 3);
		if (f_vv) {	
			cout << "tangent at point " << i << " is " << tangents[i] << endl;
			}
		}
}

void projective_space::PG_2_8_create_conic_plus_nucleus_arc_1(INT *the_arc, INT &size, INT verbose_level)
{
	INT frame_data[] = {1,0,0, 0,1,0,  0,0,1,  1,1,1 };
	INT frame[4];
	INT i, j, b, h, idx;
	INT L[3];
	INT v[3];

	if (n != 2) {
		cout << "projective_space::PG_2_8_create_conic_plus_nucleus_arc_1 n != 2" << endl;
		exit(1);
		}
	if (q != 8) {
		cout << "projective_space::PG_2_8_create_conic_plus_nucleus_arc_1 q != 8" << endl;
		exit(1);
		}
	for (i = 0; i < 4; i++) {
		frame[i] = rank_point(frame_data + i * 3);
		}

	cout << "frame: ";
	INT_vec_print(cout, frame, 4);
	cout << endl;
	
	L[0] = Line_through_two_points[frame[0] * N_points + frame[1]];
	L[1] = Line_through_two_points[frame[1] * N_points + frame[2]];
	L[2] = Line_through_two_points[frame[2] * N_points + frame[0]];
	
	cout << "l1=" << L[0] << " l2=" << L[1] << " l3=" << L[2] << endl;

	size = 0;	
	for (h = 0; h < 3; h++) {
		for (i = 0; i < r; i++) {
			b = Lines[L[h] * r + i];
			if (INT_vec_search(the_arc, size, b, idx)) {
				continue;
				}
			for (j = size; j > idx; j--) {
				the_arc[j] = the_arc[j - 1];
				}
			the_arc[idx] = b;
			size++;
			}
		}
	cout << "there are " << size << " points on the three lines: ";
	INT_vec_print(cout, the_arc, size);
	cout << endl;


	for (i = 1; i < q; i++) {
		v[0] = 1;
		v[1] = i;
		v[2] = F->mult(i, i);
		b = rank_point(v);
		if (INT_vec_search(the_arc, size, b, idx)) {
			continue;
			}
		for (j = size; j > idx; j--) {
			the_arc[j] = the_arc[j - 1];
			}
		the_arc[idx] = b;
		size++;
		
		}

	cout << "projective_space::PG_2_8_create_conic_plus_nucleus_arc_1: after adding the rest of the conic, there are " << size << " points on the arc: ";
	INT_vec_print(cout, the_arc, size);
	cout << endl;
}

void projective_space::PG_2_8_create_conic_plus_nucleus_arc_2(INT *the_arc, INT &size, INT verbose_level)
{
	INT frame_data[] = {1,0,0, 0,1,0,  0,0,1,  1,1,1 };
	INT frame[4];
	INT i, j, b, h, idx;
	INT L[3];
	INT v[3];

	if (n != 2) {
		cout << "projective_space::PG_2_8_create_conic_plus_nucleus_arc_2 n != 2" << endl;
		exit(1);
		}
	if (q != 8) {
		cout << "projective_space::PG_2_8_create_conic_plus_nucleus_arc_2 q != 8" << endl;
		exit(1);
		}
	for (i = 0; i < 4; i++) {
		frame[i] = rank_point(frame_data + i * 3);
		}

	cout << "frame: ";
	INT_vec_print(cout, frame, 4);
	cout << endl;
	
	L[0] = Line_through_two_points[frame[0] * N_points + frame[2]];
	L[1] = Line_through_two_points[frame[2] * N_points + frame[3]];
	L[2] = Line_through_two_points[frame[3] * N_points + frame[0]];
	
	cout << "l1=" << L[0] << " l2=" << L[1] << " l3=" << L[2] << endl;

	size = 0;	
	for (h = 0; h < 3; h++) {
		for (i = 0; i < r; i++) {
			b = Lines[L[h] * r + i];
			if (INT_vec_search(the_arc, size, b, idx)) {
				continue;
				}
			for (j = size; j > idx; j--) {
				the_arc[j] = the_arc[j - 1];
				}
			the_arc[idx] = b;
			size++;
			}
		}
	cout << "there are " << size << " points on the three lines: ";
	INT_vec_print(cout, the_arc, size);
	cout << endl;


	for (i = 0; i < q; i++) {
		if (i == 1) {
			v[0] = 0;
			v[1] = 1;
			v[2] = 0;
			}
		else {
			v[0] = 1;
			v[1] = i;
			v[2] = F->mult(i, i);
			}
		b = rank_point(v);
		if (INT_vec_search(the_arc, size, b, idx)) {
			continue;
			}
		for (j = size; j > idx; j--) {
			the_arc[j] = the_arc[j - 1];
			}
		the_arc[idx] = b;
		size++;
		
		}

	cout << "projective_space::PG_2_8_create_conic_plus_nucleus_arc_2: after adding the rest of the conic, there are " << size << " points on the arc: ";
	INT_vec_print(cout, the_arc, size);
	cout << endl;
}

void projective_space::create_Maruta_Hamada_arc(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT data[] = { 
		0,1,2, 0,1,3, 0,1,4, 0,1,5,
		1,0,7, 1,0,8, 1,0,9, 1,0,10, 
		1,2,0, 1,3,0, 1,4,0, 1,5,0,
		1,7,5, 1,8,4, 1,9,3, 1,10,2, 
		1,1,1, 1,1,10, 1,10,1,  1,4,4,
		1,12,0, 1,0,12 
		 };
	INT points[22];
	INT i, j, b, h, idx;
	INT L[4];
	INT v[3];

	if (n != 2) {
		cout << "projective_space::create_Maruta_Hamada_arc n != 2" << endl;
		exit(1);
		}
	if (q != 13) {
		cout << "projective_space::create_Maruta_Hamada_arc q != 13" << endl;
		exit(1);
		}
	for (i = 0; i < 22; i++) {
		points[i] = rank_point(data + i * 3);
		cout << "point " << i << " has rank " << points[i] << endl;
		}

	if (f_v) {
		cout << "projective_space::create_Maruta_Hamada_arc() points: ";
		INT_vec_print(cout, points, 22);
		cout << endl;
		}
	
	L[0] = Line_through_two_points[1 * N_points + 2];
	L[1] = Line_through_two_points[0 * N_points + 2];
	L[2] = Line_through_two_points[0 * N_points + 1];
	L[3] = Line_through_two_points[points[20] * N_points + points[21]];
	
	if (f_v) {
		cout << "L:";
		INT_vec_print(cout, L, 4);
		cout << endl;
		}

	for (h = 0; h < 4; h++) {
		cout << "h=" << h << " : L[h]=" << L[h] << " : " << endl;
		for (i = 0; i < r; i++) {
			b = Lines[L[h] * r + i];
			cout << "point " << b << " = ";
			unrank_point(v, b);
			PG_element_normalize_from_front(*F, v, 1, 3);
			INT_vec_print(cout, v, 3);
			cout << endl;
			}
		cout << endl;
		}
	size = 0;	
	for (h = 0; h < 4; h++) {
		for (i = 0; i < r; i++) {
			b = Lines[L[h] * r + i];
			if (INT_vec_search(the_arc, size, b, idx)) {
				continue;
				}
			for (j = size; j > idx; j--) {
				the_arc[j] = the_arc[j - 1];
				}
			the_arc[idx] = b;
			size++;
			}
		}
	if (f_v) {
		cout << "there are " << size << " points on the quadrilateral: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}


	// remove the first 16 points:
	for (i = 0; i < 16; i++) {
		cout << "removing point " << i << " : " << points[i] << endl;
		if (!INT_vec_search(the_arc, size, points[i], idx)) {
			cout << "error, cannot find point to be removed" << endl;
			exit(1);
			}
		for (j = idx; j < size; j++) {
			the_arc[j] = the_arc[j + 1];
			}
		size--;
		}

	// add points 16-19:
	for (i = 16; i < 20; i++) {
		if (INT_vec_search(the_arc, size, points[i], idx)) {
			cout << "error, special point already there" << endl;
			exit(1);
			}
		for (j = size; j > idx; j--) {
			the_arc[j] = the_arc[j - 1];
			}
		the_arc[idx] = points[i];
		size++;
		}
		
	if (f_v) {
		cout << "projective_space::create_Maruta_Hamada_arc: after adding the special point, there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}

}

void projective_space::create_Maruta_Hamada_arc2(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT data[] = { 
		1,6,2, 1,11,4, 1,5,5, 1,2,6, 1,10,7, 1,12,8, 1,7,10, 1,4,11, 1,8,12, 
		0,1,10, 0,1,12, 0,1,4, 1,0,1, 1,0,3, 1,0,9, 1,1,0, 1,3,0, 1,9,0,
		1,4,4, 1,4,12, 1,12,4, 1,10,10, 1,10,12, 1,12,10 
		 };
	INT points[24];
	INT i, j, a;
	INT L[9];

	if (n != 2) {
		cout << "projective_space::create_Maruta_Hamada_arc2 n != 2" << endl;
		exit(1);
		}
	if (q != 13) {
		cout << "projective_space::create_Maruta_Hamada_arc2 q != 13" << endl;
		exit(1);
		}
	for (i = 0; i < 24; i++) {
		points[i] = rank_point(data + i * 3);
		cout << "point " << i << " has rank " << points[i] << endl;
		}

	if (f_v) {
		cout << "projective_space::create_Maruta_Hamada_arc2() points: ";
		INT_vec_print(cout, points, 25);
		cout << endl;
		}
	for (i = 0; i < 9; i++) {
		L[i] = Polarity_point_to_hyperplane[points[i]];
		}
	size = 0;
	for (i = 0; i < 9; i++) {
		for (j = i + 1; j < 9; j++) {
			a = line_intersection(L[i], L[j]);
			the_arc[size++] = a;
			}
		}
	for (i = 9; i < 24; i++) {
		the_arc[size++] = points[i];
		}
	if (f_v) {
		cout << "projective_space::create_Maruta_Hamada_arc2: there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}
}


void projective_space::create_pasch_arc(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT data[] = {1,1,1, 1,0,0, 0,1,1,  0,1,0,  1,0,1 };
	INT points[5];
	INT i, j, b, h, idx;
	INT L[4];

	if (n != 2) {
		cout << "projective_space::create_pasch_arc n != 2" << endl;
		exit(1);
		}
#if 0
	if (q != 8) {
		cout << "projective_space::create_pasch_arc q != 8" << endl;
		exit(1);
		}
#endif
	for (i = 0; i < 5; i++) {
		points[i] = rank_point(data + i * 3);
		}

	if (f_v) {
		cout << "projective_space::create_pasch_arc() points: ";
		INT_vec_print(cout, points, 5);
		cout << endl;
		}
	
	L[0] = Line_through_two_points[points[0] * N_points + points[1]];
	L[1] = Line_through_two_points[points[0] * N_points + points[3]];
	L[2] = Line_through_two_points[points[2] * N_points + points[3]];
	L[3] = Line_through_two_points[points[1] * N_points + points[4]];
	
	if (f_v) {
		cout << "L:";
		INT_vec_print(cout, L, 4);
		cout << endl;
		}

	size = 0;	
	for (h = 0; h < 4; h++) {
		for (i = 0; i < r; i++) {
			b = Lines[L[h] * r + i];
			if (INT_vec_search(the_arc, size, b, idx)) {
				continue;
				}
			for (j = size; j > idx; j--) {
				the_arc[j] = the_arc[j - 1];
				}
			the_arc[idx] = b;
			size++;
			}
		}
	if (f_v) {
		cout << "there are " << size << " points on the pasch lines: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}

	
	v[0] = 1;
	v[1] = 1;
	v[2] = 0;
	b = rank_point(v);
	if (INT_vec_search(the_arc, size, b, idx)) {
		cout << "error, special point already there" << endl;
		exit(1);
		}
	for (j = size; j > idx; j--) {
		the_arc[j] = the_arc[j - 1];
		}
	the_arc[idx] = b;
	size++;
		
	if (f_v) {
		cout << "projective_space::create_pasch_arc: after adding the special point, there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}

}

void projective_space::create_Cheon_arc(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT data[] = {1,0,0, 0,1,0, 0,0,1 };
	INT points[3];
	INT i, j, a, b, c, h, idx, t;
	INT L[3];
	INT pencil[9];
	INT Pencil[21];

	if (n != 2) {
		cout << "projective_space::create_Cheon_arc n != 2" << endl;
		exit(1);
		}
#if 0
	if (q != 8) {
		cout << "projective_space::create_Cheon_arc q != 8" << endl;
		exit(1);
		}
#endif
	for (i = 0; i < 3; i++) {
		points[i] = rank_point(data + i * 3);
		}

	if (f_v) {
		cout << "points: ";
		INT_vec_print(cout, points, 5);
		cout << endl;
		}
	
	L[0] = Line_through_two_points[points[0] * N_points + points[1]];
	L[1] = Line_through_two_points[points[1] * N_points + points[2]];
	L[2] = Line_through_two_points[points[2] * N_points + points[0]];
	
	if (f_v) {
		cout << "L:";
		INT_vec_print(cout, L, 3);
		cout << endl;
		}

	size = 0;	
	for (h = 0; h < 3; h++) {
		for (i = 0; i < r; i++) {
			b = Lines[L[h] * r + i];
			if (INT_vec_search(the_arc, size, b, idx)) {
				continue;
				}
			for (j = size; j > idx; j--) {
				the_arc[j] = the_arc[j - 1];
				}
			the_arc[idx] = b;
			size++;
			}
		}

	if (f_v) {
		cout << "projective_space::create_Cheon_arc there are " << size << " points on the 3 lines: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}




	for (h = 0; h < 3; h++) {

		if (f_v) {
			cout << "h=" << h << endl;
			}

		for (i = 0; i < r; i++) {
			pencil[i] = Lines_on_point[points[h] * r + i];
			}


		j = 0;
		for (i = 0; i < r; i++) {
			b = pencil[i];
			if (b == L[0] || b == L[1] || b == L[2])
				continue;
			Pencil[h * 7 + j] = b;
			j++;
			}
		if (j != 7) {
			cout << "j=" << j << endl;
			exit(1);
			}
		}
	if (f_v) {
		cout << "Pencil:" << endl;
		print_integer_matrix_width(cout, Pencil, 3, 7, 7, 4);
		}

	for (i = 0; i < 7; i++) {
		a = Pencil[0 * 7 + i];
		for (j = 0; j < 7; j++) {
			b = Pencil[1 * 7 + j];
			if (f_v) {
				cout << "i=" << i << " a=" << a << " j=" << j << " b=" << b << endl;
				}
			c = Line_intersection[a * N_lines + b];
			if (f_v) {
				cout << "c=" << c << endl;
				}
			if (INT_vec_search(the_arc, size, c, idx)) {
				continue;
				}
			for (t = size; t > idx; t--) {
				the_arc[t] = the_arc[t - 1];
				}
			the_arc[idx] = c;
			size++;
#if 0
			if (size > 31) {
				cout << "create_Cheon_arc size=" << size << endl;
				}
#endif
			}
		}
	if (f_v) {
		cout << "there are " << size << " points on the Cheon lines: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}

	
}


void projective_space::create_regular_hyperoval(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	INT v[3];

	if (n != 2) {
		cout << "projective_space::create_regular_hyperoval n != 2" << endl;
		exit(1);
		}

	for (i = 0; i < q; i++) {
		v[0] = F->mult(i, i);
		v[1] = i;
		v[2] = 1;
		the_arc[i] = rank_point(v);		
		}
	v[0] = 1;
	v[1] = 0;
	v[2] = 0;
	the_arc[q] = rank_point(v);
	
	v[0] = 0;
	v[1] = 1;
	v[2] = 0;
	the_arc[q + 1] = rank_point(v);
	
	size = q + 2;

	if (f_v) {
		cout << "projective_space::create_regular_hyperoval: there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}
}

void projective_space::create_translation_hyperoval(INT *the_arc, INT &size, 
	INT exponent, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	INT v[3];

	if (f_v) {
		cout << "projective_space::create_translation_hyperoval" << endl;
		cout << "exponent = " << exponent << endl;
		}
	if (n != 2) {
		cout << "projective_space::create_translation_hyperoval n != 2" << endl;
		exit(1);
		}

	for (i = 0; i < q; i++) {
		v[0] = F->frobenius_power(i, exponent);
		v[1] = i;
		v[2] = 1;
		the_arc[i] = rank_point(v);		
		}
	v[0] = 1;
	v[1] = 0;
	v[2] = 0;
	the_arc[q] = rank_point(v);
	
	v[0] = 0;
	v[1] = 1;
	v[2] = 0;
	the_arc[q + 1] = rank_point(v);
	
	size = q + 2;

	if (f_v) {
		cout << "projective_space::create_translation_hyperoval: there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}
	if (f_v) {
		cout << "projective_space::create_translation_hyperoval done" << endl;
		}
}

void projective_space::create_Segre_hyperoval(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	INT v[3];

	if (n != 2) {
		cout << "projective_space::create_Segre_hyperoval n != 2" << endl;
		exit(1);
		}

	for (i = 0; i < q; i++) {
		v[0] = F->power(i, 6);
		v[1] = i;
		v[2] = 1;
		the_arc[i] = rank_point(v);		
		}
	v[0] = 1;
	v[1] = 0;
	v[2] = 0;
	the_arc[q] = rank_point(v);
	
	v[0] = 0;
	v[1] = 1;
	v[2] = 0;
	the_arc[q + 1] = rank_point(v);
	
	size = q + 2;

	if (f_v) {
		cout << "projective_space::create_Segre_hyperoval: there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}
}

void projective_space::create_Payne_hyperoval(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	INT v[3];
	longinteger_domain D;
	longinteger_object a, b, u, u2, g;
	INT exponent;
	INT one_sixth, one_half, five_sixth;

	if (f_v) {
		cout << "projective_space::create_Payne_hyperoval" << endl;
		}
	if (n != 2) {
		cout << "projective_space::create_Payne_hyperoval n != 2" << endl;
		exit(1);
		}
	exponent = q - 1;
	a.create(6);
	b.create(exponent);
	
	D.extended_gcd(a, b, g, u, u2, 0 /* verbose_level */);
	one_sixth = u.as_INT();
	while (one_sixth < 0) {
		one_sixth += exponent;
		}
	if (f_v) {
		cout << "one_sixth = " << one_sixth << endl;
		}

	a.create(2);
	D.extended_gcd(a, b, g, u, u2, 0 /* verbose_level */);
	one_half = u.as_INT();
	while (one_half < 0) {
		one_half += exponent;
		}
	if (f_v) {
		cout << "one_half = " << one_half << endl;
		}

	five_sixth = (5 * one_sixth) % exponent;
	if (f_v) {
		cout << "five_sixth = " << five_sixth << endl;
		}

	for (i = 0; i < q; i++) {
		v[0] = F->add3(
			F->power(i, one_sixth), 
			F->power(i, one_half), 
			F->power(i, five_sixth));
		v[1] = i;
		v[2] = 1;
		the_arc[i] = rank_point(v);		
		}
	v[0] = 1;
	v[1] = 0;
	v[2] = 0;
	the_arc[q] = rank_point(v);
	
	v[0] = 0;
	v[1] = 1;
	v[2] = 0;
	the_arc[q + 1] = rank_point(v);
	
	size = q + 2;

	if (f_v) {
		cout << "projective_space::create_Payne_hyperoval: there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}
}

void projective_space::create_Cherowitzo_hyperoval(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	INT v[3];
	INT h;
	INT sigma;
	INT exponent, one_half, e1, e2, e3;

	if (f_v) {
		cout << "projective_space::create_Cherowitzo_hyperoval" << endl;
		}
	if (n != 2) {
		cout << "projective_space::create_Cherowitzo_hyperoval n != 2" << endl;
		exit(1);
		}
	h = F->e;
	if (EVEN(h)) {
		cout << "projective_space::create_Cherowitzo_hyperoval field degree must be odd" << endl;
		exit(1);
		}
	if (F->p != 2) {
		cout << "projective_space::create_Cherowitzo_hyperoval needs characteristic 2" << endl;
		exit(1);
		}
	exponent = q - 1;
	one_half = (h + 1) >> 1;
	sigma = i_power_j(2, one_half);
	e1 = sigma;
	e2 = (sigma + 2) % exponent;
	e3 = (3 * sigma + 4) % exponent;

	for (i = 0; i < q; i++) {
		v[0] = F->add3(
			F->power(i, e1), 
			F->power(i, e2), 
			F->power(i, e3));
		v[1] = i;
		v[2] = 1;
		the_arc[i] = rank_point(v);		
		}
	v[0] = 1;
	v[1] = 0;
	v[2] = 0;
	the_arc[q] = rank_point(v);
	
	v[0] = 0;
	v[1] = 1;
	v[2] = 0;
	the_arc[q + 1] = rank_point(v);
	
	size = q + 2;

	if (f_v) {
		cout << "projective_space::create_Cherowitzo_hyperoval: there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}
}

void projective_space::create_OKeefe_Penttila_hyperoval_32(INT *the_arc, INT &size, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	INT v[3];

	if (f_v) {
		cout << "projective_space::create_OKeefe_Penttila_hyperoval_32" << endl;
		}
	if (n != 2) {
		cout << "projective_space::create_OKeefe_Penttila_hyperoval_32 n != 2" << endl;
		exit(1);
		}
	if (F->q != 32) {
		cout << "projective_space::create_OKeefe_Penttila_hyperoval_32 needs q=32" << endl;
		exit(1);
		}

	for (i = 0; i < q; i++) {
		v[0] = OKeefe_Penttila_32(F, i);
		v[1] = i;
		v[2] = 1;
		the_arc[i] = rank_point(v);		
		}
	v[0] = 1;
	v[1] = 0;
	v[2] = 0;
	the_arc[q] = rank_point(v);
	
	v[0] = 0;
	v[1] = 1;
	v[2] = 0;
	the_arc[q + 1] = rank_point(v);
	
	size = q + 2;

	if (f_v) {
		cout << "projective_space::create_OKeefe_Penttila_hyperoval_32: there are " << size << " points on the arc: ";
		INT_vec_print(cout, the_arc, size);
		cout << endl;
		}
}




void projective_space::line_intersection_type(INT *set, INT set_size, INT *type, INT verbose_level)
// type[N_lines]
{
	INT f_v = (verbose_level >= 1);
	INT i, j, a, b;

	if (f_v) {
		cout << "projective_space::line_intersection_type" << endl;
		}
	if (Lines_on_point == NULL) {
		line_intersection_type_basic(set, set_size, type, verbose_level);
		}
	else {
		for (i = 0; i < N_lines; i++) {
			type[i] = 0;
			}
		for (i = 0; i < set_size; i++) {
			a = set[i];
			for (j = 0; j < r; j++) {
				b = Lines_on_point[a * r + j];
				type[b]++;
				}
			}
		}
}

void projective_space::line_intersection_type_basic(INT *set, INT set_size, INT *type, INT verbose_level)
// type[N_lines]
{
	INT f_v = (verbose_level >= 1);
	INT rk, h, i, j, d;
	INT *M;

	if (f_v) {
		cout << "projective_space::line_intersection_type_basic" << endl;
		}
	d = n + 1;
	M = NEW_INT(3 * d);
	for (rk = 0; rk < N_lines; rk++) {
		type[rk] = 0;
		Grass_lines->unrank_INT(rk, 0 /* verbose_level */);
		for (h = 0; h < set_size; h++) {
			for (i = 0; i < 2; i++) {
				for (j = 0; j < d; j++) {
					M[i * d + j] = Grass_lines->M[i * d + j];
					}
				}
			unrank_point(M + 2 * d, set[h]);
			if (F->rank_of_rectangular_matrix(M, 3, d, 0 /*verbose_level*/) == 2) {
				type[rk]++;
				}
			} // next h
		} // next rk
	FREE_INT(M);
}

void projective_space::plane_intersection_type_basic(INT *set, INT set_size, 
	INT *type, INT verbose_level)
// type[N_planes]
{
	INT f_v = (verbose_level >= 1);
	INT rk, h, i, j, d, N_planes;
	INT *M;
	grassmann *G;

	if (f_v) {
		cout << "projective_space::plane_intersection_type_basic" << endl;
		}
	d = n + 1;
	M = NEW_INT(4 * d);
	G = new grassmann;

	G->init(d, 3, F, 0 /* verbose_level */);

	N_planes = nb_rk_k_subspaces_as_INT(3);
	if (f_v) {
		cout << "projective_space::plane_intersection_type_basic N_planes=" << N_planes << endl;
		}
	
	for (rk = 0; rk < N_planes; rk++) {
		if (rk && (rk % ONE_MILLION) == 0) {
			cout << "projective_space::plane_intersection_type_basic rk=" << rk << endl;
			}
		type[rk] = 0;
		G->unrank_INT(rk, 0 /* verbose_level */);
		for (h = 0; h < set_size; h++) {
			for (i = 0; i < 3; i++) {
				for (j = 0; j < d; j++) {
					M[i * d + j] = G->M[i * d + j];
					}
				}
			unrank_point(M + 3 * d, set[h]);
			if (F->rank_of_rectangular_matrix(M, 4, d, 0 /*verbose_level*/) == 3) {
				type[rk]++;
				}
			} // next h
		} // next rk
	FREE_INT(M);
	delete G;
}

void projective_space::hyperplane_intersection_type_basic(INT *set, INT set_size, INT *type, INT verbose_level)
// type[N_hyperplanes]
{
	INT f_v = (verbose_level >= 1);
	INT rk, h, i, j, d, N_hyperplanes;
	INT *M;
	grassmann *G;

	if (f_v) {
		cout << "projective_space::hyperplane_intersection_type_basic" << endl;
		}
	d = n + 1;
	M = NEW_INT(4 * d);
	G = new grassmann;

	G->init(d, d - 1, F, 0 /* verbose_level */);

	N_hyperplanes = nb_rk_k_subspaces_as_INT(d - 1);
	
	for (rk = 0; rk < N_hyperplanes; rk++) {
		type[rk] = 0;
		G->unrank_INT(rk, 0 /* verbose_level */);
		for (h = 0; h < set_size; h++) {
			for (i = 0; i < d - 1; i++) {
				for (j = 0; j < d; j++) {
					M[i * d + j] = G->M[i * d + j];
					}
				}
			unrank_point(M + (d - 1) * d, set[h]);
			if (F->rank_of_rectangular_matrix(M, d, d, 0 /*verbose_level*/) == d - 1) {
				type[rk]++;
				}
			} // next h
		} // next rk
	FREE_INT(M);
	delete G;
}




void projective_space::line_intersection_type_collected(INT *set, INT set_size, INT *type_collected, INT verbose_level)
// type[set_size + 1]
{
	INT i, a;
	INT *type;


	type = NEW_INT(N_lines);

	line_intersection_type(set, set_size, type, verbose_level);

	for (i = 0; i <= set_size; i++) {
		type_collected[i] = 0;
		}
	for (i = 0; i < N_lines; i++) {
		a = type[i];
		type_collected[a]++;
		}
	FREE_INT(type);
}

void projective_space::point_types(INT *set_of_lines, INT set_size, INT *type, INT verbose_level)
{
	INT i, j, a, b;

	for (i = 0; i < N_points; i++) {
		type[i] = 0;
		}
	for (i = 0; i < set_size; i++) {
		a = set_of_lines[i];
		for (j = 0; j < k; j++) {
			b = Lines[a * k + j];
			type[b]++;
			}
		}
}

void projective_space::find_external_lines(INT *set, INT set_size, 
	INT *external_lines, INT &nb_external_lines, INT verbose_level)
{
	INT *type;
	INT i;
	
	nb_external_lines = 0;
	type = NEW_INT(N_lines);
	line_intersection_type(set, set_size, type, verbose_level);
	for (i = 0; i < N_lines; i++) {
		if (type[i]) {
			continue;
			}
		external_lines[nb_external_lines++] = i;
		}
	FREE_INT(type);
}

void projective_space::find_tangent_lines(INT *set, INT set_size, 
	INT *tangent_lines, INT &nb_tangent_lines, INT verbose_level)
{
	INT *type;
	INT i;
	
	nb_tangent_lines = 0;
	type = NEW_INT(N_lines);
	line_intersection_type(set, set_size, type, verbose_level);
	for (i = 0; i < N_lines; i++) {
		if (type[i] != 1) {
			continue;
			}
		tangent_lines[nb_tangent_lines++] = i;
		}
	FREE_INT(type);
}

void projective_space::find_secant_lines(INT *set, INT set_size, 
	INT *secant_lines, INT &nb_secant_lines, INT verbose_level)
{
	INT *type;
	INT i;
	
	nb_secant_lines = 0;
	type = NEW_INT(N_lines);
	line_intersection_type(set, set_size, type, verbose_level);
	for (i = 0; i < N_lines; i++) {
		if (type[i] != 2) {
			continue;
			}
		secant_lines[nb_secant_lines++] = i;
		}
	FREE_INT(type);
}

void projective_space::find_k_secant_lines(INT *set, INT set_size, INT k, 
	INT *secant_lines, INT &nb_secant_lines, INT verbose_level)
{
	INT *type;
	INT i;
	
	nb_secant_lines = 0;
	type = NEW_INT(N_lines);
	line_intersection_type(set, set_size, type, verbose_level);
	for (i = 0; i < N_lines; i++) {
		if (type[i] != k) {
			continue;
			}
		secant_lines[nb_secant_lines++] = i;
		}
	FREE_INT(type);
}


void projective_space::Baer_subline(INT *pts3, INT *&pts, INT &nb_pts, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	INT *M;
	INT *Basis;
	INT *N; // local coordinates w.r.t. basis
	INT *base_cols;
	INT *z;
	INT rk;
	INT len;
	INT i, j;

	if (f_v) {
		cout << "projective_space::Baer_subline" << endl;
		}
	if (ODD(F->e)) {
		cout << "projective_space::Baer_subline field degree must be even (because we need a quadratic subfield)" << endl;
		exit(1);
		}
	len = n + 1;
	M = NEW_INT(3 * len);
	base_cols = NEW_INT(len);
	z = NEW_INT(len);
	for (j = 0; j < 3; j++) {
		unrank_point(M + j * len, pts3[j]);
		}
	if (f_vv) {
		cout << "projective_space::Baer_subline" << endl;
		cout << "M=" << endl;
		print_integer_matrix_width(cout, M, 3, len, len, F->log10_of_q);
		}
	rk = F->Gauss_simple(M, 3, len, base_cols, verbose_level - 3);
	if (f_vv) {
		cout << "projective_space::Baer_subline" << endl;
		cout << "has rank " << rk << endl;
		cout << "base_cols=";
		INT_vec_print(cout, base_cols, rk);
		cout << endl;
		cout << "basis:" << endl;
		print_integer_matrix_width(cout, M, rk, len, len, F->log10_of_q);
		}

	if (rk != 2) {
		cout << "projective_space::Baer_subline: rk should be 2 (points are not collinear)" << endl;
		exit(1);
		}
	
	Basis = NEW_INT(rk * len);
	for (j = 0; j < rk * len; j++) {
		Basis[j] = M[j];
		}
	if (f_vv) {
		cout << "projective_space::Baer_subline basis:" << endl;
		print_integer_matrix_width(cout, Basis, rk, len, len, F->log10_of_q);
		}
		
	N = NEW_INT(3 * rk);
	for (j = 0; j < 3; j++) {
		unrank_point(M + j * len, pts3[j]);
		//cout << "M + j * len:";
		//INT_vec_print(cout, M + j * len, len);
		//cout << endl;
		//cout << "basis:" << endl;
		//print_integer_matrix_width(cout, Basis, rk, 5, 5, P4->F->log10_of_q);
		
		F->reduce_mod_subspace_and_get_coefficient_vector(
			rk, len, Basis, base_cols, 
			M + j * len, N + j * rk, verbose_level - 3);
		}
	//cout << "after reduce_mod_subspace_and_get_coefficient_vector: M=" << endl;
	//print_integer_matrix_width(cout, M, 3, len, len, F->log10_of_q);
	//cout << "(should be all zeros)" << endl;
	if (f_vv) {
		cout << "projective_space::Baer_subline local coordinates in the subspace are N=" << endl;
		print_integer_matrix_width(cout, N, 3, rk, rk, F->log10_of_q);
		}
	INT *Frame;
	INT *base_cols2;
	INT rk2, a;

	Frame = NEW_INT(2 * 3);
	base_cols2 = NEW_INT(3);
	for (j = 0; j < 3; j++) {
		for (i = 0; i < 2; i++) {
			Frame[i * 3 + j] = N[j * 2 + i];
			}
		}
	if (f_vv) {
		cout << "projective_space::Baer_subline Frame=" << endl;
		print_integer_matrix_width(cout, Frame, 2, 3, 3, F->log10_of_q);
		}
	rk2 = F->Gauss_simple(Frame, 2, 3, base_cols2, verbose_level - 3);
	if (rk2 != 2) {
		cout << "projective_space::Baer_subline: rk2 should be 2" << endl;
		exit(1);
		}
	if (f_vv) {
		cout << "projective_space::Baer_subline after Gauss Frame=" << endl;
		print_integer_matrix_width(cout, Frame, 2, 3, 3, F->log10_of_q);
		cout << "projective_space::Baer_subline base_cols2=";
		INT_vec_print(cout, base_cols2, rk2);
		cout << endl;
		}
	for (i = 0; i < 2; i++) {
		a = Frame[i * 3 + 2];
		for (j = 0; j < 2; j++) {
			N[i * 2 + j] = F->mult(a, N[i * 2 + j]);
			}
		}
	if (f_vv) {
		cout << "projective_space::Baer_subline local coordinates in the subspace are N=" << endl;
		print_integer_matrix_width(cout, N, 3, rk, rk, F->log10_of_q);
		}

#if 0
	INT *Local_pts;
	INT *Local_pts_sorted;
	INT w[2];


	Local_pts = NEW_INT(nb_pts);
	Local_pts_sorted = NEW_INT(nb_pts);

	for (i = 0; i < nb_pts; i++) {
		for (j = 0; j < 2; j++) {
			w[j] = N[i * 2 + j];
			}
		PG_element_rank_modified(*F, w, 1, 2, a);
		Local_pts[i] = a;
		Local_pts_sorted[i] = a;
		}
	INT_vec_heapsort(Local_pts_sorted, nb_pts);
	if (f_vv) {
		cout << "Local_pts=" << endl;
		INT_vec_print(cout, Local_pts, nb_pts);
		cout << endl;
		cout << "Local_pts_sorted=" << endl;
		INT_vec_print(cout, Local_pts_sorted, nb_pts);
		cout << endl;
		}
#endif


	INT q0, index, t;


	q0 = i_power_j(F->p, F->e >> 1);
	index = (F->q - 1) / (q0 - 1);
	
	nb_pts = q0 + 1;
	pts = NEW_INT(nb_pts);


	if (f_vv) {
		cout << "projective_space::Baer_subline q0=" << q0 << endl;
		cout << "projective_space::Baer_subline index=" << index << endl;
		cout << "projective_space::Baer_subline nb_pts=" << nb_pts << endl;
		}

#if 0
	for (i = 0; i < 3; i++) {
		for (j = 0; j < len; j++) {
			if (i < 2) {
				z[j] = Basis[i * len + j];
				}
			else {
				z[j] = F->add(Basis[0 * len + j], Basis[1 * len + j]);
				}
			}
		pts[i] = rank_point(z);
		}
#endif
	for (t = 0; t < 3; t++) {
		if (f_vvv) {
			cout << "t=" << t << endl;
			}
		F->mult_vector_from_the_left(N + t * 2, Basis, z, 2, len);
		if (f_vvv) {
			cout << "z=w*Basis";
			INT_vec_print(cout, z, len);
			cout << endl;
			}
		a = rank_point(z);
		pts[t] = a;
		}
	for (t = 2; t < q0; t++) {
		a = F->alpha_power((t - 1) * index);
		if (f_vvv) {
			cout << "t=" << t << " a=" << a << endl;
			}
		for (j = 0; j < 2; j++) {
			w[j] = F->add(N[0 * 2 + j], F->mult(a, N[1 * 2 + j]));
			}
		if (f_vvv) {
			cout << "w=";
			INT_vec_print(cout, w, 2);
			cout << endl;
			}
		F->mult_vector_from_the_left(w, Basis, z, 2, len);
		if (f_vvv) {
			cout << "z=w*Basis";
			INT_vec_print(cout, z, len);
			cout << endl;
			}
		a = rank_point(z);
		pts[t + 1] = a;
		if (f_vvv) {
			cout << "rank=" << a << endl;
			}
#if 0
		PG_element_rank_modified(*F, w, 1, 2, a);
		pts[t] = a;
		if (!INT_vec_search(Local_pts_sorted, nb_pts, a, idx)) {
			ret = FALSE;
			if (f_vv) {
				cout << "did not find this point in the list of points, hence not contained in Baer subline" << endl;
				}
			goto done;
			}
#endif
		
		}

	if (f_vv) {
		cout << "projective_space::Baer_subline The Baer subline is";
		INT_vec_print(cout, pts, nb_pts);
		cout << endl;
		print_set(pts, nb_pts);
		}
	



	FREE_INT(N);
	FREE_INT(M);
	FREE_INT(base_cols);
	FREE_INT(Basis);
	FREE_INT(Frame);
	FREE_INT(base_cols2);
	FREE_INT(z);
}




void projective_space::print_set_numerical(INT *set, INT set_size)
{
	INT i, a;
	INT *v;
	
	v = NEW_INT(n + 1);
	for (i = 0; i < set_size; i++) {
		a = set[i];
		unrank_point(v, a);
		cout << setw(3) << i << " : " << setw(5) << a << " : ";
		INT_vec_print(cout, v, n + 1);
		cout << "=";
		PG_element_normalize_from_front(*F, v, 1, n + 1);
		INT_vec_print(cout, v, n + 1);
		cout << endl;
		}
	FREE_INT(v);
}

void projective_space::print_set(INT *set, INT set_size)
{
	INT i, a;
	INT *v;
	
	v = NEW_INT(n + 1);
	for (i = 0; i < set_size; i++) {
		a = set[i];
		unrank_point(v, a);
		cout << setw(3) << i << " : " << setw(5) << a << " : ";
		F->INT_vec_print(cout, v, n + 1);
		cout << "=";
		PG_element_normalize_from_front(*F, v, 1, n + 1);
		F->INT_vec_print(cout, v, n + 1);
		cout << endl;
		}
	FREE_INT(v);
}

void projective_space::print_line_set_numerical(INT *set, INT set_size)
{
	INT i, a;
	INT *v;
	
	v = NEW_INT(2 * (n + 1));
	for (i = 0; i < set_size; i++) {
		a = set[i];
		unrank_line(v, a);
		cout << setw(3) << i << " : " << setw(5) << a << " : ";
		INT_vec_print(cout, v, 2 * (n + 1));
		cout << endl;
		}
	FREE_INT(v);
}


INT projective_space::is_contained_in_Baer_subline(INT *pts, INT nb_pts, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *subline;
	INT sz;
	INT i, idx, a;
	INT ret = TRUE;

	if (f_v) {
		cout << "projective_space::is_contained_in_Baer_subline pts=" << endl;
		INT_vec_print(cout, pts, nb_pts);
		cout << endl;
		cout << "computing Baer subline determined by the first three points:" << endl;
		}
	Baer_subline(pts, subline, sz, verbose_level - 2);
	if (f_vv) {
		cout << "projective_space::is_contained_in_Baer_subline The Baer subline is:" << endl;
		INT_vec_print(cout, subline, sz);
		cout << endl;
		}
	INT_vec_heapsort(subline, sz);
	for (i = 0; i < nb_pts; i++) {
		a = pts[i];
		if (!INT_vec_search(subline, sz, a, idx)) {
			ret = FALSE;
			if (f_vv) {
				cout << "did not find " << i << "-th point " << a << " in the list of points, hence not contained in Baer subline" << endl;
				}
			goto done;
			}
		
		}
done:
	FREE_INT(subline);
	
	return ret;
}

INT projective_space::determine_hermitian_form_in_plane(INT *pts, INT nb_pts, INT *six_coeffs, INT verbose_level)
// there is a memory problem in this function
// detected 7/14/11
// solved June 17, 2012: 
// coords and system were not freed
// system was allocated too short
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *coords; //[nb_pts * 3];
	INT *system; //[nb_pts * 9];
	INT kernel[9 * 9];
	INT base_cols[9];
	INT i, x, y, z, xq, yq, zq, rk;
	INT Q, q, little_e;
	INT kernel_m, kernel_n;

	if (f_v) {
		cout << "projective_space::determine_hermitian_form_in_plane" << endl;
		}
	coords = NEW_INT(nb_pts * 3);
	system = NEW_INT(nb_pts * 9);
	Q = F->q;
	if (ODD(F->e)) {
		cout << "projective_space::determine_hermitian_form_in_plane field degree must be even" << endl;
		exit(1);
		}
	little_e = F->e >> 1;
	q = i_power_j(F->p, little_e);
	if (f_v) {
		cout << "projective_space::determine_hermitian_form_in_plane Q=" << Q << " q=" << q << endl;
		}
	if (n != 2) {
		cout << "projective_space::determine_hermitian_form_in_plane n != 2" << endl;
		exit(1);
		}
	for (i = 0; i < nb_pts; i++) {
		unrank_point(coords + i * 3, pts[i]);
		}
	if (f_vv) {
		cout << "projective_space::determine_hermitian_form_in_plane points:" << endl;
		print_integer_matrix_width(cout, coords, nb_pts, 3, 3, F->log10_of_q);
		}
	for (i = 0; i < nb_pts; i++) {
		x = coords[i * 3 + 0];
		y = coords[i * 3 + 1];
		z = coords[i * 3 + 2];
		xq = F->frobenius_power(x, little_e);
		yq = F->frobenius_power(y, little_e);
		zq = F->frobenius_power(z, little_e);
		system[i * 9 + 0] = F->mult(x, xq);
		system[i * 9 + 1] = F->mult(y, yq);
		system[i * 9 + 2] = F->mult(z, zq);
		system[i * 9 + 3] = F->mult(x, yq);
		system[i * 9 + 4] = F->mult(y, xq);
		system[i * 9 + 5] = F->mult(x, zq);
		system[i * 9 + 6] = F->mult(z, xq);
		system[i * 9 + 7] = F->mult(y, zq);
		system[i * 9 + 8] = F->mult(z, yq);
		}
	if (f_v) {
		cout << "projective_space::determine_hermitian_form_in_plane system:" << endl;
		print_integer_matrix_width(cout, system, nb_pts, 9, 9, F->log10_of_q);
		}



	rk = F->Gauss_simple(system, nb_pts, 9, base_cols, verbose_level - 2);
	if (f_v) {
		cout << "projective_space::determine_hermitian_form_in_plane rk=" << rk << endl;
		print_integer_matrix_width(cout, system, rk, 9, 9, F->log10_of_q);
		}
#if 0
	if (rk != 8) {
		if (f_v) {
			cout << "projective_space::determine_hermitian_form_in_plane system underdetermined" << endl;
			}
		return FALSE;
		}
#endif
	F->matrix_get_kernel(system, MINIMUM(nb_pts, 9), 9, base_cols, rk, 
		kernel_m, kernel_n, kernel);
	if (f_v) {
		cout << "projective_space::determine_hermitian_form_in_plane kernel:" << endl;
		print_integer_matrix_width(cout, kernel, kernel_m, kernel_n, kernel_n, F->log10_of_q);
		}
	six_coeffs[0] = kernel[0 * kernel_n + 0];
	six_coeffs[1] = kernel[1 * kernel_n + 0];
	six_coeffs[2] = kernel[2 * kernel_n + 0];
	six_coeffs[3] = kernel[3 * kernel_n + 0];
	six_coeffs[4] = kernel[5 * kernel_n + 0];
	six_coeffs[5] = kernel[7 * kernel_n + 0];
	if (f_v) {
		cout << "projective_space::determine_hermitian_form_in_plane six_coeffs:" << endl;
		INT_vec_print(cout, six_coeffs, 6);
		cout << endl;
		}
	FREE_INT(coords);
	FREE_INT(system);
	return TRUE;
}

void projective_space::circle_type_of_line_subset(INT *pts, INT nb_pts, INT *circle_type, INT verbose_level)
// circle_type[nb_pts]
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *subline;
	INT subset[3];
	INT idx_set[3];
	INT sz;
	INT i, idx, a, b;

	if (f_v) {
		cout << "projective_space::circle_type_of_line_subset pts=" << endl;
		INT_vec_print(cout, pts, nb_pts);
		cout << endl;
		//cout << "computing Baer subline determined by the first three points:" << endl;
		}
	for (i = 0; i < nb_pts; i++) {
		circle_type[i] = 0;
		}
	
	first_k_subset(idx_set, nb_pts, 3);
	do {
		for (i = 0; i < 3; i++) {
			subset[i] = pts[idx_set[i]];
			}
		Baer_subline(subset, subline, sz, verbose_level - 2);
		b = 0;
		INT_vec_heapsort(subline, sz);
		for (i = 0; i < nb_pts; i++) {
			a = pts[i];
			if (INT_vec_search(subline, sz, a, idx)) {
				b++;
				}
			}


		if (f_v) {
			cout << "projective_space::circle_type_of_line_subset The Baer subline determined by " << endl;
			INT_vec_print(cout, subset, 3);
			cout << " is ";
			INT_vec_print(cout, subline, sz);
			cout << " which intersects in " << b << " points" << endl;
			}



		FREE_INT(subline);
		circle_type[b]++;
		} while (next_k_subset(idx_set, nb_pts, 3));

	if (f_vv) {
		cout << "projective_space::circle_type_of_line_subset circle_type before fixing =" << endl;
		INT_vec_print(cout, circle_type, nb_pts);
		cout << endl;
		}
	for (i = 4; i < nb_pts; i++) {
		a = INT_n_choose_k(i, 3);
		if (circle_type[i] % a) {
			cout << "projective_space::circle_type_of_line_subset  circle_type[i] % a" << endl;
			exit(1);
			}
		circle_type[i] /= a;
		}
	if (f_vv) {
		cout << "projective_space::circle_type_of_line_subset circle_type after fixing =" << endl;
		INT_vec_print(cout, circle_type, nb_pts);
		cout << endl;
		}
}

void projective_space::create_unital_XXq_YZq_ZYq(INT *U, INT &sz, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	finite_field *FQ;
	INT *v;
	INT e, i, a;
	INT X, Y, Z, Xq, Yq, Zq;

	if (f_v) {
		cout << "projective_space::create_unital_XXq_YZq_ZYq" << endl;
		}
	if (n != 2) {
		cout << "projective_space::create_unital_XXq_YZq_ZYq n != 2" << endl;
		exit(1);
 		}
	if (ODD(F->e)) {
		cout << "projective_space::create_unital_XXq_YZq_ZYq ODD(F->e)" << endl;
		exit(1);
 		}
	FQ = F;
	
	v = NEW_INT(3);
	e = F->e >> 1;
	if (f_vv) {
		cout << "e=" << e << endl;
		}
	sz = 0;
	for (i = 0; i < N_points; i++) {
		unrank_point(v, i);
		if (f_vvv) {
			cout << "i=" << i << " : ";
			INT_vec_print(cout, v, 3);
			//cout << endl;
			}
		X = v[0];
		Y = v[1];
		Z = v[2];
		Xq = F->frobenius_power(X, e);
		Yq = F->frobenius_power(Y, e);
		Zq = F->frobenius_power(Z, e);
		a = F->add3(F->mult(X, Xq), F->mult(Y, Zq), F->mult(Z, Yq));
		if (f_vvv) {
			cout << " a=" << a << endl;
			}
		if (a == 0) {
			//cout << "a=0, adding i=" << i << endl;
			U[sz++] = i;
			//INT_vec_print(cout, U, sz);
			//cout << endl;
			}
		}
	if (f_vv) {
		cout << "we found " << sz << " points:" << endl;	
		INT_vec_print(cout, U, sz);
		cout << endl;
		print_set(U, sz);
		}
	FREE_INT(v);

	if (f_v) {
		cout << "projective_space::create_unital_XXq_YZq_ZYq done" << endl;
		}
}


void projective_space::intersection_of_subspace_with_point_set(
	grassmann *G, INT rk, INT *set, INT set_size, 
	INT *&intersection_set, INT &intersection_set_size, INT verbose_level)
{
	INT h;
	INT d = n + 1;
	INT k = G->k;
	INT *M;

	intersection_set = NEW_INT(set_size);
	M = NEW_INT((k + 1) * d);
	intersection_set_size = 0;

	G->unrank_INT(rk, 0 /* verbose_level */);

	for (h = 0; h < set_size; h++) {
		INT_vec_copy(G->M, M, k * d);
#if 0
		for (i = 0; i < k; i++) {
			for (j = 0; j < d; j++) {
				M[i * d + j] = G->M[i * d + j];
				}
			}
#endif
		unrank_point(M + k * d, set[h]);
		if (F->rank_of_rectangular_matrix(M, k + 1, d, 0 /*verbose_level*/) == k) {
			intersection_set[intersection_set_size++] = set[h];
			}
		} // next h

	FREE_INT(M);
}

void projective_space::intersection_of_subspace_with_point_set_rank_is_longinteger(
	grassmann *G, longinteger_object &rk, INT *set, INT set_size, 
	INT *&intersection_set, INT &intersection_set_size, INT verbose_level)
{
	INT h;
	INT d = n + 1;
	INT k = G->k;
	INT *M;

	intersection_set = NEW_INT(set_size);
	M = NEW_INT((k + 1) * d);
	intersection_set_size = 0;

	G->unrank_longinteger(rk, 0 /* verbose_level */);

	for (h = 0; h < set_size; h++) {
		INT_vec_copy(G->M, M, k * d);
#if 0
		for (i = 0; i < k; i++) {
			for (j = 0; j < d; j++) {
				M[i * d + j] = G->M[i * d + j];
				}
			}
#endif
		unrank_point(M + k * d, set[h]);
		if (F->rank_of_rectangular_matrix(M, k + 1, d, 0 /*verbose_level*/) == k) {
			intersection_set[intersection_set_size++] = set[h];
			}
		} // next h

	FREE_INT(M);
}

void projective_space::plane_intersection_invariant(grassmann *G, 
	INT *set, INT set_size, 
	INT *&intersection_type, INT &highest_intersection_number, 
	INT *&intersection_matrix, INT &nb_planes, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	longinteger_object *R;
	INT **Pts_on_plane;
	INT *nb_pts_on_plane;
	INT nb_planes_total;
	INT i, j, a, u, f, l, ii;

	if (f_v) {
		cout << "projective_space::plane_intersection_invariant" << endl;
		}
	plane_intersection_type_fast(G, 
		set, set_size, 
		R, Pts_on_plane, nb_pts_on_plane, nb_planes_total, 
		verbose_level - 1);

	classify C;
	INT f_second = FALSE;

	C.init(nb_pts_on_plane, nb_planes_total, f_second, 0);
	if (f_v) {
		cout << "projective_space::plane_intersection_invariant plane-intersection type: ";
		C.print(FALSE /* f_backwards*/);
		}

	if (f_v) {
		cout << "The plane intersection type is (";
		C.print_naked(FALSE /*f_backwards*/);
		cout << ")" << endl << endl;
		}
	f = C.type_first[C.nb_types - 1];
	highest_intersection_number = C.data_sorted[f];
	intersection_type = NEW_INT(highest_intersection_number + 1);

	INT_vec_zero(intersection_type, highest_intersection_number + 1);
	
	for (i = 0; i < C.nb_types; i++) {
		f = C.type_first[i];
		l = C.type_len[i];
		a = C.data_sorted[f];
		intersection_type[a] = l;
		}
	f = C.type_first[C.nb_types - 1];
	nb_planes = C.type_len[C.nb_types - 1];

	INT *Incma, *Incma_t, *IIt, *ItI;
	
	Incma = NEW_INT(set_size * nb_planes);
	Incma_t = NEW_INT(nb_planes * set_size);
	IIt = NEW_INT(set_size * set_size);
	ItI = NEW_INT(nb_planes * nb_planes);


	for (i = 0; i < set_size * nb_planes; i++) {
		Incma[i] = 0;
		}
	for (i = 0; i < nb_planes; i++) {
		ii = C.sorting_perm_inv[f + i];
		for (j = 0; j < nb_pts_on_plane[ii]; j++) {
			a = Pts_on_plane[ii][j];
			Incma[a * nb_planes + i] = 1;
			}
		}
	if (f_vv) {
		cout << "Incidence matrix:" << endl;
		print_integer_matrix_width(cout, Incma, set_size, nb_planes, nb_planes, 1);
		}
	for (i = 0; i < set_size; i++) {
		for (j = 0; j < set_size; j++) {
			a = 0;
			for (u = 0; u < nb_planes; u++) {
				a += Incma[i * nb_planes + u] * Incma_t[u * set_size + j];
				}
			IIt[i * set_size + j] = a;
			}
		}
	if (f_vv) {
		cout << "I * I^\\top = " << endl;
		print_integer_matrix_width(cout, IIt, set_size, set_size, set_size, 2);
		}
	for (i = 0; i < nb_planes; i++) {
		for (j = 0; j < nb_planes; j++) {
			a = 0;
			for (u = 0; u < set_size; u++) {
				a += Incma[u * nb_planes + i] * Incma[u * nb_planes + j];
				}
			ItI[i * nb_planes + j] = a;
			}
		}
	if (f_v) {
		cout << "I^\\top * I = " << endl;
		print_integer_matrix_width(cout, ItI, nb_planes, nb_planes, nb_planes, 3);
		}
	
	intersection_matrix = NEW_INT(nb_planes * nb_planes);
	INT_vec_copy(ItI, intersection_matrix, nb_planes * nb_planes);
#if 0
	for (i = 0; i < nb_planes; i++) {
		for (j = 0; j < nb_planes; j++) {
			intersection_matrix[i * nb_planes + j] = ItI[i * nb_planes + j];
			}
		}
#endif

	FREE_INT(Incma);
	FREE_INT(Incma_t);
	FREE_INT(IIt);
	FREE_INT(ItI);


	for (i = 0; i < nb_planes_total; i++) {
		FREE_INT(Pts_on_plane[i]);
		}
	FREE_PINT(Pts_on_plane);
	FREE_INT(nb_pts_on_plane);
	delete [] R;
}

void projective_space::plane_intersection_type(grassmann *G, 
	INT *set, INT set_size, 
	INT *&intersection_type, INT &highest_intersection_number, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	longinteger_object *R;
	INT **Pts_on_plane;
	INT *nb_pts_on_plane;
	INT nb_planes;
	INT i, f, l, a;

	if (f_v) {
		cout << "projective_space::plane_intersection_type" << endl;
		}
	plane_intersection_type_fast(G, 
		set, set_size, 
		R, Pts_on_plane, nb_pts_on_plane, nb_planes, 
		verbose_level - 1);

	classify C;
	INT f_second = FALSE;

	C.init(nb_pts_on_plane, nb_planes, f_second, 0);
	if (f_v) {
		cout << "projective_space::plane_intersection_type plane-intersection type: ";
		C.print(FALSE /*f_backwards*/);
		}

	if (f_v) {
		cout << "The plane intersection type is (";
		C.print_naked(FALSE /*f_backwards*/);
		cout << ")" << endl << endl;
		}

	f = C.type_first[C.nb_types - 1];
	highest_intersection_number = C.data_sorted[f];
	intersection_type = NEW_INT(highest_intersection_number + 1);
	INT_vec_zero(intersection_type, highest_intersection_number + 1);
#if 0
	for (i = 0; i <= highest_intersection_number; i++) {
		intersection_type[i] = 0;
		}
#endif
	for (i = 0; i < C.nb_types; i++) {
		f = C.type_first[i];
		l = C.type_len[i];
		a = C.data_sorted[f];
		intersection_type[a] = l;
		}

	for (i = 0; i < nb_planes; i++) {
		FREE_INT(Pts_on_plane[i]);
		}
	FREE_PINT(Pts_on_plane);
	FREE_INT(nb_pts_on_plane);
	delete [] R;

}

void projective_space::plane_intersections(grassmann *G, 
	INT *set, INT set_size, 
	longinteger_object *&R, set_of_sets &SoS, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT **Pts_on_plane;
	INT *nb_pts_on_plane;
	INT nb_planes;
	INT i;
	
	if (f_v) {
		cout << "projective_space::plane_intersections" << endl;
		}
	plane_intersection_type_fast(G, 
		set, set_size, 
		R, Pts_on_plane, nb_pts_on_plane, nb_planes, 
		verbose_level - 1);
	if (f_v) {
		cout << "projective_space::plane_intersections before Sos.init" << endl;
		}
	SoS.init(set_size, nb_planes, Pts_on_plane, nb_pts_on_plane, verbose_level - 1);
	if (f_v) {
		cout << "projective_space::plane_intersections after Sos.init" << endl;
		}
	for (i = 0; i < nb_planes; i++) {
		FREE_INT(Pts_on_plane[i]);
		}
	FREE_PINT(Pts_on_plane);
	FREE_INT(nb_pts_on_plane);
	if (f_v) {
		cout << "projective_space::plane_intersections done" << endl;
		}
}

void projective_space::plane_intersection_type_slow(grassmann *G, 
	INT *set, INT set_size, 
	longinteger_object *&R, INT **&Pts_on_plane, INT *&nb_pts_on_plane, INT &len, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	//INT f_v3 = (verbose_level >= 3);
	INT r, rk, i, u, d, N_planes, l;

	INT *Basis;
	INT *Basis_save;
	INT *Coords;

	if (f_v) {
		cout << "projective_space::plane_intersection_type_slow" << endl;
		}
	if (f_vv) {
		print_set_numerical(set, set_size);
		}
	if (!test_if_set_with_return_value(set, set_size)) {
		cout << "projective_space::plane_intersection_type_slow the input set if not a set" << endl;
		exit(1);
		}
	d = n + 1;
	N_planes = nb_rk_k_subspaces_as_INT(3);

	if (f_v) {
		cout << "N_planes=" << N_planes << endl;
		}
	// allocate data that is returned:
	R = new longinteger_object[N_planes];
	Pts_on_plane = NEW_PINT(N_planes);
	nb_pts_on_plane = NEW_INT(N_planes);

	// allocate temporary data:
	Basis = NEW_INT(4 * d);
	Basis_save = NEW_INT(4 * d);
	Coords = NEW_INT(set_size * d);
	
	for (i = 0; i < set_size; i++) {
		unrank_point(Coords + i * d, set[i]);
		}
	if (f_vv) {
		cout << "projective_space::plane_intersection_type_fast Coords:" << endl;
		INT_matrix_print(Coords, set_size, d);
		}

	l = 0;
	for (rk = 0; rk < N_planes; rk++) {

		G->unrank_INT(rk, 0 /* verbose_level */);
		for (i = 0; i < 3 * d; i++) {
			Basis_save[i] = G->M[i];
			}
		INT *pts_on_plane;
		INT nb = 0;
	
		pts_on_plane = NEW_INT(set_size);
			
		for (u = 0; u < set_size; u++) {
			for (i = 0; i < 3 * d; i++) {
				Basis[i] = Basis_save[i];
				}
			for (i = 0; i < d; i++) {
				Basis[3 * d + i] = Coords[u * d + i];
				}
			r = F->rank_of_rectangular_matrix(Basis, 4, d, 0 /* verbose_level */);
			if (r < 4) {
				pts_on_plane[nb++] = u;
				}
			}
		Pts_on_plane[l] = pts_on_plane;
		nb_pts_on_plane[l] = nb;
		R[l].create(rk);
		l++;
		} // rk
	len = l;
	
	FREE_INT(Basis);
	FREE_INT(Basis_save);
	FREE_INT(Coords);
	if (f_v) {
		cout << "projective_space::plane_intersection_type_slow done" << endl;
		}
}

void projective_space::plane_intersection_type_fast(grassmann *G, 
	INT *set, INT set_size, 
	longinteger_object *&R, INT **&Pts_on_plane, INT *&nb_pts_on_plane, INT &len, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_v3 = (verbose_level >= 3);
	INT r, rk, rr, h, i, j, a, d, N_planes, N, N2, idx, l;

	INT *Basis;
	INT *Basis_save;
	INT *f_subset_done;
	INT *rank_idx;
	INT *Coords;

	INT subset[3];
	INT subset2[3];
	INT subset3[3];
	longinteger_object plane_rk, aa;
	INT *pts_on_plane;

	if (f_v) {
		cout << "projective_space::plane_intersection_type_fast" << endl;
		}
	if (f_vv) {
		print_set_numerical(set, set_size);
		}

	if (!test_if_set_with_return_value(set, set_size)) {
		cout << "projective_space::plane_intersection_type_fast the input set if not a set" << endl;
		exit(1);
		}
	d = n + 1;
	N_planes = nb_rk_k_subspaces_as_INT(3);
	N = INT_n_choose_k(set_size, 3);

	if (f_v) {
		cout << "N_planes=" << N_planes << endl;
		cout << "N=number of 3-subsets of the set=" << N << endl;
		}
	
	// allocate data that is returned:
	R = new longinteger_object[N];
	Pts_on_plane = NEW_PINT(N);
	nb_pts_on_plane = NEW_INT(N);

	// allocate temporary data:
	Basis = NEW_INT(4 * d);
	Basis_save = NEW_INT(4 * d);
	rank_idx = NEW_INT(N);
	f_subset_done = NEW_INT(N);
	Coords = NEW_INT(set_size * d);
	
	for (i = 0; i < N; i++) {
		f_subset_done[i] = FALSE;
		}
	for (i = 0; i < set_size; i++) {
		unrank_point(Coords + i * d, set[i]);
		}
	if (f_vv) {
		cout << "projective_space::plane_intersection_type_fast Coords:" << endl;
		INT_matrix_print(Coords, set_size, d);
		}


	len = 0;
	for (rk = 0; rk < N; rk++) {
		unrank_k_subset(rk, subset, set_size, 3);
		if (f_v) {
			cout << rk << "-th subset ";
			INT_vec_print(cout, subset, 3);
			cout << endl;
			}
		if (f_subset_done[rk]) {
			if (f_v) {
				cout << "skipping" << endl;
				}
			continue;
			}
		for (j = 0; j < 3; j++) {
			a = subset[j];
			//a = set[subset[j]];
			INT_vec_copy(Coords + a * d, Basis + j * d, d);
#if 0
			for (u = 0; u < d; u++) {
				Basis[j * d + u] = Coords[a * d + u];
				}
#endif
			//unrank_point(Basis + j * d, a);
			}
		if (f_v3) {
			cout << "subset: ";
			INT_vec_print(cout, subset, 3);
			cout << " corresponds to Basis:" << endl;
			INT_matrix_print(Basis, 3, d);
			}
		r = F->rank_of_rectangular_matrix(Basis, 3, d, 0 /* verbose_level */);
		if (r < 3) {
			if (TRUE || f_v) {
				cout << "projective_space::plane_intersection_type_fast not independent, skip" << endl;
				cout << "subset: ";
				INT_vec_print(cout, subset, 3);
				cout << endl;
				}
			rank_idx[rk] = -1;
			continue;
			}
#if 0
		for (j = 0; j < 3 * d; j++) {
			G->M[j] = Basis[j];
			}
#endif
		G->rank_longinteger_here(Basis, plane_rk, 0 /* verbose_level */);
		if (f_v) {
			cout << rk << "-th subset ";
			INT_vec_print(cout, subset, 3);
			cout << " plane_rk=" << plane_rk << endl;
			}

		if (longinteger_vec_search(R, len, plane_rk, idx)) {
			//rank_idx[rk] = idx;
			// this case should never happen:
			cout << "projective_space::plane_intersection_type_fast longinteger_vec_search(R, len, plane_rk, idx) is TRUE" << endl;
			exit(1);
			}
		else {
			if (f_v3) {
				cout << "plane_rk=" << plane_rk << " was not found" << endl;
				}
			pts_on_plane = NEW_INT(set_size);
			//if (f_v3) {
				//cout << "after allocating pts_on_plane, plane_rk=" << plane_rk << endl;
				//}

			plane_rk.assign_to(aa);
			G->unrank_longinteger_here(Basis_save, aa, 0 /* verbose_level */);
#if 0
			for (j = 0; j < 3 * d; j++) {
				Basis_save[j] = G->M[j];
				}
#endif
			if (f_v3) {
				cout << "after unrank " << plane_rk << ", Basis:" << endl;
				INT_matrix_print(Basis_save, 3, d);
				}

#if 0
			{
				INT *test;
				INT r1;

				test = NEW_INT(6 * d);
				for (i = 0; i < 3 * d; i++) {
					test[i] = Basis[i];
					}
				for (i = 0; i < 3 * d; i++) {
					test[3 * d + i] = Basis_save[i];
					}
				r1 = F->rank_of_rectangular_matrix(test, 6, d, 0 /* verbose_level */);
				if (r1 != 3) {
					cout << "projective_space::plane_intersection_type_fast r1 != 3" << endl;
					cout << "r1=" << r1 << endl;
					cout << "Basis:" << endl;
					INT_matrix_print(Basis, 3, d);
					cout << "Basis_save:" << endl;
					INT_matrix_print(Basis_save, 3, d);
					exit(1);
					}
			}
#endif
					
			l = 0;
			for (h = 0; h < set_size; h++) {
				if (FALSE && f_v3) {
					cout << "testing point " << h << ":" << endl;
					cout << "plane_rk=" << plane_rk << endl;
					}
#if 0
				plane_rk.assign_to(aa);
				G->unrank_longinteger(aa, 0 /* verbose_level */);
				if (f_v3) {
					cout << "after unrank " << plane_rk << ":" << endl;
					INT_matrix_print(G->M, 3, d);
					}	
				for (j = 0; j < 3 * d; j++) {
					Basis[j] = G->M[j];
					}
#endif

				INT_vec_copy(Basis_save, Basis, 3 * d);
#if 0
				for (j = 0; j < 3 * d; j++) {
					Basis[j] = Basis_save[j];
					}
#endif
				//a = set[h];

				INT_vec_copy(Coords + h * d, Basis + 3 * d, d);
#if 0
				for (u = 0; u < d; u++) {
					Basis[3 * d + u] = Coords[h * d + u];
					}
#endif

				//unrank_point(Basis + 3 * d, set[h]);
				if (FALSE && f_v3) {
					cout << "Basis and point:" << endl;
					INT_matrix_print(Basis, 4, d);
					}	
				r = F->rank_of_rectangular_matrix(Basis, 4, d, 0 /* verbose_level */);
				if (r == 3) {
					pts_on_plane[l++] = h;
					if (f_v3) {
						cout << "point " << h << " is on the plane" << endl;
						}
					}
				else {
					if (FALSE && f_v3) {
						cout << "point " << h << " is not on the plane" << endl;
						}
					}
				}
			if (f_v) {
				cout << "We found an " << l << "-plane, its rank is " << plane_rk << endl;
				cout << "The ranks of points on that plane are : ";
				INT_vec_print(cout, pts_on_plane, l);
				cout << endl;
				}


			if (l >= 3) {
				for (j = len; j > idx; j--) {
					R[j].swap_with(R[j - 1]);
					Pts_on_plane[j] = Pts_on_plane[j - 1];
					nb_pts_on_plane[j] = nb_pts_on_plane[j - 1];
					}
				for (j = 0; j < N; j++) {
					if (f_subset_done[j] && rank_idx[j] >= idx) {
						rank_idx[j]++;
						}
					}
				plane_rk.assign_to(R[idx]);
				if (f_v3) {
					cout << "after assign_to, plane_rk=" << plane_rk << endl;
					}
				rank_idx[rk] = idx;
				len++;




				N2 = INT_n_choose_k(l, 3);
				for (i = 0; i < N2; i++) {
					unrank_k_subset(i, subset2, l, 3);
					for (h = 0; h < 3; h++) {
						subset3[h] = pts_on_plane[subset2[h]];
						}
					rr = rank_k_subset(subset3, set_size, 3);
					if (f_v) {
						cout << i << "-th subset3 ";
						INT_vec_print(cout, subset3, 3);
						cout << " rr=" << rr << endl;
						}
					if (!f_subset_done[rr]) {
						f_subset_done[rr] = TRUE;
						rank_idx[rr] = idx;
						}
					else if (rank_idx[rr] == -1) {
						rank_idx[rr] = idx;
						}
					else if (rank_idx[rr] != idx) {
						cout << "projective_space::plane_intersection_type_fast f_subset_done[rr] && rank_idx[rr] >= 0 && rank_idx[rr] != idx" << endl;
						exit(1);
						}
					}
				Pts_on_plane[idx] = pts_on_plane;
				nb_pts_on_plane[idx] = l;
				}
			else {
				// now l <= 2, we skip those planes:
				
				FREE_INT(pts_on_plane);
				f_subset_done[rk] = TRUE;
				rank_idx[rk] = -2;
				}
			} // else
		} // next rk
	
	FREE_INT(Basis);
	FREE_INT(Basis_save);
	FREE_INT(f_subset_done);
	FREE_INT(rank_idx);
	FREE_INT(Coords);
}


void projective_space::klein_correspondence(projective_space *P5, 
	INT *set_in, INT set_size, INT *set_out, INT verbose_level)
// Computes the Pluecker coordinates for a line in PG(3,q) in the following order:
// (x_1,x_2,x_3,x_4,x_5,x_6) = 
// (Pluecker_12, Pluecker_34, Pluecker_13, Pluecker_42, Pluecker_14, Pluecker_23)
// satisfying the quadratic form x_1x_2 + x_3x_4 + x_5x_6 = 0
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT d = n + 1;
	INT h;
	INT basis8[8];
	INT v6[6];
	INT *x4, *y4;
	INT a, b, c;
	INT f_elements_exponential = TRUE;
	const BYTE *symbol_for_print = "\\alpha";


	if (f_v) {
		cout << "projective_space::klein_correspondence" << endl;
		}
	for (h = 0; h < set_size; h++) {
		a = set_in[h];
		Grass_lines->unrank_INT(a, 0 /* verbose_level */);
		if (f_vv) {
			cout << setw(5) << h << " : " << setw(5) << a << " :" << endl;
			F->latex_matrix(cout, f_elements_exponential, 
				symbol_for_print, Grass_lines->M, 2, 4);
			cout << endl;
			}
		INT_vec_copy(Grass_lines->M, basis8, 8);
		if (f_vv) {
			INT_matrix_print(basis8, 2, 4);
			}
		x4 = basis8;
		y4 = basis8 + 4;
		v6[0] = F->Pluecker_12(x4, y4);
		v6[1] = F->Pluecker_34(x4, y4);
		v6[2] = F->Pluecker_13(x4, y4);
		v6[3] = F->Pluecker_42(x4, y4);
		v6[4] = F->Pluecker_14(x4, y4);
		v6[5] = F->Pluecker_23(x4, y4);
		if (f_vv) {
			cout << "v6 : ";
			INT_vec_print(cout, v6, 6);
			cout << endl;
			}
		a = F->mult(v6[0], v6[1]);
		b = F->mult(v6[2], v6[3]);
		c = F->mult(v6[4], v6[5]);
		d = F->add3(a, b, c);
		//cout << "a=" << a << " b=" << b << " c=" << c << endl;
		//cout << "d=" << d << endl;
		if (d) {
			cout << "d != 0" << endl;
			exit(1);
			}
		set_out[h] = P5->rank_point(v6);
		}
	if (f_v) {
		cout << "projective_space::klein_correspondence done" << endl;
		}
	
}

void projective_space::Pluecker_coordinates(INT line_rk, INT *v6, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT basis8[8];
	INT *x4, *y4;
	INT f_elements_exponential = FALSE;
	const BYTE *symbol_for_print = "\\alpha";
	
	if (f_v) {
		cout << "projective_space::Pluecker_coordinates" << endl;
		}
	Grass_lines->unrank_INT(line_rk, 0 /* verbose_level */);
	if (f_vv) {
		cout << setw(5) << line_rk << " :" << endl;
		F->latex_matrix(cout, f_elements_exponential, 
			symbol_for_print, Grass_lines->M, 2, 4);
		cout << endl;
		}
	INT_vec_copy(Grass_lines->M, basis8, 8);
	if (f_vv) {
		INT_matrix_print(basis8, 2, 4);
		}
	x4 = basis8;
	y4 = basis8 + 4;
	v6[0] = F->Pluecker_12(x4, y4);
	v6[1] = F->Pluecker_34(x4, y4);
	v6[2] = F->Pluecker_13(x4, y4);
	v6[3] = F->Pluecker_42(x4, y4);
	v6[4] = F->Pluecker_14(x4, y4);
	v6[5] = F->Pluecker_23(x4, y4);
	if (f_vv) {
		cout << "v6 : ";
		INT_vec_print(cout, v6, 6);
		cout << endl;
		}
	if (f_v) {
		cout << "projective_space::Pluecker_coordinates done" << endl;
		}
}

void projective_space::klein_correspondence_special_model(projective_space *P5, 
	INT *table, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT d = n + 1;
	INT h;
	INT basis8[8];
	INT x6[6];
	INT y6[6];
	INT *x4, *y4;
	INT a, b, c;
	INT half;
	INT f_elements_exponential = TRUE;
	const BYTE *symbol_for_print = "\\alpha";
	//INT *table;

	if (f_v) {
		cout << "projective_space::klein_correspondence" << endl;
		}
	half = F->inverse(F->add(1, 1));
	if (f_v) {
		cout << "half=" << half << endl;
		cout << "N_lines=" << N_lines << endl;
		}
	//table = NEW_INT(N_lines);
	for (h = 0; h < N_lines; h++) {
		Grass_lines->unrank_INT(h, 0 /* verbose_level */);
		if (f_vv) {
			cout << setw(5) << h << " :" << endl;
			F->latex_matrix(cout, f_elements_exponential, 
				symbol_for_print, Grass_lines->M, 2, 4);
			cout << endl;
			}
		INT_vec_copy(Grass_lines->M, basis8, 8);
#if 0
		for (i = 0; i < 8; i++) {
			basis8[i] = Grass_lines->M[i];
			}
#endif
		if (f_vv) {
			INT_matrix_print(basis8, 2, 4);
			}
		x4 = basis8;
		y4 = basis8 + 4;
		x6[0] = F->Pluecker_12(x4, y4);
		x6[1] = F->Pluecker_34(x4, y4);
		x6[2] = F->Pluecker_13(x4, y4);
		x6[3] = F->Pluecker_42(x4, y4);
		x6[4] = F->Pluecker_14(x4, y4);
		x6[5] = F->Pluecker_23(x4, y4);
		if (f_vv) {
			cout << "x6 : ";
			INT_vec_print(cout, x6, 6);
			cout << endl;
			}
		a = F->mult(x6[0], x6[1]);
		b = F->mult(x6[2], x6[3]);
		c = F->mult(x6[4], x6[5]);
		d = F->add3(a, b, c);
		//cout << "a=" << a << " b=" << b << " c=" << c << endl;
		//cout << "d=" << d << endl;
		if (d) {
			cout << "d != 0" << endl;
			exit(1);
			}
		y6[0] = F->negate(x6[0]);
		y6[1] = x6[1];
		y6[2] = F->mult(half, F->add(x6[2], x6[3]));
		y6[3] = F->mult(half, F->add(x6[2], F->negate(x6[3])));
		y6[4] = x6[4];
		y6[5] = x6[5];
		if (f_vv) {
			cout << "y6 : ";
			INT_vec_print(cout, y6, 6);
			cout << endl;
			}
		table[h] = P5->rank_point(y6);
		}

	cout << "lines in PG(3,q) to points in PG(5,q) in special model:" << endl;
	for (h = 0; h < N_lines; h++) {
		cout << setw(4) << h << " : " << setw(5) << table[h] << endl;
		}
	
	//FREE_INT(table);
	if (f_v) {
		cout << "projective_space::klein_correspondence_special_model done" << endl;
		}
	
}

void projective_space::cheat_sheet_points_and_lines(ostream &f, INT verbose_level)
{

	cheat_sheet_points(f, verbose_level);

}

void projective_space::cheat_sheet_points(ostream &f, INT verbose_level)
{
	INT i, d;
	INT *v;

	d = n + 1;

	v = NEW_INT(d);

	f << "PG$(" << n << ", " << q << ")$ has " << N_points << " points:\\\\" << endl;
	if (F->e == 1) {
		for (i = 0; i < N_points; i++) {
			PG_element_unrank_modified(*F, v, 1, d, i);
			f << "$P_{" << i << "}=";
			INT_vec_print(f, v, d);
			f << "$\\\\" << endl;
			}
		}
	else {
		for (i = 0; i < N_points; i++) {
			PG_element_unrank_modified(*F, v, 1, d, i);
			f << "$P_{" << i << "}=";
			INT_vec_print(f, v, d);
			f << "=";
			F->INT_vec_print_elements_exponential(f, v, d, "\\alpha");
			f << "$\\\\" << endl;
			}
		}
	f << "Normalized from the left:\\\\" << endl;
		for (i = 0; i < N_points; i++) {
			PG_element_unrank_modified(*F, v, 1, d, i);
			PG_element_normalize_from_front(*F, v, 1, d);
			f << "$P_{" << i << "}=";
			INT_vec_print(f, v, d);
			f << "$\\\\" << endl;
			}
	FREE_INT(v);
}

void projective_space::cheat_sheet_subspaces(ostream &f, INT k, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	grassmann *Gr;
	INT *v;
	INT n1, k1;
	INT nb_points;
	INT nb_k_subspaces;
	INT i, j, u;


	if (f_v) {
		cout << "projective_space::cheat_sheet_subspaces k=" << k << endl;
		}
	n1 = n + 1;
	k1 = k + 1;
	v = NEW_INT(n1);



	Gr = new grassmann;
	Gr->init(n1, k1, F, 0 /*verbose_level*/);


	nb_points = N_points;
	nb_k_subspaces = generalized_binomial(n1, k1, q);


	f << "PG$(" << n << ", " << q << ")$ has " << nb_k_subspaces << " $" << k << "$-subspaces:\\\\" << endl;
	for (u = 0; u < nb_k_subspaces; u++) {
		Gr->unrank_INT(u, 0 /* verbose_level*/);
		f << "$L_{" << u << "}=";
		f << "\\left[" << endl;
		f << "\\begin{array}{c}" << endl;
		for (i = 0; i < k1; i++) {
			for (j = 0; j < n1; j++) {
				f << Gr->M[i * n1 + j];
				}
			f << "\\\\" << endl;
			}
		f << "\\end{array}" << endl;
		f << "\\right]" << endl;

		if (n == 3 && k == 1) {
			INT v6[6];

			Pluecker_coordinates(u, v6, 0 /* verbose_level */);
			f << "Pl=(" << v6[0] << "," << v6[1] << "," << v6[2] << "," << v6[3] << "," << v6[4] << "," << v6[5] << " ";
			f << ")" << endl;

			}
		f << "$\\\\" << endl;
		}

	delete Gr;
	FREE_INT(v);

	if (f_v) {
		cout << "projective_space::cheat_sheet_subspaces done" << endl;
		}
}

void projective_space::cheat_sheet_line_intersection(ostream &f, INT verbose_level)
{
	INT i, j, a;


	f << "intersection of 2 lines:" << endl;
	f << "$$" << endl;
	f << "\\begin{array}{|r|*{" << N_points << "}{r}|}" << endl;
	f << "\\hline" << endl;
	for (j = 0; j < N_points; j++) {
		f << "& " << j << endl;
		}
	f << "\\\\" << endl;
	f << "\\hline" << endl;
	for (i = 0; i < N_points; i++) {
		f << i;
		for (j = 0; j < N_points; j++) {
			a = Line_intersection[i * N_lines + j];
			f << " & ";
			if (i != j)
				f << a;
			}
		f << "\\\\[-3pt]" << endl;
		}
	f << "\\hline" << endl;
	f << "\\end{array}" << endl;
	f << "$$" << endl;

}

void projective_space::cheat_sheet_line_through_pairs_of_points(ostream &f, INT verbose_level)
{
	INT i, j, a;



	f << "line through 2 points:" << endl;
	f << "$$" << endl;
	f << "\\begin{array}{|r|*{" << N_points << "}{r}|}" << endl;
	f << "\\hline" << endl;
	for (j = 0; j < N_points; j++) {
		f << "& " << j << endl;
		}
	f << "\\\\" << endl;
	f << "\\hline" << endl;
	for (i = 0; i < N_points; i++) {
		f << i;
		for (j = 0; j < N_points; j++) {

			a = Line_through_two_points[i * N_points + j];
			f << " & ";
			if (i != j)
				f << a;
			}
		f << "\\\\[-3pt]" << endl;
		}
	f << "\\hline" << endl;
	f << "\\end{array}" << endl;
	f << "$$" << endl;

}

void projective_space::conic_type_randomized(INT nb_times, 
	INT *set, INT set_size, 
	INT **&Pts_on_conic, INT *&nb_pts_on_conic, INT &len, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_v3 = (verbose_level >= 3);
	INT rk, h, i, j, a, d, N, l, cnt;

	INT input_pts[5];
	INT six_coeffs[6];
	INT vec[3];

	INT subset[5];
	longinteger_object conic_rk, aa;
	INT *pts_on_conic;
	INT allocation_length;

	if (f_v) {
		cout << "projective_space::conic_type_randomized" << endl;
		}
	if (n != 2) {
		cout << "projective_space::conic_type_randomized n != 2" << endl;
		exit(1);
		}
	if (f_vv) {
		print_set_numerical(set, set_size);
		}

	if (!test_if_set_with_return_value(set, set_size)) {
		cout << "projective_space::conic_type_randomized the input set if not a set" << endl;
		exit(1);
		}
	d = n + 1;
	N = INT_n_choose_k(set_size, 5);

	if (f_v) {
		cout << "set_size=" << set_size << endl;
		cout << "N=number of 5-subsets of the set=" << N << endl;
		}
	
	// allocate data that is returned:
	allocation_length = 1024;
	//R = new longinteger_object[allocation_length];
	Pts_on_conic = NEW_PINT(allocation_length);
	nb_pts_on_conic = NEW_INT(allocation_length);


	len = 0;
	for (cnt = 0; cnt < nb_times; cnt++) {

		rk = random_integer(N);
		unrank_k_subset(rk, subset, set_size, 5);
		if (cnt && ((cnt % 1000) == 0)) {
			cout << cnt << " / " << nb_times << " : ";
			INT_vec_print(cout, subset, 5);
			cout << endl;
			}

		for (i = 0; i < len; i++) {
			if (INT_vec_is_subset_of(subset, 5, Pts_on_conic[i], nb_pts_on_conic[i])) {

#if 0
				cout << "The set ";
				INT_vec_print(cout, subset, 5);
				cout << " is a subset of the " << i << "th conic ";
				INT_vec_print(cout, Pts_on_conic[i], nb_pts_on_conic[i]);
				cout << endl;
#endif

				break;
				}
			}
		if (i < len) {
			continue;
			}
		for (j = 0; j < 5; j++) {
			a = subset[j];
			input_pts[j] = set[a];
			}
		if (FALSE /* f_v3 */) {
			cout << "subset: ";
			INT_vec_print(cout, subset, 5);
			cout << "input_pts: ";
			INT_vec_print(cout, input_pts, 5);
			}

		if (!determine_conic_in_plane(input_pts, 5, six_coeffs, 0 /* verbose_level */)) {
			continue;
			}


		PG_element_normalize(*F, six_coeffs, 1, 6);
		AG_element_rank_longinteger(F->q, six_coeffs, 1, 6, conic_rk);
		if (FALSE /* f_vv */) {
			cout << rk << "-th subset ";
			INT_vec_print(cout, subset, 5);
			cout << " conic_rk=" << conic_rk << endl;
			}

		if (FALSE /* longinteger_vec_search(R, len, conic_rk, idx) */) {

#if 0
			cout << "projective_space::conic_type_randomized longinteger_vec_search(R, len, conic_rk, idx) is TRUE" << endl;
			cout << "The current set is ";
			INT_vec_print(cout, subset, 5);
			cout << endl;
			cout << "conic_rk=" << conic_rk << endl;
			cout << "The set where it should be is ";
			INT_vec_print(cout, Pts_on_conic[idx], nb_pts_on_conic[idx]);
			cout << endl;
			cout << "R[idx]=" << R[idx] << endl;
			cout << "This is the " << idx << "th conic" << endl;
			exit(1);
#endif

			}
		else {
			if (f_v3) {
				cout << "conic_rk=" << conic_rk << " was not found" << endl;
				}
			pts_on_conic = NEW_INT(set_size);
			l = 0;
			for (h = 0; h < set_size; h++) {
				if (FALSE && f_v3) {
					cout << "testing point " << h << ":" << endl;
					cout << "conic_rk=" << conic_rk << endl;
					}
				
				unrank_point(vec, set[h]);
				a = F->evaluate_conic_form(six_coeffs, vec);

				
				if (a == 0) {
					pts_on_conic[l++] = h;
					if (f_v3) {
						cout << "point " << h << " is on the conic" << endl;
						}
					}
				else {
					if (FALSE && f_v3) {
						cout << "point " << h << " is not on the conic" << endl;
						}
					}
				}
			if (FALSE /*f_v*/) {
				cout << "We found an " << l << "-conic, its rank is " << conic_rk << endl;

				
				}


			if (l >= 8) {

				if (f_v) {
					cout << "We found an " << l << "-conic, its rank is " << conic_rk << endl;
					cout << "The " << l << " points on the " << len << "th conic are: ";
					INT_vec_print(cout, pts_on_conic, l);
					cout << endl;


				
					}


#if 0
				for (j = len; j > idx; j--) {
					R[j].swap_with(R[j - 1]);
					Pts_on_conic[j] = Pts_on_conic[j - 1];
					nb_pts_on_conic[j] = nb_pts_on_conic[j - 1];
					}
				conic_rk.assign_to(R[idx]);
				Pts_on_conic[idx] = pts_on_conic;
				nb_pts_on_conic[idx] = l;
#else

				//conic_rk.assign_to(R[len]);
				Pts_on_conic[len] = pts_on_conic;
				nb_pts_on_conic[len] = l;

#endif


				len++;
				if (f_v) {
					cout << "We now have found " << len << " conics" << endl;


					classify C;
					INT f_second = FALSE;

					C.init(nb_pts_on_conic, len, f_second, 0);

					if (f_v) {
						cout << "The conic intersection type is (";
						C.print_naked(FALSE /*f_backwards*/);
						cout << ")" << endl << endl;
						}



					}

				if (len == allocation_length) {
					INT new_allocation_length = allocation_length + 1024;


					//longinteger_object *R1;
					INT **Pts_on_conic1;
					INT *nb_pts_on_conic1;
					
					//R1 = new longinteger_object[new_allocation_length];
					Pts_on_conic1 = NEW_PINT(new_allocation_length);
					nb_pts_on_conic1 = NEW_INT(new_allocation_length);
					for (i = 0; i < len; i++) {
						//R1[i] = R[i];
						Pts_on_conic1[i] = Pts_on_conic[i];
						nb_pts_on_conic1[i] = nb_pts_on_conic[i];
						}
					//delete [] R;
					FREE_PINT(Pts_on_conic);
					FREE_INT(nb_pts_on_conic);
					//R = R1;
					Pts_on_conic = Pts_on_conic1;
					nb_pts_on_conic = nb_pts_on_conic1;
					allocation_length = new_allocation_length;
					} 




				}
			else {
				// we skip this conic:
				
				FREE_INT(pts_on_conic);
				}
			} // else
		} // next rk
	
}

void projective_space::conic_intersection_type(INT f_randomized, INT nb_times, 
	INT *set, INT set_size, 
	INT *&intersection_type, INT &highest_intersection_number, 
	INT f_save_largest_sets, set_of_sets *&largest_sets, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//longinteger_object *R;
	INT **Pts_on_conic;
	INT *nb_pts_on_conic;
	INT nb_conics;
	INT i, j, idx, f, l, a, t;

	if (f_v) {
		cout << "projective_space::conic_intersection_type" << endl;
		}

	if (f_randomized) {
		if (f_v) {
			cout << "projective_space::conic_intersection_type randomized" << endl;
			}
		conic_type_randomized(nb_times, 
			set, set_size, 
			Pts_on_conic, nb_pts_on_conic, nb_conics, 
			verbose_level - 1);
		}
	else {
		if (f_v) {
			cout << "projective_space::conic_intersection_type not randomized" << endl;
			}
		conic_type(
			set, set_size, 
			Pts_on_conic, nb_pts_on_conic, nb_conics, 
			verbose_level - 1);
		}

	classify C;
	INT f_second = FALSE;

	C.init(nb_pts_on_conic, nb_conics, f_second, 0);
	if (f_v) {
		cout << "projective_space::conic_intersection_type conic-intersection type: ";
		C.print(FALSE /*f_backwards*/);
		}

	if (f_v) {
		cout << "The conic intersection type is (";
		C.print_naked(FALSE /*f_backwards*/);
		cout << ")" << endl << endl;
		}

	f = C.type_first[C.nb_types - 1];
	highest_intersection_number = C.data_sorted[f];
	intersection_type = NEW_INT(highest_intersection_number + 1);
	INT_vec_zero(intersection_type, highest_intersection_number + 1);
	for (i = 0; i < C.nb_types; i++) {
		f = C.type_first[i];
		l = C.type_len[i];
		a = C.data_sorted[f];
		intersection_type[a] = l;
		}

	if (f_save_largest_sets) {
		largest_sets = new set_of_sets;
		t = C.nb_types - 1;
		f = C.type_first[t];
		l = C.type_len[t];
		largest_sets->init_basic_constant_size(set_size, l, highest_intersection_number, verbose_level);
		for (j = 0; j < l; j++) {
			idx = C.sorting_perm_inv[f + j];
			INT_vec_copy(Pts_on_conic[idx], largest_sets->Sets[j], highest_intersection_number);
			}
		}

	for (i = 0; i < nb_conics; i++) {
		FREE_INT(Pts_on_conic[i]);
		}
	FREE_PINT(Pts_on_conic);
	FREE_INT(nb_pts_on_conic);
	//delete [] R;

}

void projective_space::conic_type(
	INT *set, INT set_size, 
	INT **&Pts_on_conic, INT *&nb_pts_on_conic, INT &len, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_v3 = (verbose_level >= 3);
	INT rk, h, i, j, a, d, N, l;

	INT input_pts[5];
	INT six_coeffs[6];
	INT vec[3];

	INT subset[5];
	longinteger_object conic_rk, aa;
	INT *pts_on_conic;
	INT allocation_length;

	if (f_v) {
		cout << "projective_space::conic_type" << endl;
		}
	if (n != 2) {
		cout << "projective_space::conic_type n != 2" << endl;
		exit(1);
		}
	if (f_vv) {
		print_set_numerical(set, set_size);
		}

	if (!test_if_set_with_return_value(set, set_size)) {
		cout << "projective_space::conic_type the input set if not a set" << endl;
		exit(1);
		}
	d = n + 1;
	N = INT_n_choose_k(set_size, 5);

	if (f_v) {
		cout << "set_size=" << set_size << endl;
		cout << "N=number of 5-subsets of the set=" << N << endl;
		}
	
	// allocate data that is returned:
	allocation_length = 1024;
	//R = new longinteger_object[allocation_length];
	Pts_on_conic = NEW_PINT(allocation_length);
	nb_pts_on_conic = NEW_INT(allocation_length);


	len = 0;
	for (rk = 0; rk < N; rk++) {

		unrank_k_subset(rk, subset, set_size, 5);
		if (f_v3 || (rk && ((rk % 1000) == 0))) {
			cout << rk << " / " << N << " : ";
			INT_vec_print(cout, subset, 5);
			cout << endl;
			}

		for (i = 0; i < len; i++) {
			if (INT_vec_is_subset_of(subset, 5, Pts_on_conic[i], nb_pts_on_conic[i])) {

#if 0
				cout << "The set ";
				INT_vec_print(cout, subset, 5);
				cout << " is a subset of the " << i << "th conic ";
				INT_vec_print(cout, Pts_on_conic[i], nb_pts_on_conic[i]);
				cout << endl;
#endif

				break;
				}
			}
		if (i < len) {
			continue;
			}
		for (j = 0; j < 5; j++) {
			a = subset[j];
			input_pts[j] = set[a];
			}
		if (f_v3) {
			cout << "subset: ";
			INT_vec_print(cout, subset, 5);
			cout << "input_pts: ";
			INT_vec_print(cout, input_pts, 5);
			}

		if (!determine_conic_in_plane(input_pts, 5, six_coeffs, verbose_level - 2)) {
			continue;
			}


		PG_element_normalize(*F, six_coeffs, 1, 6);
		AG_element_rank_longinteger(F->q, six_coeffs, 1, 6, conic_rk);
		if (FALSE /* f_vv */) {
			cout << rk << "-th subset ";
			INT_vec_print(cout, subset, 5);
			cout << " conic_rk=" << conic_rk << endl;
			}

		if (FALSE /* longinteger_vec_search(R, len, conic_rk, idx) */) {

#if 0
			cout << "projective_space::conic_type_randomized longinteger_vec_search(R, len, conic_rk, idx) is TRUE" << endl;
			cout << "The current set is ";
			INT_vec_print(cout, subset, 5);
			cout << endl;
			cout << "conic_rk=" << conic_rk << endl;
			cout << "The set where it should be is ";
			INT_vec_print(cout, Pts_on_conic[idx], nb_pts_on_conic[idx]);
			cout << endl;
			cout << "R[idx]=" << R[idx] << endl;
			cout << "This is the " << idx << "th conic" << endl;
			exit(1);
#endif

			}
		else {
			if (f_v3) {
				cout << "conic_rk=" << conic_rk << " was not found" << endl;
				}
			pts_on_conic = NEW_INT(set_size);
			l = 0;
			for (h = 0; h < set_size; h++) {
				if (FALSE && f_v3) {
					cout << "testing point " << h << ":" << endl;
					cout << "conic_rk=" << conic_rk << endl;
					}
				
				unrank_point(vec, set[h]);
				a = F->evaluate_conic_form(six_coeffs, vec);

				
				if (a == 0) {
					pts_on_conic[l++] = h;
					if (f_v3) {
						cout << "point " << h << " is on the conic" << endl;
						}
					}
				else {
					if (FALSE && f_v3) {
						cout << "point " << h << " is not on the conic" << endl;
						}
					}
				}
			if (FALSE /*f_v*/) {
				cout << "We found an " << l << "-conic, its rank is " << conic_rk << endl;

				
				}


			if (l >= 8) {

				if (f_v) {
					cout << "We found an " << l << "-conic, its rank is " << conic_rk << endl;
					cout << "The " << l << " points on the " << len << "th conic are: ";
					INT_vec_print(cout, pts_on_conic, l);
					cout << endl;


				
					}


#if 0
				for (j = len; j > idx; j--) {
					R[j].swap_with(R[j - 1]);
					Pts_on_conic[j] = Pts_on_conic[j - 1];
					nb_pts_on_conic[j] = nb_pts_on_conic[j - 1];
					}
				conic_rk.assign_to(R[idx]);
				Pts_on_conic[idx] = pts_on_conic;
				nb_pts_on_conic[idx] = l;
#else

				//conic_rk.assign_to(R[len]);
				Pts_on_conic[len] = pts_on_conic;
				nb_pts_on_conic[len] = l;

#endif


				len++;
				if (f_v) {
					cout << "We now have found " << len << " conics" << endl;


					classify C;
					INT f_second = FALSE;

					C.init(nb_pts_on_conic, len, f_second, 0);

					if (f_v) {
						cout << "The conic intersection type is (";
						C.print_naked(FALSE /*f_backwards*/);
						cout << ")" << endl << endl;
						}



					}

				if (len == allocation_length) {
					INT new_allocation_length = allocation_length + 1024;


					//longinteger_object *R1;
					INT **Pts_on_conic1;
					INT *nb_pts_on_conic1;
					
					//R1 = new longinteger_object[new_allocation_length];
					Pts_on_conic1 = NEW_PINT(new_allocation_length);
					nb_pts_on_conic1 = NEW_INT(new_allocation_length);
					for (i = 0; i < len; i++) {
						//R1[i] = R[i];
						Pts_on_conic1[i] = Pts_on_conic[i];
						nb_pts_on_conic1[i] = nb_pts_on_conic[i];
						}
					//delete [] R;
					FREE_PINT(Pts_on_conic);
					FREE_INT(nb_pts_on_conic);
					//R = R1;
					Pts_on_conic = Pts_on_conic1;
					nb_pts_on_conic = nb_pts_on_conic1;
					allocation_length = new_allocation_length;
					} 




				}
			else {
				// we skip this conic:
				
				FREE_INT(pts_on_conic);
				}
			} // else
		} // next rk
	
}

void projective_space::find_nucleus(INT *set, INT set_size, INT &nucleus, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j, a, b, l, sz, idx, t1, t2;
	INT *Lines;

	if (f_v) {
		cout << "projective_space::find_nucleus" << endl;
		}

	if (n != 2) {
		cout << "projective_space::find_nucleus n != 2" << endl;
		exit(1);
		}
	if (set_size != F->q + 1) {
		cout << "projective_space::find_nucleus set_size != F->q + 1" << endl;
		exit(1);
		}

	if (Lines_on_point == NULL) {
		init_incidence_structure(verbose_level);
		}

	Lines = NEW_INT(r);
	a = set[0];
	for (i = 0; i < r; i++) {
		Lines[i] = Lines_on_point[a * r + i];
		}
	sz = r;
	INT_vec_heapsort(Lines, r);

	for (i = 0; i < set_size - 1; i++) {
		b = set[1 + i];
		l = line_through_two_points(a, b);
		if (!INT_vec_search(Lines, sz, l, idx)) {
			cout << "projective_space::find_nucleus cannot find secant in pencil" << endl;	
			exit(1);
			}
		for (j = idx + 1; j < sz; j++) {
			Lines[j - 1] = Lines[j];
			}
		sz--;
		}
	if (sz != 1) {
		cout << "projective_space::find_nucleus sz != 1" << endl;	
		exit(1);
		}
	t1 = Lines[0];
	if (f_v) {
		cout << "projective_space::find_nucleus t1 = " << t1 << endl;
		}
	


	a = set[1];
	for (i = 0; i < r; i++) {
		Lines[i] = Lines_on_point[a * r + i];
		}
	sz = r;
	INT_vec_heapsort(Lines, r);

	for (i = 0; i < set_size - 1; i++) {
		if (i == 0) {
			b = set[0];
			}
		else {
			b = set[1 + i];
			}
		l = line_through_two_points(a, b);
		if (!INT_vec_search(Lines, sz, l, idx)) {
			cout << "projective_space::find_nucleus cannot find secant in pencil" << endl;	
			exit(1);
			}
		for (j = idx + 1; j < sz; j++) {
			Lines[j - 1] = Lines[j];
			}
		sz--;
		}
	if (sz != 1) {
		cout << "projective_space::find_nucleus sz != 1" << endl;	
		exit(1);
		}
	t2 = Lines[0];
	if (f_v) {
		cout << "projective_space::find_nucleus t2 = " << t2 << endl;
		}
	
	nucleus = line_intersection(t1, t2);
	if (f_v) {
		cout << "projective_space::find_nucleus nucleus = " << nucleus << endl;
		INT v[3];
		unrank_point(v, nucleus);
		cout << "nucleus = ";
		INT_vec_print(cout, v, 3);
		cout << endl;
		}



	if (f_v) {
		cout << "projective_space::find_nucleus done" << endl;
		}
}


