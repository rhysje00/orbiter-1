// strong_generators.C
//
// Anton Betten
// December 4, 2013

#include "galois.h"
#include "action.h"

strong_generators::strong_generators()
{
	null();
}

strong_generators::~strong_generators()
{
	freeself();
}

void strong_generators::null()
{
	A = NULL;
	tl = NULL;
	gens = NULL;
}

void strong_generators::freeself()
{
	if (tl) {
		FREE_INT(tl);
		}
	if (gens) {
		delete gens;
		}
	null();
}

void strong_generators::init(action *A, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "strong_generators::init" << endl;
		}
	strong_generators::A = A;
	if (f_v) {
		cout << "strong_generators::init done" << endl;
		}
}

void strong_generators::init_from_sims(sims *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);

	if (f_v) {
		cout << "strong_generators::init_from_sims" << endl;
		}
	A = S->A;
	tl = NEW_INT(A->base_len);
	gens = new vector_ge;
	S->extract_strong_generators_in_order(*gens, tl, 0 /*verbose_level*/);
	if (f_v) {
		cout << "strong_generators::init_from_sims done" << endl;
		}
}

void strong_generators::init_from_ascii_coding(action *A, BYTE *ascii_coding, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	longinteger_object go;
	group *G;

	if (f_v) {
		cout << "strong_generators::init_from_ascii_coding" << endl;
		}
	G = NEW_OBJECT(group);
	G->init(A);
	if (f_vv) {
		cout << "strong_generators::init_from_ascii_coding before G->init_ascii_coding_to_sims" << endl;
		}
	G->init_ascii_coding_to_sims(ascii_coding);
	if (f_vv) {
		cout << "strong_generators::init_from_ascii_coding after G->init_ascii_coding_to_sims" << endl;
		}
		

	G->S->group_order(go);

	if (f_vv) {
		cout << "strong_generators::init_from_ascii_coding Group order=" << go << endl;
		}

	init_from_sims(G->S, 0 /* verbose_level */);

	FREE_OBJECT(G);
	if (f_v) {
		cout << "strong_generators::init_from_ascii_coding done" << endl;
		}
}



void strong_generators::init_copy(strong_generators *S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT i;

	if (f_v) {
		cout << "strong_generators::init_copy" << endl;
		}
	A = S->A;
	tl = NEW_INT(A->base_len);
	//cout << "strong_generators::init_copy before INT_vec_copy" << endl;
	INT_vec_copy(S->tl, tl, A->base_len);
	gens = new vector_ge;
	gens->init(A);
	gens->allocate(S->gens->len);
	for (i = 0; i < S->gens->len; i++) {
		//cout << "strong_generators::init_copy before element_move i=" << i << endl;
		A->element_move(S->gens->ith(i), gens->ith(i), 0);
		}
	if (f_v) {
		cout << "strong_generators::init_copy done" << endl;
		}
}

void strong_generators::init_by_hdl(action *A, INT *gen_hdl, INT nb_gen, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "strong_generators::init_by_hdl" << endl;
		}
	init(A, 0);
	tl = NEW_INT(A->base_len);
	for (i = 0; i < A->base_len; i++) {
		tl[i] = 1;
		}
	gens = new vector_ge;
	gens->init(A);
	gens->allocate(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		A->element_retrieve(gen_hdl[i], gens->ith(i), 0);
		}
	if (f_v) {
		cout << "strong_generators::init_by_hdl done" << endl;
		}
}

void strong_generators::init_from_data(action *A, INT *data, 
	INT nb_elements, INT elt_size, INT *transversal_length, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "strong_generators::init_from_data" << endl;
		}
	init(A, verbose_level - 2);
	gens = new vector_ge;

	gens->init_from_data(A, data, 
		nb_elements, elt_size, verbose_level);
	
	tl = NEW_INT(A->base_len);
	INT_vec_copy(transversal_length, tl, A->base_len);

	if (f_v) {
		cout << "strong_generators::init_from_data done" << endl;
		}
}

sims *strong_generators::create_sims(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	sims *S;


	if (f_v) {
		cout << "strong_generators::create_sims" << endl;
		}
	
	if (gens == NULL) {
		cout << "strong_generators::create_sims gens == NULL" << endl;
		exit(1);
		}
	S = create_sims_from_generators_with_target_group_order_factorized(A, 
		gens, tl, A->base_len, 0 /* verbose_level */);

	if (f_v) {
		cout << "strong_generators::create_sims done" << endl;
		}
	return S;
}

sims *strong_generators::create_sims_in_different_action(action *A_given, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	sims *S;


	if (f_v) {
		cout << "strong_generators::create_sims_in_different_action" << endl;
		}
	
	if (gens == NULL) {
		cout << "strong_generators::create_sims_in_different_action gens == NULL" << endl;
		exit(1);
		}
	S = create_sims_from_generators_with_target_group_order_factorized(A_given, 
		gens, tl, A->base_len, 0 /* verbose_level */);

	if (f_v) {
		cout << "strong_generators::create_sims_in_different_action done" << endl;
		}
	return S;
}

void strong_generators::add_generators(vector_ge *coset_reps, INT group_index, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	sims *S;
	vector_ge *gens1; 
	INT *tl1;
	INT *coset_reps_vec;
	INT i;

	if (f_v) {
		cout << "strong_generators::add_generators" << endl;
		}
	if (f_vv) {
		cout << "group_index=" << group_index << endl;
		cout << "action=";
		A->print_info();
		}

	coset_reps_vec = NEW_INT(group_index * A->elt_size_in_INT);
	for (i = 0; i < group_index; i++) {
		A->element_move(coset_reps->ith(i), coset_reps_vec + i * A->elt_size_in_INT, 0);
		}

	gens1 = new vector_ge;
	tl1 = NEW_INT(A->base_len);
	
	S = create_sims(verbose_level - 1);

	S->transitive_extension_using_coset_representatives_extract_generators(
		coset_reps_vec, group_index, 
		*gens1, tl1, 
		verbose_level - 2);

	FREE_INT(coset_reps_vec);

	if (gens) {
		delete gens;
		}
	if (tl) {
		FREE_INT(tl);
		}
	gens = gens1;
	tl = tl1;

	delete S;
	if (f_v) {
		cout << "strong_generators::add_generators done" << endl;
		}
}

void strong_generators::add_single_generator(INT *Elt, INT group_index, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	sims *S;
	vector_ge *gens1; 
	INT *tl1;

	if (f_v) {
		cout << "strong_generators::add_single_generator" << endl;
		cout << "action=";
		A->print_info();
		}

	gens1 = new vector_ge;
	tl1 = NEW_INT(A->base_len);
	
	S = create_sims(verbose_level - 1);

	S->transitive_extension_using_generators(
		Elt, 1, group_index, 
		*gens1, tl1, 
		verbose_level);

	if (gens) {
		delete gens;
		}
	if (tl) {
		FREE_INT(tl);
		}
	gens = gens1;
	tl = tl1;

	delete S;
	if (f_v) {
		cout << "strong_generators::add_single_generator done" << endl;
		}
}

void strong_generators::group_order(longinteger_object &go)
{
	longinteger_domain D;

	D.multiply_up(go, tl, A->base_len);
}

INT strong_generators::group_order_as_INT()
{
	longinteger_domain D;
	longinteger_object go;

	D.multiply_up(go, tl, A->base_len);
	return go.as_INT();
}

void strong_generators::print_generators()
{
	INT i;
	longinteger_object go;

	group_order(go);
	cout << "Strong generators for a group of order " << go << " tl=";
	INT_vec_print(cout, tl, A->base_len);
	cout << endl;
	for (i = 0; i < gens->len; i++) {
		cout << "Generator " << i << " / " << gens->len << " is:" << endl;
		A->element_print(gens->ith(i), cout);
		}
}

void strong_generators::print_generators_tex()
{
	INT i;
	longinteger_object go;

	group_order(go);
	cout << "Strong generators for a group of order " << go << ":" << endl;
	cout << "$$" << endl;
	for (i = 0; i < gens->len; i++) {
		//cout << "Generator " << i << " / " << gens->len << " is:" << endl;
		A->element_print_latex(gens->ith(i), cout);
		}
	cout << "$$" << endl;
}

void strong_generators::print_generators_as_permutations()
{
	INT i;
	longinteger_object go;

	group_order(go);
	cout << "Strong generators for a group of order " << go << ":" << endl;
	for (i = 0; i < gens->len; i++) {
		cout << "Generator " << i << " / " << gens->len << " is:" << endl;
		A->element_print(gens->ith(i), cout);
		A->element_print_as_permutation(gens->ith(i), cout);
		}
}

void strong_generators::print_with_given_action(ostream &ost, action *A2)
{
	INT i;
	
	for (i = 0; i < gens->len; i++) {
		cout << "Generator " << i << " / " << gens->len << " is:" << endl;
		A2->element_print(gens->ith(i), cout);
		A2->element_print_as_permutation(gens->ith(i), cout);
		}
}

void strong_generators::compute_schreier_with_given_action(action *A_given, schreier *&Sch, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "strong_generators::compute_schreier_with_given_action" << endl;
		cout << "action=";
		A->print_info();
		}
	Sch = new schreier;

	Sch->init(A_given);
	Sch->initialize_tables();
	Sch->init_generators(*gens);
	Sch->compute_all_point_orbits(verbose_level - 2);


	if (f_v) {
		cout << "strong_generators::compute_schreier_with_given_action done, we found " << Sch->nb_orbits << " orbits" << endl;
		}
}

void strong_generators::compute_schreier_with_given_action_on_a_given_set(action *A_given, schreier *&Sch, INT *set, INT len, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "strong_generators::compute_schreier_with_given_action_on_a_given_set" << endl;
		cout << "action=";
		A->print_info();
		}
	Sch = new schreier;

	Sch->init(A_given);
	Sch->initialize_tables();
	Sch->init_generators(*gens);
	Sch->compute_all_orbits_on_invariant_subset(len, set, 0 /* verbose_level */);
	//Sch->compute_all_point_orbits(verbose_level);


	if (f_v) {
		cout << "strong_generators::compute_schreier_with_given_action_on_a_given_set done, we found " << Sch->nb_orbits << " orbits" << endl;
		}
}

void strong_generators::orbits_on_points(INT &nb_orbits, INT *&orbit_reps, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	schreier *Sch;
	INT i, f, a;

	if (f_v) {
		cout << "strong_generators::orbits_on_points" << endl;
		cout << "action=";
		A->print_info();
		}

	compute_schreier_with_given_action(A, Sch, verbose_level - 1);

#if 0
	Sch = new schreier;

	Sch->init(A);
	Sch->initialize_tables();
	Sch->init_generators(*gens);
	Sch->compute_all_point_orbits(verbose_level);
#endif

	nb_orbits = Sch->nb_orbits;
	orbit_reps = NEW_INT(nb_orbits);
	for (i = 0; i < nb_orbits; i++) {
		f = Sch->orbit_first[i];
		a = Sch->orbit[f];
		orbit_reps[i] = a;
		}

	delete Sch;

	if (f_v) {
		cout << "strong_generators::orbits_on_points done, we found " << nb_orbits << " orbits" << endl;
		}
}

void strong_generators::orbits_on_points_with_given_action(action *A_given, INT &nb_orbits, INT *&orbit_reps, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	schreier *Sch;
	INT i, f, a;

	if (f_v) {
		cout << "strong_generators::orbits_on_points_with_given_action" << endl;
		cout << "action=";
		A->print_info();
		}
	compute_schreier_with_given_action(A_given, Sch, verbose_level - 1);

#if 0
	Sch = new schreier;

	Sch->init(A_given);
	Sch->initialize_tables();
	Sch->init_generators(*gens);
	Sch->compute_all_point_orbits(verbose_level);
#endif

	nb_orbits = Sch->nb_orbits;
	orbit_reps = NEW_INT(nb_orbits);
	for (i = 0; i < nb_orbits; i++) {
		f = Sch->orbit_first[i];
		a = Sch->orbit[f];
		orbit_reps[i] = a;
		}

	delete Sch;

	if (f_v) {
		cout << "strong_generators::orbits_on_points_with_given_action done, we found " << nb_orbits << " orbits" << endl;
		}
}

schreier *strong_generators::orbits_on_points_schreier(action *A_given, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	schreier *Sch;

	if (f_v) {
		cout << "strong_generators::orbits_on_points_schreier degree = " << A_given->degree << endl;
		}
	Sch = new schreier;

	Sch->init(A_given);
	Sch->initialize_tables();
	Sch->init_generators(*gens);
	Sch->compute_all_point_orbits(verbose_level);

	if (f_v) {
		cout << "strong_generators::orbits_on_points_schreier done, we found " << Sch->nb_orbits << " orbits" << endl;
		}
	return Sch;
}

void strong_generators::orbits_light(action *A_given, 
	INT *&Orbit_reps, INT *&Orbit_lengths, INT &nb_orbits, 
	INT **&Pts_per_generator, INT *&Nb_per_generator, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);	
	INT f_vv = FALSE; //(verbose_level >= 2);	
	UBYTE *reached;
	INT Orbit_allocated;
	INT Orbit_len;
	INT *Orbit;
	INT *Q;
	INT Q_allocated;
	INT Q_len;
	INT pt, i, h, nb_gens, a, b, idx;
	INT Orbit_reps_allocated;
	INT nb_reached;
	INT *Generator_idx;

	if (f_v) {
		cout << "strong_generators::orbits_light degree = " << A_given->degree << endl;
		}

	Orbit_reps_allocated = 1024;
	Orbit_reps = NEW_INT(Orbit_reps_allocated);
	Orbit_lengths = NEW_INT(Orbit_reps_allocated);
	nb_orbits = 0;

	if (f_v) {
		cout << "strong_generators::orbits_light allocating array Generator_idx" << endl;
		}
	Generator_idx = NEW_INT(A_given->degree);
	if (f_v) {
		cout << "strong_generators::orbits_light allocating array Generator_idx done" << endl;
		}
	for (pt = 0; pt < A_given->degree; pt++) {
		Generator_idx[pt] = -1;
		}
	reached = bitvector_allocate(A_given->degree);
	nb_reached = 0;

	Orbit_allocated = 1024;
	Orbit = NEW_INT(Orbit_allocated);

	Q_allocated = 1024;
	Q = NEW_INT(Q_allocated);
	
	nb_gens = gens->len;

	if (A_given->degree > ONE_MILLION) {
		f_v = TRUE;
		}

	Nb_per_generator = NEW_INT(nb_gens);
	INT_vec_zero(Nb_per_generator, nb_gens);
	Pts_per_generator = NEW_PINT(nb_gens);

	for (pt = 0; pt < A_given->degree; pt++) {
		if (bitvector_s_i(reached, pt)) {
			continue;
			}
		if (f_vv) {
			cout << "strong_generators::orbits_light computing orbit of point " << pt << endl;
			}
		Q[0] = pt;
		Q_len = 1;
		
		Orbit[0] = pt;
		Orbit_len = 1;

		while (Q_len) {
			if (f_vv) {
				cout << "strong_generators::orbits_light considering the next element in the queue" << endl;
				}
			a = Q[0];
			for (i = 1; i < Q_len; i++) {
				Q[i - 1] = Q[i];
				}
			Q_len--;
			if (f_vv) {
				cout << "strong_generators::orbits_light looking at element " << a << endl;
				}
			for (h = 0; h < nb_gens; h++) {
				if (f_vv) {
					cout << "strong_generators::orbits_light applying generator " << h << endl;
					}
				b = A_given->element_image_of(a, gens->ith(h), FALSE);
				if (f_vv) {
					cout << "strong_generators::orbits_light under generator " << h << " it maps to " << b << endl;
					}
				if (!INT_vec_search(Orbit, Orbit_len, b, idx)) {
					if (Orbit_len == Orbit_allocated) {
						INT new_oa;
						INT *O;

						new_oa = 2 * Orbit_allocated;
						O = NEW_INT(new_oa);
						for (i = 0; i < Orbit_len; i++) {
							O[i] = Orbit[i];
							}
						FREE_INT(Orbit);
						Orbit = O;
						Orbit_allocated = new_oa;
						}
					for (i = Orbit_len; i > idx; i--) {
						Orbit[i] = Orbit[i - 1];
						}
					Orbit[idx] = b;
					Orbit_len++;
					Generator_idx[b] = h;
					Nb_per_generator[h]++;

					if (f_vv) {
						cout << "current orbit: ";
						INT_vec_print(cout, Orbit, Orbit_len);
						cout << endl;
						}

					bitvector_m_ii(reached, b, 1);
					nb_reached++;
					if (f_v && ((nb_reached & ((1 << 18) - 1)) == 0)) {
						cout << "strong_generators::orbits_light nb_reached =  " << nb_reached << " / " << A_given->degree << endl;					
						}

					if (Q_len == Q_allocated) {
						INT new_qa;
						INT *new_Q;

						new_qa = 2 * Q_allocated;
						new_Q = NEW_INT(new_qa);
						for (i = 0; i < Q_len; i++) {
							new_Q[i] = Q[i];
							}
						FREE_INT(Q);
						Q = new_Q;
						Q_allocated = new_qa;
						}

					Q[Q_len++] = b;

					if (f_vv) {
						cout << "current Queue: ";
						INT_vec_print(cout, Q, Q_len);
						cout << endl;
						}

					}
				} // next h
			} // while (Q_len)

		if (f_vv) {
			cout << "Orbit of point " << pt << " has length " << Orbit_len << endl;
			}
		if (nb_orbits == Orbit_reps_allocated) {
			INT an;
			INT *R;
			INT *L;

			an = 2 * Orbit_reps_allocated;
			R = NEW_INT(an);
			L = NEW_INT(an);
			for (i = 0; i < nb_orbits; i++) {
				R[i] = Orbit_reps[i];
				L[i] = Orbit_lengths[i];
				}
			FREE_INT(Orbit_reps);
			FREE_INT(Orbit_lengths);
			Orbit_reps = R;
			Orbit_lengths = L;
			Orbit_reps_allocated = an;
			}
		Orbit_reps[nb_orbits] = pt;
		Orbit_lengths[nb_orbits] = Orbit_len;
		nb_orbits++;
		} // for pt
	if (f_v) {
		cout << "strong_generators::orbits_light degree = " << A_given->degree << " we found " << nb_orbits << " orbits" << endl;
		cout << i << " : " << Nb_per_generator[i] << endl;
		for (i = 0; i < nb_gens; i++) {
			cout << i << " : " << Nb_per_generator[i] << endl;
			}
		}


	if (f_v) {
		cout << "strong_generators::orbits_light computing the arrays Pts_per_generator" << endl;
		}
	for (i = 0; i < nb_gens; i++) { 
		INT *v;
		INT j;

		v = NEW_INT(Nb_per_generator[i]);
		j = 0;
		for (pt = 0; pt < A_given->degree; pt++) {
			if (Generator_idx[pt] == i) {
				v[j] = pt;
				j++;
				}
			}
		if (j != Nb_per_generator[i]) {
			cout << "strong_generators::orbits_light j != Nb_per_generator[i]" << endl;
			exit(1);
			}
		Pts_per_generator[i] = v;
		}

	FREE_INT(Orbit);
	FREE_INT(Q);
	FREE_UBYTE(reached);
	FREE_INT(Generator_idx);
	//FREE_INT(Nb_per_generator);
	if (f_v) {
		cout << "strong_generators::orbits_light degree = " << A_given->degree << " we found " << nb_orbits << " orbits" << endl;
		}
}

void strong_generators::write_to_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "strong_generators::write_to_memory_object" << endl;
		}
	m->write_int(A->base_len);
	for (i = 0; i < A->base_len; i++) {
		m->write_int(tl[i]);
		}
	gens->write_to_memory_object(m, verbose_level - 1);
	if (f_v) {
		cout << "strong_generators::write_to_memory_object done" << endl;
		}
}

void strong_generators::read_from_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, l;

	if (f_v) {
		cout << "strong_generators::read_from_memory_object" << endl;
		}
	m->read_int(&l);
	cout << "strong_generators::read_from_memory_object l=" << l << endl;
	if (l != A->base_len) {
		cout << "strong_generators::read_from_memory_object l != A->base_len" << endl;
		}
	tl = NEW_INT(A->base_len);
	gens = new vector_ge;
	gens->init(A);
	for (i = 0; i < A->base_len; i++) {
		m->read_int(&tl[i]);
		}

	cout << "strong_generators::read_from_memory_object tl=";
	INT_vec_print(cout, tl, A->base_len);
	cout << endl;
	
	gens->read_from_memory_object(m, verbose_level - 1);
	if (f_v) {
		cout << "strong_generators::read_from_memory_object done" << endl;
		}
}

void strong_generators::write_to_file_binary(ofstream &fp, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "strong_generators::write_to_file_binary" << endl;
		}
	fp.write((char *) &A->base_len, sizeof(INT));
	for (i = 0; i < A->base_len; i++) {
		fp.write((char *) &tl[i], sizeof(INT));
		}
	gens->write_to_file_binary(fp, verbose_level - 1);
}

void strong_generators::read_from_file_binary(action *A, ifstream &fp, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, l;

	if (f_v) {
		cout << "strong_generators::read_from_file_binary" << endl;
		}
	init(A, 0);
	fp.read((char *) &l, sizeof(INT));
	if (l != A->base_len) {
		cout << "strong_generators::read_from_file_binary l != A->base_len" << endl;
		exit(1);
		}
	tl = NEW_INT(A->base_len);
	for (i = 0; i < A->base_len; i++) {
		fp.read((char *) &tl[i], sizeof(INT));
		}
	gens = new vector_ge;
	gens->init(A);
	gens->read_from_file_binary(fp, verbose_level - 1);
	if (f_v) {
		cout << "strong_generators::read_from_file_binary done" << endl;
		}
}

void strong_generators::write_file(BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	
	if (f_v) {
		cout << "strong_generators::write_file" << endl;
		}
	M.alloc(1024 /* length */, verbose_level - 1);
	M.used_length = 0;
	M.cur_pointer = 0;
	write_to_memory_object(&M, verbose_level - 1);
	M.write_file(fname, verbose_level - 1);
	if (f_v) {
		cout << "strong_generators::write_file done" << endl;
		}
}

void strong_generators::read_file(action *A, const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	
	if (f_v) {
		cout << "strong_generators::read_file reading file " << fname << " of size " << file_size(fname) << endl;
		}
	init(A, 0);
	M.read_file(fname, verbose_level - 1);
	if (f_v) {
		cout << "strong_generators::read_file read file " << fname << endl;
		}
	M.cur_pointer = 0;
	read_from_memory_object(&M, verbose_level - 1);
	if (f_v) {
		cout << "strong_generators::read_file done" << endl;
		}
}

void strong_generators::generators_for_shallow_schreier_tree(BYTE *label, vector_ge *chosen_gens, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	action *AR;
	sims *S;
	INT go;
	double avg;
	double log_go;
	INT cnt = 0;
	INT i;
	
	go = group_order_as_INT();
	log_go = log(go);
	if (f_v) {
		cout << "strong_generators::generators_for_shallow_schreier_tree group of order " << go << endl;
		cout << "log_go = " << log_go << endl;
		}
	S = create_sims(verbose_level - 2);
	if (f_v) {
		cout << "strong_generators::generators_for_shallow_schreier_tree created sims" << endl;
		}
	AR = new_action_by_right_multiplication(S, TRUE /* f_transfer_ownership */, verbose_level - 2);
	if (f_v) {
		cout << "strong_generators::generators_for_shallow_schreier_tree created action by right multiplication" << endl;
		}

	chosen_gens->init(A);
	chosen_gens->allocate(gens->len);
	for (i = 0; i < gens->len; i++) {
		A->element_move(gens->ith(i), chosen_gens->ith(i), 0);
		}

	while (TRUE) {

		schreier *Sch;

		Sch = new schreier;
		Sch->init(AR);
		Sch->initialize_tables();
		Sch->init_generators(*chosen_gens);
		if (f_v) {
			cout << "strong_generators::generators_for_shallow_schreier_tree before computing all orbits" << endl;
			}
		Sch->compute_all_point_orbits(verbose_level - 2);
		if (f_v) {
			cout << "strong_generators::generators_for_shallow_schreier_tree after computing all orbits" << endl;
			}
		if (Sch->nb_orbits > 1) {
			cout << "strong_generators::generators_for_shallow_schreier_tree Sch->nb_orbits > 1" << endl;
			exit(1);
			}
		BYTE label1[1000];
		INT xmax = 1000000;
		INT ymax = 1000000;
		INT f_circletext = TRUE;
		INT rad = 3000;

		sprintf(label1, "%s_%ld", label, cnt);
		Sch->draw_tree(label1, 0 /* orbit_no */, xmax, ymax, f_circletext, rad, 0 /* verbose_level */);
		INT *Depth;
		INT avgi, f, l, idx;

		Depth = NEW_INT(Sch->A->degree);
		for (i = 0; i < Sch->A->degree; i++) {
			Depth[i] = Sch->depth_in_tree(i);
			}
		classify Cl;

		Cl.init(Depth, Sch->A->degree, FALSE, 0);
		if (f_v) {
			cout << "distribution of depth in tree is: ";
			Cl.print(TRUE);
			cout << endl;
			}
		avg = Cl.average();
		if (f_v) {
			cout << "average = " << avg << endl;
			cout << "log_go = " << log_go << endl;
			}

		if (avg < log_go) {
			if (f_v) {
				cout << "strong_generators::generators_for_shallow_schreier_tree average < log_go, we are done" << endl;
				}
			break;
			}

		avgi = (INT) avg;
		if (f_v) {
			cout << "average as INT = " << avgi << endl;
			}
		for (i = 0; i < Cl.nb_types; i++) {
			f = Cl.type_first[i];
			l = Cl.type_len[i];
			if (Cl.data_sorted[f] == avgi) {
				break;
				}
			}
		if (i == Cl.nb_types) {
			cout << "strong_generators::generators_for_shallow_schreier_tree cannot find element of depth " << avgi << endl;
			exit(1);
			}
		idx = Cl.sorting_perm_inv[f];
		if (f_v) {
			cout << "strong_generators::generators_for_shallow_schreier_tree idx = " << idx << endl;
			}
		Sch->coset_rep(idx);
		chosen_gens->append(Sch->cosetrep);
		

		FREE_INT(Depth);
		delete Sch;	
		cnt++;
		}

		
	if (f_v) {
		cout << "strong_generators::generators_for_shallow_schreier_tree done" << endl;
		}
}


void strong_generators::compute_ascii_coding(BYTE *&ascii_coding, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT sz, i, j;
	BYTE *p;

	if (f_v) {
		cout << "strong_generators::compute_ascii_coding" << endl;
		}
	sz = 2 * ((2 + A->base_len + A->base_len) * sizeof(INT4) + A->coded_elt_size_in_char * gens->len) + 1;
	ascii_coding = NEW_BYTE(sz);
	p = ascii_coding;
	code_INT4(p, (INT4) A->base_len);
		// in GALOIS/util.C
	code_INT4(p, (INT4) gens->len);
	for (i = 0; i < A->base_len; i++) {
		code_INT4(p, (INT4) A->base[i]);
		}
	for (i = 0; i < A->base_len; i++) {
		code_INT4(p, (INT4) tl[i]);
		}
	for (i = 0; i < gens->len; i++) {
		A->element_pack(gens->ith(i), A->elt1, FALSE);
		for (j = 0; j < A->coded_elt_size_in_char; j++) {
			code_UBYTE(p, A->elt1[j]);
			}
		}
	*p++ = 0;
	if (p - ascii_coding != sz) {
		cout << "strong_generators::compute_ascii_coding p - ascii_coding != sz" << endl;
		exit(1);
		}
	

	if (f_v) {
		cout << "strong_generators::compute_ascii_coding done" << endl;
		}
}

void strong_generators::decode_ascii_coding(BYTE *ascii_coding, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT str_len, len, nbsg, i, j;
	BYTE *p, *p0;
	action *A_save;
	INT *base1;

	if (f_v) {
		cout << "strong_generators::decode_ascii_coding" << endl;
		}

	// clean up before we go:
	A_save = A;
	freeself();
	A = A_save;

	p = ascii_coding;
	p0 = p;
	str_len = strlen(ascii_coding);
	len = decode_INT4(p);
	nbsg = decode_INT4(p);
	if (len != A->base_len) {
		cout << "strong_generators::decode_ascii_coding len != A->base_len" << endl;
		cout << "len=" << len << " (from file)" << endl;
		cout << "A->base_len=" << A->base_len << endl;
		cout << "action A is " << A->label << endl;
		exit(1);
		}
	gens = new vector_ge;
	gens->init(A);
	gens->allocate(nbsg);
	base1 = NEW_INT(A->base_len);
	tl = NEW_INT(A->base_len);
	for (i = 0; i < A->base_len; i++) {
		base1[i] = decode_INT4(p);
		}
	for (i = 0; i < A->base_len; i++) {
		if (base1[i] != A->base[i]) {
			cout << "strong_generators::decode_ascii_coding base element " << i << " does not match current base" << endl;
			exit(1);
			}
		}
	for (i = 0; i < A->base_len; i++) {
		tl[i] = decode_INT4(p);
		}
	for (i = 0; i < nbsg; i++) {
		for (j = 0; j < A->coded_elt_size_in_char; j++) {
			decode_UBYTE(p, A->elt1[j]);
			}
		A->element_unpack(A->elt1, gens->ith(i), FALSE);
		}
	FREE_INT(base1);
	if (p - p0 != str_len) {
		cout << "strong_generators::decode_ascii_coding p - p0 != str_len" << endl;
		cout << "p - p0 = " << p - p0 << endl;
		cout << "str_len = " << str_len << endl;
		exit(1);
		}
	if (f_v) {
		cout << "strong_generators::decode_ascii_coding done" << endl;
		}
}

void strong_generators::export_permutation_group_to_magma(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "strong_generators::export_permutation_group_to_magma" << endl;
		}
	{
		ofstream fp(fname);

		fp << "G := sub< Sym(" << A->degree << ") |" << endl;
		for (i = 0; i < gens->len; i++) {
			A->element_print_as_permutation_with_offset(gens->ith(i), fp, 
				1 /* offset */, TRUE /* f_do_it_anyway_even_for_big_degree */, 
				FALSE /* f_print_cycles_of_length_one */, 0 /* verbose_level */);
			if (i < gens->len - 1) {
				fp << ", ";
				}
			}
		fp << ">;" << endl;

	}
	cout << "Written file " << fname << " of size " << file_size(fname) << endl;

	if (f_v) {
		cout << "strong_generators::export_permutation_group_to_magma done" << endl;
		}
}


void strong_generators::compute_and_print_orbits_on_a_given_set(action *A_given, INT *set, INT len, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	schreier *Sch;
	INT i, j, f, l, a;

	if (f_v) {
		cout << "strong_generators::compute_and_print_orbits_on_a_given_set" << endl;
		}
	compute_schreier_with_given_action_on_a_given_set(A_given, Sch, set, len, verbose_level - 2);

	cout << "orbits on the set: " << endl;
	for (i = 0; i < Sch->nb_orbits; i++) {
		f = Sch->orbit_first[i];
		for (j = 0; j < Sch->orbit_len[i]; j++) {
			a = Sch->orbit[f + j];
			cout << a << " ";
			}
		if (i < Sch->nb_orbits - 1) {
			cout << "| ";
			}
		}
	cout << endl;
	cout << "partition: " << len << " = ";
	for (i = 0; i < Sch->nb_orbits; i++) {
		l = Sch->orbit_len[i];
		cout << l << " ";
		if (i < Sch->nb_orbits - 1) {
			cout << "+ ";
			}
		}
	cout << endl;
	cout << "representatives for each of the " << Sch->nb_orbits << " orbits:" << endl;
	for (i = 0; i < Sch->nb_orbits; i++) {
		f = Sch->orbit_first[i];
		l = Sch->orbit_len[i];
		a = Sch->orbit[f + j];
		cout << setw(5) << a << " : " << setw(5) << l << " : ";
		A_given->print_point(a, cout);
		cout << endl;
		}
	delete Sch;

}

void strong_generators::compute_and_print_orbits(action *A_given, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	schreier *Sch;
	INT i, j, f, l, a;

	if (f_v) {
		cout << "strong_generators::compute_and_print_orbits" << endl;
		}
	compute_schreier_with_given_action(A_given, Sch, verbose_level - 2);

	cout << "orbits on the set: " << endl;
	for (i = 0; i < Sch->nb_orbits; i++) {
		f = Sch->orbit_first[i];
		l = Sch->orbit_len[i];
		if (l >= 10) {
			cout << "too long to list ";
			}
		else {
			for (j = 0; j < Sch->orbit_len[i]; j++) {
				a = Sch->orbit[f + j];
				cout << a << " ";
				}
			}
		if (i < Sch->nb_orbits - 1) {
			cout << "| ";
			}
		}
	cout << endl;
	cout << "partition: " << A_given->degree << " = ";
	for (i = 0; i < Sch->nb_orbits; i++) {
		l = Sch->orbit_len[i];
		cout << l << " ";
		if (i < Sch->nb_orbits - 1) {
			cout << "+ ";
			}
		}
	cout << endl;
	cout << "representatives for each of the " << Sch->nb_orbits << " orbits:" << endl;
	for (i = 0; i < Sch->nb_orbits; i++) {
		f = Sch->orbit_first[i];
		l = Sch->orbit_len[i];
		a = Sch->orbit[f + 0];
		cout << setw(5) << a << " : " << setw(5) << l << " : ";
		A_given->print_point(a, cout);
		cout << endl;
		}

#if 0
	set_of_sets *S;

	Sch->orbits_as_set_of_sets(S, 0 /* verbose_level */);

	const BYTE *fname = "orbits.csv";
	
	cout << "writing orbits to file " << fname << endl;
	S->save_csv(fname, 1 /* verbose_level */);
	
	delete S;
#endif
	delete Sch;

}

//#######################################################################################################
// global functions:
//#######################################################################################################

void strong_generators_write_file(const BYTE *fname, strong_generators *p, INT nb, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	INT i;
	
	if (f_v) {
		cout << "strong_generators_write_file" << endl;
		}
	M.alloc(1024 /* length */, verbose_level - 1);
	M.used_length = 0;
	M.cur_pointer = 0;

	M.write_int(nb);
	for (i = 0; i < nb; i++) {
		p[i].write_to_memory_object(&M, verbose_level - 1);
		}
	M.write_file(fname, verbose_level - 1);
	cout << "strong_generators_write_file written file " << fname << " of size " << file_size(fname) << endl;
}

void strong_generators_read_from_file(const BYTE *fname, action *A, strong_generators *&p, INT &nb, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	INT i;
	
	if (f_v) {
		cout << "strong_generators_read_from_file" << endl;
		}

	M.read_file(fname, verbose_level - 1);
	if (f_v) {
		cout << "strong_generators_read_from_file reading file " << fname << " of size " << file_size(fname) << endl;
		}
	M.cur_pointer = 0;

	M.read_int(&nb);
	if (f_v) {
		cout << "strong_generators_read_from_file reading " << nb << " stabilizers" << endl;
		}
	p = new strong_generators[nb];
	for (i = 0; i < nb; i++) {
		p[i].A = A;
		p[i].read_from_memory_object(&M, verbose_level - 2);
		}
	if (f_v) {
		cout << "strong_generators_read_from_file reading " << nb << " stabilizers done" << endl;
		}

	INT *Go;
	Go = NEW_INT(nb);
	for (i = 0; i < nb; i++) {
		Go[i] = p[i].group_order_as_INT();
		}
	classify C;

	C.init(Go, nb, FALSE, 0);
	cout << "Stabilizer orders:" << endl;
	C.print(TRUE);

#if 0
	INT f, l, j, pos;

	cout << "stabilizer order : orbit" << endl;
	for (i = C.nb_types - 1; i > 0; i--) {
		// don't do i = 0, it corresponds to to orbits with trivial stabilizer
		f = C.type_first[i];
		l = C.type_len[i];
		for (j = 0; j < l; j++) {
			pos = C.sorting_perm_inv[f + j];
			cout << Go[pos] << " : " << pos << endl;
			}
		}
	#endif
	
	FREE_INT(Go);


}


