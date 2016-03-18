// schreier.C
//
// Anton Betten
// December 9, 2003

#include "galois.h"
#include "action.h"

INT schreier::cntr_new = 0;
INT schreier::cntr_objects = 0;
INT schreier::f_debug_memory = FALSE;




void *schreier::operator new(size_t bytes)
{
	cntr_new++;
	cntr_objects++;
	if (f_debug_memory) {
		cout << "schreier::operator new bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	return malloc(bytes);
}

void *schreier::operator new[](size_t bytes)
{
	INT n;
	
	n = bytes / sizeof(schreier);
	cntr_new++;
	cntr_objects += n;
	if (f_debug_memory) {
		cout << "schreier::operator new[] n=" << n 
			<< " bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	return malloc(bytes);
}

void schreier::operator delete(void *ptr, size_t bytes)
{
	if (f_debug_memory) {
		cout << "schreier::operator delete bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects--;
	return free(ptr);
}

void schreier::operator delete[](void *ptr, size_t bytes)
{
	INT n;
	
	n = bytes / sizeof(schreier);
	if (f_debug_memory) {
		cout << "schreier::operator delete[] n=" << n 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects -= n;
	return free(ptr);
}

schreier::schreier()
{
	A = NULL;
	nb_images = 0;
	images = NULL;
	f_print_function = FALSE;
}

schreier::schreier(action *A)
{
	init(A);
};

schreier::~schreier()
{
	//cout << "in ~schreier()" << endl;
	freeself();
	//cout << "~schreier() finished" << endl;
};

void schreier::freeself()
{
	//cout << "deleting A" << endl;
	if (A) {
		//cout << "deleting orbit" << endl;
		FREE_INT(orbit);
		//cout << "deleting orbit_inv" << endl;
		FREE_INT(orbit_inv);
		//cout << "deleting prev" << endl;
		FREE_INT(prev);
		//cout << "deleting label" << endl;
		FREE_INT(label);
		//cout << "deleting orbit_no" << endl;
		FREE_INT(orbit_no);
		//cout << "deleting orbit_first" << endl;
		FREE_INT(orbit_first);
		//cout << "deleting orbit_len" << endl;
		FREE_INT(orbit_len);
		//cout << "deleting Elt1" << endl;
		FREE_INT(Elt1);
		//cout << "deleting Elt2" << endl;
		FREE_INT(Elt2);
		//cout << "deleting Elt3" << endl;
		FREE_INT(Elt3);
		//cout << "deleting schreier_gen" << endl;
		FREE_INT(schreier_gen);
		//cout << "deleting schreier_gen1" << endl;
		FREE_INT(schreier_gen1);
		//cout << "deleting cosetrep" << endl;
		FREE_INT(cosetrep);
		//cout << "deleting cosetrep_tmp" << endl;
		FREE_INT(cosetrep_tmp);
		//cout << "A = NULL" << endl;
		A = NULL;
		}
	//cout << "deleting images" << endl;
	delete_images();
}

void schreier::delete_images()
{
	INT i;
	
	if (images) {
		for (i = 0; i < nb_images; i++) {
			FREE_INT(images[i]);
			}
		FREE_PINT(images);
		images = NULL;
		nb_images = 0;
		}
}

void schreier::init_images(INT nb_images, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, j;
	
	if (f_v) {
		cout << "schreier::init_images" << endl;
		}
	if (A == NULL) {
		cout << "schreier::init_images() action is NULL" << endl;
		exit(1);
		}
	delete_images();
	schreier::nb_images = nb_images;
	images = NEW_PINT(nb_images);
	for (i = 0; i < nb_images; i++) {
		images[i] = NEW_INT(2 * A->degree);
		for (j = 0; j < 2 * A->degree; j++) {
			images[i][j] = -1;
			}
		}
	if (f_v) {
		cout << "schreier::init_images done" << endl;
		}
}

void schreier::images_append()
{
	INT **new_images = NEW_PINT(nb_images + 1);
	INT i, j;
	
	new_images[nb_images] = NEW_INT(2 * A->degree);
	for (j = 0; j < 2 * A->degree; j++) {
		new_images[nb_images][j] = -1;
		}
	for (i = 0; i < nb_images; i++) {
		new_images[i] = images[i];
		}
	FREE_PINT(images);
	images = new_images;
	nb_images++;
}

void schreier::init(action *A)
{
	schreier::A = A;
	orbit = NEW_INT(A->degree);
	orbit_inv = NEW_INT(A->degree);
	prev = NEW_INT(A->degree);
	label = NEW_INT(A->degree);
	orbit_no = NEW_INT(A->degree);
	orbit_first = NEW_INT(A->degree + 1);
	orbit_len = NEW_INT(A->degree);
	gens.init(A);
	gens_inv.init(A);
	initialize_tables();
	init2();
}

void schreier::init2()
{
	Elt1 = NEW_INT(A->elt_size_in_INT);
	Elt2 = NEW_INT(A->elt_size_in_INT);
	Elt3 = NEW_INT(A->elt_size_in_INT);
	schreier_gen = NEW_INT(A->elt_size_in_INT);
	schreier_gen1 = NEW_INT(A->elt_size_in_INT);
	cosetrep = NEW_INT(A->elt_size_in_INT);
	cosetrep_tmp = NEW_INT(A->elt_size_in_INT);
}

void schreier::initialize_tables()
{
	INT i;
	
	nb_orbits = 0;
	perm_identity(orbit, A->degree);
	perm_identity(orbit_inv, A->degree);
	orbit_first[0] = 0;
	for (i = 0; i < A->degree; i++) {
		prev[i] = -1;
		label[i] = -1;
		orbit_no[i] = -1;
		}
}

void schreier::init_single_generator(INT *elt)
{
	init_generators(1, elt);
}

void schreier::init_generators(vector_ge &generators)
{
	if (generators.len) {
		init_generators(generators.len, generators.ith(0));
		}
	else {
		init_generators(generators.len, NULL);
		}
}

void schreier::init_generators(INT nb, INT *elt)
// elt must point to nb * A->elt_size_in_INT INT's that are 
// group elements in INT format
{
	INT i;
	
	gens.allocate(nb);
	gens_inv.allocate(nb);
	for (i = 0; i < nb; i++) {
		//cout << "schreier::init_generators i = " << i << endl;
		gens.copy_in(i, elt + i * A->elt_size_in_INT);
		A->element_invert(elt + i * A->elt_size_in_INT, gens_inv.ith(i), 0);
		}
	init_images(nb, 0 /* verbose_level */);	
}

void schreier::init_generators_by_hdl(INT nb_gen, INT *gen_hdl, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "schreier::init_generators_by_hdl" << endl;
		cout << "nb_gen = " << nb_gen << endl;
		cout << "degree = " << A->degree << endl;
		}
	gens.allocate(nb_gen);
	gens_inv.allocate(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		//cout << "schreier::init_generators_by_hdl i = " << i << endl;
		A->element_retrieve(gen_hdl[i], gens.ith(i), 0);
		
		//cout << "schreier::init_generators_by_hdl generator i = " << i << ":" << endl;
		//A->element_print_quick(gens.ith(i), cout);

		A->element_invert(gens.ith(i), gens_inv.ith(i), 0);
		}
	if (f_v) {
		cout << "schreier::init_generators_by_hdl before init_images()" << endl;
		}
	init_images(nb_gen, verbose_level);	
	if (f_v) {
		cout << "schreier::init_generators_by_hdl done" << endl;
		}
}

INT schreier::get_image(INT i, INT gen_idx, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT a;
	
	if (f_v) {
		cout << "schreier::get_image i=" << i << endl;
		}
	if (images == NULL) {
		cout << "schreier::get_image() images == NULL" << endl;
		exit(1);
		}
	a = images[gen_idx][i];
	if (a == -1) {
		a = A->element_image_of(i, gens.ith(gen_idx), verbose_level - 2);
		if (f_v) {
			cout << "schreier::get_image image of i=" << i << " is " << a << endl;
			}
		images[gen_idx][i] = a;
		images[gen_idx][A->degree + a] = i;
		}
	return a;
}

void schreier::print_orbit_lengths(ostream &ost)
{
	INT i, f, l, m;
	INT *orbit_len_sorted;
	INT *sorting_perm;
	INT *sorting_perm_inv;
	INT nb_types;
	INT *type_first;
	INT *type_len;
	
	INT_vec_classify(nb_orbits, orbit_len, orbit_len_sorted, 
		sorting_perm, sorting_perm_inv, 
		nb_types, type_first, type_len);

	ost << nb_orbits << " orbits: " << endl;
	for (i = 0; i < nb_types; i++) {
		f = type_first[i];
		l = type_len[i];
		m = orbit_len_sorted[f];
		if (l > 1) {
			cout << l << " \\times ";
			}
		cout << m;
		if (i < nb_types - 1)
			cout << ", ";
		}
	ost << endl;
	FREE_INT(orbit_len_sorted);
	FREE_INT(sorting_perm);
	FREE_INT(sorting_perm_inv);
	FREE_INT(type_first);
	FREE_INT(type_len);
	
}

void schreier::print_orbit_length_distribution(ostream &ost)
{
	INT *val, *mult, len;	
	
	INT_vec_distribution(orbit_len, nb_orbits, val, mult, len);
	INT_distribution_print(ost, val, mult, len);
	ost << endl;
	
	FREE_INT(val);
	FREE_INT(mult);
}


void schreier::print_orbit_reps(ostream &ost)
{
	INT i, c, r;
	
	ost << nb_orbits << " orbits" << endl;
	ost << "orbits of a group with " << gens.len << " generators:" << endl;
	ost << "i : orbit_first[i] : orbit_len[i] : rep" << endl;
	for (i = 0; i < nb_orbits; i++) {
		ost << setw(3) << i << " : " << setw(6) << orbit_first[i] << " : " << setw(6) << orbit_len[i];
		c = orbit_first[i];
		r = orbit[c];
		ost << " : " << setw(6) << r << endl;
		//<< " : ";
		//print_orbit(ost, i);
		//ost << endl;
		}
	ost << endl;
}


void schreier::print(ostream &ost)
{
	INT i;
	
	ost << nb_orbits << " orbits" << endl;
	ost << "orbit group with " << gens.len << " generators:" << endl;
	ost << "i : orbit_first[i] : orbit_len[i]" << endl;
	for (i = 0; i < nb_orbits; i++) {
		ost << i << " : " << orbit_first[i] << " : " << orbit_len[i] << endl;
		//<< " : ";
		//print_orbit(ost, i);
		//ost << endl;
		}
	ost << endl;
}

void schreier::print_and_list_orbits(ostream &ost)
{
	INT i;
	
	ost << nb_orbits << " orbits" << endl;
	ost << "orbit group with " << gens.len << " generators:" << endl;
	ost << "i : orbit_first[i] : orbit_len[i]" << endl;
	for (i = 0; i < nb_orbits; i++) {
		ost << i << " : " << orbit_first[i] << " : " << orbit_len[i];
		ost << " : ";
		print_orbit(ost, i);
		ost << endl;
		}
	ost << endl;
}

void schreier::print_and_list_orbits_using_labels(ostream &ost, INT *labels)
{
	INT i;
	
	ost << nb_orbits << " orbits" << endl;
	ost << "orbit group with " << gens.len << " generators:" << endl;
	ost << "i : orbit_first[i] : orbit_len[i]" << endl;
	for (i = 0; i < nb_orbits; i++) {
		ost << i << " : " << orbit_first[i] << " : " << orbit_len[i];
		ost << " : ";
		print_orbit_using_labels(ost, i, labels);
		ost << endl;
		}
	ost << endl;
}

void schreier::print_tables(ostream &ost, INT f_with_cosetrep)
{
	INT i, w; //  j, k;
	
#if 0
	ost << gens.len << " generators:" << endl;
	for (i = 0; i < A->degree; i++) {
		ost << i;
		for (j = 0; j < gens.len; j++) {
			k = A->element_image_of(i, gens.ith(j), FALSE);
			ost << " : " << k;
			}
		ost << endl;
		}
	ost << endl;
#endif
	w = INT_log10(A->degree) + 1;
	ost << "i : orbit_no[i] : orbit[i] : orbit_inv[i] : prev[i] : label[i]";
	if (f_with_cosetrep)
		ost << " : coset_rep";
	ost << endl;
	for (i = 0; i < A->degree; i++) {
		coset_rep(i);
		//coset_rep_inv(i);
		ost << setw(w) << i << " : " << setw(w) << orbit_no[i] << " : " 
			<< setw(w) << orbit[i] << " : " << setw(w) << orbit_inv[i] << " : " 
			<< setw(w) << prev[i] << " : " << setw(w) << label[i];
		if (f_with_cosetrep) {
			ost << " : ";
			//A->element_print(Elt1, cout);
			A->element_print_as_permutation(cosetrep, ost);
			ost << endl;
			A->element_print_quick(cosetrep, ost);
			}
		ost << endl;
		}
	ost << endl;
}

void schreier::print_generators()
{
	INT j;
	
	cout << gens.len << " generators in action " << A->label << " of degree " << A->degree << ":" << endl;
	for (j = 0; j < gens.len; j++) {
		cout << "generator " << j << ":" << endl;
		//A->element_print(gens.ith(j), cout);
		A->element_print_quick(gens.ith(j), cout);
		A->element_print_as_permutation(gens.ith(j), cout);
		if (j < gens.len - 1) {
			cout << ", " << endl;
			}
		}
}

void schreier::print_orbit(INT orbit_no)
{
	print_orbit(cout, orbit_no);
}

void schreier::print_orbit_using_labels(INT orbit_no, INT *labels)
{
	print_orbit_using_labels(cout, orbit_no, labels);
}

void schreier::print_orbit(ostream &ost, INT orbit_no)
{
	INT i, first, len;
	INT *v;
	
	first = orbit_first[orbit_no];
	len = orbit_len[orbit_no];
	v = NEW_INT(len);
	for (i = 0; i < len; i++) {
		v[i] = orbit[first + i];
		}
	//INT_vec_print(ost, v, len);
	INT_vec_heapsort(v, len);
	INT_vec_print_fully(ost, v, len);
	
	FREE_INT(v);
}

void schreier::print_orbit_using_labels(ostream &ost, INT orbit_no, INT *labels)
{
	INT i, first, len;
	INT *v;
	
	first = orbit_first[orbit_no];
	len = orbit_len[orbit_no];
	v = NEW_INT(len);
	for (i = 0; i < len; i++) {
		v[i] = labels[orbit[first + i]];
		}
	//INT_vec_print(ost, v, len);
	INT_vec_heapsort(v, len);
	INT_vec_print_fully(ost, v, len);
	
	FREE_INT(v);
}

void schreier::print_orbit_type(INT f_backwards)
{
	classify C;

	C.init(orbit_len, nb_orbits, FALSE, 0);
	C.print_naked(f_backwards);
}

void schreier::list_all_orbits_tex(ostream &ost)
{
	INT i, j, f, l, a;

	ost << "$";
	for (i = 0; i < nb_orbits; i++) {
		f = orbit_first[i];
		l = orbit_len[i];
		for (j = 0; j < l; j++) {
			a = orbit[f + j];
			ost << a;
			if (j < l - 1) {
				ost << ", ";
				}
			}
		if (i < nb_orbits - 1) {
			ost << " \\mid ";
			}
		}
	ost << "$";
}

void schreier::print_orbit_through_labels(ostream &ost, INT orbit_no, INT *point_labels)
{
	INT i, first, len;
	INT *v;
	
	first = orbit_first[orbit_no];
	len = orbit_len[orbit_no];
	v = NEW_INT(len);
	for (i = 0; i < len; i++) {
		v[i] = point_labels[orbit[first + i]];
		}
	INT_vec_heapsort(v, len);
	INT_vec_print_fully(ost, v, len);
	FREE_INT(v);
}

void schreier::print_orbit_sorted(ostream &ost, INT orbit_no)
{
	INT i, len;
	INT *v;
	
	len = orbit_first[orbit_no + 1] - orbit_first[orbit_no];
	v = NEW_INT(len);
	for (i = 0; i < len; i++) {
		v[i] = orbit[orbit_first[orbit_no] + i];
		}
	INT_vec_sort(len, v);
	
	ost << "{ ";
	for (i = 0; i < len; i++) {
		if (f_print_function) {
			ost << v[i] << "=";
			(*print_function)(ost, v[i], print_function_data);
			}
		else {
			ost << v[i];
			}
		if (i < len - 1)
			ost << ", ";
		}
	ost << " }";
	FREE_INT(v);
}

void schreier::print_orbit(INT cur, INT last)
{
	INT i;
	
	for (i = 0; i < A->degree; i++) {
		if (i == cur) 
			cout << ">";
		if (i == last)
			cout << ">";
		cout << i << " : " << orbit[i] << " : " << orbit_inv[i] << endl;
		}
	cout << endl;
}

void schreier::swap_points(INT i, INT j)
{
	INT pi, pj;
	
	pi = orbit[i];
	pj = orbit[j];
	orbit[i] = pj;
	orbit[j] = pi;
	orbit_inv[pi] = j;
	orbit_inv[pj] = i;
}

void schreier::move_point_here(INT here, INT pt)
{
	INT a, loc;
	if (orbit[here] == pt)
		return;
	a = orbit[here];
	loc = orbit_inv[pt];
	orbit[here] = pt;
	orbit[loc] = a;
	orbit_inv[a] = loc;
	orbit_inv[pt] = here;
}

INT schreier::orbit_representative(INT pt)
{
	INT j;
	
	while (TRUE) {
		j = orbit_inv[pt];
		if (prev[j] == -1)
			return pt;
		pt = prev[j];
		}
}

INT schreier::depth_in_tree(INT j)
// j is a coset, not a point
{
	if (prev[j] == -1) {
		return 0;
		}
	else {
		return depth_in_tree(orbit_inv[prev[j]]) + 1;
		}
}

void schreier::coset_rep(INT j)
// j is a coset, not a point
// result is in cosetrep
// determines an element in the group that moves the orbit representative 
// to the j-th point in the orbit.
{
	INT *gen;
	
	if (prev[j] != -1) {
		coset_rep(orbit_inv[prev[j]]);
		gen = gens.ith(label[j]);
		A->element_mult(cosetrep, gen, cosetrep_tmp, 0);
		A->element_move(cosetrep_tmp, cosetrep, 0);
		}
	else {
		A->element_one(cosetrep, 0);
		}
}

void schreier::coset_rep_inv(INT j)
// j is a coset, not a point
// result is in cosetrep
{
	INT *gen;
	
	if (prev[j] != -1) {
		coset_rep_inv(orbit_inv[prev[j]]);
		gen = gens_inv.ith(label[j]);
		A->element_mult(gen, cosetrep, cosetrep_tmp, 0);
		A->element_move(cosetrep_tmp, cosetrep, 0);
		}
	else {
		A->element_one(cosetrep, 0);
		}
}

void schreier::get_schreier_vector(INT *&sv, INT f_trivial_group, INT f_compact)
{
	if (f_compact) {
		get_schreier_vector_compact(sv, f_trivial_group);
		}
	else {
		get_schreier_vector_ordinary(sv);
		}
}

void schreier::get_schreier_vector_compact(INT *&sv, INT f_trivial_group)
// allocated and creates array sv[size] using NEW_INT
// where size is n + 1 if  f_trivial_group is TRUE
// and size is 3 * n + 1 otherwise
// Here, n is the combined size of all orbits counted by nb_orbits
// sv[0] is equal to n
// sv + 1 is the array point_list of size [n], listing the point in increasing order
// Unless f_trivial_group, sv + 1 + n is the array prev[n] and 
// sv + 1 + 2 * n is the array label[n] 
{
	INT i, j, k, f, ff, l, p, pr, la, n = 0;
	INT *point_list;
	
	for (k = 0; k < nb_orbits; k++) {
		n += orbit_len[k];
		}
	point_list = NEW_INT(n);
	
	ff = 0;
	for (k = 0; k < nb_orbits; k++) {
		f = orbit_first[k];
		l = orbit_len[k];
		for (j = 0; j < l; j++) {
			i = f + j;
			p = orbit[i];
			point_list[ff + j] = p;
			}
		ff += l;
		}
	if (ff != n) {
		cout << "schreier::get_schreier_vector_compact ff != n" << endl;
		exit(1);
		}
	INT_vec_heapsort(point_list, n);
	
	
	if (f_trivial_group) {
		sv = NEW_INT(n + 1);
		}
	else {
		sv = NEW_INT(3 * n + 1);
		}
	sv[0] = n;
	for (i = 0; i < n; i++) {
		sv[1 + i] = point_list[i];
		}
	if (!f_trivial_group) {
		for (i = 0; i < n; i++) {
			p = point_list[i];
			j = orbit_inv[p];
			pr = prev[j];
			la = label[j];
			sv[1 + n + i] = pr;
			sv[1 + 2 * n + i] = la;
			}
		}
	FREE_INT(point_list);
}

void schreier::get_schreier_vector_ordinary(INT *&sv)
// allocates and creates array sv[2 * A->degree] using NEW_INT
// sv[i * 2 + 0] is prev[i]
// sv[i * 2 + 1] is label[i]
{
	INT i, j;
	
	sv = NEW_INT(2 * A->degree);
	for (i = 0; i < A->degree; i++) {
		j = orbit_inv[i];
		if (prev[j] != -1) {
			sv[i * 2 + 0] = prev[j];
			sv[i * 2 + 1] = label[j];
			//cout << "label[" << i << "] = " << label[j] << endl;
			}
		else {
			sv[i * 2 + 0] = -1;
			sv[i * 2 + 1] = -1;
			}
		}
}

void schreier::extend_orbit(INT *elt, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	//INT f_vvv = (verbose_level >= 3);
	INT cur, total0, total, cur_pt, gen_first, i, next_pt, next_pt_loc;
	
	if (f_v) {
		cout << "extend_orbit() extending orbit " << nb_orbits - 1 << " of length " 
			<< orbit_len[nb_orbits - 1] << endl;
		}

	gens.append(elt);
	A->element_invert(elt, A->Elt1, FALSE);
	gens_inv.append(A->Elt1);
	images_append();
	
	cur = orbit_first[nb_orbits - 1];
	total = total0 = orbit_first[nb_orbits];
	while (cur < total) {
		cur_pt = orbit[cur];
		if (FALSE) {
			cout << "schreier::extend_orbit applying generator to " << cur_pt << endl;
			}
#if 0
		if (cur < total0)
			gen_first = gens.len - 1;
		else 
			gen_first = 0;
#endif
		gen_first = 0;
		for (i = gen_first; i < gens.len; i++) {
			next_pt = get_image(cur_pt, i, 0/*verbose_level - 3*/);
				// A->element_image_of(cur_pt, gens.ith(i), FALSE);
			next_pt_loc = orbit_inv[next_pt];
			if (FALSE) {
				cout << "schreier::extend_orbit generator " << i << " maps " << cur_pt << " to " << next_pt << endl;
				}
			if (next_pt_loc < total)
				continue;
			if (FALSE) {
				cout << "schreier::extend_orbit new pt " << next_pt << " reached from " << cur_pt << " under generator " << i << endl;
				}
			swap_points(total, next_pt_loc);
			prev[total] = cur_pt;
			label[total] = i;
			orbit_no[total] = nb_orbits - 1;
			total++;
			if (FALSE) {
				cout << "cur = " << cur << endl;
				cout << "total = " << total << endl;
				print_orbit(cur, total - 1);
				}
			}
		cur++;
		}
	orbit_first[nb_orbits] = total;
	orbit_len[nb_orbits - 1] = total - orbit_first[nb_orbits - 1];
	//orbit_first[nb_orbits + 1] = A->degree;
	//orbit_len[nb_orbits] = A->degree - total;
	if (f_v) {
		cout << "schreier::extend_orbit orbit extended to length " << orbit_len[nb_orbits - 1] << endl;
		}
	if (FALSE) {
		cout << "{ ";
		for (i = orbit_first[nb_orbits - 1]; i < orbit_first[nb_orbits]; i++) {
			cout << orbit[i];
			if (i < orbit_first[nb_orbits] - 1)
				cout << ", ";
			}
		cout << " }" << endl;
		}
}

void schreier::compute_all_point_orbits(INT verbose_level)
{
	INT pt, pt_loc, cur, pt0;
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		cout << "schreier::compute_all_point_orbits" << endl;
		}
	if (A->degree > ONE_MILLION) {
		f_vv = FALSE;
		}
	initialize_tables();
	for (pt0 = 0, pt = 0; pt < A->degree; pt++) {
		pt_loc = orbit_inv[pt];
		cur = orbit_first[nb_orbits];
		if (pt_loc < cur) {
			continue;
			}
		if (f_vv) {
			cout << "schreier::compute_all_point_orbits pt = " << pt << " / " << A->degree << " nb_orbits=" << nb_orbits << " computing orbit" << endl;
			}
		if (A->degree > ONE_MILLION && (pt - pt0) > 50000) {
			cout << "schreier::compute_all_point_orbits pt = " << pt << " / " << A->degree << " nb_orbits=" << nb_orbits << " computing orbit" << endl;
			pt0 = pt;
			}
		compute_point_orbit(pt, verbose_level - 2);
		}
	if (f_v) {
		cout << "schreier::compute_all_point_orbits found " << nb_orbits << " orbits" << endl;
		classify Cl;

		Cl.init(orbit_len, nb_orbits, FALSE, 0);
		cout << "The distribution of orbit lengths is: ";
		Cl.print(FALSE);
		}
}

void schreier::compute_all_point_orbits_with_prefered_reps(
	INT *prefered_reps, INT nb_prefered_reps, INT verbose_level)
{
	INT i, pt, pt_loc, cur;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "schreier::compute_all_point_orbits_with_prefered_reps" << endl;
		}
	initialize_tables();
	for (i = 0; i < nb_prefered_reps; i++) {
		pt = prefered_reps[i];
		pt_loc = orbit_inv[pt];
		cur = orbit_first[nb_orbits];
		if (pt_loc < cur) {
			continue;
			}
		compute_point_orbit(pt, verbose_level - 1);
		}
	for (pt = 0; pt < A->degree; pt++) {
		pt_loc = orbit_inv[pt];
		cur = orbit_first[nb_orbits];
		if (pt_loc < cur) {
			continue;
			}
		compute_point_orbit(pt, verbose_level - 1);
		}
	if (f_v) {
		cout << "found " << nb_orbits << " orbit";
		if (nb_orbits != 1)
			cout << "s";
		cout << " on points" << endl;
		}
}


void schreier::compute_all_point_orbits_with_preferred_labels(INT *preferred_labels, INT verbose_level)
{
	INT pt, pt_loc, cur, a, i;
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT *labels, *perm, *perm_inv;
	
	if (f_v) {
		cout << "schreier::compute_all_point_orbits_with_preferred_labels" << endl;
		//cout << "preferred_labels :";
		//INT_vec_print(cout, preferred_labels, A->degree);
		//cout << endl;
		cout << "A->degree = " << A->degree << endl;
		}
	if (f_v) {
		cout << "schreier::compute_all_point_orbits_with_preferred_labels allocating tables" << endl;
		}
	initialize_tables();
	labels = NEW_INT(A->degree);
	perm = NEW_INT(A->degree);
	perm_inv = NEW_INT(A->degree);
	for (i = 0; i < A->degree; i++) {
		labels[i] = preferred_labels[i];
		}
	if (f_v) {
		cout << "schreier::compute_all_point_orbits_with_preferred_labels allocating tables done, sorting" << endl;
		}
	INT_vec_sorting_permutation(labels, A->degree, perm, perm_inv, TRUE /* f_increasingly */);

	if (f_v) {
		cout << "schreier::compute_all_point_orbits_with_preferred_labels sorting done" << endl;
		}
	
	for (a = 0; a < A->degree; a++) {
		pt = perm_inv[a];
		pt_loc = orbit_inv[pt];
		cur = orbit_first[nb_orbits];
		if (pt_loc < cur) {
			continue;
			}
		// now we need to make sure that the point pt is moved to position cur:
		// actually this is not needed as the function compute_point_orbit does this, too.
		swap_points(cur, pt_loc);
		
		if (f_v) {
			cout << "schreier::compute_all_point_orbits_with_preferred_labels computing orbit of point " << pt << " = " << a << " / " << A->degree << endl;
			}
		compute_point_orbit(pt, verbose_level - 2);
		if (f_v) {
			cout << "schreier::compute_all_point_orbits_with_preferred_labels computing orbit of point " << pt << " done, found an orbit of length " << orbit_len[nb_orbits - 1] << " nb_orbits = " << nb_orbits << endl;
			}
		}
	if (f_v) {
		cout << "found " << nb_orbits << " orbit";
		if (nb_orbits != 1)
			cout << "s";
		cout << " on points" << endl;
		}
	FREE_INT(labels);
	FREE_INT(perm);
	FREE_INT(perm_inv);
}

void schreier::compute_all_orbits_on_invariant_subset(INT len, INT *subset, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, f;
	
	if (f_v) {
		cout << "schreier::compute_all_orbits_on_invariant_subset() computing orbits on a set of size " << len << endl;
		}
	initialize_tables();
	for (i = 0; i < len; i++) {
		move_point_here(i, subset[i]);
		}
	while (TRUE) {
		f = orbit_first[nb_orbits];
		if (f >= len)
			break;
		compute_point_orbit(orbit[f], 0 /* verbose_level */);
		}
	if (f > len) {
		cout << "schreier::compute_all_orbits_on_invariant_subset the set is not G-invariant" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "found " << nb_orbits << " orbits" << endl;
		print_orbit_length_distribution(cout);
		}
}

INT schreier::sum_up_orbit_lengths()
{
	INT i, l, N;
	
		N = 0;
	for (i = 0; i < nb_orbits; i++) {
		l = orbit_len[i];
		N += l;
		}
	return N;
}

void schreier::compute_point_orbit(INT pt, INT verbose_level)
{
	INT pt_loc, cur, cur_pt, total, i, next_pt, next_pt_loc, total1, cur1;
	INT f_v = (verbose_level >= 1);
	INT f_vv = FALSE; // (verbose_level >= 5);
	//INT f_vvv = FALSE; //(verbose_level >= 3);
	
	if (f_v) {
		cout << "schreier::compute_point_orbit computing orbit of point " << pt << " in action " << A->label << endl;
		}
	pt_loc = orbit_inv[pt];
	cur = orbit_first[nb_orbits];
	if (pt_loc < cur) {
		cout << "schreier::compute_point_orbit() i < orbit_first[nb_orbits]" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "schreier::compute_point_orbit computing orbit of pt " << pt << endl;
		}
	if (pt_loc > orbit_first[nb_orbits]) {
		swap_points(orbit_first[nb_orbits], pt_loc);
		}
	orbit_no[orbit_first[nb_orbits]] = nb_orbits;
	total = cur + 1;
	while (cur < total) {
		cur_pt = orbit[cur];
		if (f_vv) {
			cout << "schreier::compute_point_orbit cur=" << cur << " total=" << total << " applying generators to " << cur_pt << endl;
			}
		for (i = 0; i < gens.len; i++) {
			if (f_vv) {
				cout << "schreier::compute_point_orbit applying generator " << i << " to point " << cur_pt << endl;
				}
			next_pt = get_image(cur_pt, i, 0 /*verbose_level*/);
				// A->element_image_of(cur_pt, gens.ith(i), FALSE);
			next_pt_loc = orbit_inv[next_pt];
			if (f_vv) {
				cout << "schreier::compute_point_orbit generator " << i << " maps " << cur_pt << " to " << next_pt << endl;
				}
			if (next_pt_loc < total)
				continue;
			if (f_vv) {
				cout << "schreier::compute_point_orbit new pt " << next_pt << " reached from " << cur_pt << " under generator " << i << endl;
				}
			swap_points(total, next_pt_loc);
			prev[total] = cur_pt;
			label[total] = i;
			orbit_no[total] = nb_orbits;
			total++;
			total1 = total - orbit_first[nb_orbits];
			cur1 = cur - orbit_first[nb_orbits];
			if ((total1 % 10000) == 0 || (cur1 > 0 && (cur1 % 10000) == 0)) {
				cout << "schreier::compute_point_orbit degree = " << A->degree << " length = " << total1 
					<< " processed = " << cur1 << " nb_orbits=" << nb_orbits << " cur_pt=" << cur_pt << " next_pt=" << next_pt << " orbit_first[nb_orbits]=" << orbit_first[nb_orbits] << endl;
				}
			if (FALSE) {
				cout << "cur = " << cur << endl;
				cout << "total = " << total << endl;
				print_orbit(cur, total - 1);
				}
			}
		cur++;
		}
	orbit_first[nb_orbits + 1] = total;
	orbit_len[nb_orbits] = total - orbit_first[nb_orbits];
	//orbit_first[nb_orbits + 2] = A->degree;
	//orbit_len[nb_orbits + 1] = A->degree - total;
	if (f_v) {
		cout << "found orbit of length " << orbit_len[nb_orbits] << " total length " << total << " degree=" << A->degree << endl;
		}
	if (FALSE) {
		cout << "{ ";
		for (i = orbit_first[nb_orbits]; i < orbit_first[nb_orbits + 1]; i++) {
			cout << orbit[i];
			if (i < orbit_first[nb_orbits + 1] - 1)
				cout << ", ";
			}
		cout << " }" << endl;
		}
	if (FALSE) {
		cout << "coset reps:" << endl;
		for (i = orbit_first[nb_orbits]; i < orbit_first[nb_orbits + 1]; i++) {
			cout << i << " : " << endl;
			coset_rep(i);
			A->element_print(cosetrep, cout);
			cout << "image = " << orbit[i] << " = " << A->element_image_of(pt, cosetrep, 0) << endl;
			cout << endl;
			
			}
		}
	nb_orbits++;
}

void schreier::non_trivial_random_schreier_generator(action *A_original, INT verbose_level)
// computes non trivial random Schreier generator into schreier_gen
// non-trivial is with respect to A_original
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = FALSE; //(verbose_level >= 3);
	INT f_v4 = FALSE; //(verbose_level >= 4);
	INT cnt = 0;
	
	if (f_v) {
		cout << "schreier::non_trivial_random_schreier_generator verbose_level=" << verbose_level << endl;
		}
	while (TRUE) {
		if (f_v) {
			cout << "schreier::non_trivial_random_schreier_generator calling random_schreier_generator" << endl;
			}
		random_schreier_generator(verbose_level - 1);
		cnt++;
		if (!A_original->element_is_one(schreier_gen, verbose_level - 5)) {
			if (f_vv) {
				cout << "schreier::non_trivial_random_schreier_generator found a non-trivial random Schreier generator in " << cnt << " trials" << endl;
				}
			if (f_vvv) {
				A->element_print(schreier_gen, cout);
				cout << endl;
				}
			return;
			}
		else {
			if (f_v4) {
				A->element_print(schreier_gen, cout);
				cout << endl;
				}
			if (f_vv) {
				cout << "schreier::non_trivial_random_schreier_generator the element is the identity in action " << A_original->label << ", trying again" << endl;
				}
			}
		}
}

void schreier::random_schreier_generator_ith_orbit(INT orbit_no, INT verbose_level)
{
	INT first, len, r1, r2, pt, pt2, pt2_coset;
	INT *gen;
	INT f_v = (verbose_level >= 1);
	INT f_vv = FALSE; //(verbose_level >= 2);
	INT f_vvv = FALSE; //(verbose_level >= 3);
	
	if (f_v) {
		cout << "schreier::random_schreier_generator_ith_orbit, orbit " << orbit_no << endl;
		}
	if (f_vvv) {
		cout << "generators are:" << endl;
		gens.print(cout);
		}
	first = orbit_first[orbit_no];
	len = orbit_len[orbit_no];
	pt = orbit[first];
	if (f_vv) {
		cout << "pt=" << pt << endl;
		cout << "orbit_first[orbit_no]=" << orbit_first[orbit_no] << endl;
		cout << "orbit_len[orbit_no]=" << orbit_len[orbit_no] << endl;
		cout << "gens.len=" << gens.len << endl;
		}
	
	// get a random coset:
	r1 = random_integer(orbit_len[orbit_no]);
	if (f_vv) {
		cout << "r1=" << r1 << endl;
		}
	//pt1 = orbit[r1];
	coset_rep(orbit_first[orbit_no] + r1);
	// coset rep now in cosetrep
	if (f_vvv) {
		cout << "cosetrep " << orbit_first[orbit_no] + r1 << endl;
		A->element_print_quick(cosetrep, cout);
		if (A->degree < 100) {
			A->element_print_as_permutation(cosetrep, cout);
			cout << endl;
			}
		}
		
	// get a random generator:
	r2 = random_integer(gens.len);
	if (f_vv) {
		cout << "r2=" << r2 << endl;
		}
	gen = gens.ith(r2);
	if (f_vvv) {
		cout << "generator " << r2 << endl;
		A->element_print(gen, cout);
		if (A->degree < 100) {
			A->element_print_as_permutation(gen, cout);
			cout << endl;
			}
		}
	if (f_vv) {
		cout << "random coset " << r1 << ", random generator " << r2 << endl;
		}
	
	A->element_mult(cosetrep, gen, schreier_gen1, 0);
	if (f_vvv) {
		cout << "cosetrep * generator " << endl;
		A->element_print_quick(schreier_gen1, cout);
		if (A->degree < 100) {
			A->element_print_as_permutation(schreier_gen1, cout);
			cout << endl;
			}
		}
	pt2 = A->element_image_of(pt, schreier_gen1, 0);
	if (f_vv) {
		//cout << "pt2=" << pt2 << endl;
		cout << "maps " << pt << " to " << pt2 << endl;
		}
	pt2_coset = orbit_inv[pt2];
	if (f_vv) {
		cout << "pt2_coset=" << pt2_coset << endl;
		}
	if (pt2_coset < first) {
		cout << "schreier::random_schreier_generator_ith_orbit pt2_coset < first" << endl;
		exit(1);
		}
	if (pt2_coset >= first + len) {
		cout << "schreier::random_schreier_generator_ith_orbit pt2_coset >= first + len" << endl;
		exit(1);
		}
	
	coset_rep_inv(pt2_coset);
	// coset rep now in cosetrep
	if (f_vvv) {
		cout << "cosetrep (inverse) " << pt2_coset << endl;
		A->element_print_quick(cosetrep, cout);
		if (A->degree < 100) {
			A->element_print_as_permutation(cosetrep, cout);
			cout << endl;
			}
		}
	
	A->element_mult(schreier_gen1, cosetrep, schreier_gen, 0);
	if (A->element_image_of(pt, schreier_gen, 0) != pt) {
		cout << "schreier::random_schreier_generator_ith_orbit() fatal: schreier generator does not stabilize pt" << endl;
		exit(1);
		}
	if (f_vv) {
		cout << "schreier::random_schreier_generator_ith_orbit done" << endl;
		}
	if (f_vvv) {
		A->element_print_quick(schreier_gen, cout);
		cout << endl;
		if (A->degree < 100) {
			A->element_print_as_permutation(cosetrep, cout);
			cout << endl;
			}
		}
}

void schreier::random_schreier_generator(INT verbose_level)
// computes random Schreier generator into schreier_gen
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = FALSE; // (verbose_level >= 2);
	INT r1, r2, pt, pt2, pt2_coset;
	INT *gen;
	INT pt1, pt1b;
	
	if (f_v) {
		cout << "schreier::random_schreier_generator orbit_len = " 
			<< orbit_len[0] << " nb generators = " << gens.len << " in action " << A->label << endl;
		}
	pt = orbit[0];
	if (f_vv) {
		cout << "pt=" << pt << endl;
		}
	
	// get a random coset:
	r1 = random_integer(orbit_len[0]);
	pt1 = orbit[r1];
	
	coset_rep(r1);
	// coset rep now in cosetrep
	pt1b = A->element_image_of(pt, cosetrep, 0);
	if (f_vv) {
		cout << "random coset " << r1 << endl;
		cout << "pt1=" << pt1 << endl;
		cout << "cosetrep:" << endl;
		A->element_print_quick(cosetrep, cout);
		cout << "image of pt under cosetrep = " << pt1b << endl;
		}
	if (pt1b != pt1) {
		cout << "schreier::random_schreier_generator fatal: cosetrep does not work" << endl;
		cout << "pt=" << pt << endl;
		cout << "random coset " << r1 << endl;
		cout << "pt1=" << pt1 << endl;
		cout << "cosetrep:" << endl;
		A->element_print_quick(cosetrep, cout);
		cout << "image of pt under cosetrep = " << pt1b << endl;
		A->element_image_of(pt, cosetrep, 10);	
		exit(1);
		}
	
	// get a random generator:
	r2 = random_integer(gens.len);
	gen = gens.ith(r2);
	if (f_vv) {
		cout << "random coset " << r1 << ", random generator " << r2 << endl;
		cout << "generator:" << endl;
		A->element_print_quick(gen, cout);
		cout << "image of pt1 under generator = pt2 = " << A->element_image_of(pt1, gen, 0) << endl;
		}
	
	A->element_mult(cosetrep, gen, schreier_gen1, 0);
	if (f_vv) {
		cout << "cosetrep * gen=" << endl;
		A->element_print_quick(schreier_gen1, cout);
		}
	pt2 = A->element_image_of(pt, schreier_gen1, 0);
	if (f_vv) {
		cout << "image of pt under cosetrep*gen = " << pt2 << endl;
		}
	//cout << "maps " << pt << " to " << pt2 << endl;
	pt2_coset = orbit_inv[pt2];
	
	coset_rep_inv(pt2_coset);
	// coset rep now in cosetrep
	if (f_vv) {
		cout << "cosetrep:" << endl;
		A->element_print_quick(cosetrep, cout);
		cout << "image of pt2 under cosetrep = " << A->element_image_of(pt2, cosetrep, 0) << endl;
		}
	
	A->element_mult(schreier_gen1, cosetrep, schreier_gen, 0);
	if (f_vv) {
		cout << "schreier_gen=cosetrep*gen*cosetrep:" << endl;
		A->element_print_quick(schreier_gen, cout);
		cout << "image of pt under schreier_gen = " << A->element_image_of(pt, schreier_gen, 0) << endl;
		}
	if (A->element_image_of(pt, schreier_gen, 0) != pt) {
		cout << "schreier::random_schreier_generator() fatal: schreier generator does not stabilize pt" << endl;
		exit(1);
		}
	if (FALSE) {
		cout << "random Schreier generator:" << endl;
		A->element_print(schreier_gen, cout);
		cout << endl;
		}
}

void schreier::trace_back(INT *path, INT i, INT &j)
{
	INT ii = orbit_inv[i];
	
	if (prev[ii] == -1) {
		if (path) {
			path[0] = i;
			}
		j = 1;
		}
	else {
		trace_back(path, prev[ii], j);
		if (path) {
			path[j] = i;
			}
		j++;
		}
}

void schreier::print_tree(INT orbit_no)
{
	INT *path;
	INT i, j, l;
	
	path = NEW_INT(A->degree);
	i = orbit_first[orbit_no];
	while (i < orbit_first[orbit_no + 1]) {
		trace_back(path, orbit[i], l);
		// now l is the distance from the root
		cout << l;
		for (j = 0; j < l; j++) {
			cout << " " << path[j];
			}
		cout << " 0 ";
		if (label[i] != -1) {
			cout << " $s_{" << label[i] << "}$";
			}
		cout << endl;
		i++;
		}
	FREE_INT(path);
}

void schreier::draw_tree(char *label, INT orbit_no, INT xmax, INT ymax, INT f_circletext, INT rad, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *path;
	INT *weight;
	INT *placement_x;
	INT i, j, last, max_depth = 0;


	if (f_v) {
		cout << "schreier::draw_tree" << endl;
		}
	path = NEW_INT(A->degree);
	weight = NEW_INT(A->degree);
	placement_x = NEW_INT(A->degree);
		
	i = orbit_first[orbit_no];
	last = orbit_first[orbit_no + 1];
	
	for (j = 0; j < A->degree; j++) {
		weight[j] = 0;
		placement_x[j] = 0;
		}
	subtree_calc_weight(weight, max_depth, i, last);
	if (f_vv) {
		cout << "the weights: " << endl;
		for (j = i; j < last; j++) {
			cout << j << " : " << weight[j] << " : " << endl;
			}
		cout << endl;
		cout << "max_depth = " << max_depth << endl;
		}
	subtree_place(weight, placement_x, 0, 1000000, i, last);
	if (f_vv) {
		for (j = i; j < last; j++) {
			cout << j << " : " << placement_x[j] << endl;
			}
		cout << endl;
		}
	if (orbit_len[orbit_no] > 100) {
		f_circletext = FALSE;
		}
	draw_tree2(label, xmax, ymax, f_circletext, weight, placement_x, max_depth, i, last, rad, verbose_level - 2);
	
#if 0
	{
	ofstream f("tree");
	
	f << "COORDS_RANGE 1000 600" << endl;
	f << "TREE " << A->degree << endl;
	trace_back(path, orbit[i], l);
	// now l is the distance from the root
	print_path(f, path, l);
		
	subtree_depth_first(f, path, i, last);
	f << "TREE_END" << endl;
	}
#endif

	FREE_INT(path);
	FREE_INT(weight);
	FREE_INT(placement_x);
	if (f_v) {
		cout << "schreier::draw_tree done" << endl;
		}
}

static void calc_y_coordinate(INT &y, INT l, INT max_depth)
{
	INT dy;
	
	dy = (INT)((double)1000000 / (double)max_depth);
	y = (INT)(dy * ((double)l + 0.5));
	y = 1000000 - y;
}

void schreier::draw_tree2(char *fname, INT xmax, INT ymax, INT f_circletext, 
	INT *weight, INT *placement_x, INT max_depth, INT i, INT last, INT rad, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT x_min = 0, x_max = 1000000;
	INT y_min = 0, y_max = 1000000;
	INT factor_1000 = 1000;
	BYTE fname_full[1000];
	INT f_embedded = TRUE;
	INT f_sideways = FALSE;
	
	if (f_v) {
		cout << "schreier::draw_tree2" << endl;
		}
	sprintf(fname_full, "%s.mp", fname);
	mp_graphics G(fname_full, x_min, y_min, x_max, y_max, f_embedded, f_sideways);
	G.out_xmin() = 0;
	G.out_ymin() = 0;
	G.out_xmax() = xmax;
	G.out_ymax() = ymax;
	
	G.header();
	G.begin_figure(factor_1000);
	
	INT x = 500000, y;
	calc_y_coordinate(y, 0, max_depth);
	
	if (f_circletext) {
		G.circle_text(x, y, "$\\emptyset$");
		}
	else {
		G.circle(x, y, 5);
		}

	subtree_draw_lines(G, f_circletext, x, y, weight, placement_x, max_depth, i, last, verbose_level);

	subtree_draw_vertices(G, f_circletext, x, y, weight, placement_x, max_depth, i, last, rad, verbose_level);


	G.draw_boxes_final();
	G.end_figure();
	G.footer();
	if (f_v) {
		cout << "schreier::draw_tree2 done" << endl;
		}
}

void schreier::subtree_draw_lines(mp_graphics &G, INT f_circletext, INT parent_x, INT parent_y, INT *weight, 
	INT *placement_x, INT max_depth, INT i, INT last, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT pt = orbit[i];
	INT x, y, l, ii;
	INT Px[2], Py[2];
	
	if (f_v) {
		cout << "schreier::subtree_draw_lines" << endl;
		}
	trace_back(NULL, pt, l);
	x = placement_x[pt];
	calc_y_coordinate(y, l, max_depth);

	//G.circle(x, y, 2000);
	Px[0] = parent_x;
	Py[0] = parent_y;
	Px[1] = x;
	Py[1] = y;
	G.polygon2(Px, Py, 0, 1);
	
	for (ii = i + 1; ii < last; ii++) {
		if (prev[ii] == pt) {
			subtree_draw_lines(G, f_circletext, x, y, weight, placement_x, max_depth, ii, last, verbose_level);
			}
		}

	if (f_v) {
		cout << "schreier::subtree_draw_lines done" << endl;
		}
}

void schreier::subtree_draw_vertices(mp_graphics &G, INT f_circletext, INT parent_x, INT parent_y, INT *weight, 
	INT *placement_x, INT max_depth, INT i, INT last, INT rad, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT pt = orbit[i];
	INT x, y, l, ii;
	INT Px[2], Py[2];
	char str[1000];
	
	if (f_v) {
		cout << "schreier::subtree_draw_vertices" << endl;
		}
	trace_back(NULL, pt, l);
	x = placement_x[pt];
	calc_y_coordinate(y, l, max_depth);

	G.circle(x, y, rad);
	Px[0] = parent_x;
	Py[0] = parent_y;
	Px[1] = x;
	Py[1] = y;
	//G.polygon2(Px, Py, 0, 1);
	
	for (ii = i + 1; ii < last; ii++) {
		if (prev[ii] == pt) {
			subtree_draw_vertices(G, f_circletext, x, y, weight, placement_x, max_depth, ii, last, rad, verbose_level);
			}
		}
	sprintf(str, "%ld", pt);
	if (f_circletext) {
		G.circle_text(x, y, str);
		}
	else {
		//G.aligned_text(Px, Py, 1, "tl", str);
		}
	if (f_v) {
		cout << "schreier::subtree_draw_vertices done" << endl;
		}
}

void schreier::subtree_place(INT *weight, INT *placement_x, INT left, INT right, INT i, INT last)
{
	INT pt = orbit[i];
	INT ii, l, w, w0, w1, lft, rgt, width;
	double dx;
	
	placement_x[pt] = (left + right) >> 1;
	w = weight[pt];
	width = right - left;
	dx = width / (double) (w - 1);
		// the node itself counts for the weight, so we subtract one
	w0 = 0;
	
	trace_back(NULL, pt, l);
	for (ii = i + 1; ii < last; ii++) {
		if (prev[ii] == pt) {
			w1 = weight[orbit[ii]];
			lft = left + (INT)((double)w0 * dx);
			rgt = left + (INT)((double)(w0 + w1) * dx);
			subtree_place(weight, placement_x, lft, rgt, ii, last);
			w0 += w1;
			}
		}
}

INT schreier::subtree_calc_weight(INT *weight, INT &max_depth, INT i, INT last)
{
	INT pt = orbit[i];
	INT ii, l, w = 1, w1;
	
	trace_back(NULL, pt, l);
	if (l > max_depth)
		max_depth = l;
	for (ii = i + 1; ii < last; ii++) {
		if (prev[ii] == pt) {
			w1 = subtree_calc_weight(weight, max_depth, ii, last);
			w += w1;
			}
		}
	weight[pt] = w;
	return w;
}

INT schreier::subtree_depth_first(ostream &ost, INT *path, INT i, INT last)
{
	INT pt = orbit[i];
	INT ii, l, w = 1, w1;
	
	for (ii = i + 1; ii < last; ii++) {
		if (prev[ii] == pt) {
		
			
			trace_back(path, orbit[ii], l);
			// now l is the distance from the root
			print_path(ost, path, l);

			w1 = subtree_depth_first(ost, path, ii, last);
			w += w1;
			}
		}
	return w;
}

void schreier::print_path(ostream &ost, INT *path, INT l)
{
	INT j;
	
	ost << l;
	for (j = 0; j < l; j++) {
		ost << " " << path[j];
		}
	ost << endl;
}

void schreier::intersection_vector(INT *set, INT len, INT *intersection_cnt)
{
	INT i, pt, pt_loc, o;
	
	for (i = 0; i < nb_orbits; i++) {
		intersection_cnt[i] = 0;
		}
	for (i = 0; i < len; i++) {
		pt = set[i];
		pt_loc = orbit_inv[pt];
		o = orbit_no[pt_loc];
		intersection_cnt[o]++;
		}
}

void schreier::orbits_on_invariant_subset_fast(INT len, INT *subset, INT verbose_level)
{
	INT i, p, j;
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	
	if (f_v) {
		cout << "schreier::orbits_on_invariant_subset_fast "
			"computing orbits on invariant subset of size " << len << " in action ";
		A->print_info();
		}
	
	for (i = 0; i < len; i++) {
		p = subset[i];
		j = orbit_inv[p];
		if (j >= orbit_first[nb_orbits]) {
			if (f_vvv) {
				cout << "computing orbit no " << nb_orbits << endl;
				}
			compute_point_orbit(p, 0);
			}
		}
#if 0
	if (orbit_first[nb_orbits] != len) {
		cout << "schreier::orbits_on_invariant_subset_fast orbit_first[nb_orbits] != len" << endl;
		cout << "orbit_first[nb_orbits] = " << orbit_first[nb_orbits] << endl;
		cout << "len = " << len << endl;
		cout << "subset:" << endl;
		INT_vec_print(cout, subset, len);
		cout << endl;
		print_tables(cout, FALSE);
		exit(1);
		}
#endif
	if (f_v) {
		cout << "schreier::orbits_on_invariant_subset_fast "
			"found " << nb_orbits << " orbits on the invariant subset of size " << len << endl;
		}
}

void schreier::orbits_on_invariant_subset(INT len, INT *subset, INT &nb_orbits_on_subset, INT *&orbit_perm, INT *&orbit_perm_inv)
{
	INT i, j, a, pos;
	
	compute_all_point_orbits(0);
	nb_orbits_on_subset = 0;
	orbit_perm = NEW_INT(nb_orbits);
	orbit_perm_inv = NEW_INT(nb_orbits);
	for (i = 0; i < nb_orbits; i++) {
		orbit_perm_inv[i] = -1;
		}
	for (i = 0; i < nb_orbits; i++) {
		j = orbit_first[i];
		a = orbit[j];
		for (pos = 0; pos < len; pos++) {
			if (subset[pos] == a) {
				orbit_perm[nb_orbits_on_subset] = i;
				orbit_perm_inv[i] = nb_orbits_on_subset;
				nb_orbits_on_subset++;
				break;
				}
			}
		}
	j = nb_orbits_on_subset;
	for (i = 0; i < nb_orbits; i++) {
		if (orbit_perm_inv[i] == -1) {
			orbit_perm[j] = i;
			orbit_perm_inv[i] = j;
			j++;
			}
		}	
}

void schreier::get_orbit_partition_of_points_and_lines(partitionstack &S, INT verbose_level)
{
	INT first_column_element, pos, first_column_orbit, i, j, f, l, a;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "schreier::get_orbit_partition_of_points_and_lines" << endl;
		}
	first_column_element = S.startCell[1];
	if (f_v) {
		cout << "first_column_element = " << first_column_element << endl;
		}
	pos = orbit_inv[first_column_element];
	first_column_orbit = orbit_no[pos];
	
	for (i = first_column_orbit - 1; i > 0; i--) {
		f = orbit_first[i];
		l = orbit_len[i];
		for (j = 0; j < l; j++) {
			pos = f + j;
			a = orbit[pos];
			S.subset[j] = a;
			}
		S.subset_size = l;
		S.split_cell(FALSE);
		}
	for (i = nb_orbits - 1; i > first_column_orbit; i--) {
		f = orbit_first[i];
		l = orbit_len[i];
		for (j = 0; j < l; j++) {
			pos = f + j;
			a = orbit[pos];
			S.subset[j] = a;
			}
		S.subset_size = l;
		S.split_cell(FALSE);
		}
}

void schreier::get_orbit_partition(partitionstack &S, INT verbose_level)
{
	INT pos, i, j, f, l, a;
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "schreier::get_orbit_partition" << endl;
		}
	for (i = nb_orbits - 1; i > 0; i--) {
		f = orbit_first[i];
		l = orbit_len[i];
		for (j = 0; j < l; j++) {
			pos = f + j;
			a = orbit[pos];
			S.subset[j] = a;
			}
		S.subset_size = l;
		S.split_cell(FALSE);
		}
}

#if 0
void test_schreier(INT k, INT q, INT f_semilinear, INT verbose_level)
{
	INT i, j, h;
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	INT f_basis = TRUE;
	
	action A;
	
	//matrix_group M;
	
	cout << "calling M.init_PGL " << k << "," << q << " (in standard action)" << endl;
	A.init_matrix_group(FALSE /* f_projective */, k, q, "", 
		f_semilinear, f_basis, verbose_level);
	cout << "a group of degree " << A.degree << endl;
	
	{
	INT v[10];
	for (i = 0; i < MINIMUM(10, A.degree); i++) {
		PG_element_unrank_modified(*A.G.matrix_grp->GFq, v, 1, k, i);
		cout << i << " : ";
		print_set(cout, k, v);
		cout << endl;
		}
	}

	INT *elt1;
	INT *elt2;
	INT *elt3;
	INT *elt4;
	
	elt1 = NEW_INT(A.elt_size_in_INT);
	elt2 = NEW_INT(A.elt_size_in_INT);
	elt3 = NEW_INT(A.elt_size_in_INT);
	elt4 = NEW_INT(A.elt_size_in_INT);
	


	vector_ge V(&A);
	INT l = 3;
	
	V.allocate(l);
	
	for (h = 0; h < l; h++) {
		if (h == 0) {
			i = 0;
			}
		else {
			i = random_integer(A.base_len);
			}
		j = random_integer(A.transversal_length[i]);
		A.get_transversal_rep(i, j, elt1, FALSE);
		// A.element_print(elt1, cout);
		V.copy_in(h, elt1);
		}
	V.print(cout);

	schreier O;
	
	cout << "calling O.init()" << endl;
	
	O.init(&A);

	cout << "calling O.init_generators()" << endl;
	O.init_generators(V);
		
	O.print_generators();

	cout << "calling O.compute_point_orbit()" << endl;
	O.compute_point_orbit(0, verbose_level);
	O.print(cout);
	if (f_vvv)
		O.print_tables(cout, FALSE);

	while (O.orbit_len[0] < A.transversal_length[0]) {
		
		cout << "orbit length " << O.orbit_len[0] << " / " << A.transversal_length[0] << endl;
		
		if (f_v) {
			cout << "action of degree " << A.degree << endl;
			O.print(cout);
			//cout << O.orbit_len[1] << " elements not in the orbit:" << endl;
			//O.print_orbit(1);
			}
#if 0
		INT f_circletext = TRUE;
		INT xmax = 2000;
		INT ymax = 1000;
		O.draw_tree("tree", 0, xmax, ymax, f_circletext);
#endif
		
		i = O.orbit_first[1];
		j = O.orbit[i];
		cout << "point " << j << " is not in the orbit" << endl;
		A.get_transversal_rep(0, j, elt1, FALSE);
		
		V.append(elt1);
		
		cout << "calling extend_orbit" << endl;
		
		O.extend_orbit(elt1, verbose_level - 1);
		}
	
	cout << "generators: " << endl;
	V.print(cout);
}
#endif

void schreier::point_stabilizer(action *default_action, longinteger_object &go, 
	sims *&Stab, INT orbit_no, INT verbose_level)
// this function allocates a sims structure into Stab.
{
	Stab = new sims;
	longinteger_object cur_go, target_go;
	longinteger_domain D;
	INT len, r, cnt = 0, f_added, *p_gen, drop_out_level, image;
	INT *residue;
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT f_vvv = (verbose_level >= 3);
	INT f_v4 = (verbose_level >= 4);
	//INT f_v5 = (verbose_level >= 5);
	
	
	if (f_v) {
		cout << "schreier::point_stabilizer computing stabilizer of representative of orbit " 
			<< orbit_no << " inside a group of order " << go << " in action ";
		default_action->print_info();
		cout << endl;
		}
	residue = NEW_INT(default_action->elt_size_in_INT);
	len = orbit_len[orbit_no];
	D.integral_division_by_INT(go, len, target_go, r);
	if (r) {	
		cout << "schreier::point_stabilizer orbit length does not divide group order" << endl;
		exit(1);
		}
	if (f_vv) {
		cout << "expecting group of order " << target_go << endl;
		}
	
	Stab->init(default_action);
	Stab->init_trivial_group(verbose_level - 1);
	while (TRUE) {
		Stab->group_order(cur_go);
		if (D.compare(cur_go, target_go) == 0) {
			break;
			}
		if (cnt % 2 || Stab->nb_gen[0] == 0) {
			random_schreier_generator_ith_orbit(orbit_no, 0 /* verbose_level */);
			p_gen = schreier_gen;
			if (f_vvv) {
				cout << "random Schreier generator from the orbit:" << endl;
				default_action->element_print(p_gen, cout);
				}
			}
		else {
			Stab->random_schreier_generator(0 /* verbose_level */);
			p_gen = Stab->schreier_gen;
			if (f_v4) {
				cout << "random schreier generator from sims:" << endl;
				default_action->element_print(p_gen, cout);
				}
			}



		if (Stab->strip(p_gen, residue, drop_out_level, image, 0 /*verbose_level - 3*/)) {
			if (f_vvv) {
				cout << "element strips through" << endl;
				if (f_v4) {
					cout << "residue:" << endl;
					A->element_print(residue, cout);
					cout << endl;
					}
				}
			f_added = FALSE;
			}
		else {
			f_added = TRUE;
			if (f_vvv) {
				cout << "element needs to be inserted at level = " 
					<< drop_out_level << " with image " << image << endl;
				if (FALSE) {
					A->element_print(residue, cout);
					cout  << endl;
					}
				}
			Stab->add_generator_at_level(residue, drop_out_level, verbose_level - 4);
			}
		Stab->group_order(cur_go);
		if ((f_vv && f_added) || f_vvv) {
			cout << "iteration " << cnt << " the new group order is " << cur_go 
				<< " expecting a group of order " << target_go << endl; 
			}
		cnt++;
		}
	FREE_INT(residue);
	if (f_v) {
		cout << "schreier::point_stabilizer finished" << endl;
		}
}

void schreier::get_orbit(INT orbit_idx, INT *set, INT &len, INT verbose_level)
{
	INT f, i;

	f = orbit_first[orbit_idx];
	len = orbit_len[orbit_idx];
	for (i = 0; i < len; i++) {
		set[i] = orbit[f + i];
		}
}

void schreier::compute_orbit_statistic(INT *set, INT set_size, INT *orbit_count, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, a, j, o;

	if (f_v) {
		cout << "schreier::compute_orbit_statistic" << endl;
		}
	INT_vec_zero(orbit_count, nb_orbits);
#if 0
	for (i = 0; i < nb_orbits; i++) {
		orbit_count[i] = 0;
		}
#endif
	for (i = 0; i < set_size; i++) {
		a = set[i];
		j = orbit_inv[a];
		o = orbit_no[j];
		orbit_count[o]++;
		}
	if (f_v) {
		cout << "schreier::compute_orbit_statistic done" << endl;
		}
}

void schreier::test_sv(action *A, INT *hdl_strong_generators, INT *sv, 
	INT f_trivial_group, INT f_compact, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *crep, *Elt1, *Elt2, *Elt3;
	
	crep = NEW_INT(A->elt_size_in_INT);
	Elt1 = NEW_INT(A->elt_size_in_INT);
	Elt2 = NEW_INT(A->elt_size_in_INT);
	Elt3 = NEW_INT(A->elt_size_in_INT);
	INT k, i, j, pt, pt0;
	INT f_check_image = FALSE;
	
	if (f_v) {
		cout << "testing the schreier vector" << endl;
		}
	for (k = 0; k < nb_orbits; k++) {
		for (j = 0; j < orbit_len[k]; j++) {
			i = orbit_first[k] + j;
			pt = orbit[i];
			coset_rep_inv(i);
			schreier_vector_coset_rep_inv(A, sv, hdl_strong_generators, pt, pt0, 
				crep, Elt1, Elt2, Elt3, 
				f_trivial_group, f_compact, f_check_image, verbose_level - 4);
			A->element_invert(crep, Elt1, 0);
			A->element_mult(cosetrep, Elt1, Elt2, 0);
			if (!A->element_is_one(Elt2, 0)) {
				cout << "schreier::test_sv() test fails" << endl;
				exit(1);
				}
			}
		}
	if (f_v) {
		cout << "sv test passed" << endl;
		}
	FREE_INT(crep);
	FREE_INT(Elt1);
	FREE_INT(Elt2);
	FREE_INT(Elt3);
}

void schreier::write_to_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "schreier::write_to_memory_object" << endl;
		}
	m->write_int(A->degree);
	m->write_int(nb_orbits);
	for (i = 0; i < nb_orbits; i++) {
		m->write_int(orbit_first[i]);
		m->write_int(orbit_len[i]);
		}
	for (i = 0; i < A->degree; i++) {
		m->write_int(orbit[i]);
		m->write_int(prev[i]);
		m->write_int(label[i]);
		m->write_int(orbit_no[i]);
		}
	gens.write_to_memory_object(m, verbose_level - 1);
	gens_inv.write_to_memory_object(m, verbose_level - 1);
	if (f_v) {
		cout << "schreier::write_to_memory_object done" << endl;
		}
}

void schreier::read_from_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, deg;

	if (f_v) {
		cout << "schreier::read_from_memory_object" << endl;
		}
	init2();
	m->read_int(&deg);
	m->read_int(&nb_orbits);
	if (deg != A->degree) {
		cout << "schreier::read_from_memory_object deg != A->degree" << endl;
		}
	orbit_first = NEW_INT(nb_orbits);
	orbit_len = NEW_INT(nb_orbits);
	for (i = 0; i < nb_orbits; i++) {
		m->read_int(&orbit_first[i]);
		m->read_int(&orbit_len[i]);
		}
	orbit = NEW_INT(A->degree);
	orbit_inv = NEW_INT(A->degree);
	prev = NEW_INT(A->degree);
	label = NEW_INT(A->degree);
	orbit_no = NEW_INT(A->degree);
	for (i = 0; i < A->degree; i++) {
		m->read_int(&orbit[i]);
		m->read_int(&prev[i]);
		m->read_int(&label[i]);
		m->read_int(&orbit_no[i]);
		}
	perm_inverse(orbit, orbit_inv, A->degree);
	gens.init(A);
	gens.read_from_memory_object(m, verbose_level - 1);
	gens_inv.init(A);
	gens_inv.read_from_memory_object(m, verbose_level - 1);
	if (f_v) {
		cout << "schreier::read_from_memory_object done" << endl;
		}
}

void schreier::write_file(BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	
	if (f_v) {
		cout << "schreier::write_file" << endl;
		}
	M.alloc(1024 /* length */, verbose_level - 1);
	M.used_length = 0;
	M.cur_pointer = 0;
	write_to_memory_object(&M, verbose_level - 1);
	M.write_file(fname, verbose_level - 1);
	if (f_v) {
		cout << "schreier::write_file done" << endl;
		}
}

void schreier::read_file(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	memory_object M;
	
	if (f_v) {
		cout << "schreier::read_file reading file " << fname << " of size " << file_size(fname) << endl;
		}
	M.read_file(fname, verbose_level - 1);
	if (f_v) {
		cout << "schreier::read_file read file " << fname << endl;
		}
	M.cur_pointer = 0;
	read_from_memory_object(&M, verbose_level - 1);
	if (f_v) {
		cout << "schreier::read_file done" << endl;
		}
}

void schreier::write_to_file_binary(ofstream &fp, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "schreier::write_to_file_binary" << endl;
		}
	fp.write((char *) &A->degree, sizeof(INT));
	fp.write((char *) &nb_orbits, sizeof(INT));
	for (i = 0; i < nb_orbits; i++) {
		fp.write((char *) &orbit_first[i], sizeof(INT));
		fp.write((char *) &orbit_len[i], sizeof(INT));
		}
	for (i = 0; i < A->degree; i++) {
		fp.write((char *) &orbit[i], sizeof(INT));
		fp.write((char *) &prev[i], sizeof(INT));
		fp.write((char *) &label[i], sizeof(INT));
		fp.write((char *) &orbit_no[i], sizeof(INT));
		}
	gens.write_to_file_binary(fp, verbose_level - 1);
	gens_inv.write_to_file_binary(fp, verbose_level - 1);
	if (f_v) {
		cout << "schreier::write_to_file_binary done" << endl;
		}
}

void schreier::read_from_file_binary(ifstream &fp, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, deg;

	if (f_v) {
		cout << "schreier::read_from_file_binary" << endl;
		}
	init2();
	fp.read((char *) &deg, sizeof(INT));
	fp.read((char *) &nb_orbits, sizeof(INT));
	if (deg != A->degree) {
		cout << "schreier::read_from_file_binary deg != A->degree" << endl;
		}
	orbit_first = NEW_INT(nb_orbits);
	orbit_len = NEW_INT(nb_orbits);
	for (i = 0; i < nb_orbits; i++) {
		fp.read((char *) &orbit_first[i], sizeof(INT));
		fp.read((char *) &orbit_len[i], sizeof(INT));
		}
	orbit = NEW_INT(A->degree);
	orbit_inv = NEW_INT(A->degree);
	prev = NEW_INT(A->degree);
	label = NEW_INT(A->degree);
	orbit_no = NEW_INT(A->degree);
	for (i = 0; i < A->degree; i++) {
		fp.read((char *) &orbit[i], sizeof(INT));
		fp.read((char *) &prev[i], sizeof(INT));
		fp.read((char *) &label[i], sizeof(INT));
		fp.read((char *) &orbit_no[i], sizeof(INT));
		}
	perm_inverse(orbit, orbit_inv, A->degree);
	
	gens.init(A);
	gens.read_from_file_binary(fp, verbose_level - 1);
	gens_inv.init(A);
	gens_inv.read_from_file_binary(fp, verbose_level - 1);
	if (f_v) {
		cout << "schreier::read_from_file_binary done" << endl;
		}
}


void schreier::write_file_binary(BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "schreier::write_file_binary" << endl;
		}
	{
		ofstream fp(fname, ios::binary);

		write_to_file_binary(fp, verbose_level - 1);
	}
	cout << "schreier::write_file_binary Written file " << fname << " of size " << file_size(fname) << endl;
	if (f_v) {
		cout << "schreier::write_file_binary done" << endl;
		}
}

void schreier::read_file_binary(const BYTE *fname, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "schreier::read_file_binary reading file " << fname << " of size " << file_size(fname) << endl;
		}
	cout << "schreier::read_file_binary Reading file " << fname << " of size " << file_size(fname) << endl;
	{
		ifstream fp(fname, ios::binary);

		read_from_file_binary(fp, verbose_level - 1);
	}
	if (f_v) {
		cout << "schreier::read_file_binary done" << endl;
		}
}

void schreier::orbits_as_set_of_sets(set_of_sets *&S, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *Sz;
	INT i, j, a, f, l;
	
	if (f_v) {
		cout << "schreier::orbits_as_set_of_sets" << endl;
		}
	S = new set_of_sets;
	Sz = NEW_INT(nb_orbits);
	for (i = 0; i < nb_orbits; i++) {
		l = orbit_len[i];
		Sz[i] = l;
		}
	
	S->init_basic(A->degree /* underlying_set_size */, nb_orbits, Sz, 0 /* verbose_level */);
	for (i = 0; i < nb_orbits; i++) {
		f = orbit_first[i];
		l = orbit_len[i];
		for (j = 0; j < l; j++) {
			a = orbit[f + j];
			S->Sets[i][j] = a;
			}
		}
	FREE_INT(Sz);
	if (f_v) {
		cout << "schreier::orbits_as_set_of_sets done" << endl;
		}
}



