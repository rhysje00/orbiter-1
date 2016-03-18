// generator_trace.C
//
// generator_trace.C
//
// Anton Betten
//
// moved out of generator.C: Jan 21/2010

#include "orbiter.h"

#if 0
trace_result generator::find_automorphism_by_tracing(INT len, 
	INT my_node, INT cur_node, INT cur_extension, INT cur_coset, 
	INT f_debug, INT f_implicit_fusion, INT verbose_level)
// This routine is called from upstep.
// It in turn calls oracle::find_automorphism_by_tracing_recursion
// It tries to compute an isomorphism of the set in set[0][0,...,len] 
// (i.e. of size len+1) to the 
// set in S[0,...,len] which sends set[0][len] to S[len].
// Since set[0][0,...,len] is a permutation of S[0,...,len], this isomorphism is 
// in fact an automorphism which maps S[len] to one of the points in S[0,...,len - 1].
// If this is done for all possible points in S[0,...,len - 1], 
// a transversal for H in the stabilizer of S[0,...,len] results, 
// where H is the point stabilizer of S[len] in the set-stabilizer of S[0,...,len-1], 
// (which is a subgroup of S[0,...,len]).  
{
	trace_result r;
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		print_level_extension_coset_info(len + 1, 
			my_node, cur_extension, cur_coset);
		cout << "generator::find_automorphism_by_tracing()" << endl;
		}
	nb_times_trace++;
	if ((nb_times_trace % 100000) == 0) {
		cout << "generator::find_automorphism_by_tracing() " << endl;
		cout << "nb_times_trace=" << nb_times_trace << endl;
		cout << "nb_times_trace_was_saved=" << nb_times_trace_was_saved << endl;
		}
	if (f_vv) {
		print_level_extension_coset_info(len + 1, 
			my_node, cur_extension, cur_coset);
		INT_vec_print(cout, set[0], len + 1);
		cout << endl;
		if (f_print_function) {
			(*print_function)(cout, len + 1, set[0], print_function_data);
			}
		}
	if (len >= sz) {
		cout << "generator::find_automorphism_by_tracing() len >= sz" << endl;
		exit(1);
		}
	
	INT_vec_copy(len + 1, set[0], set0);
	INT_vec_heapsort(set0, len); // INT_vec_sort(len, set0);  // we keep the last point fix
	
	r = root->find_automorphism_by_tracing_recursion(this, 
		0, 0, len, my_node, cur_extension, cur_coset, 
		f_debug, f_implicit_fusion, verbose_level - 1);
	if (f_v) {
		if (r == found_automorphism) {
			print_level_extension_coset_info(len + 1, 
				my_node, cur_extension, cur_coset);
			cout << "found an automorphism" << endl;
			}
		else if (r == not_canonical) {
			print_level_extension_coset_info(len + 1, 
				my_node, cur_extension, cur_coset);
			cout << "not canonical" << endl;
			}
		else {
			print_level_extension_coset_info(len + 1, 
				my_node, cur_extension, cur_coset);
			cout << "no result" << endl;
			}
		}
	return r;
}
#endif

#if 1
void generator::generator_apply_fusion_element_no_transporter(
	INT cur_level, INT size, INT cur_node, INT cur_ex, 
	INT *set_in, INT *set_out, 
	INT verbose_level)
// Called by upstep_work::handle_extension_fusion_type
{
	INT *Elt1;
	INT *Elt2;
	INT *set_tmp;
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "generator::generator_apply_fusion_element_no_transporter" << endl;
		}

	Elt1 = NEW_INT(A->elt_size_in_INT);
	Elt2 = NEW_INT(A->elt_size_in_INT);
	set_tmp = NEW_INT(size);
	A->element_one(Elt1, 0);
	generator_apply_fusion_element(cur_level, size, cur_node, cur_ex, 
		set_in, set_out, set_tmp, 
		Elt1, Elt2, 
		TRUE /* f_tolerant */, 
		0 /*verbose_level*/);
#if 0
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "generator::generator_apply_fusion_element_no_transporter" << endl;
		}
	A->element_retrieve(root[current_node].E[current_extension].data, Elt1, 0);
	
	A2->map_a_set(S1, S2, len, Elt1, 0);

	INT_vec_heapsort(S2, len);

	//A->element_mult(gen->transporter->ith(lvl + 1), Elt1, gen->transporter->ith(0), 0);
#endif
	FREE_INT(Elt1);
	FREE_INT(Elt2);
	FREE_INT(set_tmp);

	if (f_v) {
		cout << "generator::generator_apply_fusion_element_no_transporter done" << endl;
		}
}
#endif

#if 0
void generator::generator_apply_fusion_element(INT cur_level, INT cur_node, INT size, INT level, 
	INT current_extension, 
	INT *canonical_set, INT *tmp_set, 
	INT *Elt_transporter, INT *tmp_Elt, 
	INT verbose_level)
{
	oracle *O = &root[cur_node];
	INT i;

	A->element_retrieve(O->E[current_extension].data, Elt1, 0);
	
	A2->map_a_set(canonical_set, tmp_set, size, Elt1, 0);

	
	INT_vec_heapsort(tmp_set, level); //INT_vec_sort(level, tmp_set);

	A->element_mult(Elt_transporter, Elt1, tmp_Elt, 0);

	for (i = 0; i < size; i++) {
		canonical_set[i] = tmp_set[i];
		}
	A->element_move(tmp_Elt, Elt_transporter, 0);

}
#endif

INT generator::generator_apply_fusion_element(INT level, INT size, 
	INT current_node, INT current_extension, 
	INT *set_in, INT *set_out, INT *set_tmp, 
	INT *transporter_in, INT *transporter_out, 
	INT f_tolerant, 
	INT verbose_level)
// returns next_node
{
	INT f_v = (verbose_level >= 1);
	INT next_node;
	oracle *O;

	O = &root[current_node];

	if (f_v) {
		cout << "generator::generator_apply_fusion_element current_node=" << current_node << " current_extension=" << current_extension << endl;
		cout << "level=" << level << endl;		
		cout << "applying fusion element to the set ";
		INT_set_print(cout, set_in, size);
		cout << endl;
		}

	A2->element_retrieve(O->E[current_extension].data, Elt1, 0);
	
	if (f_v) {
		cout << "generator::generator_apply_fusion_element applying fusion element" << endl;
		A2->element_print_quick(Elt1, cout);
		cout << "in action " << A2->label << ":" << endl;
		A2->element_print_as_permutation(Elt1, cout);
		cout << "to the set ";
		INT_vec_print(cout, set_in, size);
		cout << endl;
		}
	A2->map_a_set(set_in, set_tmp, size, Elt1, 0);
	if (f_v) {
		cout << "generator::generator_apply_fusion_element the set becomes: ";
		INT_vec_print(cout, set_tmp, size);
		cout << endl;
		}

	A2->element_mult(transporter_in, Elt1, Elt2, 0);
	if (f_v) {
		INT_vec_print(cout, set_in, size);
		cout << endl;
		}
	A2->move(Elt2, transporter_out);

	if (f_on_subspaces) {
		next_node = find_node_for_subspace_by_rank(set_tmp, level + 1, verbose_level - 1);
		INT_vec_copy(set_tmp, set_out, size);
		}
	else {
		INT_vec_heapsort(set_tmp, level + 1);
		INT_vec_copy(set_tmp, set_out, size);
		if (f_v) {
			cout << "generator::generator_apply_fusion_element after sorting: ";
			}
		if (f_v) {
			cout << "generator::generator_apply_fusion_element calling find_oracle_node_for_set: ";
			INT_vec_print(cout, set_out, size);
			cout << endl;
			}

		next_node = find_oracle_node_for_set(level + 1 /*size*/, set_out, f_tolerant, 0);
		// changed A Betten 2/19/2011

		}
	if (f_v) {
		cout << "generator::generator_apply_fusion_element from ";
		INT_vec_print(cout, set_in, size);
		cout << " to ";
		INT_vec_print(cout, set_out, size);
		cout << ", which is node " << next_node << endl;
		cout << "we are done" << endl;
		}
	return next_node;
}


INT generator::trace_set_recursion(INT cur_level, INT cur_node, INT size, INT level, 
	INT *canonical_set, INT *tmp_set1, INT *tmp_set2, 
	INT *Elt_transporter, INT *tmp_Elt1, 
	INT f_implicit_fusion, 
	INT f_tolerant, 
	INT verbose_level)
// called by generator::trace_set
// returns the node in the generator that corresponds to the canonical_set
// or -1 if f_tolerant and the node could not be found
{
	INT f_v = (verbose_level >= 1);
	INT pt, pt0, current_extension, i, t, next_node;
	INT f_failure_to_find_point;
	oracle *O = &root[cur_node];
	
	if (f_v) {
		cout << "generator::trace_set_recursion cur_level = " << cur_level << " cur_node = " << cur_node << " : ";
		INT_vec_print(cout, canonical_set, size);
		cout << endl;
		}
	pt = canonical_set[cur_level];
	if (f_v) {
		cout << "tracing point " << pt << endl;
		}
	if (!O->trace_next_point_in_place(this, 
		cur_level, cur_node, size, 
		canonical_set, tmp_set1,
		Elt_transporter, tmp_Elt1, 
		f_implicit_fusion, 
		f_failure_to_find_point, 
		verbose_level - 1)) {
		
		if (f_v) {
			cout << "generator::trace_set_recursion cur_level = " << cur_level << " cur_node = " << cur_node << " : ";
			cout << "O->trace_next_point_in_place returns FALSE, sorting and restarting" << endl;
			}
		// this can only happen if f_implicit_fusion is TRUE
		// we need to sort and restart the trace:

		INT_vec_heapsort(canonical_set, cur_level + 1);
		
		
		return trace_set_recursion(0, 0, 
			size, level, 
			canonical_set, tmp_set1, tmp_set2, 
			Elt_transporter, tmp_Elt1, 
			f_implicit_fusion, 
			f_tolerant, 
			verbose_level);
		}

	if (f_failure_to_find_point) {
		cout << "generator::trace_set_recursion: f_failure_to_find_point" << endl;
		exit(1);
		}
	pt0 = canonical_set[cur_level];
	if (f_v) {
		cout << "generator::trace_set_recursion cur_level = " << cur_level << " cur_node = " << cur_node << " : ";
		INT_vec_print(cout, canonical_set, size);
		cout << " point " << pt << " has been mapped to " << pt0 << endl;
		}
	current_extension = O->find_extension_from_point(this, pt0, FALSE);
	if (current_extension < 0) {
		cout << "generator::trace_set_recursion: did not find point" << endl;
		exit(1);
		}
	t = O->E[current_extension].type;
	if (t == EXTENSION_TYPE_EXTENSION) {
		// extension node
		next_node = O->E[current_extension].data;
		if (f_v) {
			cout << "generator::trace_set_recursion cur_level = " << cur_level << " cur_node = " << cur_node << " : ";
			INT_vec_print(cout, canonical_set, size);
			cout << " point " << pt << " has been mapped to " << pt0 << " next node is node " << next_node << endl;
			}
		if (cur_level + 1 == level) {
			return next_node;
			}
		else {
			return trace_set_recursion(cur_level + 1, next_node, 
				size, level, canonical_set, tmp_set1, tmp_set2,  
				Elt_transporter, tmp_Elt1, 
				f_implicit_fusion, 
				f_tolerant, 
				verbose_level);
			}
		}
	else if (t == EXTENSION_TYPE_FUSION) {
		// fusion node

#if 0
		generator_apply_fusion_element(cur_level, cur_node, size, level, current_extension, 
			canonical_set, tmp_set, 
			Elt_transporter, tmp_Elt, 
			verbose_level - 1);
#endif
		if (f_v) {
			cout << "generator::trace_set_recursion cur_level = " << cur_level << " cur_node = " << cur_node << " : ";
			cout << "before generator_apply_fusion_element" << endl;
			}
		next_node = generator_apply_fusion_element(cur_level, size, 
			cur_node, current_extension, 
			canonical_set, tmp_set1, tmp_set2, 
			Elt_transporter, tmp_Elt1, 
			f_tolerant, 
			verbose_level);
		if (f_v) {
			cout << "generator::trace_set_recursion cur_level = " << cur_level << " cur_node = " << cur_node << " : ";
			cout << "after generator_apply_fusion_element" << endl;
			}
		if (f_v) {
			cout << "generator::trace_set_recursion cur_level = " << cur_level 
				<< " cur_node = " << cur_node << " : " 
				<< " current_extension = " << current_extension 
				<< " : fusion from ";
			INT_vec_print(cout, canonical_set, size);
			cout << " to ";
			INT_vec_print(cout, tmp_set1, size);
			cout << " : we continue with node " << next_node << endl; 
			cout << endl;
			}

		if (next_node == -1) { // can only happen if f_tolerant is TRUE
			if (f_v) {
				cout << "generator::trace_set_recursion cur_level = " << cur_level 
					<< " cur_node = " << cur_node << " : " 
					<< " current_extension = " << current_extension 
					<< " : fusion from ";
				INT_vec_print(cout, canonical_set, size);
				cout << " to ";
				INT_vec_print(cout, tmp_set1, size);
				cout << "we stop tracing" << endl;
				}
			return -1;
			}
		A->element_move(tmp_Elt1, Elt_transporter, 0);
		for (i = 0; i < size; i++) {
			canonical_set[i] = tmp_set1[i];
			}

		if (cur_level + 1 == level) {
			return next_node;
			}
		else {
			return trace_set_recursion(cur_level + 1, next_node, 
				size, level, canonical_set, tmp_set1, tmp_set2,  
				Elt_transporter, tmp_Elt1, 
				f_implicit_fusion, 
				f_tolerant, 
				verbose_level);
			}
#if 0
		// we need to restart the trace:
		return trace_set_recursion(0, 0, 
			size, level, 
			canonical_set, tmp_set1, tmp_set2, 
			Elt_transporter, tmp_Elt1, 
			f_implicit_fusion, verbose_level);
#endif
		}
	cout << "generator::trace_set_recursion unknown type " << t << endl;
	exit(1);
}

INT generator::trace_set(INT *set, INT size, INT level, 
	INT *canonical_set, INT *Elt_transporter, 
	INT f_implicit_fusion, INT verbose_level)
// called by map_set_to_set_BLT in orbits.C
// returns the case number of the canonical set
{
	INT i, n, case_nb;
	INT f_v = (verbose_level >= 1);
	INT *tmp_set1, *tmp_set2;
	INT *tmp_Elt;

	tmp_set1 = NEW_INT(size);
	tmp_set2 = NEW_INT(size);
	tmp_Elt = NEW_INT(A->elt_size_in_INT);

	if (f_v) {
		cout << "generator::trace_set" << endl;
		cout << "tracing set ";
		INT_vec_print(cout, set, size);	
		cout << endl;
		cout << "level=" << level << endl;
		cout << "f_implicit_fusion=" << f_implicit_fusion << endl;
		}
	
	for (i = 0; i < size; i++) {
		canonical_set[i] = set[i];
		}
	A->element_one(Elt_transporter, 0);
	n = trace_set_recursion(0 /* cur_level */, 0 /* cur_node */,  size, level, 
		canonical_set, tmp_set1, tmp_set2, 
		Elt_transporter, tmp_Elt, 
		f_implicit_fusion, 
		FALSE /*f_tolerant*/, 
		verbose_level);
	case_nb = n - first_oracle_node_at_level[level];
	if (case_nb < 0) {
		cout << "generator::trace_set, case_nb < 0, case_nb = " << case_nb << endl;
		exit(1);
		}
	FREE_INT(tmp_set1);
	FREE_INT(tmp_set2);
	FREE_INT(tmp_Elt);
	return case_nb;
}

INT generator::find_node_for_subspace_by_rank(INT *set, INT len, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *v;
	INT *basis;
	INT *base_cols;
	INT rk, node, i, j, pt;

	if (f_v) {
		cout << "generator::find_node_for_subspace_by_rank for set ";
		INT_vec_print(cout, set, len);
		cout << endl;
		}
	v = tmp_find_node_for_subspace_by_rank1;
	basis = tmp_find_node_for_subspace_by_rank2;
	base_cols = tmp_find_node_for_subspace_by_rank3;
	//v = NEW_INT(vector_space_dimension);
	//basis = NEW_INT(len * vector_space_dimension);
	//base_cols = NEW_INT(vector_space_dimension);
	for (i = 0; i < len; i++) {
		(*unrank_point_func)(basis + i * vector_space_dimension, set[i], 				rank_point_data);		
		}
	rk = F->Gauss_simple(basis, len, vector_space_dimension, base_cols, 0 /* verbose_level */);
	if (rk != len) {
		cout << "generator::find_node_for_subspace_by_rank rk != len" << endl;
		exit(1);
		}
	node = 0;
	for (i = 0; i < len; i++) {
		oracle *O;

		O = &root[node];
		for (j = 0; j < O->nb_extensions; j++) {
			if (O->E[j].type != EXTENSION_TYPE_EXTENSION) {
				continue;
				}
			pt = O->E[j].pt;
			(*unrank_point_func)(v, pt, rank_point_data);		
			if (!F->is_contained_in_subspace(len, vector_space_dimension, basis, base_cols, 
				v, verbose_level)) {
				continue;
				}
			if (f_vv) {
				cout << "generator::find_node_for_subspace_by_rank at node " << node << " extension " << j << " with point " << pt << " to node " << O->E[j].data << endl;
				}
			node = O->E[j].data;
			set[i] = pt;
			break;
			}
		if (j == O->nb_extensions) {
			cout << "generator::find_node_for_subspace_by_rank at node " << node << " fatal, could not find extension" << endl;
			exit(1);
			}
		}
	if (f_v) {
		cout << "generator::find_node_for_subspace_by_rank the canonical set is ";
		INT_vec_print(cout, set, len);
		cout << " at node " << node << endl;
		}
	
	//FREE_INT(v);
	//FREE_INT(basis);
	//FREE_INT(base_cols);
	return node;
}

