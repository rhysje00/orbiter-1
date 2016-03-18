// linear_group.C
//
// Anton Betten
// December 24, 2015

#include "galois.h"
#include "action.h"



linear_group::linear_group()
{
	null();
}

linear_group::~linear_group()
{
	freeself();
}

void linear_group::null()
{
	description = NULL;
	initial_strong_gens = NULL;
	A_linear = NULL;
	Mtx = NULL;
	f_has_strong_generators = FALSE;
	Strong_gens = NULL;
}

void linear_group::freeself()
{
	null();
}

void linear_group::init(linear_group_description *description, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "linear_group::init" << endl;
		}
	linear_group::description = description;
	n = description->n;
	F = description->F;
	input_q = F->q;
	f_semilinear = description->f_semilinear;


	if (f_v) {
		cout << "linear_group::init initializing projective group" << endl;
		}


	
	initial_strong_gens = new strong_generators;

	if (f_v) {
		cout << "linear_group::init before initial_strong_gens->init_linear_group_from_scratch done" << endl;
		}
	
	initial_strong_gens->init_linear_group_from_scratch(A_linear, 
		F, n, 
		description->f_projective, description->f_general, description->f_affine, 
		description->f_semilinear, description->f_special, 
		verbose_level);


	if (f_v) {
		cout << "linear_group::init initializing initial_strong_gens done" << endl;
		}

#if 0
		create_linear_group(S, A, 
			F, linear_m, 
			linear_f_projective, linear_f_general, linear_f_affine, 
			linear_f_semilinear, linear_f_special, 
			verbose_level - 1);
			// ACTION/action_global.C


	A_linear = new action;
	A_linear->init_projective_group(n, 
		F, 
		f_semilinear, 
		TRUE /*f_basis */, 
		0 /*verbose_level - 1*/);
	if (f_v) {
		cout << "linear_group::init initializing projective group done" << endl;
		}
#endif

	Mtx = A_linear->G.matrix_grp;

	if (description->f_PGL2OnConic) {
		init_PGL2q_OnConic(prefix, verbose_level);
		}
	else if (description->f_wedge_action) {
		init_wedge_action(prefix, verbose_level);
		}
	else if (description->f_monomial_group) {
		init_monomial_group(prefix, verbose_level);
		}
	else if (description->f_null_polarity_group) {
		init_null_polarity_group(prefix, verbose_level);
		}
	else if (description->f_singer_group) {
		init_singer_group(prefix, verbose_level);
		}
	else if (description->f_subfield_structure_action) {
		init_subfield_structure_action(prefix, description->s, verbose_level);
		}
	else if (description->f_subgroup_from_file) {
		init_subgroup_from_file(prefix, 
			description->subgroup_fname, description->subgroup_label, 
			verbose_level);
		}
	else {
		Strong_gens = initial_strong_gens;
		sprintf(prefix, "PGL_%ld_%ld", n, q);
		}
}

void linear_group::init_PGL2q_OnConic(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "linear_group::init_PGL2q_OnConic initializing action of PGL(2,q) on conic" << endl;
		}
	if (!A_linear->f_has_sims) {
		cout << "linear_group::init_PGL2q_OnConic A_linear does not have sims, so we create it" << endl;
		A_linear->create_sims(verbose_level);
		}
	if (!A_linear->f_has_strong_generators) {
		cout << "linear_group::init_PGL2q_OnConic A_linear does not have strong generators" << endl;
		//A_linear->create_sims(verbose_level);
		exit(1);
		}
	A2 = new action;
	A2->induced_action_by_representation_on_conic(A_linear, 
		FALSE /* f_induce_action */, NULL, 
		verbose_level);

	vector_space_dimension = A2->G.Rep->dimension;
	q = input_q;
	Strong_gens = initial_strong_gens; //A_linear->Strong_gens;
	f_has_strong_generators = FALSE;

	if (f_v) {
		cout << "linear_group::init_PGL2q_OnConic vector_space_dimension=" << vector_space_dimension << endl;
		}
	if (f_v) {
		cout << "linear_group::init_PGL2q_OnConic created action of PGL2_on conic:" << endl;
		A2->print_info();
		}
	sprintf(prefix, "PGL2_OnConic_%ld_%ld", n, q);
	if (f_v) {
		cout << "linear_group::init_PGL2q_OnConic created group " << prefix << endl;
		}
}

void linear_group::init_wedge_action(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "linear_group::init_wedge_action initializing wedge action" << endl;
		}
	if (!A_linear->f_has_sims) {
		cout << "linear_group::init_wedge_action A_linear does not have sims, so we create it" << endl;
		A_linear->create_sims(verbose_level);
		}
	if (!A_linear->f_has_strong_generators) {
		cout << "linear_group::init_wedge_action A_linear does not have strong generators" << endl;
		//>create_sims(verbose_level);
		exit(1);
		}
	A2 = new action;
	action_on_wedge_product *AW;

	


	if (f_v) {
		cout << "linear_group::init_wedge_action before induced_wedge_action:" << endl;
		}
	AW = new action_on_wedge_product;

	AW->init(*A_linear, verbose_level);
	
	vector_space_dimension = AW->wedge_dimension;
	q = input_q;
	Strong_gens = initial_strong_gens; //A_linear->Strong_gens;
	f_has_strong_generators = FALSE;


	if (f_v) {
		cout << "linear_group::init_wedge_action vector_space_dimension=" << vector_space_dimension << endl;
		}
		
	A2->induced_action_on_wedge_product(A_linear, 
		AW, 
		FALSE /* f_induce_action */, NULL, 
		verbose_level);
	if (f_v) {
		cout << "linear_group::init_wedge_action created wedge action:" << endl;
		A2->print_info();
		}
	sprintf(prefix, "Wedge_%ld_%ld", n, q);
	if (f_v) {
		cout << "linear_group::init_wedge_action created group " << prefix << endl;
		}
}

void linear_group::init_monomial_group(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);


	if (f_v) {
		cout << "linear_group::init_monomial_group initializing monomial group" << endl;
		}
		
	vector_space_dimension = n;
	q = input_q;
	
	Strong_gens = new strong_generators;
	Strong_gens->generators_for_the_monomial_group(A_linear, 
		Mtx, verbose_level - 1);
	f_has_strong_generators = TRUE;
	
	A2 = A_linear;

	sprintf(prefix, "Monomial_%ld_%ld", n, q);
	if (f_v) {
		cout << "linear_group::init_monomial_group created group " << prefix << endl;
		}

	if (f_v) {
		cout << "linear_group::init_monomial_group done, prefix = " << prefix << endl;
		}
}

void linear_group::init_singer_group(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);


	if (f_v) {
		cout << "linear_group::init_singer_group initializing singer group" << endl;
		}

	vector_space_dimension = n;
	q = input_q;
	
	Strong_gens = new strong_generators;
	Strong_gens->generators_for_the_singer_cycle(A_linear, Mtx, verbose_level - 1);
	f_has_strong_generators = TRUE;
	

	A2 = A_linear;

	sprintf(prefix, "Singer_%ld_%ld", n, q);
	if (f_v) {
		cout << "linear_group::init_singer_group created group " << prefix << endl;
		}

	if (f_v) {
		cout << "linear_group::init_singer_group done, prefix = " << prefix << endl;
		}
}

void linear_group::init_null_polarity_group(BYTE *prefix, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);


	if (f_v) {
		cout << "linear_group::init_null_polarity_group initializing null polarity group" << endl;
		}

	vector_space_dimension = n;
	q = input_q;
	
	Strong_gens = new strong_generators;
	Strong_gens->generators_for_the_null_polarity_group(A_linear, Mtx, verbose_level - 1);
	f_has_strong_generators = TRUE;
	
	A2 = A_linear;

	sprintf(prefix, "NullPolarity_%ld_%ld", n, q);
	if (f_v) {
		cout << "linear_group::init_null_polarity_group created group " << prefix << endl;
		}

	if (f_v) {
		cout << "linear_group::init_null_polarity_group done, prefix = " << prefix << endl;
		}
}

void linear_group::init_subfield_structure_action(BYTE *prefix, INT s, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		cout << "linear_group::init_subfield_structure_action" << endl;
		cout << "s=" << s << endl;
		}

	if (f_v) {
		cout << "linear_group::init_subfield_structure_action before field_reduction" << endl;
		}

	vector_space_dimension = n;
	q = input_q;
	
	Strong_gens = new strong_generators;
	Strong_gens->field_reduction(A_linear, n, s, F, verbose_level - 1);
	//lift_generators_to_subfield_structure(A_linear, P->n + 1, s, P->F, SGens, verbose_level - 1);
	f_has_strong_generators = TRUE;

	A2 = A_linear;

	sprintf(prefix, "Subfield_%ld_%ld_%ld", n, q, s);
	if (f_v) {
		cout << "linear_group::init_subfield_structure_action created group " << prefix << endl;
		}
	
	if (f_v) {
		cout << "linear_group::init_subfield_structure_action done" << endl;
		}
}

void linear_group::init_subgroup_from_file(BYTE *prefix, 
	const BYTE *subgroup_fname, const BYTE *subgroup_label, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	
	if (f_v) {
		cout << "linear_group::init_subgroup_from_file" << endl;
		cout << "fname=" << subgroup_fname << endl;
		cout << "label=" << subgroup_label << endl;
		}

	if (f_v) {
		cout << "linear_group::init_subgroup_from_file before field_reduction" << endl;
		}

	vector_space_dimension = n;
	q = input_q;
	
	Strong_gens = new strong_generators;
	if (f_v) {
		cout << "linear_group::init_subgroup_from_file reading generators from file " << subgroup_fname << endl;
		}

	Strong_gens->read_file(A_linear, subgroup_fname, verbose_level - 1);

	if (f_v) {
		cout << "linear_group::init_subgroup_from_file read generators from file" << endl;
		}

	f_has_strong_generators = TRUE;

	A2 = A_linear;

	sprintf(prefix, "Subgroup_%s_%ld_%ld", subgroup_label, n, q);
	if (f_v) {
		cout << "linear_group::init_subgroup_from_file created group " << prefix << endl;
		}
	
	if (f_v) {
		cout << "linear_group::init_subgroup_from_file done" << endl;
		}
}



