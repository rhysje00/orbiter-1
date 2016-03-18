// linear_group_description.C
//
// Anton Betten
// December 25, 2015

#include "galois.h"
#include "action.h"



linear_group_description::linear_group_description()
{
	null();
}

linear_group_description::~linear_group_description()
{
	freeself();
}

void linear_group_description::null()
{
	f_projective = FALSE;
	f_general = FALSE;
	f_affine = FALSE;
	n = 0;
	F = NULL;
	f_semilinear = FALSE;
	f_special = FALSE;
	
	f_wedge_action = FALSE;
	f_PGL2OnConic = FALSE;
	f_monomial_group = FALSE;
	f_null_polarity_group = FALSE;
	f_singer_group = FALSE;
	f_subfield_structure_action = FALSE;
	f_subgroup_from_file = FALSE;
}

void linear_group_description::freeself()
{
	null();
}

void linear_group_description::read_arguments(int argc, const char **argv, 
	INT verbose_level)
{
	INT i;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			continue;
			}

		// the general linear groups:
		// GL, GGL, SL, SSL
		if (strcmp(argv[i], "-GL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = TRUE;
			f_affine = FALSE;
			f_semilinear = FALSE;
			f_special = FALSE;
			cout << "-GL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-GGL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = TRUE;
			f_affine = FALSE;
			f_semilinear = TRUE;
			f_special = FALSE;
			cout << "-GGL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-SL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = TRUE;
			f_affine = FALSE;
			f_semilinear = FALSE;
			f_special = TRUE;
			cout << "-SL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-SSL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = TRUE;
			f_affine = FALSE;
			f_semilinear = TRUE;
			f_special = TRUE;
			cout << "-SSL " << n << " " << input_q << endl;
			}


		// the projective linear groups:
		// PGL, PGGL, PSL, PSSL
		else if (strcmp(argv[i], "-PGL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = TRUE;
			f_general = FALSE;
			f_affine = FALSE;
			f_semilinear = FALSE;
			f_special = FALSE;
			cout << "-PGL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-PGGL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = TRUE;
			f_general = FALSE;
			f_affine = FALSE;
			f_semilinear = TRUE;
			f_special = FALSE;
			cout << "-PGGL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-PSL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = TRUE;
			f_general = FALSE;
			f_affine = FALSE;
			f_semilinear = FALSE;
			f_special = TRUE;
			cout << "-PSL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-PSSL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = TRUE;
			f_general = FALSE;
			f_affine = FALSE;
			f_semilinear = TRUE;
			f_special = TRUE;
			cout << "-PSSL " << n << " " << input_q << endl;
			}



		// the affine groups:
		// AGL, AGGL, ASL, ASSL
		else if (strcmp(argv[i], "-AGL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = FALSE;
			f_affine = TRUE;
			f_semilinear = FALSE;
			f_special = FALSE;
			cout << "-AGL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-AGGL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = FALSE;
			f_affine = TRUE;
			f_semilinear = TRUE;
			f_special = FALSE;
			cout << "-AGGL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-ASL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = FALSE;
			f_affine = TRUE;
			f_semilinear = FALSE;
			f_special = TRUE;
			cout << "-ASL " << n << " " << input_q << endl;
			}
		else if (strcmp(argv[i], "-ASSL") == 0) {
			n = atoi(argv[++i]);
			input_q = atoi(argv[++i]);
			f_projective = FALSE;
			f_general = FALSE;
			f_affine = TRUE;
			f_semilinear = TRUE;
			f_special = TRUE;
			cout << "-ASSL " << n << " " << input_q << endl;
			}

		else if (strcmp(argv[i], "-wedge") == 0) {
			f_wedge_action = TRUE;
			cout << "-wedge" << endl;
			}
		else if (strcmp(argv[i], "-PGL2OnConic") == 0) {
			f_PGL2OnConic = TRUE;
			cout << "-PGL2OnConic" << endl;
			}
		else if (strcmp(argv[i], "-monomial") == 0) {
			f_monomial_group = TRUE;
			cout << "-monomial " << endl;
			}
		else if (strcmp(argv[i], "-null_polarity_group") == 0) {
			f_null_polarity_group = TRUE;
			cout << "-null_polarity_group" << endl;
			}
		else if (strcmp(argv[i], "-singer") == 0) {
			f_singer_group = TRUE;
			cout << "-singer" << endl;
			}
		else if (strcmp(argv[i], "-subfield_structure_action") == 0) {
			f_subfield_structure_action = TRUE;
			s = atoi(argv[++i]);
			cout << "-subfield_structure_action " << s << endl;
			}
		else if (strcmp(argv[i], "-subgroup_from_file") == 0) {
			f_subgroup_from_file = TRUE;
			subgroup_fname = argv[++i];
			subgroup_label = argv[++i];
			cout << "-subgroup_from_file " << subgroup_fname << " " << subgroup_label << endl;
			}
		} // next i
}


