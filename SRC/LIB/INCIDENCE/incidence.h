// incidence.h
//
// Anton Betten
//
// started:  August 25, 2007

#include <iostream>
#include <map>
#include <vector>
#include <deque>

typedef class dynamic_memory dynamic_memory;
typedef class partition_backtrack partition_backtrack;
typedef class point_line point_line;
typedef struct plane_data PLANE_DATA;
typedef class tdo_scheme tdo_scheme;
typedef class tdo_data tdo_data;
typedef struct solution_file_data solution_file_data;
typedef class geo_parameter geo_parameter;



// #################################################################################
// dynamic_memory.C:
// #################################################################################



class dynamic_memory {
public:
	INT *ptr;
	INT size_allocated;
	INT size_used;
	dynamic_memory();
	~dynamic_memory();
	void allocate();
	void reallocate(INT required_size, INT f_copy_over);
};



// #################################################################################
// inc_gen_global.C:
// #################################################################################

int ordered_pair_rank(int i, int j, int n);
void ordered_pair_unrank(int &i, int &j, int n, int rk);
int ijk_rank(int i, int j, int k, int n);
void ijk_unrank(int &i, int &j, int &k, int n, int rk);
int largest_binomial2_below(int a2);
int largest_binomial3_below(int a3);
int binomial2(int a);
int binomial3(int a);
int minus_one_if_positive(int i);
void int_vec_bubblesort_increasing(int len, int *p);
int int_vec_search(int *v, int len, int a, int &idx);
void int_vec_print(int *v, int len);
int integer_vec_compare(int *p, int *q, int len);
int int_ij2k(int i, int j, int n);
void int_k2ij(int k, int & i, int & j, int n);


// #################################################################################
// plane.C:
// #################################################################################




struct plane_data {
	INT *points_on_lines; // [nb_pts * (plane_order + 1)]
	INT *line_through_two_points; // [nb_pts * nb_pts]
};


class point_line {
	
public:
	//partition_backtrack *PB;
	partitionstack *P;
	
	INT m, n;
	int *a; // the same as in PB
#if 0
	INT f_joining;
	INT f_point_pair_joining_allocated;
	INT m2; // m choose 2
	INT *point_pair_to_idx; // [m * m]
	INT *idx_to_point_i; // [m choose 2]
	INT *idx_to_point_j; // [m choose 2]
	INT max_point_pair_joining;
	INT *nb_point_pair_joining; // [m choose 2]
	INT *point_pair_joining; // [(m choose 2) * max_point_pair_joining]
	
	INT f_block_pair_joining_allocated;
	INT n2; // n choose 2
	INT *block_pair_to_idx; // [n * n]
	INT *idx_to_block_i; // [n choose 2]
	INT *idx_to_block_j; // [n choose 2]
	INT max_block_pair_joining;
	INT *nb_block_pair_joining; // [n choose 2]z
	INT *block_pair_joining; // [(n choose 2) * max_block_pair_joining]
#endif

	// plane_data:
	INT f_projective_plane;
	INT plane_order; // order = prime ^ exponent
	INT plane_prime;
	INT plane_exponent;
	INT nb_pts;
	INT f_plane_data_computed; 
		// indicats whether or not plane and dual_plane 
		// have been computed by init_plane_data()
	
	PLANE_DATA plane;
	PLANE_DATA dual_plane;

	// data for the coordinatization:
	INT line_x_eq_y;
	INT line_infty;
	INT line_x_eq_0;
	INT line_y_eq_0;
	
	INT quad_I, quad_O, quad_X, quad_Y, quad_C;
	INT *pt_labels;  // [m]
	INT *points;  // [m]
		// pt_labels and points are mutually inverse permutations of {0,1,...,m-1}
		// the affine point (x,y) is labeled as x * plane_order + y

	INT *pts_on_line_x_eq_y;  // [plane_order + 1];
	INT *pts_on_line_x_eq_y_labels;  // [plane_order + 1];
	INT *lines_through_X;  // [plane_order + 1];
	INT *lines_through_Y;  // [plane_order + 1];
	INT *pts_on_line;  // [plane_order + 1];
	INT *MOLS;  // [(plane_order + 1) * plane_order * plane_order]
	INT *field_element; // [plane_order]
	INT *field_element_inv; // [plane_order]


	INT is_desarguesian_plane(INT f_v, INT f_vv);
	INT identify_field_not_of_prime_order(INT f_v, INT f_vv);
	void init_projective_plane(INT order, INT f_v);
	void free_projective_plane();
	void plane_report(ostream &ost);
	INT plane_line_through_two_points(INT pt1, INT pt2);
	INT plane_line_intersection(INT line1, INT line2);
	void plane_get_points_on_line(INT line, INT *pts);
	void plane_get_lines_through_point(INT pt, INT *lines);
	INT plane_points_collinear(INT pt1, INT pt2, INT pt3);
	INT plane_lines_concurrent(INT line1, INT line2, INT line3);
	INT plane_first_quadrangle(INT &pt1, INT &pt2, INT &pt3, INT &pt4);
	INT plane_next_quadrangle(INT &pt1, INT &pt2, INT &pt3, INT &pt4);
	INT plane_quadrangle_first_i(INT *pt, INT i);
	INT plane_quadrangle_next_i(INT *pt, INT i);
	void coordinatize_plane(INT O, INT I, INT X, INT Y, INT *MOLS, INT f_v);
	// needs pt_labels, points, pts_on_line_x_eq_y, pts_on_line_x_eq_y_labels, 
	// lines_through_X, lines_through_Y, pts_on_line, MOLS to be allocated
	INT &MOLSsxb(INT s, INT x, INT b);
	INT &MOLSaddition(INT a, INT b);
	INT &MOLSmultiplication(INT a, INT b);
	INT ternary_field_is_linear(INT *MOLS, INT f_v);
	void print_MOLS(ostream &ost);

	INT is_projective_plane(partitionstack &P, INT &order, INT f_v, INT f_vv);
		// if it is a projective plane, the order is returned.
		// otherwise, 0 is returned.
	INT count_RC(partitionstack &P, INT row_cell, INT col_cell);
	INT count_CR(partitionstack &P, INT col_cell, INT row_cell);
	INT count_RC_representative(partitionstack &P, 
		INT row_cell, INT row_cell_pt, INT col_cell);
	INT count_CR_representative(partitionstack &P, 
		INT col_cell, INT col_cell_pt, INT row_cell);
	INT count_pairs_RRC(partitionstack &P, INT row_cell1, INT row_cell2, INT col_cell);
	INT count_pairs_CCR(partitionstack &P, INT col_cell1, INT col_cell2, INT row_cell);
	INT count_pairs_RRC_representative(partitionstack &P, INT row_cell1, INT row_cell_pt, INT row_cell2, INT col_cell);
		// returns the number of joinings from a point of row_cell1 to elements of row_cell2 within col_cell
		// if that number exists, -1 otherwise
	INT count_pairs_CCR_representative(partitionstack &P, INT col_cell1, INT col_cell_pt, INT col_cell2, INT row_cell);
		// returns the number of joinings from a point of col_cell1 to elements of col_cell2 within row_cell
		// if that number exists, -1 otherwise

};

void get_MOLm(INT *MOLS, INT order, INT m, INT *&M);

// #################################################################################
// tdo_scheme.C:
// #################################################################################


#define MAX_SOLUTION_FILE 100

#define NUMBER_OF_SCHEMES 5
#define ROW 0
#define COL 1
#define LAMBDA 2
#define EXTRA_ROW 3
#define EXTRA_COL 4

struct solution_file_data {
	INT nb_solution_files;
	INT system_no[MAX_SOLUTION_FILE];
	BYTE *solution_file[MAX_SOLUTION_FILE];
};


// #################################################################################
// tdo_data.C: TDO parameter refinement
// #################################################################################


class tdo_data {
	public:
	INT *types_first;
	INT *types_len;
	INT *only_one_type;
	INT nb_only_one_type;
	INT *multiple_types;
	INT nb_multiple_types;
	INT *types_first2;
	diophant *D1;
	diophant *D2;

	tdo_data();
	~tdo_data();
	void free();
	void allocate(int R);
	INT solve_first_system(INT verbose_level, 
		INT *&line_types, INT &nb_line_types, INT &line_types_allocated);
	void solve_second_system_omit(int verbose_level, 
		INT *classes_len, 
		INT *&line_types, INT &nb_line_types, 
		INT *&distributions, INT &nb_distributions,
		INT omit);
	void solve_second_system_with_help(int verbose_level, 
		INT f_use_mckay_solver, INT f_once, 
		INT *classes_len, int f_scale, int scaling, 
		INT *&line_types, INT &nb_line_types, 
		INT *&distributions, INT &nb_distributions,
		INT cnt_second_system, solution_file_data *Sol);
	void solve_second_system_from_file(int verbose_level, 
		INT *classes_len, int f_scale, int scaling, 
		INT *&line_types, INT &nb_line_types, 
		INT *&distributions, INT &nb_distributions, BYTE *solution_file_name);
	void solve_second_system(int verbose_level, INT f_use_mckay_solver, INT f_once, 
		INT *classes_len, int f_scale, int scaling, 
		INT *&line_types, INT &nb_line_types, 
		INT *&distributions, INT &nb_distributions);
};


class tdo_scheme {

public:

	// the following is needed by the TDO process:
	// allocated in init_partition_stack
	// freed in exit_partition_stack
		 
	//partition_backtrack PB;

	partitionstack *P;

	int part_length;
	int *part;
	int nb_entries;
	int *entries;
	int row_level;
	int col_level;
	int lambda_level;
	int extra_row_level;
	int extra_col_level;
		
	int mn; // m + n
	int m; // # of rows
	int n; // # of columns
		
	int level[NUMBER_OF_SCHEMES];
	INT *row_classes[NUMBER_OF_SCHEMES], nb_row_classes[NUMBER_OF_SCHEMES];
	INT *col_classes[NUMBER_OF_SCHEMES], nb_col_classes[NUMBER_OF_SCHEMES];
	INT *row_class_index[NUMBER_OF_SCHEMES];
	INT *col_class_index[NUMBER_OF_SCHEMES];
	INT *row_classes_first[NUMBER_OF_SCHEMES];
	INT *row_classes_len[NUMBER_OF_SCHEMES];
	INT *row_class_no[NUMBER_OF_SCHEMES];
	INT *col_classes_first[NUMBER_OF_SCHEMES];
	INT *col_classes_len[NUMBER_OF_SCHEMES];
	INT *col_class_no[NUMBER_OF_SCHEMES];
		
	int *the_row_scheme;
	int *the_col_scheme;
	int *the_extra_row_scheme;
	int *the_extra_col_scheme;
	INT *the_row_scheme_cur; // [m * nb_col_classes[ROW]]
	INT *the_col_scheme_cur; // [n * nb_row_classes[COL]]
	INT *the_extra_row_scheme_cur; // [m * nb_col_classes[EXTRA_ROW]]
	INT *the_extra_col_scheme_cur; // [n * nb_row_classes[EXTRA_COL]]
		
	// end of TDO process data

	tdo_scheme();
	~tdo_scheme();
	
	void init_part_and_entries(int *part, int *entries, INT verbose_level);
	void init_part_and_entries_INT(INT *part, INT *entries, INT verbose_level);
	void init_TDO(int *Part, int *Entries, 
		int Row_level, int Col_level, int Extra_row_level, int Extra_col_level, 
		int Lambda_level, int verbose_level);
	void exit_TDO();
	void init_partition_stack(int verbose_level);
	void exit_partition_stack();
	void get_partition(int h, int l, INT verbose_level);
	void free_partition(int h);
	void complete_partition_info(int h, INT verbose_level);
	void get_row_or_col_scheme(int h, int l, INT verbose_level);
	void get_column_split_partition(int verbose_level, partitionstack &P);
	void get_row_split_partition(int verbose_level, partitionstack &P);
	void print_all_schemes();
	void print_scheme(int h, int f_v);
	void print_scheme_tex(ostream &ost, int h);
	void print_scheme_tex_fancy(ostream &ost, int h, int f_label, BYTE *label);
	void compute_whether_first_inc_must_be_moved(INT *f_first_inc_must_be_moved, INT verbose_level);
	INT count_nb_inc_from_row_scheme(int verbose_level);
	INT count_nb_inc_from_extra_row_scheme(int verbose_level);


// #################################################################################
// geometric_tests.C: 
// #################################################################################

	INT geometric_test_for_row_scheme(partitionstack &P, 
		INT *point_types, INT nb_point_types, INT point_type_len, 
		INT *distributions, INT nb_distributions, 
		INT f_omit1, INT omit1, INT verbose_level);
	INT geometric_test_for_row_scheme_level_s(partitionstack &P, INT s, 
		INT *point_types, INT nb_point_types, INT point_type_len, 
		INT *distribution, 
		INT *non_zero_blocks, INT nb_non_zero_blocks, 
		INT f_omit1, INT omit1, 
		INT verbose_level);


// #################################################################################
// refine_rows.C: parameter refinement (TDO)
// #################################################################################

	INT refine_rows(int verbose_level, 
		INT f_use_mckay, INT f_once, 
		partitionstack &P, 
		INT *&point_types, INT &nb_point_types, INT &point_type_len, 
		INT *&distributions, INT &nb_distributions, 
		INT &cnt_second_system, solution_file_data *Sol, 
		INT f_omit1, INT omit1, INT f_omit2, INT omit2, 
		INT f_use_packing_numbers, INT f_dual_is_linear_space, INT f_do_the_geometric_test);
	INT refine_rows_easy(int verbose_level, 
		INT *&point_types, INT &nb_point_types, INT &point_type_len, 
		INT *&distributions, INT &nb_distributions, 
		INT &cnt_second_system);
	INT refine_rows_hard(partitionstack &P, int verbose_level, 
		INT f_use_mckay, INT f_once, 
		INT *&point_types, INT &nb_point_types, INT &point_type_len,  
		INT *&distributions, INT &nb_distributions, 
		INT &cnt_second_system, 
		INT f_omit1, INT omit1, INT f_omit, INT omit, 
		INT f_use_packing_numbers, INT f_dual_is_linear_space);
	void row_refinement_L1_L2(partitionstack &P, INT f_omit, INT omit, 
		INT &L1, INT &L2, INT verbose_level);
	INT tdo_rows_setup_first_system(INT verbose_level, 
		tdo_data &T, INT r, partitionstack &P, 
		INT f_omit, INT omit, 
		INT *&point_types, INT &nb_point_types);
	INT tdo_rows_setup_second_system(INT verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit, INT f_use_packing_numbers, INT f_dual_is_linear_space, 
		INT *&point_types, INT &nb_point_types);
	INT tdo_rows_setup_second_system_eqns_joining(INT verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit, INT f_dual_is_linear_space, 
		INT *point_types, INT nb_point_types, 
		INT eqn_offset);
	INT tdo_rows_setup_second_system_eqns_counting(INT verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit, 
		INT *point_types, INT nb_point_types, 
		INT eqn_offset);
	INT tdo_rows_setup_second_system_eqns_packing(INT verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit, 
		INT *point_types, INT nb_point_types,
		INT eqn_start, INT &nb_eqns_used);

// #################################################################################
// refine_columns.C: parameter refinement (TDO)
// #################################################################################

	INT refine_columns(int verbose_level, INT f_once, partitionstack &P, 
		INT *&line_types, INT &nb_line_types, INT &line_type_len, 
		INT *&distributions, INT &nb_distributions, 
		INT &cnt_second_system, solution_file_data *Sol, 
		INT f_omit1, INT omit1, INT f_omit, INT omit, 
		INT f_D1_upper_bound_x0, INT D1_upper_bound_x0, 
		INT f_use_mckay_solver, 
		INT f_use_packing_numbers);
	INT refine_cols_hard(partitionstack &P, int verbose_level, INT f_once, 
		INT *&line_types, INT &nb_line_types, INT &line_type_len,  
		INT *&distributions, INT &nb_distributions, 
		INT &cnt_second_system, solution_file_data *Sol, 
		INT f_omit1, INT omit1, INT f_omit, INT omit, 
		INT f_D1_upper_bound_x0, INT D1_upper_bound_x0, 
		INT f_use_mckay_solver, 
		INT f_use_packing_numbers);
	void column_refinement_L1_L2(partitionstack &P, INT f_omit, INT omit, 
		INT &L1, INT &L2, INT verbose_level);
	INT tdo_columns_setup_first_system(INT verbose_level, 
		tdo_data &T, INT r, partitionstack &P, 
		INT f_omit, INT omit, 
		INT *&line_types, INT &nb_line_types);
	INT tdo_columns_setup_second_system(int verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit,  
		INT f_use_packing_numbers, 
		INT *&line_types, INT &nb_line_types);
	INT tdo_columns_setup_second_system_eqns_joining(int verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit, 
		INT *line_types, INT nb_line_types,
		INT eqn_start);
	void tdo_columns_setup_second_system_eqns_counting(int verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit, 
		INT *line_types, INT nb_line_types,
		INT eqn_start);
	INT tdo_columns_setup_second_system_eqns_upper_bound(int verbose_level, 
		tdo_data &T, partitionstack &P, 
		INT f_omit, INT omit, 
		INT *line_types, INT nb_line_types,
		INT eqn_start, INT &nb_eqns_used);


// #################################################################################
// refine_3design.C: TDO parameter refinement for 3-designs
// #################################################################################


	INT td3_refine_rows(int verbose_level, INT f_once, 
		int lambda3, int block_size, 
		INT *&point_types, INT &nb_point_types, INT &point_type_len,  
		INT *&distributions, INT &nb_distributions);
	INT td3_rows_setup_first_system(int verbose_level, 
		int lambda3, int block_size, int lambda2, 
		tdo_data &T, int r, partitionstack &P, 
		int &nb_vars,int &nb_eqns, 
		INT *&point_types, INT &nb_point_types);
	INT td3_rows_setup_second_system(int verbose_level, 
		int lambda3, int block_size, int lambda2, 
		tdo_data &T, 
		INT nb_vars, INT &Nb_vars, INT &Nb_eqns, 
		INT *&point_types, INT &nb_point_types);
	INT td3_rows_counting_flags(int verbose_level, 
		int lambda3, int block_size, int lambda2, int &S,
		tdo_data &T, 
		INT nb_vars, INT Nb_vars, 
		INT *&point_types, INT &nb_point_types, INT eqn_offset);
	INT td3_refine_columns(int verbose_level, INT f_once, 
		int lambda3, int block_size, int f_scale, int scaling, 
		INT *&line_types, INT &nb_line_types, INT &line_type_len,  
		INT *&distributions, INT &nb_distributions);
	INT td3_columns_setup_first_system(int verbose_level, 
		int lambda3, int block_size, int lambda2, 
		tdo_data &T, int r, partitionstack &P, 
		int &nb_vars,int &nb_eqns, 
		INT *&line_types, INT &nb_line_types);
	INT td3_columns_setup_second_system(int verbose_level, 
		int lambda3, int block_size, int lambda2, int f_scale, int scaling, 
		tdo_data &T, 
		INT nb_vars, INT &Nb_vars, INT &Nb_eqns, 
		INT *&line_types, INT &nb_line_types);
	INT td3_columns_triples_same_class(int verbose_level, 
		int lambda3, int block_size, 
		tdo_data &T, 
		INT nb_vars, INT Nb_vars, 
		INT *&line_types, INT &nb_line_types, INT eqn_offset);
	INT td3_columns_pairs_same_class(int verbose_level, 
		int lambda3, int block_size, int lambda2, 
		tdo_data &T, 
		INT nb_vars, INT Nb_vars, 
		INT *&line_types, INT &nb_line_types, INT eqn_offset);
	INT td3_columns_counting_flags(int verbose_level, 
		int lambda3, int block_size, int lambda2, int &S,
		tdo_data &T, 
		INT nb_vars, INT Nb_vars, 
		INT *&line_types, INT &nb_line_types, INT eqn_offset);
	INT td3_columns_lambda2_joining_pairs_from_different_classes(int verbose_level, 
		int lambda3, int block_size, int lambda2, 
		tdo_data &T, 
		INT nb_vars, INT Nb_vars, 
		INT *&line_types, INT &nb_line_types, INT eqn_offset);
	INT td3_columns_lambda3_joining_triples_2_1(int verbose_level, 
		int lambda3, int block_size, int lambda2, 
		tdo_data &T, 
		INT nb_vars, INT Nb_vars, 
		INT *&line_types, INT &nb_line_types, INT eqn_offset);
	INT td3_columns_lambda3_joining_triples_1_1_1(int verbose_level, 
		int lambda3, int block_size, int lambda2, 
		tdo_data &T, 
		INT nb_vars, INT Nb_vars, 
		INT *&line_types, INT &nb_line_types, INT eqn_offset);


};


// #################################################################################
// packing.C: packing numbers and maxfit numbers
// #################################################################################

INT &TDO_upper_bound(INT i, INT j);
INT &TDO_upper_bound_internal(INT i, INT j);
INT &TDO_upper_bound_source(INT i, INT j);
INT braun_test_single_type(INT v, INT k, INT ak);
INT braun_test_upper_bound(INT v, INT k);
void TDO_refine_init_upper_bounds(INT v_max);
void TDO_refine_extend_upper_bounds(INT new_v_max);
INT braun_test_on_line_type(INT v, INT *type);
INT &maxfit(INT i, INT j);
INT &maxfit_internal(INT i, INT j);
void maxfit_table_init(INT v_max);
void maxfit_table_reallocate(INT v_max);
void maxfit_table_compute();
INT packing_number_via_maxfit(INT n, INT k);

// #################################################################################
// mckay.C: solver for systems of diophantine equations
// #################################################################################

namespace mckay {
// we use the MCKAY algorithm for now...


#include <stdio.h>
#include <math.h>

/* bigger gets more diagnostic output */
//#define VERBOSE 0

#define MCKAY_DEBUG
#define INTERVAL_IN_SECONDS 1

typedef struct {int var,coeff;} term;
typedef vector<term> equation;

class tMCKAY {
public:
  tMCKAY();
	void Init(diophant *lgs, const char *label, int aEqnAnz, int aVarAnz);
	void possolve(vector<int> &lo, vector<int> &hi, 
		vector<equation> &eqn, vector<int> &lorhs, vector<int> &hirhs, 
		vector<int> &neqn, int numeqn, int numvar, int verbose_level);

	INT nb_calls_to_solve;
	INT first_moved;
	INT second_moved;
	const char *problem_label;

protected:
	bool subtract(INT eqn1, equation &e1, int l1, int lors1, int hirs1, INT eqn2, equation &e2, int *pl2, int *plors2, int *phirs2, INT verbose_level);
	void pruneqn(vector<int> &lorhs, vector<int> &hirhs, vector<equation> &eqn, vector<int> &neqn, int numeqn, INT verbose_level);
	void varprune(vector<int> &lo, vector<int> &hi, vector<int> &lorhs, vector<int> &hirhs, vector<equation> &eqn, vector<int> &neqn, int numeqn, INT verbose_level);
	void puteqns(vector<int> &lo, vector<int> &hi, int numvar, vector<int> &lorhs, vector<int> &hirhs, vector<equation> &eqn, vector<int> &neqn, int numeqn);
	int divideeqns(vector<int> &lorhs, vector<int> &hirhs, vector<equation> &eqn, vector<int> &neqn, int numeqn);
	int gcd(int n1,int n2);
	void solve(int level, 
		vector<int> &alo, vector<int> &ahi, 
		vector<bool> &aactive, int numvar, 
		vector<int> &lorhs, vector<int> &hirhs, 
		vector<equation> &eqn, vector<int> &neqn, int numeqn, int verbose_level);
	INT restrict_variables(int level, 
		vector<int> &lo, vector<int> &hi, 
		vector<bool> &active, int numvar, 
		vector<int> &lorhs, vector<int> &hirhs, 
		vector<equation> &eqn, vector<int> &neqn, 
		int numeqn, INT &f_restriction_made, int verbose_level);
	void log_12l(INT current_node, int level);

	int _eqnanz;
	int _varanz;
	vector<bool> unitcoeffs;
	vector<bool> active;
	int rekurs;
	bool _break;

	diophant *D;
	//tLGS *_lgs;

#ifdef MCKAY_DEBUG
	vector<int> range,split,branch;
	INT ticks0;
#endif

};
}


// #################################################################################
// geo_parameter.C:
// #################################################################################


#define MODE_UNDEFINED 0
#define MODE_SINGLE 1
#define MODE_STACK 2

#define UNKNOWNTYPE 0
#define POINTTACTICAL 1
#define BLOCKTACTICAL 2
#define POINTANDBLOCKTACTICAL 3

#define FUSE_TYPE_NONE 0
#define FUSE_TYPE_SIMPLE 1
#define FUSE_TYPE_DOUBLE 2
//#define FUSE_TYPE_MULTI 3
//#define FUSE_TYPE_TDO 4

class geo_parameter {
public:
	INT decomposition_type;
	INT fuse_type;
	INT v, b;
	
	INT mode;
	BYTE label[1000];
	
	// for MODE_SINGLE
	INT nb_V, nb_B;
	INT *V, *B;
	INT *scheme;
	INT *fuse;
	
	// for MODE_STACK
	INT nb_parts, nb_entries;
	
	INT *part;
	INT *entries;
	INT part_nb_alloc;
	INT entries_nb_alloc;
	
	
	//vector<int> part;
	//vector<int> entries;
	
	INT lambda_level;
	INT row_level, col_level;
	INT extra_row_level, extra_col_level;
	
geo_parameter();
~geo_parameter();
void append_to_part(INT a);
void append_to_entries(INT a1, INT a2, INT a3, INT a4);
void write(ofstream &aStream, BYTE *label);
void write_mode_single(ofstream &aStream, BYTE *label);
void write_mode_stack(ofstream &aStream, BYTE *label);
void convert_single_to_stack(INT verbose_level);
INT partition_number_row(INT row_idx);
INT partition_number_col(INT col_idx);
INT input(ifstream &aStream);
INT input_mode_single(ifstream &aStream);
INT input_mode_stack(ifstream &aStream, INT verbose_level);
void init_tdo_scheme(tdo_scheme &G, INT verbose_level);
void print_schemes(tdo_scheme &G);
void print_schemes_tex(tdo_scheme &G);
void print_scheme_tex(ostream &ost, tdo_scheme &G, INT h);
void print_C_source();
void convert_single_to_stack_fuse_simple_pt(INT verbose_level);
void convert_single_to_stack_fuse_simple_bt(INT verbose_level);
void convert_single_to_stack_fuse_double_pt(INT verbose_level);
void cut_off_two_lines(geo_parameter &GP2, 
	int *&part_relabel, int *&part_length, 
	int verbose_level);
void cut_off(geo_parameter &GP2, int w, 
	int *&part_relabel, int *&part_length, 
	int verbose_level);
void copy(geo_parameter &GP2);
void print_schemes();
};

void INT_vec_classify(INT *v, INT len, INT *class_first, INT *class_len, INT &nb_classes);
INT tdo_scheme_get_row_class_length_fused(tdo_scheme &G, INT h, INT class_first, INT class_len);
INT tdo_scheme_get_col_class_length_fused(tdo_scheme &G, INT h, INT class_first, INT class_len);



// incidence_global.C:

INT diophant_solve_first_mckay(diophant *Dio, INT f_once, INT verbose_level);
INT diophant_solve_all_mckay(diophant *Dio, INT &nb_backtrack_nodes, INT verbose_level);
INT diophant_solve_once_mckay(diophant *Dio, INT verbose_level);
INT diophant_solve_next_mckay(diophant *Dio, INT verbose_level);
void diophant_solve_mckay(diophant *Dio, const BYTE *label, INT maxresults, INT &nb_backtrack_nodes, INT &nb_sol, INT verbose_level);
void diophant_solve_mckay_override_minrhs_in_inequalities(diophant *Dio, const BYTE *label, 
	INT maxresults, INT &nb_backtrack_nodes, 
	INT minrhs, INT &nb_sol, INT verbose_level);


