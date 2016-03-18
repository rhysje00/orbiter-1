// galois.h
//
// Anton Betten
//
// started as orbiter:  October 23, 2002
// 2nd version started:  December 7, 2003
// galois started:  August 12, 2005

// History:
//
// added class unipoly_domain: November 16, 2002
// added class finite_field: October 23, 2002
// added class longinteger: October 26, 2002
// added class mp: March 6, 2003
// added partitionstack: July 3 2007
// added class orthogonal: July 9 2007
// added class vector_hashing: October 14, 2008
// added file tensor: Dec 25 2008
// added class grassmann: June 5, 2009
// added class unusual: June 10 2009
// added file memory: June 25, 2009
// added file geometry: July 9, 2009
// added class classify: Oct 31, 2009
// added class grassmann_embedded: Jan 24, 2010
// added class hermitian: March 19, 2010
// added class incidence_structure: June 20, 2010
// added class finite_ring: June 21, 2010
// added class hjelmslev: June 22, 2010
// added class fancy_set: June 29, 2010
// added class norm_tables: Sept 23, 2010 (started 11/28/2008)
// added struct grid_frame: Sept 8, 2011
// added class data_file: Oct 13, 2011
// added class subfield_structure: November 14, 2011
// added class clique_finder: December 13, 2011
// added class colored_graph: October 28, 2012
// added class rainbow_cliques: October 28, 2012
// added class set_of_sets: November 30, 2012
// added class decomposition: December 1, 2012
// added file dlx.C: April 7, 2013
// added class spreadsheet: March 15, 2013
// added class andre_construction andre_construction: June 2, 2013
// added class andre_construction_point_element: June 2, 2013
// added class andre_construction_line_element: June 2, 2013
// added class INT_matrix: October 23, 2013
// added class gl_classes: October 23, 2013
// added class layered_graph: January 6, 2014
// added class graph_layer: January 6, 2014
// added class graph_node: January 6, 2014
// added class INT_vector: August 12, 2014
// added class projective_space (moved here from ACTION): December 31, 2014
// added class buekenhout_metz (moved here from TOP_LEVEL): December 31, 2014
// added class a_domain March 14, 2015
// added class diophant (moved here from INCIDENCE) April 16, 2015
// added class null_polarity_generator December 11, 2015
// added class layered_graph_draw_options December 15, 2015
// added class klein_correspondence January 1, 2016
// added class file_output January 8, 2016

#include <iostream>
#include <fstream>
//#include <sstream>
#include <iomanip>

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <map>
#include <vector>
#include <deque>


using namespace std;

#define SYSTEMUNIX
#undef SYSTEMWINDOWS


#define MEMORY_DEBUG

#define MAGIC_SYNC 762873656L



// define exactly one of the following to match your system:
#undef INT_HAS_2_BYTES
#define INT_HAS_4_BYTES
#undef INT_HAS_8_BYTES

#ifdef INT_HAS_2_BYTES
typedef short INT2;
typedef long INT4;
typedef long INT8; // that's a lie!
typedef unsigned short UINT2;
typedef unsigned long UINT4;
typedef unsigned long UINT8; // that's a lie!
#endif
#ifdef INT_HAS_4_BYTES
typedef short INT2;
typedef int INT4;
typedef long INT8;
typedef unsigned short UINT2;
typedef unsigned int UINT4;
typedef unsigned long UINT8;
#endif
#ifdef INT_HAS_8_BYTES
typedef short INT2; // that's a lie!
typedef short int INT4;
typedef int INT8;
typedef unsigned short UINT2; // that's a lie!
typedef unsigned short int UINT4;
typedef unsigned int UINT8;
#endif


typedef long int INT;
typedef INT *PINT;
typedef INT **PPINT;
typedef unsigned long UINT;
typedef UINT *PUINT;
typedef long LONG;
typedef LONG *PLONG;
typedef unsigned long ULONG;
typedef ULONG *PULONG;
typedef short SHORT;
typedef SHORT *PSHORT;
typedef char BYTE;
typedef BYTE *PBYTE;
typedef unsigned char UBYTE;
typedef UBYTE *PUBYTE;
typedef char SCHAR;
typedef SCHAR *PSCHAR;
typedef float FLOAT;
typedef FLOAT *PFLOAT;
typedef BYTE TSTRING;
typedef int *pint;
typedef void *pvoid;



#define PAGE_LENGTH_LOG 20
#define MAX_PAGE_SIZE_IN_BYTES (5 * 1L << 20)
#define BUFSIZE 100000
#undef DEBUG_PAGE_STORAGE


#define MINIMUM(x, y)   ( ((x) < (y)) ?  (x) : (y) )
#define MAXIMUM(x, y)   ( ((x) > (y)) ?  (x) : (y) )
#define MIN(x, y)   ( ((x) < (y)) ?  (x) : (y) )
#define MAX(x, y)   ( ((x) > (y)) ?  (x) : (y) )
#define ABS(x)      ( ((x) <  0 ) ? (-(x)) : (x) )
#define EVEN(x)     ( ((x) % 2) == 0 )
#define ODD(x)      ( ((x) % 2) == 1 )
#define DOUBLYEVEN(x)     ( ((x) % 4) == 0 )
#define SINGLYEVEN(x)     ( ((x) % 4) == 2 )
#define ONE_BYTE_INT(a) (((a) > -126) && ((a) < 127))
#define ONE_MILLION 1000000
#define ONE_HUNDRED_THOUSAND 100000


#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846264
#endif

typedef class finite_field finite_field;
typedef class longinteger_object longinteger_object;
typedef longinteger_object *plonginteger_object;
typedef class longinteger_domain longinteger_domain;
typedef class rank_checker rank_checker;
typedef class classify classify;
typedef void *unipoly_object;
typedef class unipoly_domain unipoly_domain;
typedef class mp_graphics mp_graphics;
typedef class partitionstack partitionstack;
typedef class orthogonal orthogonal;
typedef class vector_hashing vector_hashing;
typedef class unusual_model unusual_model;
typedef class grassmann grassmann;
typedef class grassmann_embedded grassmann_embedded;
typedef class hermitian hermitian;
typedef class incidence_structure incidence_structure;
typedef class finite_ring finite_ring;
typedef class hjelmslev hjelmslev;
typedef class norm_tables norm_tables;
typedef struct coordinate_frame coordinate_frame;
typedef class data_file data_file;
typedef class subfield_structure subfield_structure;
typedef class clique_finder clique_finder;
typedef class colored_graph colored_graph;
typedef class rainbow_cliques rainbow_cliques;
typedef class set_of_sets set_of_sets;
typedef class decomposition decomposition;
typedef class brick_domain brick_domain;
typedef class spreadsheet spreadsheet;
typedef class andre_construction andre_construction;
typedef class andre_construction_point_element andre_construction_point_element;
typedef class andre_construction_line_element andre_construction_line_element;
typedef class memory_object memory_object;
typedef class tree_node tree_node;
typedef tree_node *ptree_node;
typedef class tree tree;
typedef class INT_matrix INT_matrix;
typedef class gl_classes gl_classes;
typedef class gl_class_rep gl_class_rep;
typedef class matrix_block_data matrix_block_data;
typedef class layered_graph layered_graph;
typedef class graph_layer graph_layer;
typedef class graph_node graph_node;
typedef class INT_vector INT_vector;
typedef class projective_space projective_space;
typedef class buekenhout_metz buekenhout_metz;
typedef class a_domain a_domain;
typedef class diophant dophant;
typedef class null_polarity_generator null_polarity_generator;
typedef class layered_graph_draw_options layered_graph_draw_options;
typedef class klein_correspondence klein_correspondence;
typedef class file_output file_output;


#ifdef MEMORY_DEBUG
#define NEW_int(n) allocate_int(n, __FILE__, __LINE__)
#define NEW_pint(n) allocate_pint(n, __FILE__, __LINE__)
#define NEW_INT(n) allocate_INT(n, __FILE__, __LINE__)
#define NEW_PINT(n) allocate_PINT(n, __FILE__, __LINE__)
#define NEW_PPINT(n) allocate_PPINT(n, __FILE__, __LINE__)
#define NEW_BYTE(n) allocate_BYTE(n, __FILE__, __LINE__)
#define NEW_UBYTE(n) allocate_UBYTE(n, __FILE__, __LINE__)
#define NEW_PBYTE(n) allocate_PBYTE(n, __FILE__, __LINE__)
#define NEW_pvoid(n) allocate_pvoid(n, __FILE__, __LINE__)
//#define NEW_CLASS(type, n) (type *)allocate_OBJECT(new type[n], n, sizeof(type), __FILE__, __LINE__)
#define NEW_OBJECT(type) (type *)allocate_OBJECT(new type, sizeof(type), __FILE__, __LINE__)
#define NEW_OBJECTS(type, n) (type *)allocate_OBJECTS(new type[n], n, sizeof(type), __FILE__, __LINE__)
#define FREE_int(p) free_int(p, __FILE__, __LINE__)
#define FREE_pint(p) free_pint(p, __FILE__, __LINE__)
#define FREE_INT(p) free_INT(p, __FILE__, __LINE__)
#define FREE_PINT(p) free_PINT(p, __FILE__, __LINE__)
#define FREE_PPINT(p) free_PPINT(p, __FILE__, __LINE__)
#define FREE_BYTE(p) free_BYTE(p, __FILE__, __LINE__)
#define FREE_UBYTE(p) free_UBYTE(p, __FILE__, __LINE__)
#define FREE_PBYTE(p) free_PBYTE(p, __FILE__, __LINE__)
#define FREE_pvoid(p) free_pvoid(p, __FILE__, __LINE__)
//#define FREE_CLASS(p) free_OBJECT(p, __FILE__, __LINE__); delete [] p
#define FREE_OBJECT(p) free_OBJECT(p, __FILE__, __LINE__); delete p
#define FREE_OBJECTS(p) free_OBJECTS(p, __FILE__, __LINE__); delete [] p
#else
#define NEW_int(n) new int[n]
#define NEW_pint(n) new pint[n]
#define NEW_INT(n) new INT[n]
#define NEW_PINT(n) new PINT[n]
#define NEW_PPINT(n) new PPINT[n]
#define NEW_BYTE(n) new BYTE[n]
#define NEW_UBYTE(n) new UBYTE[n]
#define NEW_PBYTE(n) new PBYTE[n]
#define NEW_pvoid(n) new pvoid[n]
#define NEW_CLASS(n, type) new type[n]
#define FREE_int(p) delete [] p
#define FREE_pint(p) delete [] p
#define FREE_INT(p) delete [] p
#define FREE_PINT(p) delete [] p
#define FREE_PPINT(p) delete [] p
#define FREE_BYTE(p) delete [] p
#define FREE_UBYTE(p) delete [] p
#define FREE_PBYTE(p) delete [] p
#define FREE_pvoid(p) delete [] p
//#define FREE_CLASS(p) delete [] p
#define FREE_OBJECT(p) delete [] p
#endif



enum diophant_equation_type {
	t_EQ, 
	t_LE,
	t_ZOR
}; 

typedef enum diophant_equation_type diophant_equation_type;



// ####################################################################################
// fancy_set.C:
// ####################################################################################

class fancy_set {
	
	public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;

	INT n;
	INT k;
	INT *set;
	INT *set_inv;

	fancy_set();
	~fancy_set();
	void null();
	void freeself();
	void init(INT n, INT verbose_level);
	void init_with_set(INT n, INT k, INT *subset, INT verbose_level);
	void print();
	void println();
	void swap(INT pos, INT a);
	INT is_contained(INT a);
	void copy_to(fancy_set *to);
	void add_element(INT elt);
	void add_elements(INT *elts, INT nb);
	void delete_elements(INT *elts, INT nb);
	void delete_element(INT elt);
	void select_subset(INT *elts, INT nb);
	void intersect_with(INT *elts, INT nb);
	void subtract_set(fancy_set *set_to_subtract);
	void sort();
	INT compare_lexicographically(fancy_set *second_set);
	void complement(fancy_set *compl_set);
	INT is_subset(fancy_set *set2);
	INT is_equal(fancy_set *set2);

};

// ####################################################################################
// finite_ring.C:
// ####################################################################################

class finite_ring {

	INT *add_table; // [q * q]
	INT *mult_table; // [q * q]

	INT *f_is_unit_table;
	INT *negate_table;
	INT *inv_table;

	public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;

	INT q;
	INT p;
	INT e;

	finite_field *Fp;


	finite_ring();
	~finite_ring();
	void null();
	void freeself();
	void init(INT q, INT verbose_level);
	INT zero();
	INT one();
	INT is_zero(INT i);
	INT is_one(INT i);
	INT is_unit(INT i);
	INT add(INT i, INT j);
	INT mult(INT i, INT j);
	INT negate(INT i);
	INT inverse(INT i);
	INT Gauss_INT(INT *A, INT f_special, INT f_complete, INT *base_cols, 
		INT f_P, INT *P, INT m, INT n, INT Pn, INT verbose_level);
		// returns the rank which is the number of entries in base_cols
		// A is a m x n matrix,
		// P is a m x Pn matrix (if f_P is TRUE)
};

// ####################################################################################
// finite_field.C:
// ####################################################################################

class finite_field {

private:
	INT f_has_table;
	INT *add_table; // [q * q]
	INT *mult_table; // [q * q]
		// add_table and mult_table are needed in mindist

	INT *negate_table;
	INT *inv_table;
	INT *frobenius_table; // x \mapsto x^p
	INT *absolute_trace_table;
	INT *log_alpha_table;
	INT *alpha_power_table;
	INT *v1, *v2, *v3; // vectors of length e.
	BYTE *symbol_for_print;

public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;
	
	const BYTE *override_poly;
	BYTE *polynomial; // the actual polynomial we consider as integer (in text form)
	INT q, p, e;
	INT alpha; // primitive element
	INT log10_of_q; // needed for printing purposes
	INT f_print_as_exponentials;
	
	finite_field();
	void null();
	~finite_field();
	void init(INT q, INT verbose_level);
	void init_symbol_for_print(const BYTE *symbol);
	void init_override_polynomial(INT q, const BYTE *poly, INT verbose_level);
	void print_minimum_polynomial(INT p, const BYTE *polynomial);
	INT compute_subfield_polynomial(INT order_subfield, INT verbose_level);
	void compute_subfields(INT verbose_level);
	void create_alpha_table(INT verbose_level);
	void create_alpha_table_extension_field(INT verbose_level);
	void create_alpha_table_prime_field(INT verbose_level);
	void create_tables_prime_field(INT verbose_level);
	void create_tables_extension_field(INT verbose_level);
	void print(INT f_add_mult_table);
	void print_add_mult_tables();
	void print_tables();
	void print_tables_extension_field(const BYTE *poly);
	void display_T2(ostream &ost);
	void display_T3(ostream &ost);
	void display_N2(ostream &ost);
	void display_N3(ostream &ost);
	void print_integer_matrix_zech(ostream &ost, INT *p, INT m, INT n);

	INT *private_add_table();
	INT *private_mult_table();
	INT zero();
	INT one();
	INT is_zero(INT i);
	INT is_one(INT i);
	INT mult(INT i, INT j);
	INT mult3(INT a1, INT a2, INT a3);
	INT product3(INT a1, INT a2, INT a3);
	INT product4(INT a1, INT a2, INT a3, INT a4);
	INT product5(INT a1, INT a2, INT a3, INT a4, INT a5);
	INT square(INT a);
	INT twice(INT a);
	INT four_times(INT a);
	INT add(INT i, INT j);
	INT add3(INT i1, INT i2, INT i3);
	INT add4(INT i1, INT i2, INT i3, INT i4);
	INT add5(INT i1, INT i2, INT i3, INT i4, INT i5);
	INT add6(INT i1, INT i2, INT i3, INT i4, INT i5, INT i6);
	INT add7(INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7);
	INT add8(INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8);
	INT negate(INT i);
	INT inverse(INT i);
	INT power(INT a, INT n); // computes a^n
	INT frobenius_power(INT a, INT i); // computes a^{p^i}
	INT absolute_trace(INT i);
	INT absolute_norm(INT i);
	INT alpha_power(INT i);
	INT log_alpha(INT i);
	INT primitive_root();
	INT N2(INT a);
	INT N3(INT a);
	INT T2(INT a);
	INT T3(INT a);
	INT bar(INT a);
	void abc2xy(INT a, INT b, INT c, INT &x, INT &y, INT verbose_level);
	// given a, b, c, determine x and y such that 
	// c = a * x^2 + b * y^2
	// such elements x and y exist for any choice of a, b, c.
	INT retract(finite_field &subfield, INT index, INT a, INT verbose_level);
	void retract_INT_vec(finite_field &subfield, INT index, INT *v_in, INT *v_out, INT len, INT verbose_level);
	INT embed(finite_field &subfield, INT index, INT b, INT verbose_level);
	void subfield_embedding_2dimensional(finite_field &subfield, 
		INT *&components, INT *&embedding, INT *&pair_embedding, INT verbose_level);
	void print_embedding(finite_field &subfield, 
		INT *components, INT *embedding, INT *pair_embedding);
		// we think of F as two dimensional vector space 
		// over f with basis (1,alpha)
		// for i,j \in f, with x = i + j * alpha \in F, we have 
		// pair_embedding[i * q + j] = x;
		// also, 
		// components[x * 2 + 0] = i;
		// components[x * 2 + 1] = j;
		// also, for i \in f, embedding[i] is the element 
		// in F that corresponds to i 
		// components[Q * 2]
		// embedding[q]
		// pair_embedding[q * q]
	void print_embedding_tex(finite_field &subfield, 
		INT *components, INT *embedding, INT *pair_embedding);
	void print_indicator_square_nonsquare(INT a);
	void print_element(ostream &ost, INT a);
	void print_element_with_symbol(ostream &ost, INT a, INT f_exponential, INT width, const BYTE *symbol);
	void INT_vec_print(ostream &ost, INT *v, INT len);
	void INT_vec_print_elements_exponential(ostream &ost, INT *v, INT len, const BYTE *symbol_for_print);
	void latex_addition_table(ostream &f, INT f_elements_exponential, const BYTE *symbol_for_print);
	void latex_multiplication_table(ostream &f, INT f_elements_exponential, const BYTE *symbol_for_print);
	void latex_matrix(ostream &f, INT f_elements_exponential, const BYTE *symbol_for_print, INT *M, INT m, INT n);
	void power_table(INT t, INT *power_table, INT len);
	INT evaluate_conic_form(INT *six_coeffs, INT *v3);
	INT evaluate_quadric_form_in_PG_three(INT *ten_coeffs, INT *v4);
	INT Pluecker_12(INT *x4, INT *y4);
	INT Pluecker_21(INT *x4, INT *y4);
	INT Pluecker_13(INT *x4, INT *y4);
	INT Pluecker_31(INT *x4, INT *y4);
	INT Pluecker_14(INT *x4, INT *y4);
	INT Pluecker_41(INT *x4, INT *y4);
	INT Pluecker_23(INT *x4, INT *y4);
	INT Pluecker_32(INT *x4, INT *y4);
	INT Pluecker_24(INT *x4, INT *y4);
	INT Pluecker_42(INT *x4, INT *y4);
	INT Pluecker_34(INT *x4, INT *y4);
	INT Pluecker_43(INT *x4, INT *y4);
	INT Pluecker_ij(INT i, INT j, INT *x4, INT *y4);
	INT evaluate_symplectic_form(INT len, INT *x, INT *y);
	INT is_totally_isotropic_wrt_symplectic_form(INT k, INT n, INT *Basis);
	void cheat_sheet(ostream &f, INT verbose_level);
	void cheat_sheet_top(ostream &f, INT nb_cols);
	void cheat_sheet_bottom(ostream &f);

	// ####################################################################################
	// finite_field_linear_algebra.C:
	// ####################################################################################

	void copy_matrix(INT *A, INT *B, INT ma, INT na);
	void reverse_matrix(INT *A, INT *B, INT ma, INT na);
	void identity_matrix(INT *A, INT n);
	INT is_identity_matrix(INT *A, INT n);
	INT is_diagonal_matrix(INT *A, INT n);
	INT is_scalar_multiple_of_identity_matrix(INT *A, INT n, INT &scalar);
	void diagonal_matrix(INT *A, INT n, INT alpha);
	void mult_matrix(INT *A, INT *B, INT *C, INT ma, INT na, INT nb);
	void mult_vector_from_the_left(INT *v, INT *A, INT *vA, INT m, INT n);
		// v[m], A[m][n], vA[n]
	void mult_vector_from_the_right(INT *A, INT *v, INT *Av, INT m, INT n);
		// A[m][n], v[n], Av[m]
	void mult_matrix_matrix_verbose(INT *A, INT *B, INT *C, INT m, INT n, INT o, INT verbose_level);
	void mult_matrix_matrix(INT *A, INT *B, INT *C, INT m, INT n, INT o);
		// multiplies C := A * B, where A is m x n and B is n x o, so that C is m by o
		// C must already be allocated
	void semilinear_matrix_mult(INT *A, INT *B, INT *AB, INT n);
	void matrix_mult_affine(INT *A, INT *B, INT *AB, INT n, INT verbose_level);
	void semilinear_matrix_mult_affine(INT *A, INT *B, INT *AB, INT n);
	INT matrix_determinant(INT *A, INT n, INT verbose_level);
	void matrix_inverse(INT *A, INT *Ainv, INT n, INT verbose_level);
	void matrix_invert(INT *A, INT *Tmp, INT *Tmp_basecols, INT *Ainv, INT n, INT verbose_level);
		// Tmp points to n * n + 1 INT's
		// Tmp_basecols points to n INT's
	void semilinear_matrix_invert(INT *A, INT *Tmp, INT *Tmp_basecols, INT *Ainv, INT n, INT verbose_level);
		// Tmp points to n * n + 1 INT's
		// Tmp_basecols points to n INT's
	void semilinear_matrix_invert_affine(INT *A, INT *Tmp, INT *Tmp_basecols, INT *Ainv, INT n, INT verbose_level);
		// Tmp points to n * n + 1 INT's
		// Tmp_basecols points to n INT's
	void matrix_invert_affine(INT *A, INT *Tmp, INT *Tmp_basecols, INT *Ainv, INT n, INT verbose_level);
		// Tmp points to n * n + 1 INT's
		// Tmp_basecols points to n INT's
	void projective_action_from_the_right(INT f_semilinear, INT *v, INT *A, INT *vA, INT n, INT verbose_level);
	void general_linear_action_from_the_right(INT f_semilinear, INT *v, INT *A, INT *vA, INT n, INT verbose_level);
	void semilinear_action_from_the_right(INT *v, INT *A, INT *vA, INT n);
		// vA = (v * A)^{p^f} 
	void semilinear_action_from_the_left(INT *A, INT *v, INT *Av, INT n);
		// Av = A * v^{p^f}
	void affine_action_from_the_right(INT f_semilinear, INT *v, INT *A, INT *vA, INT n);
		// vA = (v * A)^{p^f} + b
	void zero_vector(INT *A, INT m);
	void all_one_vector(INT *A, INT m);
	void support(INT *A, INT m, INT *&support, INT &size);
	void characteristic_vector(INT *A, INT m, INT *set, INT size);
	INT is_zero_vector(INT *A, INT m);
	void add_vector(INT *A, INT *B, INT *C, INT m);
	void negate_vector_in_place(INT *A, INT m);
	void scalar_multiply_vector_in_place(INT c, INT *A, INT m);
	void vector_frobenius_power_in_place(INT *A, INT m, INT f);
	INT dot_product(INT len, INT *v, INT *w);
	void transpose_matrix(INT *A, INT *At, INT ma, INT na);
	void transpose_matrix_in_place(INT *A, INT m);
	void invert_matrix(INT *A, INT *A_inv, INT n);
	void transform_form_matrix(INT *A, INT *Gram, INT *new_Gram, INT d);
		// computes new_Gram = A * Gram * A^\top
	INT rank_of_matrix(INT *A, INT m, INT verbose_level);
	INT rank_of_rectangular_matrix(INT *A, INT m, INT n, INT verbose_level);
	INT rank_and_basecols(INT *A, INT m, INT *base_cols, INT verbose_level);
	void Gauss_step(INT *v1, INT *v2, INT len, INT idx, INT verbose_level);
		// afterwards: v2[idx] = 0 and v1,v2 span the same space as before
		// v1 is not changed if v1[idx] is nonzero
	void Gauss_step_make_pivot_one(INT *v1, INT *v2, 
		INT len, INT idx, INT verbose_level);
		// afterwards: v2[idx] = 0 and v1,v2 span the same space as before
		// v1[idx] is zero
	INT base_cols_and_embedding(INT m, INT n, INT *A, 
		INT *base_cols, INT *embedding, INT verbose_level);
		// returns the rank rk of the matrix.
		// It also computes base_cols[rk] and embedding[m - rk]
		// It leaves A unchanged
	INT Gauss_easy(INT *A, INT m, INT n);
		// returns the rank
	INT Gauss_simple(INT *A, INT m, INT n, INT *base_cols, INT verbose_level);
		// returns the rank which is the number of entries in base_cols
	INT Gauss_INT(INT *A, INT f_special, INT f_complete, INT *base_cols, 
		INT f_P, INT *P, INT m, INT n, INT Pn, INT verbose_level);
		// returns the rank which is the number of entries in base_cols
		// A is m x n,
		// P is m x Pn (provided f_P is TRUE)
	INT Gauss_INT_with_pivot_strategy(INT *A, 
		INT f_special, INT f_complete, INT *pivot_perm, 
		INT m, INT n, 
		INT (*find_pivot_function)(INT *A, INT m, INT n, INT r, INT *pivot_perm, void *data),
		void *find_pivot_data,  
		INT verbose_level);
		// returns the rank which is the number of entries in pivots
		// A is a m x n matrix
	void Gauss_INT_with_given_pivots(INT *A, 
		INT f_special, INT f_complete, INT *pivots, INT nb_pivots, 
		INT m, INT n, 
		INT verbose_level);
		// A is a m x n matrix
	void kernel_columns(INT n, INT nb_base_cols, INT *base_cols, INT *kernel_cols);
	void matrix_get_kernel_as_INT_matrix(INT *M, INT m, INT n, INT *base_cols, INT nb_base_cols, 
		INT_matrix *kernel);
	void matrix_get_kernel(INT *M, INT m, INT n, INT *base_cols, INT nb_base_cols, 
		INT &kernel_m, INT &kernel_n, INT *kernel);
		// kernel must point to the appropriate amount of memory! 
		// (at least n * (n - nb_base_cols) INT's)
	INT perp(INT n, INT k, INT *A, INT *Gram);
	INT perp_standard(INT n, INT k, INT *A, INT verbose_level);
	INT perp_standard_with_temporary_data(INT n, INT k, INT *A, 
		INT *B, INT *K, INT *base_cols, 
		INT verbose_level);
	INT intersect_subspaces(INT n, INT k1, INT *A, INT k2, INT *B, 
		INT &k3, INT *intersection, INT verbose_level);
	INT n_choose_k_mod_p(INT n, INT k, INT verbose_level);
	void Dickson_polynomial(INT *map, INT *coeffs);
		// compute the coefficients of a degree q-1 polynomial which interpolates a given map
		// from F_q to F_q
	void projective_action_on_columns_from_the_left(INT *A, INT *M, INT m, INT n, INT *perm, INT verbose_level);
	void builtin_transversal_rep_GLnq(INT *A, INT n, INT f_semilinear, INT i, INT j, INT verbose_level);
	void affine_translation(INT n, INT coordinate_idx, INT field_base_idx, INT *perm);
		// perm points to q^n INT's
		// field_base_idx is the base element whose translation we compute, 0 \le field_base_idx < e
		// coordinate_idx is the coordinate in which we shift, 0 \le coordinate_idx < n
	void affine_multiplication(INT n, INT multiplication_order, INT *perm);
		// perm points to q^n INT's
		// compute the diagonal multiplication by alpha, i.e. 
		// the multiplication by alpha of each component
	void affine_frobenius(INT n, INT k, INT *perm);
		// perm points to q^n INT's
		// compute the diagonal action of the Frobenius automorphism to the power k, i.e., 
		// raises each component to the p^k-th power
	INT all_affine_translations_nb_gens(INT n);
	void all_affine_translations(INT n, INT *gens);
	void affine_generators(INT n, INT f_translations, 
		INT f_semilinear, INT frobenius_power, 
		INT f_multiplication, INT multiplication_order, 
		INT &nb_gens, INT &degree, INT *&gens, 
		INT &base_len, INT *&the_base);
	INT evaluate_bilinear_form(INT n, INT *v1, INT *v2, INT *Gram);
	INT evaluate_standard_hyperbolic_bilinear_form(INT n, INT *v1, INT *v2);
	INT evaluate_quadratic_form(INT n, INT nb_terms, 
		INT *i, INT *j, INT *coeff, INT *x);
	void find_singular_vector_brute_force(INT n, INT form_nb_terms, 
		INT *form_i, INT *form_j, INT *form_coeff, INT *Gram, 
		INT *vec, INT verbose_level);
	void find_singular_vector(INT n, INT form_nb_terms, 
		INT *form_i, INT *form_j, INT *form_coeff, INT *Gram, 
		INT *vec, INT verbose_level);
	void complete_hyperbolic_pair(INT n, INT form_nb_terms, 
		INT *form_i, INT *form_j, INT *form_coeff, INT *Gram, 
		INT *vec1, INT *vec2, INT verbose_level);
	void find_hyperbolic_pair(INT n, INT form_nb_terms, 
		INT *form_i, INT *form_j, INT *form_coeff, INT *Gram, 
		INT *vec1, INT *vec2, INT verbose_level);
	void restrict_quadratic_form_list_coding(INT k, INT n, INT *basis, 
		INT form_nb_terms, INT *form_i, INT *form_j, INT *form_coeff, 
		INT &restricted_form_nb_terms, 
		INT *&restricted_form_i, INT *&restricted_form_j, 
		INT *&restricted_form_coeff, 
		INT verbose_level);
	void restrict_quadratic_form(INT k, INT n, INT *basis, INT *C, INT *D, INT verbose_level);
	INT compare_subspaces_ranked(INT *set1, INT *set2, INT size, 
		INT vector_space_dimension, INT verbose_level);
		// Compares the span of two sets of vectors.
		// returns 0 if equal, 1 if not
		// (this is so that it matches to the result of a compare function)
	INT compare_subspaces_ranked_with_unrank_function(
		INT *set1, INT *set2, INT size, 
		INT vector_space_dimension, 
		void (*unrank_point_func)(INT *v, INT rk, void *data), 
		void *rank_point_data, 
		INT verbose_level);
	INT Gauss_canonical_form_ranked(INT *set1, INT *set2, INT size, 
		INT vector_space_dimension, INT verbose_level);
		// Computes the Gauss canonical form for the generating set in set1.
		// The result is written to set2.
		// Returns the rank of the span of the elements in set1.
	INT lexleast_canonical_form_ranked(INT *set1, INT *set2, INT size, 
		INT vector_space_dimension, INT verbose_level);
		// Computes the lexleast generating set the subspace spanned by the elements in set1.
		// The result is written to set2.
		// Returns the rank of the span of the elements in set1.
	void reduce_mod_subspace_and_get_coefficient_vector(
		INT k, INT len, INT *basis, INT *base_cols, 
		INT *v, INT *coefficients, INT verbose_level);
	void reduce_mod_subspace(INT k, INT len, INT *basis, INT *base_cols, 
		INT *v, INT verbose_level);
	INT is_contained_in_subspace(INT k, INT len, INT *basis, INT *base_cols, 
		INT *v, INT verbose_level);
	void code_projective_weight_enumerator(INT n, INT k, 
		INT *code, // [k * n]
		INT *weight_enumerator, // [n + 1]
		INT verbose_level);
	void code_weight_enumerator(INT n, INT k, 
		INT *code, // [k * n]
		INT *weight_enumerator, // [n + 1]
		INT verbose_level);
	void code_weight_enumerator_fast(INT n, INT k, 
		INT *code, // [k * n]
		INT *weight_enumerator, // [n + 1]
		INT verbose_level);
	void code_projective_weights(INT n, INT k, 
		INT *code, // [k * n]
		INT *&weights, // will be allocated [N] where N = theta_{k-1}
		INT verbose_level);
	INT is_subspace(INT d, INT dim_U, INT *Basis_U, INT dim_V, INT *Basis_V, INT verbose_level);
	void Kronecker_product(INT *A, INT *B, 
		INT n, INT *AB);
	INT dependency(INT d, INT *v, INT *A, INT m, INT *rho, INT verbose_level);
		// Lueneburg~\cite{Lueneburg87a} p. 104.
		// A is a matrix of size d + 1 times d
		// v[d]
		// rho is a column permutation of degree d
	void order_ideal_generator(INT d, INT idx, INT *mue, INT &mue_deg, 
		INT *A, INT *Frobenius, 
		INT verbose_level);
		// Lueneburg~\cite{Lueneburg87a} p. 105.
		// Frobenius is a matrix of size d x d
		// A is (d + 1) x d
		// mue[d + 1]
	void span_cyclic_module(INT *A, INT *v, INT n, INT *Mtx, INT verbose_level);
	void random_invertible_matrix(INT *M, INT k, INT verbose_level);
	void make_all_irreducible_polynomials_of_degree_d(INT d, INT &nb, INT *&Table, INT verbose_level);
	INT count_all_irreducible_polynomials_of_degree_d(INT d, INT verbose_level);
	void choose_vector_in_here_but_not_in_here_column_spaces(INT_matrix *V, INT_matrix *W, INT *v, INT verbose_level);
	void choose_vector_in_here_but_not_in_here_or_here_column_spaces(INT_matrix *V, INT_matrix *W1, INT_matrix *W2, INT *v, INT verbose_level);
	INT choose_vector_in_here_but_not_in_here_or_here_column_spaces_coset(INT &coset, 
		INT_matrix *V, INT_matrix *W1, INT_matrix *W2, INT *v, INT verbose_level);
	void vector_add_apply(INT *v, INT *w, INT c, INT n);
	void vector_add_apply_with_stride(INT *v, INT *w, INT stride, INT c, INT n);
	INT test_if_commute(INT *A, INT *B, INT k, INT verbose_level);
	void unrank_point_in_PG(INT *v, INT len, INT rk);
		// len is the length of the vector, not the projective dimension
	INT rank_point_in_PG(INT *v, INT len);
	INT nb_points_in_PG(INT n);
		// n is projective dimension

	// ####################################################################################
	// finite_field_representations.C:
	// ####################################################################################

	void representing_matrix8_R(INT *A, INT q, INT a, INT b, INT c, INT d);
	void representing_matrix9_R(INT *A, INT q, INT a, INT b, INT c, INT d);
	void representing_matrix9_U(INT *A, INT a, INT b, INT c, INT d, INT beta);
	void representing_matrix8_U(INT *A, INT a, INT b, INT c, INT d, INT beta);
	void representing_matrix8_V(INT *A, INT beta);
	void representing_matrix9b(INT *A, INT beta);
	void representing_matrix8a(INT *A, INT a, INT b, INT c, INT d, INT beta);
	void representing_matrix8b(INT *A, INT beta);
	INT Term1(INT a1, INT e1);
	INT Term2(INT a1, INT a2, INT e1, INT e2);
	INT Term3(INT a1, INT a2, INT a3, INT e1, INT e2, INT e3);
	INT Term4(INT a1, INT a2, INT a3, INT a4, INT e1, INT e2, INT e3, INT e4);
	INT Term5(INT a1, INT a2, INT a3, INT a4, INT a5, INT e1, INT e2, INT e3, INT e4, INT e5);
	INT term1(INT a1, INT e1);
	INT term2(INT a1, INT a2, INT e1, INT e2);
	INT term3(INT a1, INT a2, INT a3, INT e1, INT e2, INT e3);
	INT term4(INT a1, INT a2, INT a3, INT a4, INT e1, INT e2, INT e3, INT e4);
	INT term5(INT a1, INT a2, INT a3, INT a4, INT a5, INT e1, INT e2, INT e3, INT e4, INT e5);
	//INT product2(INT a1, INT a2);
	INT m_term(INT q, INT a1, INT a2, INT a3);
	INT beta_trinomial(INT q, INT beta, INT a1, INT a2, INT a3);
	INT T3product2(INT a1, INT a2);
};

// ####################################################################################
// finite_field_tables.C:
// ####################################################################################

extern INT finitefield_primes[];
extern INT finitefield_nb_primes;
extern INT finitefield_largest_degree_irreducible_polynomial[];
extern const BYTE *finitefield_primitive_polynomial[][100];
const BYTE *get_primitive_polynomial(INT p, INT e, INT verbose_level);

// ####################################################################################
// subfield_structure.C:
// ####################################################################################


class subfield_structure {
public:

	finite_field *FQ;
	finite_field *Fq;
	INT Q;
	INT q;
	INT s; // subfield index: q^s = Q
	INT *Basis; // [s], entries are elements in FQ

	INT *embedding; 
			// [Q], entries are elements in FQ, 
			// indexed by elements in AG(s,q)
	INT *embedding_inv;
			// [Q], entries are ranks of elements in AG(s,q), 
			// indexed by elements in FQ
			// the inverse of embedding

	INT *components; // [Q * s], entries are elements in Fq
			// the vectors corresponding to the AG(s,q) ranks in embedding_inv[]

	INT *FQ_embedding; 
		// [q] entries are elements in FQ corresponding to 
		// the elements in Fq
	INT *Fq_element;
		// [Q], entries are the elements in Fq corresponding to a given FQ element
		// or -1 if the FQ element does not belong to Fq.
	INT *v; // [s]
	
	subfield_structure();
	~subfield_structure();
	void null();
	void freeself();
	void init(finite_field *FQ, finite_field *Fq, INT verbose_level);
	void init_with_given_basis(finite_field *FQ, finite_field *Fq, INT *given_basis, INT verbose_level);
	void print_embedding();
	INT evaluate_over_FQ(INT *v);
	INT evaluate_over_Fq(INT *v);
	void lift_matrix(INT *MQ, INT m, INT *Mq, INT verbose_level);
	void retract_matrix(INT *Mq, INT n, INT *MQ, INT m, INT verbose_level);

};




// ####################################################################################
// tensor.C:
// ####################################################################################

void twisted_tensor_product_codes(
	INT *&H_subfield, INT &m, INT &n, 
	finite_field *F, finite_field *f, 
	INT f_construction_A, INT f_hyperoval, 
	INT f_construction_B, INT verbose_level);
void create_matrix_M(
	INT *&M, 
	finite_field *F, finite_field *f,
	INT &m, INT &n, INT &beta, INT &r, INT *exponents, 
	INT f_construction_A, INT f_hyperoval, INT f_construction_B, 
	INT f_elements_exponential, const BYTE *symbol_for_print, 
	INT verbose_level);
// INT exponents[9]
void create_matrix_H_subfield(finite_field *F, finite_field*f, 
	INT *H_subfield, INT *C, INT *C_inv, INT *M, INT m, INT n, INT beta, INT beta_q, 
	INT f_elements_exponential, const BYTE *symbol_for_print, const BYTE *symbol_for_print_subfield, 
	INT f_construction_A, INT f_hyperoval, INT f_construction_B, 
	INT verbose_level);
void tt_field_reduction(finite_field &F, finite_field &f, INT m, INT n, INT *M, INT *MM, INT verbose_level);


void make_tensor_code_9dimensional_as_point_set(finite_field *F, 
	INT *&the_set, INT &length, 
	INT verbose_level);
void make_tensor_code_9_dimensional(INT q, 
	const BYTE *override_poly_Q, const BYTE *override_poly, 
	INT f_hyperoval, 
	INT *&code, INT &length, 
	INT verbose_level);

// ####################################################################################
// norm_tables.C:
// ####################################################################################

class norm_tables {
public:
	INT *norm_table;
	INT *norm_table_sorted;
	INT *sorting_perm, *sorting_perm_inv;
	INT nb_types;
	INT *type_first, *type_len;
	INT *the_type;

	norm_tables();
	~norm_tables();
	void init(unusual_model &U, INT verbose_level);
	INT choose_an_element_of_given_norm(INT norm, INT verbose_level);
	
};

// ####################################################################################
// longinteger_object.C:
// ####################################################################################

extern INT longinteger_f_print_scientific;

class longinteger_object {

private:
	char sgn; // TRUE if negative
	int l;
	char *r;
	
public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;
	
	longinteger_object();
	~longinteger_object();
	void freeself();
	
	char &sign() { return sgn; };
	int &len() { return l; };
	char *&rep() { return r; };
	void create(INT i);
	void create_product(INT nb_factors, INT *factors);
	void create_from_base_b_representation(INT b, INT *rep, INT len);
	void create_from_base_10_string(const BYTE *str, INT verbose_level);
	INT as_INT();
	void as_longinteger(longinteger_object &a);
	void assign_to(longinteger_object &b);
	void swap_with(longinteger_object &b);
	ostream& print(ostream& ost);
	ostream& print_not_scientific(ostream& ost);
	INT output_width();
	void print_width(ostream& ost, INT width);
	void print_to_string(BYTE *str);
	void normalize();
	void negate();
	int is_zero();
	void zero();
	int is_one();
	int is_mone();
	int is_one_or_minus_one();
	void one();
	void increment();
	void decrement();
	void add_INT(INT a);
	void create_i_power_j(INT i, INT j);
	INT compare_with_INT(INT a);
};

ostream& operator<<(ostream& ost, longinteger_object& p);

// ####################################################################################
// longinteger_domain.C:
// ####################################################################################

class longinteger_domain {

public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;
	
	INT compare(longinteger_object &a, longinteger_object &b);
	INT compare_unsigned(longinteger_object &a, longinteger_object &b);
	// returns -1 if a < b, 0 if a = b, and 1 if a > b, treating a and b as unsigned.
	void subtract_signless(longinteger_object &a, longinteger_object &b, longinteger_object &c);
	// c = a - b, assuming a > b
	void subtract_signless_in_place(longinteger_object &a, longinteger_object &b);
	// a := a - b, assuming a > b
	void add(longinteger_object &a, longinteger_object &b, longinteger_object &c);
	void add_in_place(longinteger_object &a, longinteger_object &b);
	void mult(longinteger_object &a, longinteger_object &b, longinteger_object &c);
	void mult_integer_in_place(longinteger_object &a, INT b);
	void mult_mod(longinteger_object &a, 
		longinteger_object &b, longinteger_object &c, 
		longinteger_object &m, INT verbose_level);
	void multiply_up(longinteger_object &a, INT *x, INT len);
	INT quotient_as_INT(longinteger_object &a, longinteger_object &b);
	void integral_division_exact(longinteger_object &a, longinteger_object &b, longinteger_object &a_over_b);
	void integral_division(
		longinteger_object &a, longinteger_object &b, 
		longinteger_object &q, longinteger_object &r, INT verbose_level);
	void integral_division_by_INT(longinteger_object &a, 
		INT b, longinteger_object &q, INT &r);
	void extended_gcd(longinteger_object &a, longinteger_object &b, 
		longinteger_object &g, longinteger_object &u, longinteger_object &v, INT verbose_level);
	INT logarithm_base_b(longinteger_object &a, INT b);
	void base_b_representation(longinteger_object &a, INT b, INT *&rep, INT &len);
	void power_int(longinteger_object &a, INT n);
	void power_int_mod(longinteger_object &a, INT n, longinteger_object &m);
	void power_longint_mod(longinteger_object &a, 
		longinteger_object &n, longinteger_object &m, INT verbose_level);
	void create_qnm1(longinteger_object &a, INT q, INT n);
	void binomial(longinteger_object &a, INT n, INT k);
	void size_of_conjugacy_class_in_sym_n(longinteger_object &a, INT n, INT *part);
	void q_binomial(longinteger_object &a, 
		INT n, INT k, INT q, INT verbose_level);
	void q_binomial_no_table(longinteger_object &a, 
		INT n, INT k, INT q, INT verbose_level);
	void krawtchouk(longinteger_object &a, INT n, INT q, INT k, INT x);
	INT is_even(longinteger_object &a);
	INT is_odd(longinteger_object &a);
	INT remainder_mod_INT(longinteger_object &a, INT p);
	INT multiplicity_of_p(longinteger_object &a, longinteger_object &residue, INT p);
	INT smallest_primedivisor(longinteger_object &a, INT p_min, INT verbose_level);
	void factor_into_longintegers(longinteger_object &a, 
		INT &nb_primes, longinteger_object *&primes, 
		INT *&exponents, INT verbose_level);
	void factor(longinteger_object &a, INT &nb_primes, 
		INT *&primes, INT *&exponents, 
		INT verbose_level);
	INT jacobi(longinteger_object &a, longinteger_object &m, INT verbose_level);
	void random_number_less_than_n(longinteger_object &n, longinteger_object &r);
	void find_probable_prime_above(
		longinteger_object &a, 
		INT nb_solovay_strassen_tests, INT f_miller_rabin_test, 
		INT verbose_level);
	INT solovay_strassen_is_prime(
		longinteger_object &n, INT nb_tests, INT verbose_level);
	INT solovay_strassen_is_prime_single_test(
		longinteger_object &n, INT verbose_level);
	INT miller_rabin_test(
		longinteger_object &n, INT verbose_level);
	void get_k_bit_random_pseudoprime(
		longinteger_object &n, INT k, 
		INT nb_tests_solovay_strassen, 
		INT f_miller_rabin_test, INT verbose_level);
	void RSA_setup(longinteger_object &n, 
		longinteger_object &p, longinteger_object &q, 
		longinteger_object &a, longinteger_object &b, 
		INT nb_bits, 
		INT nb_tests_solovay_strassen, INT f_miller_rabin_test, 
		INT verbose_level);
	void matrix_product(longinteger_object *A, longinteger_object *B, longinteger_object *&C, INT Am, INT An, INT Bn);
	void matrix_entries_integral_division_exact(longinteger_object *A, longinteger_object &b, INT Am, INT An);
	void matrix_print_GAP(ostream &ost, longinteger_object *A, INT Am, INT An);
	void matrix_print_tex(ostream &ost, longinteger_object *A, INT Am, INT An);
	void power_mod(char *aa, char *bb, char *nn, 
		longinteger_object &result, INT verbose_level);
	void factorial(longinteger_object &result, INT n);
};

void test_longinteger();
void test_longinteger2();
void test_longinteger3();
void test_longinteger4();
void test_longinteger5();
void test_longinteger6();
void test_longinteger7();
void test_longinteger8();
void mac_williams_equations(longinteger_object *&M, INT n, INT k, INT q);
void determine_weight_enumerator();
void longinteger_collect_setup(INT &nb_agos, longinteger_object *&agos, INT *&multiplicities);
void longinteger_collect_free(INT &nb_agos, longinteger_object *&agos, INT *&multiplicities);
void longinteger_collect_add(INT &nb_agos, longinteger_object *&agos, INT *&multiplicities, longinteger_object &ago);
void longinteger_collect_print(ostream &ost, INT &nb_agos, longinteger_object *&agos, INT *&multiplicities);
void longinteger_free_global_data();
void longinteger_print_digits(BYTE *rep, INT len);

// ####################################################################################
// unipoly.C:
// ####################################################################################

class unipoly_domain {
public:
	finite_field *gfq;
	INT f_factorring;
	INT factor_degree;
	INT *factor_coeffs;
	unipoly_object factor_poly;

	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;
	

	unipoly_domain(finite_field *GFq);
	unipoly_domain(finite_field *GFq, unipoly_object m);
	~unipoly_domain();
	INT &s_i(unipoly_object p, INT i) { INT *rep = (INT *) p; return rep[i + 1]; };
	void create_object_of_degree(unipoly_object &p, INT d);
	void create_object_of_degree_with_coefficients(unipoly_object &p, INT d, INT *coeff);
	void create_object_by_rank(unipoly_object &p, INT rk);
	void create_object_by_rank_longinteger(unipoly_object &p, 
		longinteger_object &rank, INT verbose_level);
	void create_object_by_rank_string(unipoly_object &p, const BYTE *rk, INT verbose_level);
	void create_Dickson_polynomial(unipoly_object &p, INT *map);
	void delete_object(unipoly_object &p);
	void unrank(unipoly_object p, INT rk);
	void unrank_longinteger(unipoly_object p, longinteger_object &rank);
	INT rank(unipoly_object p);
	void rank_longinteger(unipoly_object p, longinteger_object &rank);
	INT degree(unipoly_object p);
	ostream& print_object(unipoly_object p, ostream& ost);
	void assign(unipoly_object a, unipoly_object &b);
	void one(unipoly_object p);
	void m_one(unipoly_object p);
	void zero(unipoly_object p);
	INT is_one(unipoly_object p);
	INT is_zero(unipoly_object p);
	void negate(unipoly_object a);
	void make_monic(unipoly_object &a);
	void add(unipoly_object a, unipoly_object b, unipoly_object &c);
	void mult(unipoly_object a, unipoly_object b, unipoly_object &c);
	void mult_easy(unipoly_object a, unipoly_object b, unipoly_object &c);
	void mult_mod(unipoly_object a, unipoly_object b, unipoly_object &c, 
		INT factor_polynomial_degree, INT *factor_polynomial_coefficents_negated, 
		INT verbose_level);
	void Frobenius_matrix(INT *&Frob, unipoly_object factor_polynomial, INT verbose_level);
	void Berlekamp_matrix(INT *&B, unipoly_object factor_polynomial, INT verbose_level);
	void integral_division_exact(unipoly_object a, unipoly_object b, unipoly_object &q, INT verbose_level);
	void integral_division(unipoly_object a, unipoly_object b, 
		unipoly_object &q, unipoly_object &r, INT verbose_level);
	void derive(unipoly_object a, unipoly_object &b);
	INT compare_euclidean(unipoly_object m, unipoly_object n);
	void greatest_common_divisor(unipoly_object m, unipoly_object n, 
		unipoly_object &g, INT verbose_level);
	void extended_gcd(unipoly_object m, unipoly_object n, 
		unipoly_object &u, unipoly_object &v, 
		unipoly_object &g, INT verbose_level);
	INT is_squarefree(unipoly_object p, INT verbose_level);
	void compute_normal_basis(INT d, INT *Normal_basis, INT *Frobenius, INT verbose_level);
	void order_ideal_generator(INT d, INT idx, unipoly_object &mue, 
		INT *A, INT *Frobenius, 
		INT verbose_level);
		// Lueneburg~\cite{Lueneburg87a} p. 105.
		// Frobenius is a matrix of size d x d
		// A is a matrix of size (d + 1) x d
	void matrix_apply(unipoly_object &p, INT *Mtx, INT n, INT verbose_level);
		// The matrix is applied on the left
	void substitute_matrix_in_polynomial(unipoly_object &p, INT *Mtx_in, INT *Mtx_out, INT k, INT verbose_level);
		// The matrix is substituted into the polynomial
	INT substitute_scalar_in_polynomial(unipoly_object &p, INT scalar, INT verbose_level);
		// The scalar 'scalar' is substituted into the polynomial
	void module_structure_apply(INT *v, INT *Mtx, INT n, unipoly_object p, INT verbose_level);
	void take_away_all_factors_from_b(unipoly_object a, 
		unipoly_object b, unipoly_object &a_without_b, INT verbose_level);
		// Computes the polynomial $r$ with
		//\begin{enumerate}
		//\item
		//$r$ divides $a$
		//\item
		//$gcd(r,b) = 1$ and
		//\item
		//each irreducible polynomial dividing $a/r$ divides $b$.
		//Lueneburg~\cite{Lueneburg87a}, p. 37.
		//\end{enumerate}
	INT is_irreducible(unipoly_object a, INT verbose_level);
	void singer_candidate(unipoly_object &m, INT p, INT d, INT b, INT a);
	INT is_primitive(unipoly_object &m, 
		longinteger_object &qm1, 
		INT nb_primes, longinteger_object *primes, 
		INT verbose_level);
	void get_a_primitive_polynomial(unipoly_object &m, 
		INT f, INT verbose_level);
	void get_an_irreducible_polynomial(unipoly_object &m, 
		INT f, INT verbose_level);
	void power_INT(unipoly_object &a, INT n, INT verbose_level);
	void power_longinteger(unipoly_object &a, longinteger_object &n);
	void power_coefficients(unipoly_object &a, INT n);
	void minimum_polynomial(unipoly_object &a, 
		INT alpha, INT p, INT verbose_level);
	INT minimum_polynomial_factorring(INT alpha, INT p, INT verbose_level);
	void minimum_polynomial_factorring_longinteger(
		longinteger_object &alpha, longinteger_object &rk_minpoly, 
		INT p, INT verbose_level);
	void BCH_generator_polynomial(unipoly_object &g, INT n, 
		INT designed_distance, INT &bose_distance, 
		INT &transversal_length, INT *&transversal, 
		longinteger_object *&rank_of_irreducibles, 
		INT verbose_level);
	//void BCH_generator_polynomial_general(unipoly_object &g, INT n, 
		//INT designed_distance, INT start, INT &bose_distance, 
		//INT &transversal_length, INT *&transversal, longinteger_object *&rank_of_rreducibles, 
		//INT f_v, INT f_vv);
	void compute_generator_matrix(unipoly_object a, INT *&genma, INT n, INT &k, INT verbose_level);
	void print_vector_of_polynomials(unipoly_object *sigma, INT deg);
	void minimum_polynomial_extension_field(unipoly_object &g, unipoly_object m, 
		unipoly_object &minpol, INT d, INT *Frobenius, INT verbose_level);
		// Lueneburg~\cite{Lueneburg87a}, p. 112.
	void characteristic_polynomial(INT *Mtx, INT k, unipoly_object &char_poly, INT verbose_level);
	void print_matrix(unipoly_object *M, INT k);
	void determinant(unipoly_object *M, INT k, unipoly_object &p, INT verbose_level);
	void deletion_matrix(unipoly_object *M, INT k, INT delete_row, INT delete_column, unipoly_object *&N, INT verbose_level);

};


// ####################################################################################
// rank_checker.C:
// ####################################################################################



class rank_checker {

public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;
	
	finite_field *GFq;
	INT m, n, d;
	
	INT *M1; // [m * n]
	INT *M2; // [m * n]
	INT *base_cols; // [n]
	INT *set; // [n] used in check_mindist

	rank_checker();
	~rank_checker();
	void init(finite_field *GFq, INT m, INT n, INT d);
	INT check_rank(INT len, INT *S, INT verbose_level);
	INT check_rank_matrix_input(INT len, INT *S, INT dim_S, INT verbose_level);
	INT check_rank_last_two_are_fixed(INT len, INT *S, INT verbose_level);
	INT compute_rank(INT len, INT *S, INT f_projective, INT verbose_level);
	INT compute_rank_row_vectors(INT len, INT *S, INT f_projective, INT verbose_level);
};




// ####################################################################################
// classify.C:
// ####################################################################################

class classify {

public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;
	
	INT data_length;
	
	INT *data;
	INT *data_sorted;
	INT *sorting_perm; // computed using INT_vec_sorting_permutation
	INT *sorting_perm_inv; // perm_inv[i] is the index in data of the element in data_sorted[i]
	INT nb_types;
	INT *type_first;
	INT *type_len;
	
	INT f_second;
	INT *second_data_sorted;
	INT *second_sorting_perm;
	INT *second_sorting_perm_inv;
	INT second_nb_types;
	INT *second_type_first;
	INT *second_type_len;
	
	classify();
	~classify();
	void init(INT *data, INT data_length, INT f_second, INT verbose_level);
	INT class_of(INT pt_idx);
	void print(INT f_backwards);
	void print_file(ostream &ost, INT f_backwards);
	void print_naked(INT f_backwards);
	void print_naked_tex(ostream &ost, INT f_backwards);
	double average();
	double average_of_non_zero_values();
};



// ####################################################################################
// mp_graphics.C:
// ####################################################################################


struct grid_frame {
	INT f_matrix_notation;
	double origin_x;
	double origin_y;
	INT m; // number of rows in the grid
	INT n; // number of columns in the grid
	double dx;
	double dy;
};

class mp_graphics {

	char fname_base[1000];
	char fname_mp[1000];
	char fname_log[1000];
	char fname_tikz[1000];
	ofstream fp_mp;
	ofstream fp_log;
	ofstream fp_tikz;
	INT f_file_open;
	
	
	// coordinate systems:
	
	INT user[4]; // llx/lly/urx/ury 
	INT dev[4]; // llx/lly/urx/ury 

	INT x_min, x_max, y_min, y_max, f_min_max_set;

	INT txt_halign; // 0=left aligned, 1=centered, 2=right aligned; default=0
	INT txt_valign; // 0=bottom, 1=middle, 2=top; default=0
	INT txt_boxed; //  default=0
	INT txt_overwrite; // default=0
	INT txt_rotate; // default = 0 (in degree)

	INT line_beg_style; // default=0
	INT line_end_style; // 0=nothing, 1=arrow; default=0

	INT line_thickness; // 1,2,3

	INT fill_interior; // in 1/100th, 0= none (used for pie-drawing); default=0
	INT fill_color; // 0 = white, 1 = black; default=0
	INT fill_shape; // 0 =  .., 1 = -- ; default=1
	INT fill_outline; // default=0
	INT fill_nofill; // default=0


	INT line_dashing;
		// 0 = no dashing, otherwise scaling factor 1/100th evenly
		// default=0
	
	INT cur_path;

	INT f_embedded; // have a header so that the file can be compiled (for tikz)
	INT f_sideways;

public:
	// for tikz:
	double tikz_global_scale; // .45 works
	double tikz_global_line_width; // 1.5 works


	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;

	mp_graphics();
	mp_graphics(char *file_name, INT xmin, INT ymin, INT xmax, INT ymax, 
		INT f_embedded, INT f_sideways);
	~mp_graphics();
	void default_values();
	void init(char *file_name, INT xmin, INT ymin, INT xmax, INT ymax, INT f_embedded, INT f_sideways);
	void exit(ostream &ost, INT verbose_level);
	void setup(const char *fname_base, 
		INT in_xmin, INT in_ymin, INT in_xmax, INT in_ymax, 
		INT xmax, INT ymax, INT f_embedded, INT f_sideways, 
		double scale, double line_width);
	void set_parameters(double scale, double line_width);
	void set_scale(double scale);
	void frame(double move_out);
	void frame_constant_aspect_ratio(double move_out);
	void finish(ostream &ost, INT verbose_level);

	INT& out_xmin();
	INT& out_ymin();
	INT& out_xmax();
	INT& out_ymax();

	void user2dev(INT &x, INT &y);
	void dev2user(INT &x, INT &y);
	void user2dev_dist(INT &x, INT &y);
	void user2dev_dist_x(INT &x);
	void user2dev_dist_y(INT &y);
	void dev2user_dist(INT &x, INT &y);

	void nice_circle(INT x, INT y, INT rad);
	void grid_polygon2(grid_frame *F, INT x0, INT y0, INT x1, INT y1);
	void grid_polygon4(grid_frame *F, INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3);
	void grid_polygon5(grid_frame *F, INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3, INT x4, INT y4);
	void polygon(INT *Px, INT *Py, INT n);
	void polygon2(INT *Px, INT *Py, INT i1, INT i2);
	void polygon3(INT *Px, INT *Py, INT i1, INT i2, INT i3);
	void polygon4(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4);
	void polygon5(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5);
	void polygon6(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6);
	void polygon7(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7);
	void polygon8(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8);
	void polygon9(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9);
	void polygon10(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10);
	void polygon11(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10, INT i11);
	void polygon_idx(INT *Px, INT *Py, INT *Idx, INT n);
	void bezier(INT *Px, INT *Py, INT n);
	void bezier2(INT *Px, INT *Py, INT i1, INT i2);
	void bezier3(INT *Px, INT *Py, INT i1, INT i2, INT i3);
	void bezier4(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4);
	void bezier5(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5);
	void bezier6(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6);
	void bezier7(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7);
	void bezier_idx(INT *Px, INT *Py, INT *Idx, INT n);
	void grid_fill_polygon4(grid_frame *F, 
		INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3);
	void grid_fill_polygon5(grid_frame *F, 
		INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3, INT x4, INT y4);
	void fill_polygon3(INT *Px, INT *Py, INT i1, INT i2, INT i3);
	void fill_polygon4(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4);
	void fill_polygon5(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5);
	void fill_polygon6(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6);
	void fill_polygon7(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7);
	void fill_polygon8(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8);
	void fill_polygon9(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9);
	void fill_polygon10(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10);
	void fill_polygon11(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10, INT i11);
	void polygon2_arrow_halfway(INT *Px, INT *Py, INT i1, INT i2);
	void polygon2_arrow_halfway_and_label(INT *Px, INT *Py, INT i1, INT i2, 
		const BYTE *alignment, const BYTE *txt);
	void grid_aligned_text(grid_frame *F, INT x, INT y, const char *alignment, const char *p);
	void aligned_text(INT x, INT y, const char *alignment, const char *p);
	void aligned_text_array(INT *Px, INT *Py, INT idx, const char *alignment, const char *p);
	void aligned_text_with_offset(INT x, INT y, INT xoffset, INT yoffset, 
		const char *alignment, const char *p);

	void st_alignment(INT txt_halign, INT txt_valign);
	void get_alignment(BYTE *align);
	void sl_udsty(INT line_dashing);
	void sl_ends(INT line_beg_style, INT line_end_style);
	void sl_thickness(INT line_thickness);
		// added Oct 22 2012
	void sf_interior(INT fill_interior);
	void sf_color(INT fill_color);
	void sf_shape(INT fill_shape);
	void sf_outline(INT fill_outline);
	void sf_nofill(INT fill_nofill);
	void st_boxed(INT txt_boxed);
	void st_overwrite(INT txt_overwrite);
	void st_rotate(INT txt_rotate);
	void coords_min_max(INT x, INT y);
	INT get_label(INT x, INT y);
	void draw_boxes_final();

	// output commands:
	void header();
	void footer();
	void begin_figure(INT factor_1000);
	void end_figure();
	void output_xy_metapost(INT x, INT y);
	void output_x_metapost(INT x);
	void output_y_metapost(INT y);
	void output_xy_tikz(INT x, INT y);
	void output_x_tikz(INT x);
	void output_y_tikz(INT y);

	void comment(const BYTE *p);
	void text(INT x, INT y, const char *p);
	void circle(INT x, INT y, INT rad);
	void circle_text(INT x, INT y, const char *text);
	void output_circle_text_mp(INT x, INT y, INT idx, const char *text);
	void output_circle_text_tikz(INT x, INT y, INT idx, const char *text);
	void polygon_or_bezier_idx(INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle);
	void fill_idx(INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle);
	void fill_idx_mp(INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle);
	void fill_idx_tikz(ofstream &fp, INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle);


};


// ####################################################################################
// draw.C
// ####################################################################################

void transform_llur(INT *in, INT *out, INT &x, INT &y);
void transform_dist(INT *in, INT *out, INT &x, INT &y);
void transform_dist_x(INT *in, INT *out, INT &x);
void transform_dist_y(INT *in, INT *out, INT &y);
void transform_llur_double(double *in, double *out, double &x, double &y);
void draw(BYTE *fname);
void on_circle_int(INT *Px, INT *Py, INT idx, INT angle_in_degree, INT rad);
void on_circle_double(double *Px, double *Py, INT idx, double angle_in_degree, double rad);
void polygon3D(mp_graphics &G, INT *Px, INT *Py, INT dim, INT x0, INT y0, INT z0, INT x1, INT y1, INT z1);
void integer_4pts(mp_graphics &G, INT *Px, INT *Py, INT p1, INT p2, INT p3, INT p4, 
	const BYTE *align, INT a);
void text_4pts(mp_graphics &G, INT *Px, INT *Py, INT p1, INT p2, INT p3, INT p4, 
	const BYTE *align, const BYTE *str);
void affine_pt1(INT *Px, INT *Py, INT p0, INT p1, INT p2, double f1, INT p3);
void affine_pt2(INT *Px, INT *Py, INT p0, INT p1, INT p1b, double f1, INT p2, INT p2b, double f2, INT p3);
INT C3D(INT i, INT j, INT k);
INT C2D(INT i, INT j);
double cos_grad(double phi);
double sin_grad(double phi);
double tan_grad(double phi);
double atan_grad(double x);
void adjust_coordinates_double(double *Px, double *Py, INT *Qx, INT *Qy, 
	INT N, double xmin, double ymin, double xmax, double ymax, INT verbose_level);
void Intersection_of_lines(double *X, double *Y, double *a, double *b, double *c, INT l1, INT l2, INT pt);
void intersection_of_lines(double a1, double b1, double c1, double a2, double b2, double c2, 
	double &x, double &y);
void Line_through_points(double *X, double *Y, double *a, double *b, double *c, 
	INT pt1, INT pt2, INT line_idx);
void line_through_points(double pt1_x, double pt1_y, 
	double pt2_x, double pt2_y, double &a, double &b, double &c);
void intersect_circle_line_through(double rad, double x0, double y0, 
	double pt1_x, double pt1_y, 
	double pt2_x, double pt2_y, 
	double &x1, double &y1, double &x2, double &y2);
void intersect_circle_line(double rad, double x0, double y0, double a, double b, double c, 
	double &x1, double &y1, double &x2, double &y2);
void affine_combination(double *X, double *Y, 
	INT pt0, INT pt1, INT pt2, double alpha, INT new_pt);
void draw_graph(mp_graphics *G, INT x, INT y, INT dx, INT dy, INT nb_V, INT *Edges, INT nb_E);
void draw_graph_with_distinguished_edge(mp_graphics *G, INT x, INT y, 
	INT dx, INT dy, INT nb_V, INT *Edges, INT nb_E, 
	INT distinguished_edge, INT verbose_level);
void draw_graph_on_multiple_circles(mp_graphics *G, INT x, INT y, INT dx, INT dy, INT nb_V, INT *Edges, INT nb_E, INT nb_circles);
void draw_graph_on_2D_grid(mp_graphics *G, INT x, INT y, INT dx, INT dy, INT rad, INT nb_V, 
	INT *Edges, INT nb_E, INT *coords_2D, INT *Base, 
	INT f_point_labels, INT point_label_offset, INT f_directed);
void draw_tournament(mp_graphics *G, INT x, INT y, INT dx, INT dy, INT nb_V, INT *Edges, INT nb_E, INT verbose_level);
void draw_bitmatrix(const BYTE *fname_base, INT f_dots, 
	INT f_partition, INT nb_row_parts, INT *row_part_first, INT nb_col_parts, INT *col_part_first, 
	INT f_row_grid, INT f_col_grid, 
	INT f_bitmatrix, UBYTE *D, INT *M, 
	INT m, INT n, INT xmax_in, INT ymax_in, INT xmax, INT ymax, 
	INT f_has_labels, INT *labels);
void draw_bitmatrix2(mp_graphics &G, INT f_dots, 
	INT f_partition, INT nb_row_parts, INT *row_part_first, INT nb_col_parts, INT *col_part_first, 
	INT f_row_grid, INT f_col_grid, 
	INT f_bitmatrix, UBYTE *D, INT *M, 
	INT m, INT n, INT xmax, INT ymax, 
	INT f_has_labels, INT *labels);


// ####################################################################################
// partitionstack.C
// ####################################################################################


ostream& operator<<(ostream& ost, partitionstack& p);

class partitionstack {
	public:

	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;


	// data structure for the partition stack,
	// following Leon:
		INT n;
		INT ht;
		INT ht0;

		INT *pointList, *invPointList;
		INT *cellNumber;

		INT *startCell;
		INT *cellSize;
		INT *parent;


	// for matrix canonization:
	// INT first_column_element;

	// subset to be chosen by classify_by_..._extract_subset():
	// used as input for split_cell()
		//
	// used if SPLIT_MULTIPLY is defined:
		INT nb_subsets;
		INT *subset_first;
		INT *subset_length;
		INT *subsets;
		//
	// used if SPLIT_MULTIPLY is not defined:
		INT *subset;
		INT subset_size;

		partitionstack();
		~partitionstack();
		void allocate(INT n, INT verbose_level);
		void free();
		INT parent_at_height(INT h, INT cell);
		INT is_discrete();
		INT smallest_non_discrete_cell();
		INT biggest_non_discrete_cell();
		INT smallest_non_discrete_cell_rows_preferred();
		INT biggest_non_discrete_cell_rows_preferred();
		INT nb_partition_classes(INT from, INT len);
		INT is_subset_of_cell(INT *set, INT size, INT &cell_idx);
		void sort_cells();
		void sort_cell(INT cell);
		void reverse_cell(INT cell);
		void check();
		void print_raw();
		void print_class(ostream& ost, INT idx);
		void print_classes_tex(ostream& ost);
		void print_class_tex(ostream& ost, INT idx);
		void print_class_point_or_line(ostream& ost, INT idx);
		void print_classes(ostream& ost);
		void print_classes_points_and_lines(ostream& ost);
		ostream& print(ostream& ost);
		void print_cell(INT i);
		void print_cell_latex(ostream &ost, INT i);
		void print_subset();
		void write_cell_to_file(INT i, BYTE *fname, INT verbose_level);
		void write_cell_to_file_points_or_lines(INT i, BYTE *fname, INT verbose_level);
		void refine_arbitrary_set(INT size, INT *set, INT verbose_level);
		void split_cell(INT verbose_level);
		void split_multiple_cells(INT *set, INT set_size, INT f_front, INT verbose_level);
		void split_line_cell_front_or_back(INT *set, INT set_size, INT f_front, INT verbose_level);
		void split_cell_front_or_back(INT *set, INT set_size, INT f_front, INT verbose_level);
		void split_cell(INT *set, INT set_size, INT verbose_level);
		void join_cell();
		void reduce_height(INT ht0);
		void isolate_point(INT pt);
		void subset_continguous(INT from, INT len);
		INT is_row_class(INT c);
		INT is_col_class(INT c);
		void allocate_and_get_decomposition(
			INT *&row_classes, INT *&row_class_inv, INT &nb_row_classes,
			INT *&col_classes, INT *&col_class_inv, INT &nb_col_classes, 
			INT verbose_level);
		void get_row_and_col_permutation(
			INT *row_classes, INT nb_row_classes,
			INT *col_classes, INT nb_col_classes, 
			INT *row_perm, INT *row_perm_inv, 
			INT *col_perm, INT *col_perm_inv);
		void get_row_and_col_classes(INT *row_classes, INT &nb_row_classes,
			INT *col_classes, INT &nb_col_classes, INT verbose_level);
		void initial_matrix_decomposition(INT nbrows, INT nbcols,
			INT *V, INT nb_V, INT *B, INT nb_B, INT verbose_level);
		INT is_descendant_of(INT cell, INT ancestor_cell, INT verbose_level);
		INT is_descendant_of_at_level(INT cell, INT ancestor_cell, INT level, INT verbose_level);
		INT cellSizeAtLevel(INT cell, INT level);

	// TDO for orthogonal:
	INT compute_TDO(orthogonal &O, INT ht0, 
		INT marker1, INT marker2, INT depth, INT verbose_level);
	void get_and_print_row_decomposition_scheme(orthogonal &O, 
		INT marker1, INT marker2);
	void get_and_print_col_decomposition_scheme(orthogonal &O, 
		INT marker1, INT marker2);
	void get_and_print_decomposition_schemes(orthogonal &O, INT marker1, INT marker2);
	void print_decomposition_tex(ostream &ost, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes);
	void print_decomposition_scheme(ostream &ost, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes, 
		INT *scheme, INT marker1, INT marker2);
	void print_decomposition_scheme_tex(ostream &ost, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes, 
		INT *scheme);
	void print_tactical_decomposition_scheme_tex_internal(
		ostream &ost, INT f_enter_math_mode, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes, 
		INT *row_scheme, INT *col_scheme, INT f_print_subscripts);
	void print_tactical_decomposition_scheme_tex(ostream &ost, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes, 
		INT *row_scheme, INT *col_scheme, INT f_print_subscripts);
	void print_row_tactical_decomposition_scheme_tex(
		ostream &ost, INT f_enter_math_mode, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes, 
		INT *row_scheme, INT f_print_subscripts);
	void print_column_tactical_decomposition_scheme_tex(
		ostream &ost, INT f_enter_math_mode, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes, 
		INT *col_scheme, INT f_print_subscripts);
	void print_non_tactical_decomposition_scheme_tex(
		ostream &ost, INT f_enter_math_mode, 
		INT *row_classes, INT nb_row_classes,
		INT *col_classes, INT nb_col_classes, 
		INT f_print_subscripts);
	void row_scheme_to_col_scheme(orthogonal &O, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT *row_scheme, INT *col_scheme, INT verbose_level);
	void get_row_decomposition_scheme(orthogonal &O, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT *row_scheme, INT verbose_level);
	void get_col_decomposition_scheme(orthogonal &O, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT *col_scheme, INT verbose_level);
	INT refine_column_partition(orthogonal &O, INT ht0, INT verbose_level);
	INT refine_row_partition(orthogonal &O, INT ht0, INT verbose_level);
	INT hash_column_refinement_info(INT ht0, INT *data, INT depth, INT hash0);
	INT hash_row_refinement_info(INT ht0, INT *data, INT depth, INT hash0);
	void print_column_refinement_info(INT ht0, INT *data, INT depth);
	void print_row_refinement_info(INT ht0, INT *data, INT depth);
	void radix_sort(INT left, INT right, INT *C, 
		INT length, INT radix, INT verbose_level);
	void radix_sort_bits(INT left, INT right, 
		INT *C, INT length, INT radix, INT mask, INT verbose_level);
	void swap_ij(INT *perm, INT *perm_inv, INT i, INT j);
	INT my_log2(INT m);
	void split_by_orbit_partition(INT nb_orbits, 
		INT *orbit_first, INT *orbit_len, INT *orbit,
		INT offset, 
		INT verbose_level);
};

// ####################################################################################
// incidence_structure.C
// ####################################################################################

#define INCIDENCE_STRUCTURE_REALIZATION_BY_MATRIX 1
#define INCIDENCE_STRUCTURE_REALIZATION_BY_ORTHOGONAL 2
#define INCIDENCE_STRUCTURE_REALIZATION_BY_HJELMSLEV 3

class incidence_structure {
	public:

	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;

	BYTE label[1000];


	INT nb_rows;
	INT nb_cols;

	
	INT f_rowsums_constant;
	INT f_colsums_constant;
	INT r;
	INT k;
	INT *nb_lines_on_point;
	INT *nb_points_on_line;
	INT max_r;
	INT min_r;
	INT max_k;
	INT min_k;
	INT *lines_on_point; // [nb_rows * max_r]
	INT *points_on_line; // [nb_cols * max_k]

	INT realization_type;
		// INCIDENCE_STRUCTURE_REALIZATION_BY_MATRIX
		// INCIDENCE_STRUCTURE_REALIZATION_BY_ORTHOGONAL

	INT *M;
	orthogonal *O;
	hjelmslev *H;
	//projective_space *PG;
	
	
	incidence_structure();
	~incidence_structure();
	void null();
	void freeself();
	void check_point_pairs(INT verbose_level);
	INT lines_through_two_points(INT *lines, INT p1, INT p2, INT verbose_level);
	void init_hjelmslev(hjelmslev *H, INT verbose_level);
	void init_orthogonal(orthogonal *O, INT verbose_level);
	void init_by_incidences(INT m, INT n, INT nb_inc, INT *X, INT verbose_level);
	void init_by_R_and_X(INT m, INT n, INT *R, INT *X, INT max_r, INT verbose_level);
	void init_by_matrix(INT m, INT n, INT *M, INT verbose_level);
	void init_by_matrix_as_bitvector(INT m, INT n, UBYTE *M_bitvec, INT verbose_level);
	void init_by_matrix2(INT verbose_level);
	INT nb_points();
	INT nb_lines();
	INT get_ij(INT i, INT j);
	INT get_lines_on_point(INT *data, INT i);
	INT get_points_on_line(INT *data, INT j);
	INT get_nb_inc();
	void save_inc_file(BYTE *fname);
	void save_row_by_row_file(BYTE *fname);
	void print(ostream &ost);
	void compute_TDO_safe_first(partitionstack &PStack, 
		INT depth, INT &step, INT &f_refine, INT &f_refine_prev, INT verbose_level);
	INT compute_TDO_safe_next(partitionstack &PStack, 
		INT depth, INT &step, INT &f_refine, INT &f_refine_prev, INT verbose_level);
		// returns TRUE when we are done, FALSE otherwise
	void compute_TDO_safe(partitionstack &PStack, 
		INT depth, INT verbose_level);
	INT compute_TDO(partitionstack &PStack, INT ht0, INT depth, INT verbose_level);
	INT compute_TDO_step(partitionstack &PStack, INT ht0, INT verbose_level);
	void get_partition(partitionstack &PStack, 
		INT *row_classes, INT *row_class_idx, INT &nb_row_classes, 
		INT *col_classes, INT *col_class_idx, INT &nb_col_classes);
	INT refine_column_partition_safe(partitionstack &PStack, INT verbose_level);
	INT refine_row_partition_safe(partitionstack &PStack, INT verbose_level);
	INT refine_column_partition(partitionstack &PStack, INT ht0, INT verbose_level);
	INT refine_row_partition(partitionstack &PStack, INT ht0, INT verbose_level);
	void print_row_tactical_decomposition_scheme_incidences_tex(
		partitionstack &PStack, 
		ostream &ost, INT f_enter_math_mode, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT f_local_coordinates, INT verbose_level);
	void print_col_tactical_decomposition_scheme_incidences_tex(
		partitionstack &PStack, 
		ostream &ost, INT f_enter_math_mode, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT f_local_coordinates, INT verbose_level);
	void get_incidences_by_row_scheme(partitionstack &PStack, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT row_class_idx, INT col_class_idx, 
		INT rij, INT *&incidences, INT verbose_level);
	void get_incidences_by_col_scheme(partitionstack &PStack, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT row_class_idx, INT col_class_idx, 
		INT kij, INT *&incidences, INT verbose_level);
	void get_row_decomposition_scheme(partitionstack &PStack, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT *row_scheme, INT verbose_level);
	void get_row_decomposition_scheme_if_possible(partitionstack &PStack, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT *row_scheme, INT verbose_level);
	void get_col_decomposition_scheme(partitionstack &PStack, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT *col_scheme, INT verbose_level);
	
	void row_scheme_to_col_scheme(partitionstack &PStack, 
		INT *row_classes, INT *row_class_inv, INT nb_row_classes,
		INT *col_classes, INT *col_class_inv, INT nb_col_classes, 
		INT *row_scheme, INT *col_scheme, INT verbose_level);
	void get_and_print_row_decomposition_scheme(partitionstack &PStack, 
		INT f_list_incidences, INT f_local_coordinates);
	void get_and_print_col_decomposition_scheme(
		partitionstack &PStack, 
		INT f_list_incidences, INT f_local_coordinates);
	void get_and_print_decomposition_schemes(partitionstack &PStack);
	void get_and_print_decomposition_schemes_tex(partitionstack &PStack);
	void get_and_print_tactical_decomposition_scheme_tex(
		ostream &ost, INT f_enter_math, partitionstack &PStack);
	void get_scheme(
		INT *&row_classes, INT *&row_class_inv, INT &nb_row_classes,
		INT *&col_classes, INT *&col_class_inv, INT &nb_col_classes,
		INT *&scheme, INT f_row_scheme, partitionstack &PStack);
	void free_scheme(
		INT *row_classes, INT *row_class_inv, 
		INT *col_classes, INT *col_class_inv, 
		INT *scheme);
	void get_and_print_row_tactical_decomposition_scheme_tex(
		ostream &ost, INT f_enter_math, partitionstack &PStack);
	void get_and_print_column_tactical_decomposition_scheme_tex(
		ostream &ost, INT f_enter_math, partitionstack &PStack);
	void print_non_tactical_decomposition_scheme_tex(
		ostream &ost, INT f_enter_math, partitionstack &PStack);
	void print_line(ostream &ost, partitionstack &P, 
		INT row_cell, INT i, INT *col_classes, INT nb_col_classes, 
		INT width, INT f_labeled);
	void print_column_labels(ostream &ost, partitionstack &P, 
		INT *col_classes, INT nb_col_classes, INT width);
	void print_hline(ostream &ost, partitionstack &P, 
		INT *col_classes, INT nb_col_classes, INT width, INT f_labeled);
	void print_partitioned(ostream &ost, partitionstack &P, INT f_labeled);
	void point_collinearity_graph(INT *Adj, INT verbose_level);
		// G[nb_points() * nb_points()]
	void line_intersection_graph(INT *Adj, INT verbose_level);
		// G[nb_lines() * nb_lines()]
	void latex_it(ostream &ost, partitionstack &P);
	void rearrange(INT *&Vi, INT &nb_V, 
		INT *&Bj, INT &nb_B, INT *&R, INT *&X, partitionstack &P);
	void decomposition_print_tex(ostream &ost, 
		partitionstack &PStack, INT f_row_tactical, INT f_col_tactical, 
		INT f_detailed, INT f_local_coordinates, INT verbose_level);
	void do_tdo_high_level(partitionstack &S, 
		INT f_tdo_steps, INT f_tdo_depth, INT tdo_depth, 
		INT f_write_tdo_files, INT f_pic, 
		INT f_include_tdo_scheme, INT f_include_tdo_extra, INT f_write_tdo_class_files, 
		INT verbose_level);
	void compute_tdo(partitionstack &S, 
		INT f_write_tdo_files, 
		INT f_pic, 
		INT f_include_tdo_scheme, 
		INT verbose_level);
	void compute_tdo_stepwise(partitionstack &S, 
		INT TDO_depth, 
		INT f_write_tdo_files, 
		INT f_pic, 
		INT f_include_tdo_scheme, 
		INT f_include_extra, 
		INT verbose_level);
	void init_partitionstack_trivial(partitionstack *S, 
		INT verbose_level);
	void init_partitionstack(partitionstack *S, 
		INT f_row_part, INT nb_row_parts, INT *row_parts,
		INT f_col_part, INT nb_col_parts, INT *col_parts,
		INT nb_distinguished_point_sets, INT **distinguished_point_sets, INT *distinguished_point_set_size, 
		INT nb_distinguished_line_sets, INT **distinguished_line_sets, INT *distinguished_line_set_size, 
		INT verbose_level);
	void shrink_aut_generators(
		INT nb_distinguished_point_sets, 
		INT nb_distinguished_line_sets, 
		INT Aut_counter, INT *Aut, INT *Base, INT Base_length, 
		INT verbose_level);
	void print_aut_generators(INT Aut_counter, INT *Aut, 
		INT Base_length, INT *Base, INT *Transversal_length);
	void compute_extended_collinearity_graph(
		INT *&Adj, INT &v, INT *&partition, 
		INT f_row_part, INT nb_row_parts, INT *row_parts,
		INT f_col_part, INT nb_col_parts, INT *col_parts,
		INT nb_distinguished_point_sets, INT **distinguished_point_sets, INT *distinguished_point_set_size, 
		INT nb_distinguished_line_sets, INT **distinguished_line_sets, INT *distinguished_line_set_size, 
		INT verbose_level);
		// side effect: the distinguished sets will be sorted afterwards
	void compute_extended_matrix(
		INT *&M, INT &nb_rows, INT &nb_cols, INT &total, INT *&partition, 
		INT f_row_part, INT nb_row_parts, INT *row_parts,
		INT f_col_part, INT nb_col_parts, INT *col_parts,
		INT nb_distinguished_point_sets, INT **distinguished_point_sets, INT *distinguished_point_set_size, 
		INT nb_distinguished_line_sets, INT **distinguished_line_sets, INT *distinguished_line_set_size, 
		INT verbose_level);
};



// two functions from DISCRETA1:

void incma_latex_picture(ostream &fp, 
	INT width, INT width_10, 
	INT f_outline_thin, const BYTE *unit_length, 
	const BYTE *thick_lines, const BYTE *thin_lines, const BYTE *geo_line_width, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X, 
	INT f_labelling_points, const BYTE **point_labels, 
	INT f_labelling_blocks, const BYTE **block_labels);
// width for one box in 0.1mm 
// width_10 is 1 10th of width
// example: width = 40, width_10 = 4 */
void incma_latex(ostream &fp, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X);
void incma_latex_override_unit_length(const BYTE *override_unit_length);
void incma_latex_override_unit_length_drop();






// ####################################################################################
// orthogonal_points.C:
// ####################################################################################



INT count_Sbar(INT n, INT q);
INT count_S(INT n, INT q);
INT count_N1(INT n, INT q);
INT count_T1(INT epsilon, INT n, INT q);
INT count_T2(INT n, INT q);
INT nb_pts_Qepsilon(INT epsilon, INT k, INT q);
// number of singular points on Q^epsilon(k,q)
INT dimension_given_Witt_index(INT epsilon, INT n);
INT Witt_index(INT epsilon, INT k);
INT nb_pts_Q(INT k, INT q);
// number of singular points on Q(k,q)
INT nb_pts_Qplus(INT k, INT q);
// number of singular points on Q^+(k,q)
INT nb_pts_Qminus(INT k, INT q);
// number of singular points on Q^-(k,q)
INT evaluate_quadratic_form(finite_field &GFq, INT *v, INT stride, INT epsilon, INT k, INT form_c1, INT form_c2, INT form_c3);
void Q_epsilon_unrank(finite_field &GFq, INT *v, INT stride, INT epsilon, INT k, INT c1, INT c2, INT c3, INT a);
INT Q_epsilon_rank(finite_field &GFq, INT *v, INT stride, INT epsilon, INT k, INT c1, INT c2, INT c3);
void init_hash_table_parabolic(finite_field &GFq, INT k, INT verbose_level);
void Q_unrank(finite_field &GFq, INT *v, INT stride, INT k, INT a);
INT Q_rank(finite_field &GFq, INT *v, INT stride, INT k);
void Q_unrank_directly(finite_field &GFq, INT *v, INT stride, INT k, INT a);
// k = projective dimension, must be even
INT Q_rank_directly(finite_field &GFq, INT *v, INT stride, INT k);
// k = projective dimension, must be even
void Qplus_unrank(finite_field &GFq, INT *v, INT stride, INT k, INT a);
// k = projective dimension, must be odd
INT Qplus_rank(finite_field &GFq, INT *v, INT stride, INT k);
// k = projective dimension, must be odd
void Qminus_unrank(finite_field &GFq, INT *v, INT stride, INT k, INT a, INT c1, INT c2, INT c3);
// k = projective dimension, must be odd
// the form is 
// \sum_{i=0}^n x_{2i}x_{2i+1} + c1 x_{2n}^2 + c2 x_{2n} x_{2n+1} + c3 x_{2n+1}^2
INT Qminus_rank(finite_field &GFq, INT *v, INT stride, INT k, INT c1, INT c2, INT c3);
// k = projective dimension, must be odd
// the form is 
// \sum_{i=0}^n x_{2i}x_{2i+1} + c1 x_{2n}^2 + c2 x_{2n} x_{2n+1} + c3 x_{2n+1}^2
INT nb_pts_S(INT n, INT q);
INT nb_pts_N(INT n, INT q);
INT nb_pts_N1(INT n, INT q);
INT nb_pts_Sbar(INT n, INT q);
INT nb_pts_Nbar(INT n, INT q);
void S_unrank(finite_field &GFq, INT *v, INT stride, INT n, INT a);
void N_unrank(finite_field &GFq, INT *v, INT stride, INT n, INT a);
void N1_unrank(finite_field &GFq, INT *v, INT stride, INT n, INT a);
void Sbar_unrank(finite_field &GFq, INT *v, INT stride, INT n, INT a);
void Nbar_unrank(finite_field &GFq, INT *v, INT stride, INT n, INT a);
void S_rank(finite_field &GFq, INT *v, INT stride, INT n, INT &a);
void N_rank(finite_field &GFq, INT *v, INT stride, INT n, INT &a);
void N1_rank(finite_field &GFq, INT *v, INT stride, INT n, INT &a);
void Sbar_rank(finite_field &GFq, INT *v, INT stride, INT n, INT &a);
void Nbar_rank(finite_field &GFq, INT *v, INT stride, INT n, INT &a);
INT evaluate_hyperbolic_quadratic_form(finite_field &GFq, INT *v, INT stride, INT n);
INT evaluate_hyperbolic_bilinear_form(finite_field &GFq, INT *u, INT *v, INT n);
INT primitive_element(finite_field &GFq);
void order_POmega_epsilon(INT epsilon, INT m, INT q, longinteger_object &o, INT verbose_level);
void order_PO_epsilon(INT f_semilinear, INT epsilon, INT k, INT q, longinteger_object &go, INT verbose_level);
// k is projective dimension
void order_PO(INT epsilon, INT m, INT q, longinteger_object &o, INT verbose_level);
void order_Pomega(INT epsilon, INT k, INT q, longinteger_object &go, INT verbose_level);
void order_PO_plus(INT m, INT q, longinteger_object &o, INT verbose_level);
void order_PO_minus(INT m, INT q, longinteger_object &o, INT verbose_level);
// m = Witt index, the dimension is n = 2m+2
void order_PO_parabolic(INT m, INT q, longinteger_object &o, INT verbose_level);
void order_Pomega_plus(INT m, INT q, longinteger_object &o, INT verbose_level);
// m = Witt index, the dimension is n = 2m
void order_Pomega_minus(INT m, INT q, longinteger_object &o, INT verbose_level);
// m = half the dimension, the dimension is n = 2m, the Witt index is m - 1
void order_Pomega_parabolic(INT m, INT q, longinteger_object &o, INT verbose_level);
// m = Witt index, the dimension is n = 2m + 1
INT index_POmega_in_PO(INT epsilon, INT m, INT q, INT verbose_level);
void Gram_matrix(finite_field &GFq, INT epsilon, INT k, INT form_c1, INT form_c2, INT form_c3, INT *&Gram);
INT evaluate_bilinear_form(finite_field &GFq, INT *u, INT *v, INT d, INT *Gram);
void choose_anisotropic_form(finite_field &GFq, INT &c1, INT &c2, INT &c3, INT verbose_level);
void Siegel_Transformation(finite_field &GFq, INT epsilon, INT k, 
	INT form_c1, INT form_c2, INT form_c3, INT *M, INT *v, INT *u, INT verbose_level);
void test_Orthogonal(INT epsilon, INT k, INT q);
void test_orthogonal(INT n, INT q);
void orthogonal_Siegel_map_between_singular_points(INT *T, 
	INT rk_from, INT rk_to, INT root, 
	finite_field &GFq, INT epsilon, INT algebraic_dimension, 
	INT form_c1, INT form_c2, INT form_c3, INT *Gram_matrix, 
	INT verbose_level);
// root is not perp to from and to.
INT orthogonal_find_root(INT rk2, 
	finite_field &GFq, INT epsilon, INT algebraic_dimension, 
	INT form_c1, INT form_c2, INT form_c3, INT *Gram_matrix, 
	INT verbose_level);
void orthogonal_points_free_global_data();

// ####################################################################################
// orthogonal.C:
// ####################################################################################

class orthogonal {

public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;

	INT epsilon;
	INT n; // the algebraic dimension
	INT m; // Witt index
	INT q;
	INT f_even;
	INT form_c1, form_c2, form_c3;
	INT *Gram_matrix;
	INT *T1, *T2, *T3; // [n * n]
	INT pt_P, pt_Q;
	INT nb_points;
	INT nb_lines;
	
	INT T1_m;
	INT T1_mm1;
	INT T1_mm2;
	INT T2_m;
	INT T2_mm1;
	INT T2_mm2;
	INT N1_m;
	INT N1_mm1;
	INT N1_mm2;
	INT S_m;
	INT S_mm1;
	INT S_mm2;
	INT Sbar_m;
	INT Sbar_mm1;
	INT Sbar_mm2;
	
	INT alpha; // number of points in the subspace
	INT beta; // number of points in the subspace of the subspace
	INT gamma; // = alpha * beta / (q + 1);
	INT subspace_point_type;
	INT subspace_line_type;
	
	INT nb_point_classes, nb_line_classes;
	INT *A, *B, *P, *L;

	// for hyperbolic:
	INT p1, p2, p3, p4, p5, p6;
	INT l1, l2, l3, l4, l5, l6, l7;
	INT a11, a12, a22, a23, a26, a32, a34, a37, a41, a43, a44, a45, a46, a47, a56, a67;
	INT b11, b12, b22, b23, b26, b32, b34, b37, b41, b43, b44, b45, b46, b47, b56, b67;
	// additionally, for parabolic:
	INT p7, l8;
	INT a21, a36, a57, a22a, a33, a22b, a32b, a42b, a51, a53, a54, a55, a66, a77;
	INT b21, b36, b57, b22a, b33, b22b, b32b, b42b, b51, b53, b54, b55, b66, b77;
	INT a12b, a52a;
	INT b12b, b52a;
	INT delta, omega, lambda, mu, nu, zeta;
	// parabolic q odd requires square / nonsquare tables
	INT *minus_squares; // [(q-1)/2]
	INT *minus_squares_without; // [(q-1)/2 - 1]
	INT *minus_nonsquares; // [(q-1)/2]
	INT *f_is_minus_square; // [q]
	INT *index_minus_square; // [q]
	INT *index_minus_square_without; // [q]
	INT *index_minus_nonsquare; // [q]
	
	INT *v1, *v2, *v3, *v4, *v5, *v_tmp;
	INT *v_tmp2; // for use in parabolic_type_and_index_to_point_rk
	INT *v_neighbor5; 
	
	INT *find_root_x, *find_root_y, *find_root_z;
	INT *line1, *line2, *line3;
	finite_field *F;
	
	// stuff for rank_point
	INT *rk_pt_v;
	
	// stuff for Siegel_transformation
	INT *Sv1, *Sv2, *Sv3, *Sv4;
	INT *Gram2;
	INT *ST_N1, *ST_N2, *ST_w;
	INT *STr_B, *STr_Bv, *STr_w, *STr_z, *STr_x;
	
	// for determine_line
	INT *determine_line_v1, *determine_line_v2, *determine_line_v3;
	
	// for lines_on_point
	INT *lines_on_point_coords1; // [alpha * n]
	INT *lines_on_point_coords2; // [alpha * n]

	orthogonal *subspace;

	// in orthogonal.C:
	void unrank_point(INT *v, INT stride, INT rk, INT verbose_level);
	INT rank_point(INT *v, INT stride, INT verbose_level);
	void unrank_line(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line(INT p1, INT p2, INT verbose_level);
	INT line_type_given_point_types(INT pt1, INT pt2, INT pt1_type, INT pt2_type);
	INT type_and_index_to_point_rk(INT type, INT index, INT verbose_level);
	void point_rk_to_type_and_index(INT rk, INT &type, INT &index, INT verbose_level);
	void canonical_points_of_line(INT line_type, INT pt1, INT pt2, 
		INT &cpt1, INT &cpt2, INT verbose_level);
	INT evaluate_quadratic_form(INT *v, INT stride);
	INT evaluate_bilinear_form(INT *u, INT *v, INT stride);
	INT evaluate_bilinear_form_by_rank(INT i, INT j);
	INT find_root(INT rk2, INT verbose_level);
	void points_on_line_by_line_rank(INT line_rk, INT *line, INT verbose_level);
	void points_on_line(INT pi, INT pj, INT *line, INT verbose_level);
	void points_on_line_by_coordinates(INT pi, INT pj, INT *pt_coords, INT verbose_level);
	void lines_on_point(INT pt, INT *line_pencil_point_ranks, INT verbose_level);
	void lines_on_point_by_line_rank(INT pt, INT *line_pencil_line_ranks, INT verbose_level);
	void list_points_by_type(INT verbose_level);
	void list_points_of_given_type(INT t, INT verbose_level);
	void list_all_points_vs_points(INT verbose_level);
	void list_points_vs_points(INT t1, INT t2, INT verbose_level);
	void test_Siegel(INT index, INT verbose_level);
	void make_initial_partition(partitionstack &S, INT verbose_level);
	void point_to_line_map(INT size, INT *point_ranks, INT *&line_vector, INT verbose_level);
	void move_points_by_ranks_in_place(INT pt_from, INT pt_to, 
		INT nb, INT *ranks, INT verbose_level);
	void move_points_by_ranks(INT pt_from, INT pt_to, 
		INT nb, INT *input_ranks, INT *output_ranks, INT verbose_level);
	void move_points(INT pt_from, INT pt_to, 
		INT nb, INT *input_coords, INT *output_coords, INT verbose_level);
	INT BLT_test_full(INT size, INT *set, INT verbose_level);
	INT BLT_test(INT size, INT *set, INT verbose_level);
	INT collinearity_test(INT size, INT *set, INT verbose_level);
	
	// orthogonal_init.C:
	orthogonal();
	~orthogonal();
	void init(INT epsilon, INT n, finite_field *F, INT verbose_level);
	void init_parabolic(INT verbose_level);
	void init_parabolic_even(INT verbose_level);
	void init_parabolic_odd(INT verbose_level);
	void print_minus_square_tables();
	void init_hyperbolic(INT verbose_level);
	void print_schemes();
	void fill(INT *M, INT i, INT j, INT a);
	
	
	// orthogonal_hyperbolic.C:
	INT hyperbolic_type_and_index_to_point_rk(INT type, INT index);
	void hyperbolic_point_rk_to_type_and_index(INT rk, INT &type, INT &index);
	void hyperbolic_unrank_line(INT &p1, INT &p2, INT rk, INT verbose_level);
	INT hyperbolic_rank_line(INT p1, INT p2, INT verbose_level);
	void unrank_line_L1(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line_L1(INT p1, INT p2, INT verbose_level);
	void unrank_line_L2(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line_L2(INT p1, INT p2, INT verbose_level);
	void unrank_line_L3(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line_L3(INT p1, INT p2, INT verbose_level);
	void unrank_line_L4(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line_L4(INT p1, INT p2, INT verbose_level);
	void unrank_line_L5(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line_L5(INT p1, INT p2, INT verbose_level);
	void unrank_line_L6(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line_L6(INT p1, INT p2, INT verbose_level);
	void unrank_line_L7(INT &p1, INT &p2, INT index, INT verbose_level);
	INT rank_line_L7(INT p1, INT p2, INT verbose_level);
	void hyperbolic_canonical_points_of_line(INT line_type, INT pt1, INT pt2, INT &cpt1, INT &cpt2, INT verbose_level);
	void canonical_points_L1(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void canonical_points_L2(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void canonical_points_L3(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void canonical_points_L4(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void canonical_points_L5(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void canonical_points_L6(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void canonical_points_L7(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	INT hyperbolic_line_type_given_point_types(INT pt1, INT pt2, INT pt1_type, INT pt2_type);
	INT hyperbolic_decide_P1(INT pt1, INT pt2);
	INT hyperbolic_decide_P2(INT pt1, INT pt2);
	INT hyperbolic_decide_P3(INT pt1, INT pt2);
	INT find_root_hyperbolic(INT rk2, INT m, INT verbose_level);
	// m = Witt index
	void find_root_hyperbolic_xyz(INT rk2, INT m, INT *x, INT *y, INT *z, INT verbose_level);
	INT evaluate_hyperbolic_quadratic_form(INT *v, INT stride, INT m);
	INT evaluate_hyperbolic_bilinear_form(INT *u, INT *v, INT stride, INT m);

	// orthogonal_parabolic.C:
	INT parabolic_type_and_index_to_point_rk(INT type, INT index, INT verbose_level);
	INT parabolic_even_type_and_index_to_point_rk(INT type, INT index, INT verbose_level);
	void parabolic_even_type1_index_to_point(INT index, INT *v);
	void parabolic_even_type2_index_to_point(INT index, INT *v);
	INT parabolic_odd_type_and_index_to_point_rk(INT type, INT index, INT verbose_level);
	void parabolic_odd_type1_index_to_point(INT index, INT *v, INT verbose_level);
	void parabolic_odd_type2_index_to_point(INT index, INT *v, INT verbose_level);
	void parabolic_point_rk_to_type_and_index(INT rk, INT &type, INT &index, INT verbose_level);
	void parabolic_even_point_rk_to_type_and_index(INT rk, INT &type, INT &index, INT verbose_level);
	void parabolic_even_point_to_type_and_index(INT *v, INT &type, INT &index, INT verbose_level);
	void parabolic_odd_point_rk_to_type_and_index(INT rk, INT &type, INT &index, INT verbose_level);
	void parabolic_odd_point_to_type_and_index(INT *v, INT &type, INT &index, INT verbose_level);

	void parabolic_neighbor51_odd_unrank(INT index, INT *v, INT verbose_level);
	INT parabolic_neighbor51_odd_rank(INT *v, INT verbose_level);
	void parabolic_neighbor52_odd_unrank(INT index, INT *v, INT verbose_level);
	INT parabolic_neighbor52_odd_rank(INT *v, INT verbose_level);
	void parabolic_neighbor52_even_unrank(INT index, INT *v, INT verbose_level);
	INT parabolic_neighbor52_even_rank(INT *v, INT verbose_level);
	void parabolic_neighbor34_unrank(INT index, INT *v, INT verbose_level);
	INT parabolic_neighbor34_rank(INT *v, INT verbose_level);
	void parabolic_neighbor53_unrank(INT index, INT *v, INT verbose_level);
	INT parabolic_neighbor53_rank(INT *v, INT verbose_level);
	void parabolic_neighbor54_unrank(INT index, INT *v, INT verbose_level);
	INT parabolic_neighbor54_rank(INT *v, INT verbose_level);
	

	void parabolic_unrank_line(INT &p1, INT &p2, INT rk, INT verbose_level);
	INT parabolic_rank_line(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L1_even(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L1_even(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L1_odd(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L1_odd(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L2_even(INT &p1, INT &p2, INT index, INT verbose_level);
	void parabolic_unrank_line_L2_odd(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L2_even(INT p1, INT p2, INT verbose_level);
	INT parabolic_rank_line_L2_odd(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L3(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L3(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L4(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L4(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L5(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L5(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L6(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L6(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L7(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L7(INT p1, INT p2, INT verbose_level);
	void parabolic_unrank_line_L8(INT &p1, INT &p2, INT index, INT verbose_level);
	INT parabolic_rank_line_L8(INT p1, INT p2, INT verbose_level);
	INT parabolic_line_type_given_point_types(INT pt1, INT pt2, 
		INT pt1_type, INT pt2_type, INT verbose_level);
	INT parabolic_decide_P11_odd(INT pt1, INT pt2);
	INT parabolic_decide_P22_even(INT pt1, INT pt2);
	INT parabolic_decide_P22_odd(INT pt1, INT pt2);
	INT parabolic_decide_P33(INT pt1, INT pt2);
	INT parabolic_decide_P35(INT pt1, INT pt2);
	INT parabolic_decide_P45(INT pt1, INT pt2);
	INT parabolic_decide_P44(INT pt1, INT pt2);
	void find_root_parabolic_xyz(INT rk2, INT *x, INT *y, INT *z, INT verbose_level);
	INT find_root_parabolic(INT rk2, INT verbose_level);
	void Siegel_move_forward_by_index(INT rk1, INT rk2, INT *v, INT *w, INT verbose_level);
	void Siegel_move_backward_by_index(INT rk1, INT rk2, INT *w, INT *v, INT verbose_level);
	void Siegel_move_forward(INT *v1, INT *v2, INT *v3, INT *v4, INT verbose_level);
	void Siegel_move_backward(INT *v1, INT *v2, INT *v3, INT *v4, INT verbose_level);
	void parabolic_canonical_points_of_line(INT line_type, INT pt1, INT pt2, 
		INT &cpt1, INT &cpt2, INT verbose_level);
	void parabolic_canonical_points_L1_even(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void parabolic_canonical_points_separate_P5(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void parabolic_canonical_points_L3(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void parabolic_canonical_points_L7(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	void parabolic_canonical_points_L8(INT pt1, INT pt2, INT &cpt1, INT &cpt2);
	INT evaluate_parabolic_bilinear_form(INT *u, INT *v, INT stride, INT m);
	void parabolic_point_normalize(INT *v, INT stride, INT n);
	void parabolic_normalize_point_wrt_subspace(INT *v, INT stride);
	void parabolic_point_properties(INT *v, INT stride, INT n, 
		INT &f_start_with_one, INT &value_middle, INT &value_end, 
		INT verbose_level);
	INT parabolic_is_middle_dependent(INT *vec1, INT *vec2);

	

	// orthogonal_util.C:
	INT test_if_minimal_on_line(INT *v1, INT *v2, INT *v3);
	void find_minimal_point_on_line(INT *v1, INT *v2, INT *v3);
	void zero_vector(INT *u, INT stride, INT len);
	INT is_zero_vector(INT *u, INT stride, INT len);
	void change_form_value(INT *u, INT stride, INT m, INT multiplyer);
	void scalar_multiply_vector(INT *u, INT stride, INT len, INT multiplyer);
	INT last_non_zero_entry(INT *u, INT stride, INT len);
	void Siegel_map_between_singular_points(INT *T, 
		INT rk_from, INT rk_to, INT root, INT verbose_level);
	void Siegel_map_between_singular_points_hyperbolic(INT *T, 
		INT rk_from, INT rk_to, INT root, INT m, INT verbose_level);
	void Siegel_Transformation(INT *T, 
		INT rk_from, INT rk_to, INT root, 
		INT verbose_level);
		// root is not perp to from and to.
	void Siegel_Transformation2(INT *T, 
		INT rk_from, INT rk_to, INT root, 
		INT *B, INT *Bv, INT *w, INT *z, INT *x,
		INT verbose_level);
	void Siegel_Transformation3(INT *T, 
		INT *from, INT *to, INT *root, 
		INT *B, INT *Bv, INT *w, INT *z, INT *x,
		INT verbose_level);
	void random_generator_for_orthogonal_group(
		INT f_action_is_semilinear, 
		INT f_siegel, 
		INT f_reflection, 
		INT f_similarity,
		INT f_semisimilarity, 
		INT *Mtx, INT verbose_level);
	void create_random_Siegel_transformation(INT *Mtx, INT verbose_level);
		// Only makes a n x n matrix. Does not put a semilinear component.
	void create_random_semisimilarity(INT *Mtx, INT verbose_level);
	void create_random_similarity(INT *Mtx, INT verbose_level);
		// Only makes a d x d matrix. Does not put a semilinear component.
	void create_random_orthogonal_reflection(INT *Mtx, INT verbose_level);
		// Only makes a d x d matrix. Does not put a semilinear component.
	void make_orthogonal_reflection(INT *M, INT *z, INT verbose_level);
	void make_Siegel_Transformation(INT *M, INT *v, INT *u, 
		INT n, INT *Gram, INT verbose_level);
		// if u is singular and v \in \la u \ra^\perp, then
		// \pho_{u,v}(x) := x + \beta(x,v) u - \beta(x,u) v - Q(v) \beta(x,u) u
		// is called the Siegel transform (see Taylor p. 148)
		// Here Q is the quadratic form and \beta is the corresponding bilinear form
	void unrank_S(INT *v, INT stride, INT m, INT rk);
	INT rank_S(INT *v, INT stride, INT m);
	void unrank_N(INT *v, INT stride, INT m, INT rk);
	INT rank_N(INT *v, INT stride, INT m);
	void unrank_N1(INT *v, INT stride, INT m, INT rk);
	INT rank_N1(INT *v, INT stride, INT m);
	void unrank_Sbar(INT *v, INT stride, INT m, INT rk);
	INT rank_Sbar(INT *v, INT stride, INT m);
	void unrank_Nbar(INT *v, INT stride, INT m, INT rk);
	INT rank_Nbar(INT *v, INT stride, INT m);
	void normalize_point(INT *v, INT stride);
	INT triple_is_collinear(INT pt1, INT pt2, INT pt3);
	INT is_minus_square(INT i);
	INT is_ending_dependent(INT *vec1, INT *vec2);
	void Gauss_step(INT *v1, INT *v2, INT len, INT idx);
		// afterwards: v2[idx] = 0 and v2,v1 span the same space as before
};

// ####################################################################################
// vector_hashing.C:
// ####################################################################################

class vector_hashing {

public:
	INT data_size;
	INT N;
	INT bit_length;
	INT *vector_data;
	INT *H;
	INT *H_sorted;
	INT *perm;
	INT *perm_inv;
	INT nb_types;
	INT *type_first;
	INT *type_len;
	INT *type_value;


	vector_hashing();
	~vector_hashing();
	void allocate(INT data_size, INT N, INT bit_length);
	void compute_tables(INT verbose_level);
	void print();
	INT rank(INT *data);
	void unrank(INT rk, INT *data);
};


// ####################################################################################
// geometry.C:
// ####################################################################################

INT hyperoval_nb_reps(INT q);
INT *hyperoval_representative(INT q, INT i);
// i starts from 0
void hyperoval_gens(INT q, INT i, INT *&data, INT &nb_gens, INT &data_size, const BYTE *&stab_order);
INT DH_nb_reps(INT k, INT n);
INT *DH_representative(INT k, INT n, INT i);
// i starts from 0
void DH_stab_gens(INT k, INT n, INT i, INT *&data, INT &nb_gens, INT &data_size, const BYTE *&stab_order);
INT TP_nb_reps(INT q, INT k);
INT *TP_representative(INT q, INT k, INT i);
// i starts from 0
void TP_stab_gens(INT q, INT k, INT i, INT *&data, INT &nb_gens, INT &data_size, const BYTE *&stab_order);
// i starts from 0
INT BLT_nb_reps(INT q, INT k);
INT *BLT_representative(INT q, INT no);
// i starts from 0
void BLT_stab_gens(INT q, INT k, INT no, INT *&data, INT &nb_gens, INT &data_size, const BYTE *&stab_order);



//INT *BLT_representative(INT q, INT i);
// i starts from 1
const BYTE *override_polynomial_subfield(INT q);
const BYTE *override_polynomial_extension_field(INT q);
void create_Fisher_BLT_set(INT *Fisher_BLT, INT q, const BYTE *poly_q, const BYTE *poly_Q, INT verbose_level);
void create_Linear_BLT_set(INT *BLT, INT q, const BYTE *poly_q, const BYTE *poly_Q, INT verbose_level);
void create_Mondello_BLT_set(INT *BLT, INT q, const BYTE *poly_q, const BYTE *poly_Q, INT verbose_level);
void print_quadratic_form_list_coded(INT form_nb_terms, 
	INT *form_i, INT *form_j, INT *form_coeff);
void make_Gram_matrix_from_list_coded_quadratic_form(INT n, finite_field &F, 
	INT nb_terms, INT *form_i, INT *form_j, INT *form_coeff, INT *Gram);
void add_term(INT n, finite_field &F, INT &nb_terms, INT *form_i, INT *form_j, INT *form_coeff, INT *Gram, 
	INT i, INT j, INT coeff);
void create_BLT_point(finite_field *F, INT *v5, INT a, INT b, INT c, INT verbose_level);
// creates the point (-b/2,-c,a,-(b^2/4-ac),1) 
// check if it satisfies x_0^2 + x_1x_2 + x_3x_4:
// b^2/4 + (-c)*a + -(b^2/4-ac)
// = b^2/4 -ac -b^2/4 + ac = 0
void create_FTWKB_BLT_set(orthogonal *O, INT *set, INT verbose_level);
void create_K1_BLT_set(orthogonal *O, INT *set, INT verbose_level);
void create_K2_BLT_set(orthogonal *O, INT *set, INT verbose_level);
void create_LP_37_72_BLT_set(orthogonal *O, INT *set, INT verbose_level);
void create_LP_37_4a_BLT_set(orthogonal *O, INT *set, INT verbose_level);
void create_LP_37_4b_BLT_set(orthogonal *O, INT *set, INT verbose_level);


void GlynnI_hyperoval(finite_field *F, INT *&Pts, INT &nb_pts, INT verbose_level);
void GlynnII_hyperoval(finite_field *F, INT *&Pts, INT &nb_pts, INT verbose_level);
void Segre_hyperoval(finite_field *F, INT *&Pts, INT &nb_pts, INT verbose_level);
void Adelaide_hyperoval(subfield_structure *S, INT *&Pts, INT &nb_pts, INT verbose_level);
void Subiaco_oval(finite_field *F, INT *&Pts, INT &nb_pts, INT f_short, INT verbose_level);
void Subiaco_hyperoval(finite_field *F, INT *&Pts, INT &nb_pts, INT verbose_level);
INT OKeefe_Penttila_32(finite_field *F, INT t);
INT Subiaco64_1(finite_field *F, INT t);
// needs the field generated by beta with beta^6=beta+1
INT Subiaco64_2(finite_field *F, INT t);
// needs the field generated by beta with beta^6=beta+1
INT Adelaide64(finite_field *F, INT t);
// needs the field generated by beta with beta^6=beta+1
void LunelliSce(finite_field *Fq, INT *pts18, INT verbose_level);
INT LunelliSce_evaluate_cubic1(finite_field *F, INT *v);
INT LunelliSce_evaluate_cubic2(finite_field *F, INT *v);

void plane_invariant(INT q, orthogonal *O, unusual_model *U, 
	INT size, INT *set, 
	INT &nb_planes, INT *&intersection_matrix, 
	INT &Block_size, INT *&Blocks, 
	INT verbose_level);
// using hash values


void create_Law_71_BLT_set(orthogonal *O, INT *set, INT verbose_level);


// ####################################################################################
// unusual.C:
// ####################################################################################


class unusual_model {
public:
	finite_field F, f;
	INT q;
	INT qq;
	INT alpha;
	INT T_alpha, N_alpha;
	INT nb_terms, *form_i, *form_j, *form_coeff, *Gram;
	INT r_nb_terms, *r_form_i, *r_form_j, *r_form_coeff, *r_Gram;
	INT rr_nb_terms, *rr_form_i, *rr_form_j, *rr_form_coeff, *rr_Gram;
	INT hyperbolic_basis[4 * 4];
	INT hyperbolic_basis_inverse[4 * 4];
	INT basis[4 * 4];
	INT basis_subspace[2 * 2];
	INT *M;
	INT *components, *embedding, *pair_embedding;
		// data computed by F.subfield_embedding_2dimensional
	
	unusual_model();
	~unusual_model();
	void setup_sum_of_squares(INT q, const BYTE *poly_q, const BYTE *poly_Q, INT verbose_level);
	void setup(INT q, const BYTE *poly_q, const BYTE *poly_Q, INT verbose_level);
	void setup2(INT q, const BYTE *poly_q, const BYTE *poly_Q, INT f_sum_of_squares, INT verbose_level);
	void convert_to_ranks(INT n, INT *unusual_coordinates, INT *ranks, INT verbose_level);
	void convert_from_ranks(INT n, INT *ranks, INT *unusual_coordinates, INT verbose_level);
	INT convert_to_rank(INT *unusual_coordinates, INT verbose_level);
	void convert_from_rank(INT rank, INT *unusual_coordinates, INT verbose_level);
	void convert_to_usual(INT n, INT *unusual_coordinates, 
		INT *usual_coordinates, INT verbose_level);
	void create_Fisher_BLT_set(INT *Fisher_BLT, INT verbose_level);
	void convert_from_usual(INT n, INT *usual_coordinates, 
		INT *unusual_coordinates, INT verbose_level);
	void create_Linear_BLT_set(INT *BLT, INT verbose_level);
	void create_Mondello_BLT_set(INT *BLT, INT verbose_level);
	INT N2(INT a);
	INT T2(INT a);
	INT quadratic_form(INT a, INT b, INT c, INT verbose_level);
	INT bilinear_form(INT a1, INT b1, INT c1, INT a2, INT b2, INT c2, INT verbose_level);
	void print_coordinates_detailed_set(INT *set, INT len);
	void print_coordinates_detailed(INT pt, INT cnt);
	INT build_candidate_set(orthogonal &O, INT q, 
		INT gamma, INT delta, INT m, INT *Set, 
		INT f_second_half, INT verbose_level);
	INT build_candidate_set_with_offset(orthogonal &O, INT q, 
		INT gamma, INT delta, INT offset, INT m, INT *Set, 
		INT f_second_half, INT verbose_level);
	INT build_candidate_set_with_or_without_test(orthogonal &O, INT q, 
		INT gamma, INT delta, INT offset, INT m, INT *Set, 
		INT f_second_half, INT f_test, INT verbose_level);
	INT create_orbit_of_psi(orthogonal &O, INT q, 
		INT gamma, INT delta, INT m, INT *Set, 
		INT f_test, INT verbose_level);
	void transform_matrix_unusual_to_usual(orthogonal *O, 
		INT *M4, INT *M5, INT verbose_level);
	void transform_matrix_usual_to_unusual(orthogonal *O, 
		INT *M5, INT *M4, INT verbose_level);

	void parse_4by4_matrix(INT *M4, 
		INT &a, INT &b, INT &c, INT &d, 
		INT &f_semi1, INT &f_semi2, INT &f_semi3, INT &f_semi4);
	void create_4by4_matrix(INT *M4, 
		INT a, INT b, INT c, INT d, 
		INT f_semi1, INT f_semi2, INT f_semi3, INT f_semi4, 
		INT verbose_level);
	void print_2x2(INT *v, INT *f_semi);
	void print_M5(orthogonal *O, INT *M5);
};

// ####################################################################################
// grassmann.C:
// ####################################################################################

class grassmann {
public:
	INT n, k, q;
	finite_field *F;
	INT *base_cols;
	INT *coset;
	INT *M;
	grassmann *G;

	grassmann();
	~grassmann();
	void init(INT n, INT k, finite_field *F, INT verbose_level);
	INT nb_of_subspaces(INT verbose_level);
	INT nb_points_covered(INT verbose_level);
	void points_covered(INT *the_points, INT verbose_level);
	void unrank_INT_here(INT *Mtx, INT rk, INT verbose_level);
	INT rank_INT_here(INT *Mtx, INT verbose_level);
	void unrank_INT(INT rk, INT verbose_level);
	INT rank_INT(INT verbose_level);
	void unrank_longinteger_here(INT *Mtx, longinteger_object &rk, INT verbose_level);
	void rank_longinteger_here(INT *Mtx, longinteger_object &rk, INT verbose_level);
	void unrank_longinteger(longinteger_object &rk, INT verbose_level);
	void rank_longinteger(longinteger_object &r, INT verbose_level);
	void print();
	INT dimension_of_join(INT rk1, INT rk2, INT verbose_level);
	void unrank_INT_here_and_extend_basis(INT *Mtx, INT rk, INT verbose_level);
		// Mtx must be n x n
	void line_regulus_in_PG_3_q(INT *&regulus, INT &regulus_size, INT verbose_level);
};


// ####################################################################################
// grassmann_embedded.C:
// ####################################################################################

class grassmann_embedded {
public:
	INT big_n, n, k, q;
	finite_field *F;
	grassmann *G; // only a reference, not freed
	INT *M; // [n * big_n] the original matrix
	INT *M_Gauss; // [n * big_n] the echeolon form (RREF)
	INT *transform; // [n * n] the transformation matrix
	INT *base_cols; // [n]
	INT *embedding; // [big_n - n]
	INT *Tmp1; // [big_n]
	INT *Tmp2; // [big_n]
	INT *Tmp3; // [big_n]
	INT *tmp_M1; // [n * n]
	INT *tmp_M2; // [n * n]
	INT degree; // q_binomial n choose k


	grassmann_embedded();
	~grassmann_embedded();
	void init(INT big_n, INT n, grassmann *G, INT *M, INT verbose_level);
		// M is n x big_n
	void unrank_INT(INT *subspace_basis, INT rk, INT verbose_level);
		// subspace_basis is k x big_n
	INT rank_INT(INT *subspace_basis, INT verbose_level);
		// subspace_basis is k x big_n
};

// ####################################################################################
// hjelmslev.C:
// ####################################################################################

class hjelmslev {
public:
	INT n, k, q;
	INT n_choose_k_p;
	finite_ring *R; // do not free
	grassmann *G;
	INT *v;
	INT *Mtx;
	INT *base_cols;

	hjelmslev();
	~hjelmslev();
	void null();
	void freeself();
	void init(finite_ring *R, INT n, INT k, INT verbose_level);
	INT number_of_submodules();
	void unrank_INT(INT *M, INT rk, INT verbose_level);
	INT rank_INT(INT *M, INT verbose_level);
};

// ####################################################################################
// hermitian.C:
// ####################################################################################

class hermitian {

public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;

	finite_field *F; // only a reference, not to be freed
	INT Q;
	INT q;
	INT k; // nb_vars

	INT *cnt_N; // [k + 1]
	INT *cnt_N1; // [k + 1]
	INT *cnt_S; // [k + 1]
	INT *cnt_Sbar; // [k + 1]
	
	INT *norm_one_elements; // [q + 1]
	INT *index_of_norm_one_element; // [Q]
	INT alpha; // a primitive element for GF(Q), namely F->p
	INT beta; // alpha^(q+1), a primitive element for GF(q)
	INT *log_beta; // [Q]
	INT *beta_power; // [q - 1]
	
	hermitian();
	~hermitian();
	void null();
	void init(finite_field *F, INT nb_vars, INT verbose_level);
	INT nb_points();
	void unrank_point(INT *v, INT rk);
	INT rank_point(INT *v);
	void list_of_points_embedded_in_PG(INT *&Pts, INT &nb_pts, INT verbose_level);
	void list_all_N(INT verbose_level);
	void list_all_N1(INT verbose_level);
	void list_all_S(INT verbose_level);
	void list_all_Sbar(INT verbose_level);
	INT evaluate_hermitian_form(INT *v, INT len);
	void N_unrank(INT *v, INT len, INT rk, INT verbose_level);
	INT N_rank(INT *v, INT len, INT verbose_level);
	void N1_unrank(INT *v, INT len, INT rk, INT verbose_level);
	INT N1_rank(INT *v, INT len, INT verbose_level);
	void S_unrank(INT *v, INT len, INT rk, INT verbose_level);
	INT S_rank(INT *v, INT len, INT verbose_level);
	void Sbar_unrank(INT *v, INT len, INT rk, INT verbose_level);
	INT Sbar_rank(INT *v, INT len, INT verbose_level);
};




// ####################################################################################
// projective.C:
// ####################################################################################

INT nb_PG_elements(INT n, INT q);
// $\frac{q^{n+1} - 1}{q-1} = \sum_{i=0}^{n} q^i $
INT nb_PG_elements_not_in_subspace(INT n, INT m, INT q);
INT nb_AG_elements(INT n, INT q);
void all_PG_elements_in_subspace(finite_field *F, INT *genma, INT k, INT n, INT *&point_list, INT &nb_points, INT verbose_level);
void all_PG_elements_in_subspace_array_is_given(finite_field *F, INT *genma, INT k, INT n, INT *point_list, INT &nb_points, INT verbose_level);
void display_all_PG_elements(INT n, finite_field &GFq);
void display_all_PG_elements_not_in_subspace(INT n, INT m, finite_field &GFq);
void display_all_AG_elements(INT n, finite_field &GFq);
void PG_element_apply_frobenius(INT n, finite_field &GFq, INT *v, INT f);
void PG_element_normalize(finite_field &GFq, INT *v, INT stride, INT len);
void PG_element_normalize_from_front(finite_field &GFq, INT *v, INT stride, INT len);
void PG_element_rank_modified(finite_field &GFq, INT *v, INT stride, INT len, INT &a);
void PG_element_unrank_modified(finite_field &GFq, INT *v, INT stride, INT len, INT a);
void PG_element_rank_modified_not_in_subspace(finite_field &GFq, INT *v, INT stride, INT len, INT m, INT &a);
void PG_element_unrank_modified_not_in_subspace(finite_field &GFq, INT *v, INT stride, INT len, INT m, INT a);
void AG_element_rank(INT q, INT *v, INT stride, INT len, INT &a);
void AG_element_unrank(INT q, INT *v, INT stride, INT len, INT a);
void AG_element_rank_longinteger(INT q, INT *v, INT stride, INT len, longinteger_object &a);
void AG_element_unrank_longinteger(INT q, INT *v, INT stride, INT len, longinteger_object &a);
INT PG_element_modified_is_in_subspace(INT n, INT m, INT *v);
void PG_element_modified_not_in_subspace_perm(INT n, INT m, 
	finite_field &GFq, INT *orbit, INT *orbit_inv, INT verbose_level);
INT PG2_line_on_point_unrank(finite_field &GFq, INT *v1, INT rk);
void PG2_line_on_point_unrank_second_point(finite_field &GFq, INT *v1, INT *v2, INT rk);
INT PG2_line_rank(finite_field &GFq, INT *v1, INT *v2, INT stride);
void PG2_line_unrank(finite_field &GFq, INT *v1, INT *v2, INT stride, INT line_rk);
void diagonal_orbit_perm(INT n, finite_field &GFq, INT *orbit, INT *orbit_inv, INT verbose_level);
void frobenius_orbit_perm(INT n, finite_field &GFq, INT *orbit, INT *orbit_inv, INT verbose_level);
void test_PG(INT n, INT q);
void line_through_two_points(finite_field &GFq, INT len, INT pt1, INT pt2, INT *line);
void print_set_in_affine_plane(finite_field &GFq, INT len, INT *S);
INT consecutive_ones_property_in_affine_plane(ostream &ost, finite_field &GFq, INT len, INT *S);
void oval_polynomial(finite_field &GFq, INT *S, unipoly_domain &D, unipoly_object &poly, 
	INT verbose_level);
#if 0
void regular_oval_point_list(finite_field &GFq, INT &len, INT *&pts);
void segre_oval_point_list(finite_field &GFq, INT &len, INT *&pts);
void translation_oval_point_list(finite_field &GFq, INT &len, INT *&pts);
void penttila_okeefe_oval_point_list(finite_field &GFq, INT &len, INT *&pts);
void oval_32_5_point_list(finite_field &GFq, INT &len, INT *&pts);
void oval_32_10_point_list(finite_field &GFq, INT &len, INT *&pts);
INT regular_oval_point(finite_field &GFq, INT i);
INT segre_oval_point(finite_field &GFq, INT i);
#endif
INT line_intersection_with_oval(finite_field &GFq, INT *f_oval_point, INT line_rk, 
	INT verbose_level);
INT get_base_line(finite_field &GFq, INT plane1, INT plane2, INT verbose_level);
void translation_in_AG(finite_field &GFq, INT n, INT i, INT a, INT *perm, INT *v, INT verbose_level);
// v[n] needs to be allocated 
// p[q^n] needs to be allocated
void frobenius_in_AG(finite_field &GFq, INT n, INT *perm, INT *v, INT verbose_level);
// v[n] needs to be allocated 
// p[q^n] needs to be allocated
void frobenius_in_PG(finite_field &GFq, INT n, INT *perm, INT *v, INT verbose_level);
// v[n + 1] needs to be allocated 
// p[q^n+...+q+1] needs to be allocated
void AG_representation_of_matrix(finite_field &GFq, INT n, INT f_from_the_right, 
	INT *M, INT *v, INT *w, INT *perm, INT verbose_level);
// perm[q^n] needs to be already allocated
void AG_representation_one_dimensional(finite_field &GFq, 
	INT a, INT *perm, INT verbose_level);
// perm[q] needs to be already allocated
INT nb_generators_affine_translations(finite_field &GFq, INT n);
void generators_affine_translations(finite_field &GFq, INT n, INT *perms, INT verbose_level);
// primes[n * d] needs to be allocated, where d = q^n
void generators_AGL1xAGL1_subdirect1(finite_field &GFq1, finite_field &GFq2, INT u, INT v, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_AGL1q(finite_field &GFq, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_AGL1q_subgroup(finite_field &GFq, INT index_in_multiplicative_group, 
	INT &nb_perms, INT *&perms, INT verbose_level);
void generators_AGL1_x_AGL1(finite_field &GFq1, finite_field &GFq2, INT &deg, 
	INT &nb_perms, INT *&perms, INT verbose_level);
void generators_AGL1_x_AGL1_extension(finite_field &GFq1, finite_field &GFq2, INT u, INT v, 
	INT &deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_AGL1_x_AGL1_extended_once(finite_field &F1, finite_field &F2, INT u, INT v, 
	INT &deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_AGL1_x_AGL1_extended_twice(finite_field &F1, finite_field &F2, INT u1, INT v1, 
	INT u2, INT v2, INT &deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_symmetric_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_cyclic_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_dihedral_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_dihedral_involution(INT deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_identity_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level);
void generators_direct_product(INT deg1, INT nb_perms1, INT *perms1, 
	INT deg2, INT nb_perms2, INT *perms2, 
	INT &deg3, INT &nb_perms3, INT *&perms3, 
	INT verbose_levels);
void generators_concatenate(INT deg1, INT nb_perms1, INT *perms1, 
	INT deg2, INT nb_perms2, INT *perms2, 
	INT &deg3, INT &nb_perms3, INT *&perms3, 
	INT verbose_level);
void O4_isomorphism_4to2(finite_field *F, INT *At, INT *As, INT &f_switch, INT *B, INT verbose_level);
void O4_isomorphism_2to4(finite_field *F, INT *At, INT *As, INT f_switch, INT *B);
void O4_grid_coordinates_rank(finite_field &F, INT x1, INT x2, INT x3, INT x4, INT &grid_x, INT &grid_y, INT verbose_level);
void O4_grid_coordinates_unrank(finite_field &F, INT &x1, INT &x2, INT &x3, INT &x4, INT grid_x, INT grid_y, INT verbose_level);
void O4_find_tangent_plane(finite_field &F, INT pt_x1, INT pt_x2, INT pt_x3, INT pt_x4, 
	INT *tangent_plane, INT verbose_level);
INT PHG_element_normalize(finite_ring &R, INT *v, INT stride, INT len);
// last unit element made one
INT PHG_element_normalize_from_front(finite_ring &R, INT *v, INT stride, INT len);
// first non unit element made one
INT PHG_element_rank(finite_ring &R, INT *v, INT stride, INT len);
void PHG_element_unrank(finite_ring &R, INT *v, INT stride, INT len, INT rk);
INT nb_PHG_elements(INT n, finite_ring &R);
void display_all_PHG_elements(INT n, INT q);
INT matrix_group_base_len_projective_group(INT n, INT q, INT f_semilinear, INT verbose_level);
INT matrix_group_base_len_affine_group(INT n, INT q, INT f_semilinear, INT verbose_level);
INT matrix_group_base_len_general_linear_group(INT n, INT q, INT f_semilinear, INT verbose_level);
void projective_matrix_group_base_and_orbits(INT n, 
	finite_field *F, INT f_semilinear, 
	INT base_len, INT degree, 
	INT *base, INT *transversal_length, 
	INT **orbit, INT **orbit_inv, 
	INT verbose_level);
void affine_matrix_group_base_and_transversal_length(INT n, 
	finite_field *F, INT f_semilinear, 
	INT base_len, INT degree, 
	INT *base, INT *transversal_length, 
	INT verbose_level);
void general_linear_matrix_group_base_and_transversal_length(INT n, 
	finite_field *F, INT f_semilinear, 
	INT base_len, INT degree, 
	INT *base, INT *transversal_length, 
	INT verbose_level);
void strong_generators_for_projective_linear_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT *&data, INT &size, INT &nb_gens, 
	INT verbose_level);
void strong_generators_for_affine_linear_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT *&data, INT &size, INT &nb_gens, 
	INT verbose_level);
void strong_generators_for_general_linear_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT *&data, INT &size, INT &nb_gens, 
	INT verbose_level);



// ####################################################################################
// combinatorics.C:
// ####################################################################################

INT INT_factorial(INT a);
INT Kung_mue_i(INT *part, INT i, INT m);
void partition_dual(INT *part, INT *dual_part, INT n, INT verbose_level);
void make_all_partitions_of_n(INT n, INT *&Table, INT &nb, INT verbose_level);
INT count_all_partitions_of_n(INT n);
INT partition_first(INT *v, INT n);
INT partition_next(INT *v, INT n);
void partition_print(ostream &ost, INT *v, INT n);
INT INT_vec_is_regular_word(INT *v, INT len, INT q);
// Returns TRUE if the word v of length len is regular, i.~e. 
// lies in an orbit of length $len$ under the action of the cyclic group 
// $C_{len}$ acting on the coordinates. 
// Lueneburg~\cite{Lueneburg87a} p. 118.
// v is a vector over $\{0, 1, \ldots , q-1\}$
INT INT_vec_first_regular_word(INT *v, INT len, INT Q, INT q);
INT INT_vec_next_regular_word(INT *v, INT len, INT Q, INT q);
INT is_subset_of(INT *A, INT sz_A, INT *B, INT sz_B);
INT set_find(INT *elts, INT size, INT a);
void set_complement(INT *elts, INT &size, INT *elts_complement, INT &size_complement, INT n);
void set_add_elements(INT *elts, INT &size, INT *elts_to_add, INT nb_elts_to_add);
void set_add_element(INT *elts, INT &size, INT a);
void set_delete_elements(INT *elts, INT &size, INT *elts_to_delete, INT nb_elts_to_delete);
void set_delete_element(INT *elts, INT &size, INT a);
INT compare_lexicographically(INT a_len, INT *a, INT b_len, INT *b);
INT INT_n_choose_k(INT n, INT k);
void make_t_k_incidence_matrix(INT v, INT t, INT k, INT &m, INT &n, INT *&M, 
	INT verbose_level);
void print_k_subsets_by_rank(ostream &ost, INT v, INT k);
INT f_is_subset_of(INT v, INT t, INT k, INT rk_t_subset, INT rk_k_subset);
INT rank_k_subset(INT *set, INT n, INT k);
void unrank_k_subset(INT rk, INT *set, INT n, INT k);
INT first_k_subset(INT *set, INT n, INT k);
INT next_k_subset(INT *set, INT n, INT k);
INT next_k_subset_at_level(INT *set, INT n, INT k, INT backtrack_level);
void subset_permute_up_front(INT n, INT k, INT *set, INT *k_subset_idx, INT *permuted_set);
INT ij2k(INT i, INT j, INT n);
void k2ij(INT k, INT & i, INT & j, INT n);
INT ijk2h(INT i, INT j, INT k, INT n);
void h2ijk(INT h, INT &i, INT &j, INT &k, INT n);
void random_permutation(INT *random_permutation, INT n);
void perm_move(INT *from, INT *to, INT n);
void perm_identity(INT *a, INT n);
void perm_mult(INT *a, INT *b, INT *c, INT n);
void perm_conjugate(INT *a, INT *b, INT *c, INT n);
// c := a^b = b^-1 * a * b
void perm_inverse(INT *a, INT *b, INT n);
// b := a^-1
void perm_raise(INT *a, INT *b, INT e, INT n);
// b := a^e (e >= 0)
void perm_direct_product(INT n1, INT n2, INT *perm1, INT *perm2, INT *perm3);
void perm_print_list(ostream &ost, INT *a, INT n);
void perm_print_list_offset(ostream &ost, INT *a, INT n, INT offset);
void perm_print_product_action(ostream &ost, INT *a, INT m_plus_n, INT m, INT offset, INT f_cycle_length);
void perm_print(ostream &ost, INT *a, INT n);
void perm_print_with_cycle_length(ostream &ost, INT *a, INT n);
void perm_print_counting_from_one(ostream &ost, INT *a, INT n);
void perm_print_offset(ostream &ost, INT *a, INT n, INT offset, INT f_cycle_length, 
	INT f_max_cycle_length, INT max_cycle_length, INT f_orbit_structure);
INT perm_order(INT *a, INT n);
INT perm_signum(INT *perm, INT n);
void first_lehmercode(INT n, INT *v);
INT next_lehmercode(INT n, INT *v);
void lehmercode_to_permutation(INT n, INT *code, INT *perm);
INT disjoint_binary_representation(INT u, INT v);
INT hall_test(INT *A, INT n, INT kmax, INT *memo, INT verbose_level);
INT philip_hall_test(INT *A, INT n, INT k, INT *memo, INT verbose_level);
// memo points to free memory of n INT's
INT philip_hall_test_dual(INT *A, INT n, INT k, INT *memo, INT verbose_level);
// memo points to free memory of n INT's
void print_01_matrix_with_stars(ostream &ost, INT *A, INT m, INT n);
void print_INT_matrix(ostream &ost, INT *A, INT m, INT n);
INT create_roots_H4(finite_field *F, INT *roots);
INT generalized_binomial(INT n, INT k, INT q);
void print_tableau(INT *Tableau, INT l1, INT l2, INT *row_parts, INT *col_parts);

// ####################################################################################
// number_theory.C:
// ####################################################################################


INT INT_abs(INT a);
INT irem(INT a, INT m);
INT gcd_INT(INT m, INT n);
INT i_power_j(INT i, INT j);
INT order_mod_p(INT a, INT p);
INT INT_log2(INT n);
INT INT_log10(INT n);
INT INT_logq(INT n, INT q);
// returns the number of digits in base q representation
INT is_strict_prime_power(INT q);
// assuming that q is a prime power, this fuction tests 
// whether or not q is a strict prime power
INT is_prime(INT p);
INT is_prime_power(INT q, INT &p, INT &h);
INT smallest_primedivisor(INT n);
//Computes the smallest prime dividing $n$. 
//The algorithm is based on Lueneburg~\cite{Lueneburg87a}.
INT sp_ge(INT n, INT p_min);
INT factor_INT(INT a, INT *&primes, INT *&exponents);
void factor_prime_power(INT q, INT &p, INT &e);
INT primitive_root(INT p, INT verbose_level);
INT Legendre(INT a, INT p, INT verbose_level);
INT Jacobi(INT a, INT m, INT verbose_level);
INT ny2(INT x, INT &x1);
INT ny_p(INT n, INT p);
INT sqrt_mod_simple(INT a, INT p);
void print_factorization(INT nb_primes, INT *primes, INT *exponents);
void print_longfactorization(INT nb_primes, longinteger_object *primes, INT *exponents);
INT euler_function(INT n);
void INT_add_fractions(INT at, INT ab, INT bt, INT bb, INT &ct, INT &cb, INT verbose_level);
void INT_mult_fractions(INT at, INT ab, INT bt, INT bb, INT &ct, INT &cb, INT verbose_level);

// ####################################################################################
// util.C:
// ####################################################################################


INT INT_vec_count_number_of_nonzero_entries(INT *v, INT len);
INT INT_vec_find_first_nonzero_entry(INT *v, INT len);
void INT_vec_zero(INT *v, INT len);
void INT_vec_mone(INT *v, INT len);
void INT_vec_copy(INT *from, INT *to, INT len);
UBYTE *bitvector_allocate(INT length);
void bitvector_m_ii(UBYTE *bitvec, INT i, INT a);
INT bitvector_s_i(UBYTE *bitvec, INT i);
// returns 0 or 1
INT INT_vec_hash(INT *data, INT len);
INT INT_vec_hash_after_sorting(INT *data, INT len);
const BYTE *plus_minus_string(INT epsilon);
const BYTE *plus_minus_letter(INT epsilon);
void INT_vec_complement(INT *v, INT n, INT k);
// computes the complement to v + k (v must be allocated to n lements)
void INT_vec_complement(INT *v, INT *w, INT n, INT k);
// computes the complement of v[k] w[n - k] 
void INT_vec_init5(INT *v, INT a0, INT a1, INT a2, INT a3, INT a4);
void dump_memory_chain(void *allocated_objects);
void print_vector(ostream &ost, INT *v, int size);
INT INT_vec_minimum(INT *v, INT len);
INT INT_vec_maximum(INT *v, INT len);
//void INT_vec_copy(INT len, INT *from, INT *to);
INT INT_vec_first_difference(INT *p, INT *q, INT len);
void itoa(char *p, INT len_of_p, INT i);
void BYTE_swap(BYTE *p, BYTE *q, INT len);
void INT_vec_distribution_compute_and_print(ostream &ost, INT *v, INT v_len);
void INT_vec_distribution(INT *v, INT len_v, INT *&val, INT *&mult, INT &len);
void INT_distribution_print(ostream &ost, INT *val, INT *mult, INT len);
void INT_swap(INT& x, INT& y);
void INT_set_print(INT *v, INT len);
void INT_set_print(ostream &ost, INT *v, INT len);
void INT_vec_print(ostream &ost, INT *v, INT len);
void INT_vec_print_fully(ostream &ost, INT *v, INT len);
void integer_vec_print(ostream &ost, int *v, int len);
void print_integer_matrix(ostream &ost, INT *p, INT m, INT n);
void print_integer_matrix_width(ostream &ost, INT *p, INT m, INT n, INT dim_n, INT w);
void print_01_matrix_tex(ostream &ost, INT *p, INT m, INT n);
void print_integer_matrix_tex(ostream &ost, INT *p, INT m, INT n);
void print_integer_matrix_tex_block_by_block(ostream &ost, INT *p, INT m, INT n, INT block_width);
void print_big_integer_matrix_tex(ostream &ost, INT *p, INT m, INT n);
void INT_matrix_make_block_matrix_2x2(INT *Mtx, INT k, INT *A, INT *B, INT *C, INT *D);
// makes the 2k x 2k block matrix 
// (A B)
// (C D)
void INT_matrix_delete_column_in_place(INT *Mtx, INT k, INT n, INT pivot);
// afterwards, the matrix is k x (n - 1)
void INT_matrix_print(INT *p, INT m, INT n);
void INT_matrix_print(INT *p, INT m, INT n, INT w);
void INT_matrix_print_tex(ostream &ost, INT *p, INT m, INT n);
void UBYTE_print_bitwise(ostream &ost, UBYTE u);
void UBYTE_move(UBYTE *p, UBYTE *q, INT len);
void INT_submatrix_all_rows(INT *A, INT m, INT n, INT nb_cols, INT *cols, INT *B);
void INT_submatrix_all_cols(INT *A, INT m, INT n, INT nb_rows, INT *rows, INT *B);
void INT_submatrix(INT *A, INT m, INT n, INT nb_rows, INT *rows, INT nb_cols, INT *cols, INT *B);
void INT_matrix_transpose(INT n, INT *A);
void INT_matrix_transpose(INT *M, INT m, INT n, INT *Mt);
// Mt must point to the right amount of memory (n * m INT's)
void INT_matrix_shorten_rows(INT *&p, INT m, INT n);
void PINT_matrix_shorten_rows(PINT *&p, INT m, INT n);
void runtime(long *l);
INT os_ticks();
INT os_ticks_system();
INT os_ticks_per_second();
void os_ticks_to_dhms(INT ticks, INT tps, INT &d, INT &h, INT &m, INT &s);
void time_check_delta(ostream &ost, INT dt);
void print_elapsed_time(ostream &ost, INT d, INT h, INT m, INT s);
void time_check(ostream &ost, INT t0);
INT file_size(const BYTE *name);
void delete_file(const BYTE *fname);
void fwrite_INT4(FILE *fp, INT a);
INT fread_INT4(FILE *fp);
void fwrite_UBYTEs(FILE *fp, UBYTE *p, INT len);
void fread_UBYTEs(FILE *fp, UBYTE *p, INT len);
void latex_head_easy(ostream& ost);
void latex_head_easy_sideways(ostream& ost);
void latex_head(ostream& ost, INT f_book, INT f_title, 
	const BYTE *title, const BYTE *author, 
	INT f_toc, INT f_landscape, INT f_12pt, 
	INT f_enlarged_page, INT f_pagenumbers);
void latex_foot(ostream& ost);
void seed_random_generator_with_system_time();
void seed_random_generator(INT seed);
INT random_integer(INT p);
void print_set(ostream &ost, INT size, INT *set);
void block_swap_bytes(SCHAR *ptr, INT size, INT no);
void code_INT4(char *&p, INT4 i);
INT4 decode_INT4(char *&p);
void code_UBYTE(char *&p, UBYTE a);
void decode_UBYTE(char *&p, UBYTE &a);
void print_incidence_structure(ostream &ost, INT m, INT n, INT len, INT *S);
void scan_permutation_from_string(const char *s, INT *&perm, INT &degree, INT verbose_level);
void scan_permutation_from_stream(istream & is, INT *&perm, INT &degree, INT verbose_level);
char get_character(istream & is, INT verbose_level);
void replace_extension_with(char *p, const char *new_ext);
void chop_off_extension_if_present(char *p, const char *ext);
void get_fname_base(const char *p, BYTE *fname_base);
void get_extension_if_present(const char *p, char *ext);
INT s_scan_int(BYTE **s, INT *i);
INT s_scan_token(BYTE **s, BYTE *str);
INT s_scan_token_arbitrary(BYTE **s, BYTE *str);
INT s_scan_str(BYTE **s, BYTE *str);
INT s_scan_token_comma_separated(BYTE **s, BYTE *str);
INT hashing(INT hash0, INT a);
INT hashing_fixed_width(INT hash0, INT a, INT bit_length);
INT INT_vec_hash(INT *v, INT len, INT bit_length);
void parse_sets(INT nb_cases, BYTE **data, INT f_casenumbers, 
	INT *&Set_sizes, INT **&Sets, BYTE **&Ago_ascii, BYTE **&Aut_ascii, 
	INT *&Casenumbers, 
	INT verbose_level);
void parse_sets_and_check_sizes_easy(INT len, INT nb_cases, 
	BYTE **data, INT **&sets);
void parse_line(BYTE *line, INT &len, INT *&set, BYTE *ago_ascii, BYTE *aut_ascii);
INT count_number_of_orbits_in_file(const BYTE *fname, INT verbose_level);
INT count_number_of_lines_in_file(const BYTE *fname, INT verbose_level);
INT try_to_read_file(const BYTE *fname, INT &nb_cases, BYTE **&data, INT verbose_level);
void read_and_parse_data_file(const BYTE *fname, INT &nb_cases, BYTE **&data, INT **&sets, INT *&set_sizes, INT verbose_level);
void free_data_fancy(INT nb_cases, 
	INT *Set_sizes, INT **Sets, 
	BYTE **Ago_ascii, BYTE **Aut_ascii, 
	INT *Casenumbers);
void read_and_parse_data_file_fancy(const BYTE *fname, 
	INT f_casenumbers, 
	INT &nb_cases, 
	INT *&Set_sizes, INT **&Sets, BYTE **&Ago_ascii, BYTE **&Aut_ascii, 
	INT *&Casenumbers, 
	INT verbose_level);
void read_set_from_file(const BYTE *fname, INT *&the_set, INT &set_size, INT verbose_level);
void write_set_to_file(const BYTE *fname, INT *the_set, INT set_size, INT verbose_level);
void read_set_from_file_INT4(const BYTE *fname, INT *&the_set, INT &set_size, INT verbose_level);
void write_set_to_file_as_INT4(const BYTE *fname, INT *the_set, INT set_size, INT verbose_level);
void read_k_th_set_from_file(const BYTE *fname, INT k, INT *&the_set, INT &set_size, INT verbose_level);
void write_incidence_matrix_to_file(BYTE *fname, INT *Inc, INT m, INT n, INT verbose_level);
void read_incidence_matrix_from_inc_file(INT *&M, INT &m, INT &n, 
	BYTE *inc_file_name, INT inc_file_idx, INT verbose_level);
INT inc_file_get_number_of_geometries(
	BYTE *inc_file_name, INT verbose_level);

void print_line_of_number_signs();
void print_repeated_character(ostream &ost, BYTE c, INT n);
void print_pointer_hex(ostream &ost, void *p);
void print_hex_digit(ostream &ost, INT digit);
void count_number_of_solutions_in_file_by_case(const BYTE *fname, 
	INT *&nb_solutions, INT *&case_nb, INT &nb_cases, 
	INT verbose_level);
void read_solutions_from_file_by_case(const BYTE *fname, 
	INT *nb_solutions, INT *case_nb, INT nb_cases, 
	INT **&Solutions, INT solution_size, 
	INT verbose_level);
void copy_file_to_ostream(ostream &ost, BYTE *fname);
void INT_vec_write_csv(INT *v, INT len, const BYTE *fname, const BYTE *label);
void INT_vecs_write_csv(INT *v1, INT *v2, INT len, const BYTE *fname, const BYTE *label1, const BYTE *label2);
void INT_vec_array_write_csv(INT nb_vecs, INT **Vec, INT len, const BYTE *fname, const BYTE **column_label);
void INT_matrix_write_csv(const BYTE *fname, INT *M, INT m, INT n);
void INT_matrix_write_csv_with_labels(const BYTE *fname, INT *M, INT m, INT n, const BYTE **column_label);
void INT_matrix_read_csv(const BYTE *fname, INT *&M, INT &m, INT &n, INT verbose_level);
void INT_matrix_write_text(const BYTE *fname, INT *M, INT m, INT n);
void INT_matrix_read_text(const BYTE *fname, INT *&M, INT &m, INT &n);
INT compare_sets(INT *set1, INT *set2, INT sz1, INT sz2);
INT test_if_sets_are_disjoint(INT *set1, INT *set2, INT sz1, INT sz2);
void make_graph_of_disjoint_sets_from_rows_of_matrix(INT *M, INT m, INT n, INT *&Adj, INT verbose_level);
void write_exact_cover_problem_to_file(INT *Inc, INT nb_rows, INT nb_cols, const BYTE *fname);
void read_solution_file(BYTE *fname, 
	INT *Inc, INT nb_rows, INT nb_cols, 
	INT *&Solutions, INT &sol_length, INT &nb_sol, 
	INT verbose_level);
// sol_length must be constant

// ####################################################################################
// memory.C:
// ####################################################################################

extern int f_memory_debug;
extern int f_memory_debug_verbose;
extern INT memory_count_allocate;
void start_memory_debug();
void stop_memory_debug();
void memory_watch_list_add_pointer(void *p);
void memory_watch_list_delete_pointer(INT idx);
void add_to_registry(void *p, int pointer_type, int size, int size_of, const char *file, int line);
INT delete_from_registry(void *p);
void memory_watch_list_dump();
void registry_dump();
void registry_dump_sorted();
void registry_dump_sorted_by_size();
INT registry_entry_size(INT i);
void registry_print_entry(INT i);
INT registry_sizeof(INT t);
void registry_print_type(INT t);
int registry_search(int len, void *p, int &idx);
int memory_watch_list_search(int len, void *p, int &idx);
int *allocate_int(int n, const char *file, int line);
void free_int(int *p, const char *file, int line);
pint *allocate_pint(int n, const char *file, int line);
void free_pint(pint *p, const char *file, int line);
INT *allocate_INT(int n, const char *file, int line);
void free_INT(INT *p, const char *file, int line);
INT **allocate_PINT(int n, const char *file, int line);
void free_PINT(INT **p, const char *file, int line);
INT ***allocate_PPINT(int n, const char *file, int line);
void free_PPINT(INT ***p, const char *file, int line);
BYTE *allocate_BYTE(int n, const char *file, int line);
void free_BYTE(BYTE *p, const char *file, int line);
UBYTE *allocate_UBYTE(int n, const char *file, int line);
void free_UBYTE(UBYTE *p, const char *file, int line);
BYTE **allocate_PBYTE(int n, const char *file, int line);
void free_PBYTE(BYTE **p, const char *file, int line);
void **allocate_pvoid(int n, const char *file, int line);
void free_pvoid(void **p, const char *file, int line);
void *allocate_OBJECT(void *p, int size_of, const char *file, int line);
void free_OBJECT(void *p, const char *file, int line);
void *allocate_OBJECTS(void *p, int n, int size_of, const char *file, int line);
void free_OBJECTS(void *p, const char *file, int line);


// ####################################################################################
// sorting.C:
// ####################################################################################

INT INT_vec_is_subset_of(INT *set, INT sz, INT *big_set, INT big_set_sz);
void INT_vec_swap_points(INT *list, INT *list_inv, INT idx1, INT idx2);
INT INT_vec_is_sorted(INT *v, INT len);
void INT_vec_sort_and_remove_duplicates(INT *v, INT &len);
INT INT_vec_sort_and_test_if_contained(INT *v1, INT len1, INT *v2, INT len2);
void INT_vec_insert_and_reallocate_if_necessary(INT *&vec, INT &used_length, INT &alloc_length, INT a, INT verbose_level);
void INT_vec_append_and_reallocate_if_necessary(INT *&vec, INT &used_length, INT &alloc_length, INT a, INT verbose_level);
INT INT_vec_is_zero(INT *v, INT len);
void test_if_set(INT *set, INT set_size);
INT test_if_set_with_return_value(INT *set, INT set_size);
void rearrange_subset(INT n, INT k, INT *set, 
	INT *subset, INT *rearranged_set, INT verbose_level);
INT INT_vec_search_linear(INT *v, INT len, INT a, INT &idx);
void INT_vec_intersect(INT *v1, INT len1, INT *v2, INT len2, INT *&v3, INT &len3);
void INT_vec_intersect_sorted_vectors(INT *v1, INT len1, INT *v2, INT len2, INT *v3, INT &len3);
void INT_vec_sorting_permutation(INT *v, INT len, INT *perm, INT *perm_inv, INT f_increasingly);
// perm and perm_inv must be allocated to len elements
INT INT_compare_increasingly(void *a, void *b, void *data);
INT INT_compare_decreasingly(void *a, void *b, void *data);
void INT_vec_quicksort(INT *v, INT (*compare_func)(INT a, INT b), INT left, INT right);
INT compare_increasingly_INT(INT a, INT b);
INT compare_decreasingly_INT(INT a, INT b);
void INT_vec_quicksort_increasingly(INT *v, INT len);
void INT_vec_quicksort_decreasingly(INT *v, INT len);
void quicksort_array(INT len, void **v, 
	INT (*compare_func)(void *a, void *b, void *data), void *data);
void quicksort_array_with_perm(INT len, void **v, INT *perm, 
	INT (*compare_func)(void *a, void *b, void *data), void *data);
void INT_vec_sort(INT len, INT *p);
int int_vec_compare(int *p, int *q, int len);
INT INT_vec_compare(INT *p, INT *q, INT len);
INT INT_vec_compare_stride(INT *p, INT *q, INT len, INT stride);
INT vec_search(void **v, INT (*compare_func)(void *a, void *b, void *data), void *data_for_compare, 
	INT len, void *a, INT &idx, INT verbose_level);
INT vec_search_general(void *vec, 
	INT (*compare_func)(void *vec, void *a, INT b, void *data_for_compare), void *data_for_compare, 
	INT len, void *a, INT &idx, INT verbose_level);
INT INT_vec_search(INT *v, INT len, INT a, INT &idx);
// This function finds the last occurence of the element a.
// If a is not found, it returns in idx the position where it should be inserted if 
// the vector is assumed to be in increasing order.
INT INT_vec_search_first_occurence(INT *v, INT len, INT a, INT &idx);
// This function finds the first occurence of the element a.
INT longinteger_vec_search(longinteger_object *v, INT len, 
	longinteger_object &a, INT &idx);
void INT_vec_classify(INT length, INT *the_vec, INT *&the_vec_sorted, 
	INT *&sorting_perm, INT *&sorting_perm_inv, 
	INT &nb_types, INT *&type_first, INT *&type_len);
void INT_vec_classify_with_arrays(INT length, INT *the_vec, INT *the_vec_sorted, 
	INT *sorting_perm, INT *sorting_perm_inv, 
	INT &nb_types, INT *type_first, INT *type_len);
void INT_vec_sorted_collect_types(INT length, INT *the_vec_sorted, 
	INT &nb_types, INT *type_first, INT *type_len);
void INT_vec_print_classified(ostream &ost, INT *vec, INT len);
void INT_vec_print_types(ostream &ost, INT f_backwards, INT *the_vec_sorted, 
	INT nb_types, INT *type_first, INT *type_len);
void INT_vec_print_types_naked(ostream &ost, INT f_backwards, INT *the_vec_sorted, 
	INT nb_types, INT *type_first, INT *type_len);
void INT_vec_print_types_naked_tex(ostream &ost, INT f_backwards, INT *the_vec_sorted, 
	INT nb_types, INT *type_first, INT *type_len);
void Heapsort(void *v, INT len, INT entry_size_in_bytes, 
	INT (*compare_func)(void *v1, void *v2));
void Heapsort_general(void *data, INT len, 
	INT (*compare_func)(void *data, INT i, INT j), 
	void (*swap_func)(void *data, INT i, INT j));
void INT_vec_heapsort(INT *v, INT len);
void INT_vec_heapsort_with_log(INT *v, INT *w, INT len);
void heapsort_make_heap(INT *v, INT len);
void heapsort_make_heap_with_log(INT *v, INT *w, INT len);
void Heapsort_make_heap(void *v, INT len, INT entry_size_in_bytes, 
	INT (*compare_func)(void *v1, void *v2));
void Heapsort_general_make_heap(void *data, INT len, 
	INT (*compare_func)(void *data, INT i, INT j), 
	void (*swap_func)(void *data, INT i, INT j));
void heapsort_sift_down(INT *v, INT start, INT end);
void heapsort_sift_down_with_log(INT *v, INT *w, INT start, INT end);
void Heapsort_sift_down(void *v, INT start, INT end, INT entry_size_in_bytes, 
	INT (*compare_func)(void *v1, void *v2));
void Heapsort_general_sift_down(void *data, INT start, INT end, 
	INT (*compare_func)(void *data, INT i, INT j), 
	void (*swap_func)(void *data, INT i, INT j));
void heapsort_swap(INT *v, INT i, INT j);
void Heapsort_swap(void *v, INT i, INT j, INT entry_size_in_bytes);
INT is_all_digits(BYTE *p);


// ####################################################################################
// data_file.C:
// ####################################################################################

class data_file {
	
	public:

	BYTE fname[1000];
	INT nb_cases;
	INT *set_sizes;
	INT **sets;
	INT *casenumbers;
	BYTE **Ago_ascii;
	BYTE **Aut_ascii;

	INT f_has_candidates;
	INT *nb_candidates;
	INT **candidates;

	data_file();
	~data_file();
	void null();
	void freeself();
	void read(const BYTE *fname, INT f_casenumbers, INT verbose_level);
	void read_candidates(const BYTE *candidates_fname, INT verbose_level);
};


// ##################################################################################################
// clique_finder.C
// ##################################################################################################


class clique_finder {
public:
	BYTE label[1000];
	INT n; // number of points
	
	INT print_interval;
	
	INT f_write_tree;
	INT f_decision_nodes_only;
	BYTE fname_tree[1000];
	ofstream *fp_tree;

	
	INT f_maxdepth;
	INT maxdepth;
	
	INT *point_labels;
	INT *point_is_suspicous;
	
	INT target_depth;
	INT verbose_level;
	
	//INT *adjacency;
	//bitmatrix *adjacency;
	INT bitmatrix_m;
	INT bitmatrix_n;
	INT bitmatrix_N;
	UBYTE *bitmatrix_adjacency;

	INT *pt_list;
	INT *pt_list_inv;
	INT *nb_points;
	INT *candidates; // [max_depth * n]
	INT *nb_candidates; // [max_depth]
	INT *current_choice; // [max_depth]
	INT *level_counter; // [max_depth] (added Nov 8, 2014)
	INT *f_level_mod; // [max_depth] (added Nov 8, 2014)
	INT *level_r; // [max_depth] (added Nov 8, 2014)
	INT *level_m; // [max_depth] (added Nov 8, 2014)

	INT *current_clique; // [max_depth]

	UINT counter; // number of backtrack nodes
	UINT decision_step_counter; // number of backtrack nodes that are decision nodes

	// solution storage:
	INT f_store_solutions;
	deque<vector<int> > solutions;
	INT nb_sol;


	// callbacks:
	void (*call_back_clique_found)(clique_finder *CF, INT verbose_level);
	
	// added May 26, 2009:
	void (*call_back_add_point)(clique_finder *CF, 
		INT current_clique_size, INT *current_clique, 
		INT pt, INT verbose_level);
	void (*call_back_delete_point)(clique_finder *CF, 
		INT current_clique_size, INT *current_clique, 
		INT pt, INT verbose_level);
	INT (*call_back_find_candidates)(clique_finder *CF, 
		INT current_clique_size, INT *current_clique, 
		INT nb_pts, INT &reduced_nb_pts, 
		INT *pt_list, INT *pt_list_inv, 
		INT *candidates, INT verbose_level);
		// Jan 2, 2012: added reduced_nb_pts and pt_list_inv

	INT (*call_back_is_adjacent)(clique_finder *CF, 
		INT pt1, INT pt2, INT verbose_level);
	// added Oct 2011:
	void (*call_back_after_reduction)(clique_finder *CF, 
		INT depth, INT nb_points, INT verbose_level);

	// added Nov 2014:
	INT f_has_print_current_choice_function;
	void (*call_back_print_current_choice)(clique_finder *CF, 
		INT depth, void *user_data, INT verbose_level);
	void *print_current_choice_data;
	
	void *call_back_clique_found_data;
	
	
	void open_tree_file(const BYTE *fname_base, INT f_decision_nodes_only);
	void close_tree_file();
	void init(const BYTE *label, INT n, 
		INT target_depth, 
		INT f_has_adj_list, INT *adj_list_coded, 
		INT f_has_bitvector, UBYTE *bitvector_adjacency, 
		INT print_interval, 
		INT f_maxdepth, INT maxdepth, 
		INT f_store_solutions, 
		INT verbose_level);
	void init_restrictions(INT *restrictions, INT verbose_level);
	clique_finder();
	~clique_finder();
	void null();
	void free();
	void init_point_labels(INT *pt_labels);
	void init_suspicous_points(INT nb, INT *point_list);
	void print_suspicous_points();
	void print_set(INT size, INT *set);
	void print_suspicous_point_subset(INT size, INT *set);
	void log_position_and_choice(INT depth, INT counter_save, INT counter);
	void log_position(INT depth, INT counter_save, INT counter);
	void log_choice(INT depth);
	void swap_point(INT idx1, INT idx2);
	void degree_of_point_statistic(INT depth, INT nb_points, INT verbose_level);
	INT degree_of_point(INT depth, INT i, INT nb_points);
	//INT degree_of_point_verbose(INT i, INT nb_points);
	INT is_suspicous(INT i);
	INT point_label(INT i);
	INT is_adjacent(INT depth, INT i, INT j);
	INT is_viable(INT depth, INT pt);
	void write_entry_to_tree_file(INT depth, INT verbose_level);
	void m_iji(INT i, INT j, INT a);
	INT s_ij(INT i, INT j);
	void backtrack_search(INT depth, INT verbose_level);
	INT solve_decision_problem(INT depth, INT verbose_level);
		// returns TRUE if we found a solution
	void get_solutions(INT *&Sol, INT &nb_sol, INT &clique_sz, INT verbose_level);
};

void all_cliques_of_given_size(INT *Adj, INT nb_pts, INT clique_sz, INT *&Sol, INT &nb_sol, INT verbose_level);



// ##################################################################################################
// colored_graph.C
// ##################################################################################################


class colored_graph {
public:

	BYTE fname_base[1000];
	
	INT nb_points;
	INT nb_colors;
	
	INT bitvector_length;
	INT L;
	
	INT *points;
	INT *point_color;
	

	INT user_data_size;
	INT *user_data;

	INT f_ownership_of_bitvec;
	UBYTE *bitvector_adjacency;

	INT f_has_list_of_edges;
	INT nb_edges;
	INT *list_of_edges;

	colored_graph();
	~colored_graph();
	void null();
	void freeself();
	void compute_edges(INT verbose_level);
	INT is_adjacent(INT i, INT j);
	void set_adjacency(INT i, INT j, INT a);
	void print();
	void init(INT nb_points, INT nb_colors, 
		INT *colors, UBYTE *bitvec, INT f_ownership_of_bitvec, 
		INT verbose_level);
	void init_no_colors(INT nb_points, UBYTE *bitvec, INT f_ownership_of_bitvec, 
		INT verbose_level);
	void init_adjacency(INT nb_points, INT nb_colors, 
		INT *colors, INT *Adj, INT verbose_level);
	void init_adjacency_no_colors(INT nb_points, INT *Adj, INT verbose_level);
	void init_user_data(INT *data, INT data_size, INT verbose_level);
	void save(const BYTE *fname, INT verbose_level);
	void load(const BYTE *fname, INT verbose_level);
	void all_cliques_of_size_k_ignore_colors(INT target_depth, 
		INT &nb_sol, INT &decision_step_counter, INT verbose_level);
	void all_cliques_of_size_k_ignore_colors_and_write_solutions_to_file(INT target_depth, 
		const BYTE *fname, 
		INT &nb_sol, INT &decision_step_counter, 
		INT verbose_level);
	void all_rainbow_cliques(ofstream *fp, INT f_output_solution_raw, 
		INT f_maxdepth, INT maxdepth, 
		INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
		INT print_interval, 
		INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
		INT verbose_level);
	void all_rainbow_cliques_with_additional_test_function(ofstream *fp, INT f_output_solution_raw, 
		INT f_maxdepth, INT maxdepth, 
		INT f_restrictions, INT *restrictions, 
		INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
		INT print_interval, 
		INT f_has_additional_test_function,
		void (*call_back_additional_test_function)(rainbow_cliques *R, void *user_data, 
			INT current_clique_size, INT *current_clique, 
			INT nb_pts, INT &reduced_nb_pts, 
			INT *pt_list, INT *pt_list_inv, 
			INT verbose_level), 
		INT f_has_print_current_choice_function,
		void (*call_back_print_current_choice)(clique_finder *CF, 
			INT depth, void *user_data, INT verbose_level), 
		void *user_data, 
		INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
		INT verbose_level);
	void draw_on_circle(char *fname, 
		INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out,
		INT f_labels, INT f_embedded, INT f_sideways, 
		double tikz_global_scale, double tikz_global_line_width);
	void draw_on_circle_2(mp_graphics &G, INT f_labels);
	void draw(const BYTE *fname, 
		INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out,
		INT verbose_level);
	void draw_partitioned(const BYTE *fname, 
		INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out,
		INT verbose_level);
	colored_graph *compute_neighborhood_subgraph(INT pt, 
		fancy_set *&vertex_subset, fancy_set *&color_subset, INT verbose_level);
	colored_graph *compute_neighborhood_subgraph_with_additional_test_function(INT pt, 
		fancy_set *&vertex_subset, fancy_set *&color_subset, 
		INT (*test_function)(colored_graph *CG, INT test_point, INT pt, void *test_function_data, INT verbose_level),
		void *test_function_data, 
		INT verbose_level);
	void export_to_magma(const BYTE *fname, INT verbose_level);
	void export_to_file(const BYTE *fname, INT verbose_level);
	void export_to_file_matlab(const BYTE *fname, INT verbose_level);
	void early_test_func_for_clique_search(INT *S, INT len, 
		INT *candidates, INT nb_candidates, 
		INT *good_candidates, INT &nb_good_candidates, 
		INT verbose_level);
	void early_test_func_for_coclique_search(INT *S, INT len, 
		INT *candidates, INT nb_candidates, 
		INT *good_candidates, INT &nb_good_candidates, 
		INT verbose_level);
	void early_test_func_for_path_and_cycle_search(INT *S, INT len, 
		INT *candidates, INT nb_candidates, 
		INT *good_candidates, INT &nb_good_candidates, 
		INT verbose_level);
	INT is_cycle(INT nb_e, INT *edges, INT verbose_level);
	void draw_it(const BYTE *fname_base, INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out);

};

void colored_graph_all_cliques(const BYTE *fname, INT f_output_solution_raw, 
	INT f_draw, INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out, 
	INT f_output_fname, const BYTE *output_fname, 
	INT f_maxdepth, INT maxdepth, 
	INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
	INT print_interval, 
	INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
	INT verbose_level);
void colored_graph_all_cliques_list_of_cases(INT *list_of_cases, INT nb_cases, INT f_output_solution_raw, 
	INT f_draw, INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out, 
	const BYTE *fname_template, 
	const BYTE *fname_sol, const BYTE *fname_stats, 
	INT f_maxdepth, INT maxdepth, 
	INT f_prefix, const BYTE *prefix, 
	INT print_interval, INT verbose_level);
void call_back_clique_found_using_file_output(clique_finder *CF, INT verbose_level);

// ##################################################################################################
// rainbow_cliques.C
// ##################################################################################################


class rainbow_cliques {
public:

	rainbow_cliques();
	~rainbow_cliques();
	void null();
	void freeself();

	ofstream *fp_sol;
	INT f_output_solution_raw;
	
	colored_graph *graph;
	clique_finder *CF;
	INT *f_color_satisfied;
	INT *color_chosen_at_depth;
	INT *color_frequency;
	INT target_depth;

	// added November 5, 2014:
	INT f_has_additional_test_function;
	void (*call_back_additional_test_function)(rainbow_cliques *R, void *user_data, 
		INT current_clique_size, INT *current_clique, 
		INT nb_pts, INT &reduced_nb_pts, 
		INT *pt_list, INT *pt_list_inv, 
		INT verbose_level);
	void *user_data;


	void search(colored_graph *graph, ofstream *fp_sol, INT f_output_solution_raw, 
		INT f_maxdepth, INT maxdepth, 
		INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
		INT print_interval, 
		INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
		INT verbose_level);
	void search_with_additional_test_function(colored_graph *graph, ofstream *fp_sol, INT f_output_solution_raw, 
		INT f_maxdepth, INT maxdepth, 
		INT f_restrictions, INT *restrictions, 
		INT f_tree, INT f_decision_nodes_only, const BYTE *fname_tree,  
		INT print_interval, 
		INT f_has_additional_test_function,
		void (*call_back_additional_test_function)(rainbow_cliques *R, void *user_data, 
			INT current_clique_size, INT *current_clique, 
			INT nb_pts, INT &reduced_nb_pts, 
			INT *pt_list, INT *pt_list_inv, 
			INT verbose_level), 
		INT f_has_print_current_choice_function,
		void (*call_back_print_current_choice)(clique_finder *CF, 
			INT depth, void *user_data, INT verbose_level), 
		void *user_data, 
		INT &search_steps, INT &decision_steps, INT &nb_sol, INT &dt, 
		INT verbose_level);
	INT find_candidates(
		INT current_clique_size, INT *current_clique, 
		INT nb_pts, INT &reduced_nb_pts, 
		INT *pt_list, INT *pt_list_inv, 
		INT *candidates, INT verbose_level);
	void clique_found(INT *current_clique, INT verbose_level);
	void clique_found_record_in_original_labels(INT *current_clique, INT verbose_level);

};

void call_back_colored_graph_clique_found(clique_finder *CF, INT verbose_level);
void call_back_colored_graph_add_point(clique_finder *CF, 
	INT current_clique_size, INT *current_clique, 
	INT pt, INT verbose_level);
void call_back_colored_graph_delete_point(clique_finder *CF, 
	INT current_clique_size, INT *current_clique, 
	INT pt, INT verbose_level);
INT call_back_colored_graph_find_candidates(clique_finder *CF, 
	INT current_clique_size, INT *current_clique, 
	INT nb_pts, INT &reduced_nb_pts, 
	INT *pt_list, INT *pt_list_inv, 
	INT *candidates, INT verbose_level);


// ####################################################################################
// super_fast_hash.C:
// ####################################################################################

uint32_t SuperFastHash (const char * data, int len);


// ####################################################################################
// plot.C:
// ####################################################################################

void draw_density(BYTE *prefix, INT *the_set, INT set_size,
	INT f_title, const BYTE *title, INT out_of, 
	const BYTE *label_x, 
	INT f_circle, INT circle_at, INT circle_rad, 
	INT f_mu, INT f_sigma, INT nb_standard_deviations, 
	INT f_v_grid, INT v_grid, INT f_h_grid, INT h_grid, 
	INT xmax, INT ymax, INT offset_x, INT f_switch_x, INT no, INT f_embedded, 
	INT verbose_level);
void draw_density_multiple_curves(BYTE *prefix,
	INT **Data, INT *Data_size, INT nb_data_sets, 
	INT f_title, const BYTE *title, INT out_of, 
	const BYTE *label_x, 
	INT f_v_grid, INT v_grid, INT f_h_grid, INT h_grid, 
	INT xmax, INT ymax, INT offset_x, INT f_switch_x, 
	INT f_v_logarithmic, double log_base, INT no, INT f_embedded, 
	INT verbose_level);
void draw_density2(mp_graphics &G, INT no, 
	INT *outline_value, INT *outline_number, INT outline_sz, 
	INT min_value, INT max_value, INT offset_x, INT f_switch_x, 
	INT f_title, const BYTE *title, 
	const BYTE *label_x, 
	INT f_circle, INT circle_at, INT circle_rad, 
	INT f_mu, INT f_sigma, INT nb_standard_deviations, 
	INT f_v_grid, INT v_grid, INT f_h_grid, INT h_grid);
void draw_density2_multiple_curves(mp_graphics &G, INT no, 
	INT **outline_value, INT **outline_number, INT *outline_sz, INT nb_curves, 
	INT min_x, INT max_x, INT min_y, INT max_y, 
	INT offset_x, INT f_switch_x, 
	INT f_title, const BYTE *title, 
	const BYTE *label_x, 
	INT f_v_grid, INT v_grid, INT f_h_grid, INT h_grid, 
	INT f_v_logarithmic, double log_base);
void read_numbers_from_file(const BYTE *fname, 
	INT *&the_set, INT &set_size, INT verbose_level);
void get_coord(INT *Px, INT *Py, INT idx, INT x, INT y, INT min_x, INT min_y, INT max_x, INT max_y, INT f_switch_x);
void get_coord_log(INT *Px, INT *Py, INT idx, INT x, INT y, INT min_x, INT min_y, INT max_x, INT max_y, double log_base, INT f_switch_x);
void y_to_pt_on_curve(INT y_in, INT &x, INT &y,  
	INT *outline_value, INT *outline_number, INT outline_sz);



// ####################################################################################
// set_of_sets.C:
// ####################################################################################

class set_of_sets {

public:
	
	INT underlying_set_size;
	INT nb_sets;
	INT **Sets;
	INT *Set_size;


	set_of_sets();
	~set_of_sets();
	void null();
	void freeself();
	void init_simple(INT underlying_set_size, INT nb_sets, INT verbose_level);
	void init(INT underlying_set_size, INT nb_sets, INT **Pts, INT *Sz, INT verbose_level);
	void init_basic(INT underlying_set_size, INT nb_sets, INT *Sz, INT verbose_level);
	void init_basic_constant_size(INT underlying_set_size, 
		INT nb_sets, INT constant_size, INT verbose_level);
	void init_set(INT idx_of_set, INT *set, INT sz, INT verbose_level);
		// Stores a copy of the given set.
	void print();
	void print_table();
	void remove_sets_of_given_size(INT k, set_of_sets &S, INT *&Idx, INT verbose_level);
	void extract_largest_sets(set_of_sets &S, INT *&Idx, INT verbose_level);
	void intersection_matrix(
		INT *&intersection_type, INT &highest_intersection_number, 
		INT *&intersection_matrix, INT &nb_big_sets, 
		INT verbose_level);
	void compute_incidence_matrix(INT *&Inc, INT &m, INT &n, INT verbose_level);
	void compute_and_print_tdo_row_scheme(ofstream &file, INT verbose_level);
	void compute_and_print_tdo_col_scheme(ofstream &file, INT verbose_level);
	void init_decomposition(decomposition *&D, INT verbose_level);
	void compute_tdo_decomposition(decomposition &D, INT verbose_level);
	INT is_member(INT i, INT a, INT verbose_level);
	void sort_all(INT verbose_level);
	void all_pairwise_intersections(set_of_sets *&Intersections, INT verbose_level);
	void pairwise_intersection_matrix(INT *&M, INT verbose_level);
	void all_triple_intersections(set_of_sets *&Intersections, INT verbose_level);
	INT has_constant_size_property();
	INT largest_set_size();
	void save_csv(const BYTE *fname, INT f_make_heading, INT verbose_level);
	void sort_big(INT verbose_level);
	void compute_orbits(INT &nb_orbits, INT *&orbit, INT *&orbit_inv, 
		INT *&orbit_first, INT *&orbit_len, 
		void (*compute_image_function)(set_of_sets *S, void *compute_image_data, INT elt_idx, INT gen_idx, INT &idx_of_image, INT verbose_level), 
		void *compute_image_data, 
		INT nb_gens, 
		INT verbose_level);
};

INT set_of_sets_compare_func(void *data, INT i, INT j);
void set_of_sets_swap_func(void *data, INT i, INT j);

// ####################################################################################
// decomposition.C:
// ####################################################################################

class decomposition {

public:
	
	INT nb_points;
	INT nb_blocks;
	INT *Inc;
	incidence_structure *I;
	partitionstack *Stack;

	INT f_has_decomposition;
	INT *row_classes;
	INT *row_class_inv;
	INT nb_row_classes;
	INT *col_classes;
	INT *col_class_inv;
	INT nb_col_classes;
	INT f_has_row_scheme;
	INT *row_scheme;
	INT f_has_col_scheme;
	INT *col_scheme;
	


	decomposition();
	~decomposition();
	void null();
	void freeself();
	void init_inc_and_stack(incidence_structure *Inc, partitionstack *Stack, INT verbose_level);
	void init_incidence_matrix(INT m, INT n, INT *M, INT verbose_level);
		// copies the incidence matrix
	void setup_default_partition(INT verbose_level);
	void compute_TDO(INT max_depth, INT verbose_level);
	void print_row_decomposition_tex(ostream &ost, 
		INT f_enter_math, INT f_print_subscripts, INT verbose_level);
	void print_column_decomposition_tex(ostream &ost, 
		INT f_enter_math, INT f_print_subscripts, INT verbose_level);
	void get_row_scheme(INT verbose_level);
	void get_col_scheme(INT verbose_level);
	
};

// ####################################################################################
// brick_domain.C:
// ####################################################################################

class brick_domain {

public:
	finite_field *F;
	INT q;
	INT nb_bricks;

	brick_domain();
	~brick_domain();
	void null();
	void freeself();
	void init(finite_field *F, INT verbose_level);
	void unrank(INT rk, INT &f_vertical, INT &x0, INT &y0, INT verbose_level);
	INT rank(INT f_vertical, INT x0, INT y0, INT verbose_level);
	void unrank_coordinates(INT rk, INT &x1, INT &y1, INT &x2, INT &y2, INT verbose_level);
	INT rank_coordinates(INT x1, INT y1, INT x2, INT y2, INT verbose_level);
};

void brick_test(INT q, INT verbose_level);


// ####################################################################################
// spreadsheet.C:
// ####################################################################################

class spreadsheet {

public:

	BYTE **tokens;
	INT nb_tokens;

	INT *line_start, *line_size;
	INT nb_lines;

	INT nb_rows, nb_cols;
	INT *Table;


	spreadsheet();
	~spreadsheet();
	void null();
	void freeself();
	void init_set_of_sets(set_of_sets *S, INT f_make_heading);
	void init_INT_matrix(INT nb_rows, INT nb_cols, INT *A);
	void add_token(BYTE *label);
	void save(const BYTE *fname, INT verbose_level);
	void read_spreadsheet(const BYTE *fname, INT verbose_level);
	void print_table(ostream &ost, INT f_enclose_in_parentheses);
	void print_table_row(INT row, INT f_enclose_in_parentheses, ostream &ost);
	void print_table_row_detailed(INT row, ostream &ost);
	void print_table_with_row_selection(INT *f_selected, ostream &ost);
	void print_table_sorted(ostream &ost, const BYTE *sort_by);
	void add_column_with_constant_value(BYTE *label, BYTE *value);
	void reallocate_table();
	void reallocate_table_add_row();
	INT find_by_column(const BYTE *join_by);
	void tokenize(const BYTE *fname, 
		BYTE **&tokens, INT &nb_tokens, INT verbose_level);
	void remove_quotes(INT verbose_level);
	void remove_rows(const BYTE *drop_column, const BYTE *drop_label, INT verbose_level);
	void remove_rows_where_field_is_empty(const BYTE *drop_column, INT verbose_level);
	void find_rows(INT verbose_level);
	BYTE *get_string(INT i, INT j);
	INT get_INT(INT i, INT j);
	void join_with(spreadsheet *S2, INT by1, INT by2, INT verbose_level);
	void patch_with(spreadsheet *S2, BYTE *join_by);


};

INT my_atoi(BYTE *str);
INT compare_strings(void *a, void *b, void *data);

// ####################################################################################
// dlx.C:
// ####################################################################################

void install_callback_solution_found(
	void (*callback_solution_found)(INT *solution, INT len, INT nb_sol, void *data),
	void *callback_solution_found_data);
void de_install_callback_solution_found();
void DlxTest();
void DlxTransposeAppendAndSolve(INT *Data, INT nb_rows, INT nb_cols, 
	INT &nb_sol, INT &nb_backtrack, 
	INT f_write_file, const BYTE *solution_fname, 
	INT f_write_tree_file, const BYTE *tree_fname, 
	INT verbose_level);
void DlxTransposeAndSolveRHS(INT *Data, INT nb_rows, INT nb_cols, 
	INT *RHS, INT f_has_type, diophant_equation_type *type, 
	INT &nb_sol, INT &nb_backtrack, 
	INT f_write_file, const BYTE *solution_fname, 
	INT f_write_tree_file, const BYTE *tree_fname, 
	INT verbose_level);
void DlxAppendRowAndSolve(INT *Data, INT nb_rows, INT nb_cols, 
	INT &nb_sol, INT &nb_backtrack, 
	INT f_write_file, const BYTE *solution_fname, 
	INT f_write_tree_file, const BYTE *tree_fname, 
	INT verbose_level);
void DlxAppendRowAndSolveRHS(INT *Data, INT nb_rows, INT nb_cols, 
	INT *RHS, INT f_has_type, diophant_equation_type *type, 
	INT &nb_sol, INT &nb_backtrack, 
	INT f_write_file, const BYTE *solution_fname, 
	INT f_write_tree_file, const BYTE *tree_fname, 
	INT verbose_level);
void DlxSolve(INT *Data, INT nb_rows, INT nb_cols, 
	INT &nb_sol, INT &nb_backtrack, 
	INT f_write_file, const BYTE *solution_fname, 
	INT f_write_tree_file, const BYTE *tree_fname, 
	INT verbose_level);
void DlxSolve_with_RHS(INT *Data, INT nb_rows, INT nb_cols, 
	INT *RHS, INT f_has_type, diophant_equation_type *type, 
	INT &nb_sol, INT &nb_backtrack, 
	INT f_write_file, const BYTE *solution_fname, 
	INT f_write_tree_file, const BYTE *tree_fname, 
	INT verbose_level);
void DlxSearchRHS(INT k, INT verbose_level);



// ####################################################################################
// andre_construction.C:
// ####################################################################################

class andre_construction {
public:
	INT order; // = q^k
	INT spread_size; // order + 1
	INT n; // = 2 * k
	INT k;
	INT q;
	INT N; // order^2 + order + 1

	
	grassmann *Grass;
	finite_field *F;

	INT *spread_elements_numeric; // [spread_size]
	INT *spread_elements_numeric_sorted; // [spread_size]

	INT *spread_elements_perm;
	INT *spread_elements_perm_inv;

	INT *spread_elements_genma; // [spread_size * k * n]
	INT *pivot; //[spread_size * k]
	INT *non_pivot; //[spread_size * (n - k)]
	

	andre_construction();
	~andre_construction();
	void null();
	void freeself();
	void init(finite_field *F, INT k, INT *spread_elements_numeric, 
		INT verbose_level);
	void points_on_line(andre_construction_line_element *Line, 
		INT *pts_on_line, INT verbose_level);
	
};




// ####################################################################################
// andre_construction_point_element.C:
// ####################################################################################

class andre_construction_point_element {
public:
	andre_construction *Andre;
	INT k, n, q, spread_size;
	finite_field *F;
	INT point_rank;
	INT f_is_at_infinity;
	INT at_infinity_idx;
	INT affine_numeric;
	INT *coordinates; // [n]

	andre_construction_point_element();
	~andre_construction_point_element();
	void null();
	void freeself();
	void init(andre_construction *Andre, INT verbose_level);
	void unrank(INT point_rank, INT verbose_level);
	INT rank(INT verbose_level);
};


// ####################################################################################
// andre_construction_line_element.C:
// ####################################################################################


class andre_construction_line_element {
public:
	andre_construction *Andre;
	INT k, n, q, spread_size;
	finite_field *F;
	INT line_rank;
	INT f_is_at_infinity;
	INT affine_numeric;
	INT parallel_class_idx;
	INT coset_idx;
	INT *pivots; // [k]
	INT *non_pivots; // [n - k]
	INT *coset; // [n - k]
	INT *coordinates; // [(k + 1) * n], last row is special vector

	andre_construction_line_element();
	~andre_construction_line_element();
	void null();
	void freeself();
	void init(andre_construction *Andre, INT verbose_level);
	void unrank(INT line_rank, INT verbose_level);
	INT rank(INT verbose_level);
	INT make_affine_point(INT idx, INT verbose_level);
		// 0 \le idx \le order
};


// ####################################################################################
// memory_object.C:
// ####################################################################################




class memory_object {
public:
	memory_object();
	~memory_object();
	void null();
	void freeself();

	BYTE *char_pointer;
	INT alloc_length;
	INT used_length;
	INT cur_pointer;


	char & s_i(INT i) { return char_pointer[i]; };
	void init(INT length, char *d, INT verbose_level);
	void alloc(INT length, INT verbose_level);
	void append(INT length, char *d, INT verbose_level);
	void realloc(INT new_length, INT verbose_level);
	void write_char(char c);
	void read_char(char *c);
	void write_string(const BYTE *p);
	void read_string(BYTE *&p);
	void write_double(double f);
	void read_double(double *f);
	void write_int64(INT i);
	void read_int64(INT *i);
	void write_int(INT i);
	void read_int(INT *i);
	void read_file(const BYTE *fname, INT verbose_level);
	void write_file(const BYTE *fname, INT verbose_level);
	INT multiplicity_of_character(BYTE c);
	void compress(INT verbose_level);
	void decompress(INT verbose_level);
};


// ####################################################################################
// tree_node.C:
// ####################################################################################

class tree_node {

public:
	tree_node *parent;
	INT depth;
	INT f_value;
	INT value;
	
	INT f_int_data;
	INT int_data;
	BYTE *char_data;
	INT nb_children;
	tree_node **children;

	INT weight;
	INT placement_x;
	INT placement_y;
	INT width;

	INT DFS_rank;

	tree_node();
	~tree_node();
	void init(INT depth, tree_node *parent, INT f_value, INT value, INT f_i_data, INT i_data, BYTE *c_data, INT verbose_level);
	void print_path();
	void print_depth_first();
	void compute_DFS_rank(INT &rk);
	void get_coordinates(INT &idx, INT *coord_xy);
	void get_coordinates_and_width(INT &idx, INT *coord_xyw);
	void calc_weight();
	void place_xy(INT left, INT right, INT ymax, INT max_depth);
	void place_on_circle(INT xmax, INT ymax, INT max_depth);
	void add_node(INT l, INT depth, INT *path, INT i_data, BYTE *c_data, 
		INT verbose_level);
	INT find_child(INT val);
	void get_values(INT *v);
	void draw_edges(mp_graphics &G, INT rad, INT f_circle, INT f_circletext, INT f_i, 
		INT f_has_parent, INT parent_x, INT parent_y, INT max_depth, INT f_edge_labels, 
		INT f_has_draw_vertex_callback, 
		void (*draw_vertex_callback)(tree *T, mp_graphics *G, INT *v, INT layer, tree_node *N, INT x, INT y, INT dx, INT dy),
		tree *T
		);
	void draw_vertices(mp_graphics &G, INT rad, INT f_circle, INT f_circletext, INT f_i, 
		INT f_has_parent, INT parent_x, INT parent_y, INT max_depth, INT f_edge_labels, 
		INT f_has_draw_vertex_callback, 
		void (*draw_vertex_callback)(tree *T, mp_graphics *G, INT *v, INT layer, tree_node *N, INT x, INT y, INT dx, INT dy),
		tree *T
		);
	void draw_sideways(mp_graphics &G, INT f_circletext, INT f_i, 
		INT f_has_parent, INT parent_x, INT parent_y, INT max_depth, INT f_edge_labels);
};

INT tree_node_calc_y_coordinate(INT ymax, INT l, INT max_depth);

// ####################################################################################
// tree.C:
// ####################################################################################

class tree {

public:

	tree_node *root;
	
	INT nb_nodes;
	INT max_depth;
	
	INT *path;

	INT f_count_leaves;
	INT leaf_count;

	tree();
	~tree();
	void init(const BYTE *fname, INT xmax, INT ymax, INT verbose_level);
	void draw(char *fname, INT xmax_in, INT ymax_in, INT xmax, INT ymax, INT rad, 
		INT f_circle, INT f_circletext, INT f_i, INT f_edge_labels, 
		INT f_has_draw_vertex_callback, 
		void (*draw_vertex_callback)(tree *T, mp_graphics *G, INT *v, INT layer, tree_node *N, INT x, INT y, INT dx, INT dy), 
		INT f_embedded, INT f_sideways, INT f_on_circle, 
		double tikz_global_scale, double tikz_global_line_width
		);
	void circle_center_and_radii(INT xmax, INT ymax, INT max_depth, INT &x0, INT &y0, INT *&rad);
	void compute_DFS_ranks(INT &nb_nodes, INT verbose_level);
};


// ####################################################################################
// galois_global.C:
// ####################################################################################

void test_unipoly();
void test_unipoly2();
BYTE *search_for_primitive_polynomial_of_given_degree(INT p, INT degree, INT verbose_level);
void search_for_primitive_polynomials(INT p_min, INT p_max, INT n_min, INT n_max, INT verbose_level);
void make_linear_irreducible_polynomials(INT q, INT &nb, INT *&table, INT verbose_level);
void gl_random_matrix(INT k, INT q, INT verbose_level);
void save_colored_graph(const BYTE *fname, INT nb_vertices, INT nb_colors, 
	INT *vertex_labels, INT *vertex_colors, 
	INT *data, INT data_sz, 
	UBYTE *bitvector_adjacency, INT bitvector_length,
	INT verbose_level);
void load_colored_graph(const BYTE *fname, INT &nb_vertices, INT &nb_colors, 
	INT *&vertex_labels, INT *&vertex_colors, 
	INT *&user_data, INT &user_data_size, 
	UBYTE *&bitvector_adjacency, INT &bitvector_length,
	INT verbose_level);
INT is_diagonal_matrix(INT *A, INT n);
INT is_association_scheme(INT *color_graph, INT n, INT *&Pijk, 
	INT *&colors, INT &nb_colors, INT verbose_level);
void print_Pijk(INT *Pijk, INT nb_colors);
void write_colored_graph(ofstream &ost, BYTE *label, 
	INT point_offset, 
	INT nb_points, 
	INT f_has_adjacency_matrix, INT *Adj, 
	INT f_has_adjacency_list, INT *adj_list, 
	INT f_has_bitvector, UBYTE *bitvector_adjacency, 
	INT f_has_is_adjacent_callback, 
	INT (*is_adjacent_callback)(INT i, INT j, void *data), 
	void *is_adjacent_callback_data, 
	INT f_colors, INT nb_colors, INT *point_color, 
	INT f_point_labels, INT *point_label);
int str2int(string &str);
void print_longinteger_after_multiplying(ostream &ost, INT *factors, INT len);
void andre_preimage(projective_space *P2, projective_space *P4, 
	INT *set2, INT sz2, INT *set4, INT &sz4, INT verbose_level);
void determine_conic(INT q, const BYTE *override_poly, INT *input_pts, INT nb_pts, INT verbose_level);



// ####################################################################################
// INT_vector.C:
// ####################################################################################

class INT_vector {
public:

	INT *M;
	INT m;
	INT alloc_length;

	INT_vector();
	~INT_vector();
	void null();
	void freeself();
	void allocate(INT len);
	void allocate_and_init(INT len, INT *V);
	void init_permutation_from_string(const char *s);
	void read_ascii_file(const BYTE *fname);
	void read_binary_file_INT4(const BYTE *fname);
	INT &s_i(INT i);
	INT &length();
	void print(ostream &ost);
	void zero();
	INT search(INT a, INT &idx);
	void sort();
	void make_space();
	void append(INT a);
	void insert_at(INT a, INT idx);
	void insert_if_not_yet_there(INT a);
	void sort_and_remove_duplicates();
	void write_to_ascii_file(const BYTE *fname);
	void write_to_binary_file_INT4(const BYTE *fname);
	void write_to_csv_file(const BYTE *fname, const BYTE *label);
	INT hash();
	INT minimum();
	INT maximum();

	

};

// ####################################################################################
// INT_matrix.C:
// ####################################################################################

class INT_matrix {
public:

	INT *M;
	INT m;
	INT n;

	INT_matrix();
	~INT_matrix();
	void null();
	void freeself();
	void allocate(INT m, INT n);
	void allocate_and_init(INT m, INT n, INT *Mtx);
	INT &s_ij(INT i, INT j);
	INT &s_m();
	INT &s_n();
	void print();

	

};

// ####################################################################################
// gl_classes.C:
// ####################################################################################

class gl_classes {
public:
	INT k;
	INT q;
	finite_field *F;
	INT nb_irred;
	INT *Nb_irred;
	INT *First_irred;
	INT *Nb_part;
	INT **Tables;
	INT **Partitions;
	INT *Degree;

	gl_classes();
	~gl_classes();
	void null();
	void freeself();
	void init(INT k, finite_field *F, INT verbose_level);
	void print_polynomials(ofstream &ost);
	INT select_polynomial_first(INT *Select, INT verbose_level);
	INT select_polynomial_next(INT *Select, INT verbose_level);
	INT select_partition_first(INT *Select, INT *Select_partition, INT verbose_level);
	INT select_partition_next(INT *Select, INT *Select_partition, INT verbose_level);
	INT first(INT *Select, INT *Select_partition, INT verbose_level);
	INT next(INT *Select, INT *Select_partition, INT verbose_level);
	void print_matrix_and_centralizer_order_latex(ofstream &ost, gl_class_rep *R);
	void make_matrix_from_class_rep(INT *Mtx, gl_class_rep *R, INT verbose_level);
	void make_matrix(INT *Mtx, INT *Select, INT *Select_Partition, INT verbose_level);
	void centralizer_order_Kung_basic(INT nb_irreds, 
		INT *poly_degree, INT *poly_mult, INT *partition_idx, 
		longinteger_object &co, 
		INT verbose_level);
	void centralizer_order_Kung(INT *Select_polynomial, INT *Select_partition, longinteger_object &co, 
		INT verbose_level);
		// Computes the centralizer order of a matrix in GL(k,q) 
		// according to Kung's formula~\cite{Kung81}.
	void make_classes(gl_class_rep *&R, INT &nb_classes, INT f_no_eigenvalue_one, INT verbose_level);
	void identify_matrix(INT *Mtx, gl_class_rep *R, INT *Basis, INT verbose_level);
	void identify2(INT *Mtx, unipoly_object &poly, INT *Mult, INT *Select_partition, INT *Basis, INT verbose_level);
	void compute_data_on_blocks(INT *Mtx, INT *Irreds, INT nb_irreds, 
		INT *Degree, INT *Mult, matrix_block_data *Data,
		INT verbose_level);
	void compute_generalized_kernels(matrix_block_data *Data, INT *M2, INT d, INT b0, INT m, INT *poly_coeffs, INT verbose_level);
	INT identify_partition(INT *part, INT m, INT verbose_level);
	void choose_basis_for_rational_normal_form(INT *Mtx, matrix_block_data *Data, INT nb_irreds, 
		INT *Basis, 
		INT verbose_level);
	void choose_basis_for_rational_normal_form_block(INT *Mtx, matrix_block_data *Data, 
		INT *Basis, INT &b, 
		INT verbose_level);
	void generators_for_centralizer(INT *Mtx, gl_class_rep *R, 
		INT *Basis, INT **&Gens, INT &nb_gens, INT &nb_alloc, 
		INT verbose_level);
	void centralizer_generators(INT *Mtx, unipoly_object &poly, INT *Mult, INT *Select_partition, 
		INT *Basis, INT **&Gens, INT &nb_gens, INT &nb_alloc,  
		INT verbose_level);
	void centralizer_generators_block(INT *Mtx, matrix_block_data *Data, INT nb_irreds, INT h, 
		INT **&Gens, INT &nb_gens, INT &nb_alloc,  
		INT verbose_level);
	INT choose_basis_for_rational_normal_form_coset(INT level1, INT level2, INT &coset, 
		INT *Mtx, matrix_block_data *Data, INT &b, INT *Basis, 
		INT verbose_level);
	void factor_polynomial(unipoly_object &char_poly, INT *Mult, INT verbose_level);
	INT find_class_rep(gl_class_rep *Reps, INT nb_reps, gl_class_rep *R, INT verbose_level);

};


class gl_class_rep {
public:
	INT_matrix type_coding;
	longinteger_object centralizer_order;
	longinteger_object class_length;

	gl_class_rep();
	~gl_class_rep();
	void init(INT nb_irred, INT *Select_polynomial, INT *Select_partition, INT verbose_level);
	void compute_vector_coding(gl_classes *C, INT &nb_irred, INT *&Poly_degree, INT *&Poly_mult, INT *&Partition_idx, INT verbose_level);
	void centralizer_order_Kung(gl_classes *C, longinteger_object &co, INT verbose_level);
};


class matrix_block_data {
public:
	INT d;
	INT m;
	INT *poly_coeffs;
	INT b0;
	INT b1;
	
	INT_matrix *K;
	INT cnt;
	INT *dual_part;
	INT *part;
	INT height;
	INT part_idx;
	
	matrix_block_data();
	~matrix_block_data();
	void null();
	void freeself();
	void allocate(INT k);
};


// ##################################################################################################
// layered_graph.C
// ##################################################################################################

class layered_graph {
public:
	INT nb_layers;
	INT nb_nodes_total;
	INT id_of_first_node;
	graph_layer *L;
	BYTE fname_base[1000];
	INT data1;

	layered_graph();
	~layered_graph();
	void null();
	void freeself();
	void init(INT nb_layers, INT *Nb_nodes_layer, const BYTE *fname_base, INT verbose_level);
	void place(INT verbose_level);
	void place_with_y_stretch(double y_stretch, INT verbose_level);
	void place_with_grouping(INT **Group_sizes, INT *Nb_groups, INT verbose_level);
	void add_edge(INT l1, INT n1, INT l2, INT n2, INT verbose_level);
	void add_text(INT l, INT n, const BYTE *text, INT verbose_level);
	void add_data1(INT data, INT verbose_level);
	void add_node_vec_data(INT l, INT n, INT *v, INT len, INT verbose_level);
	void set_distinguished_element_index(INT l, INT n, INT index, INT verbose_level);
	void add_node_data1(INT l, INT n, INT data, INT verbose_level);
	void add_node_data2(INT l, INT n, INT data, INT verbose_level);
	void add_node_data3(INT l, INT n, INT data, INT verbose_level);
#if 0
	void draw(const char *fname, INT xmax, INT ymax, INT x_max, INT y_max, INT rad, 
		INT f_circle, INT f_corners, INT f_nodes_empty, 
		INT f_select_layers, INT nb_layer_select, INT *layer_select, 
		INT f_has_draw_begining_callback, 
		void (*draw_begining_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT f_rotated, INT dx, INT dy), 
		INT f_has_draw_ending_callback, 
		void (*draw_ending_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT f_rotated, INT dx, INT dy), 
		INT f_has_draw_vertex_callback, 
		void (*draw_vertex_callback)(layered_graph *LG, mp_graphics *G, INT layer, INT node, INT x, INT y, INT dx, INT dy), 
		INT f_show_level_info, 
		INT f_embedded, INT f_sideways, 
		INT f_label_edges, 
		INT f_rotated, 
		double global_scale, double global_line_width);
#endif
	void draw_with_options(const char *fname, layered_graph_draw_options *O, INT verbose_level);
	void coordinates_direct(double x_in, double y_in, INT x_max, INT y_max, INT f_rotated, INT &x, INT &y);
	void coordinates(INT id, INT x_max, INT y_max, INT f_rotated, INT &x, INT &y);
	void find_node_by_id(INT id, INT &l, INT &n);
	void write_file(BYTE *fname, INT verbose_level);
	void read_file(const BYTE *fname, INT verbose_level);
	void write_memory_object(memory_object *m, INT verbose_level);
	void read_memory_object(memory_object *m, INT verbose_level);
};

// ##################################################################################################
// layered_graph_draw_options.C
// ##################################################################################################

class layered_graph_draw_options {
public:

	INT xmax;
	INT ymax;
	INT x_max;
	INT y_max;
	INT rad;
	
	INT f_circle;
	INT f_corners;
	INT f_nodes_empty;
	INT f_select_layers;
	INT nb_layer_select;
	INT *layer_select;


	INT f_has_draw_begining_callback;
	void (*draw_begining_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT f_rotated, INT dx, INT dy);
	INT f_has_draw_ending_callback;
	void (*draw_ending_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT f_rotated, INT dx, INT dy);
	INT f_has_draw_vertex_callback;
	void (*draw_vertex_callback)(layered_graph *LG, mp_graphics *G, INT layer, INT node, INT x, INT y, INT dx, INT dy);
	
	INT f_show_level_info;
	INT f_embedded;
	INT f_sideways;
	INT f_label_edges;
	INT f_rotated;
	
	double global_scale;
	double global_line_width;

	layered_graph_draw_options();
	~layered_graph_draw_options();
	void init(
		INT xmax, INT ymax, INT x_max, INT y_max, INT rad, 
		INT f_circle, INT f_corners, INT f_nodes_empty, 
		INT f_select_layers, INT nb_layer_select, INT *layer_select, 
		INT f_has_draw_begining_callback, 
		void (*draw_begining_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT f_rotated, INT dx, INT dy), 
		INT f_has_draw_ending_callback, 
		void (*draw_ending_callback)(layered_graph *LG, mp_graphics *G, INT x_max, INT y_max, INT f_rotated, INT dx, INT dy), 
		INT f_has_draw_vertex_callback, 
		void (*draw_vertex_callback)(layered_graph *LG, mp_graphics *G, INT layer, INT node, INT x, INT y, INT dx, INT dy), 
		INT f_show_level_info, 
		INT f_embedded, INT f_sideways, 
		INT f_label_edges, 
		INT f_rotated, 
		double global_scale, double global_line_width);
};


// ##################################################################################################
// graph_layer.C
// ##################################################################################################

class graph_layer {
public:
	INT id_of_first_node;
	INT nb_nodes;
	graph_node *Nodes;
	double y_coordinate;

	graph_layer();
	~graph_layer();
	void null();
	void freeself();
	void init(INT nb_nodes, INT id_of_first_node, INT verbose_level);
	void place(INT verbose_level);
	void place_with_grouping(INT *group_size, INT nb_groups, INT verbose_level);
	void write_memory_object(memory_object *m, INT verbose_level);
	void read_memory_object(memory_object *m, INT verbose_level);
};

// ##################################################################################################
// graph_node.C
// ##################################################################################################

class graph_node {
public:
	BYTE *label;
	INT id;

	INT f_has_data1;
	INT data1;

	INT f_has_data2;
	INT data2;

	INT f_has_data3;
	INT data3;

	INT f_has_vec_data;
	INT *vec_data;
	INT vec_data_len;

	INT f_has_distinguished_element; // refers to vec_data
	INT distinguished_element_index;
		
	INT layer;
	INT neighbor_list_allocated;
	INT nb_neighbors;
	INT *neighbor_list;
	double x_coordinate;

	graph_node();
	~graph_node();
	void null();
	void freeself();
	void add_neighbor(INT l, INT n, INT id);
	void add_text(const BYTE *text);
	void add_vec_data(INT *v, INT len);
	void set_distinguished_element(INT idx);
	void add_data1(INT data);
	void add_data2(INT data);
	void add_data3(INT data);
	void write_memory_object(memory_object *m, INT verbose_level);
	void read_memory_object(memory_object *m, INT verbose_level);
};


// #################################################################################
// nauty_interface.C:
// #################################################################################

void nauty_interface_graph_bitvec(INT v, UBYTE *bitvector_adjacency, 
	INT *labeling, INT *partition, 
	INT *Aut, INT &Aut_counter, 
	INT *Base, INT &Base_length, 
	INT *Transversal_length, INT &Ago, INT verbose_level);
void nauty_interface_graph_INT(INT v, INT *Adj, 
	INT *labeling, INT *partition, 
	INT *Aut, INT &Aut_counter, 
	INT *Base, INT &Base_length, 
	INT *Transversal_length, INT &Ago, INT verbose_level);
void nauty_interface_INT(INT v, INT b, INT *X, INT nb_inc, 
	INT *labeling, INT *partition, 
	INT *Aut, INT &Aut_counter, 
	INT *Base, INT &Base_length, 
	INT *Transversal_length, INT &Ago);
void nauty_interface(int v, int b, INT *X, INT nb_inc, 
	int *labeling, int *partition, 
	INT *Aut, int &Aut_counter, 
	int *Base, int &Base_length, 
	int *Transversal_length, int &Ago);
void nauty_interface_matrix(int *M, int v, int b, 
	int *labeling, int *partition, 
	INT *Aut, int &Aut_counter, 
	int *Base, int &Base_length, 
	int *Transversal_length, int &Ago);
void nauty_interface_matrix_INT(INT *M, INT v, INT b, 
	INT *labeling, INT *partition, 
	INT *Aut, INT &Aut_counter, 
	INT *Base, INT &Base_length, 
	INT *Transversal_length, INT &Ago, INT verbose_level);

// ####################################################################################
// projective_space.C:
// ####################################################################################

class projective_space {

public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;




	grassmann *Grass_lines;
	finite_field *F;
	longinteger_object *Go;

	INT n; // projective dimension
	INT q;
	INT N_points, N_lines;
	INT r; // number of lines on a point
	INT k; // number of points on a line


	UBYTE *incidence_bitvec; // N_points * N_lines bits
	INT *Lines; // [N_lines * k]
	INT *Lines_on_point; // [N_points * r]
	INT *Line_through_two_points; // [N_points * N_points]
	INT *Line_intersection;	// [N_lines * N_lines]

	// only if n = 2:
	INT *Polarity_point_to_hyperplane; // [N_points]
	INT *Polarity_hyperplane_to_point; // [N_points]

	INT *v; // [n + 1]
	INT *w; // [n + 1]

	projective_space();
	~projective_space();
	void null();
	void freeself();
	void init(INT n, finite_field *F, 
		INT f_init_incidence_structure, 
		INT verbose_level);
	void init_incidence_structure(INT verbose_level);
	void make_incidence_matrix(INT &m, INT &n, INT *&Inc, INT verbose_level);
	INT is_incident(INT pt, INT line);
	void incidence_m_ii(INT pt, INT line, INT a);
	void make_incidence_structure_and_partition(incidence_structure *&Inc, 
		partitionstack *&Stack, INT verbose_level);
	INT nb_rk_k_subspaces_as_INT(INT k);
	void print_all_points();
	INT rank_point(INT *v);
	void unrank_point(INT *v, INT rk);
	INT rank_line(INT *basis);
	void unrank_line(INT *basis, INT rk);
	INT line_through_two_points(INT p1, INT p2);
	INT test_if_lines_are_disjoint(INT l1, INT l2);
	INT test_if_lines_are_disjoint_from_scratch(INT l1, INT l2);
	INT line_intersection(INT l1, INT l2);
		// works only for projective planes, i.e., n = 2
	
	INT arc_test(INT *input_pts, INT nb_pts, INT verbose_level);
	INT determine_conic_in_plane(INT *input_pts, INT nb_pts, INT *six_coeffs, INT verbose_level);
		// returns FALSE is the rank of the coefficient matrix is not 5. TRUE otherwise.

	void determine_quadric_in_solid(INT *nine_pts_or_more, INT nb_pts, INT *ten_coeffs, INT verbose_level);
	void conic_points_brute_force(INT *six_coeffs, 
		INT *points, INT &nb_points, INT verbose_level);
	void quadric_points_brute_force(INT *ten_coeffs, 
		INT *points, INT &nb_points, INT verbose_level);
	void conic_points(INT *five_pts, INT *six_coeffs, 
		INT *points, INT &nb_points, INT verbose_level);
	void find_tangent_lines_to_conic(INT *six_coeffs, 
		INT *points, INT nb_points, 
		INT *tangents, INT verbose_level);
	void PG_2_8_create_conic_plus_nucleus_arc_1(INT *the_arc, INT &size, INT verbose_level);
	void PG_2_8_create_conic_plus_nucleus_arc_2(INT *the_arc, INT &size, INT verbose_level);
	void create_Maruta_Hamada_arc(INT *the_arc, INT &size, INT verbose_level);
	void create_Maruta_Hamada_arc2(INT *the_arc, INT &size, INT verbose_level);
	void create_pasch_arc(INT *the_arc, INT &size, INT verbose_level);
	void create_Cheon_arc(INT *the_arc, INT &size, INT verbose_level);
	void create_regular_hyperoval(INT *the_arc, INT &size, INT verbose_level);
	void create_translation_hyperoval(INT *the_arc, INT &size, 
		INT exponent, INT verbose_level);
	void create_Segre_hyperoval(INT *the_arc, INT &size, INT verbose_level);
	void create_Payne_hyperoval(INT *the_arc, INT &size, INT verbose_level);
	void create_Cherowitzo_hyperoval(INT *the_arc, INT &size, INT verbose_level);
	void create_OKeefe_Penttila_hyperoval_32(INT *the_arc, INT &size, INT verbose_level);
	void line_intersection_type(INT *set, INT set_size, INT *type, INT verbose_level);
	void line_intersection_type_basic(INT *set, INT set_size, INT *type, INT verbose_level);
		// type[N_lines]
	void plane_intersection_type_basic(INT *set, INT set_size, INT *type, INT verbose_level);
		// type[N_planes]
	void hyperplane_intersection_type_basic(INT *set, INT set_size, INT *type, INT verbose_level);
		// type[N_hyperplanes]
	void line_intersection_type_collected(INT *set, INT set_size, 
		INT *type_collected, INT verbose_level);
		// type[set_size + 1]
	void point_types(INT *set_of_lines, INT set_size, INT *type, INT verbose_level);
	void find_external_lines(INT *set, INT set_size, 
		INT *external_lines, INT &nb_external_lines, INT verbose_level);
	void find_tangent_lines(INT *set, INT set_size, 
		INT *tangent_lines, INT &nb_tangent_lines, INT verbose_level);
	void find_secant_lines(INT *set, INT set_size, 
		INT *secant_lines, INT &nb_secant_lines, INT verbose_level);
	void find_k_secant_lines(INT *set, INT set_size, INT k, 
		INT *secant_lines, INT &nb_secant_lines, INT verbose_level);
	void Baer_subline(INT *pts3, INT *&pts, INT &nb_pts, INT verbose_level);
	INT is_contained_in_Baer_subline(INT *pts, INT nb_pts, INT verbose_level);
	void print_set_numerical(INT *set, INT set_size);
	void print_set(INT *set, INT set_size);
	void print_line_set_numerical(INT *set, INT set_size);
	INT determine_hermitian_form_in_plane(INT *pts, INT nb_pts, INT *six_coeffs, INT verbose_level);
	void circle_type_of_line_subset(INT *pts, INT nb_pts, INT *circle_type, INT verbose_level);
		// circle_type[nb_pts]
	void create_unital_XXq_YZq_ZYq(INT *U, INT &sz, INT verbose_level);
	void intersection_of_subspace_with_point_set(
		grassmann *G, INT rk, INT *set, INT set_size, 
		INT *&intersection_set, INT &intersection_set_size, INT verbose_level);
	void intersection_of_subspace_with_point_set_rank_is_longinteger(
		grassmann *G, longinteger_object &rk, INT *set, INT set_size, 
		INT *&intersection_set, INT &intersection_set_size, INT verbose_level);
	void plane_intersection_invariant(grassmann *G, 
		INT *set, INT set_size, 
		INT *&intersection_type, INT &highest_intersection_number, 
		INT *&intersection_matrix, INT &nb_planes, 
		INT verbose_level);
	void plane_intersection_type(grassmann *G, 
		INT *set, INT set_size, 
		INT *&intersection_type, INT &highest_intersection_number, 
		INT verbose_level);
	void plane_intersections(grassmann *G, 
		INT *set, INT set_size, 
		longinteger_object *&R, set_of_sets &SoS, 
		INT verbose_level);
	void plane_intersection_type_slow(grassmann *G, 
		INT *set, INT set_size, 
		longinteger_object *&R, INT **&Pts_on_plane, INT *&nb_pts_on_plane, INT &len, 
		INT verbose_level);
	void plane_intersection_type_fast(grassmann *G, 
		INT *set, INT set_size, 
		longinteger_object *&R, INT **&Pts_on_plane, INT *&nb_pts_on_plane, INT &len, 
		INT verbose_level);
	void klein_correspondence(projective_space *P5, 
		INT *set_in, INT set_size, INT *set_out, INT verbose_level);
	// Computes the Pluecker coordinates for a line in PG(3,q) in the following order:
	// (x_1,x_2,x_3,x_4,x_5,x_6) = 
	// (Pluecker_12, Pluecker_34, Pluecker_13, Pluecker_42, Pluecker_14, Pluecker_23)
	// satisfying the quadratic form x_1x_2 + x_3x_4 + x_5x_6 = 0
	void Pluecker_coordinates(INT line_rk, INT *v6, INT verbose_level);
	void klein_correspondence_special_model(projective_space *P5, 
		INT *table, INT verbose_level);
	void cheat_sheet_points_and_lines(ostream &f, INT verbose_level);
	void cheat_sheet_points(ostream &f, INT verbose_level);
	void cheat_sheet_subspaces(ostream &f, INT k, INT verbose_level);
	void cheat_sheet_line_intersection(ostream &f, INT verbose_level);
	void cheat_sheet_line_through_pairs_of_points(ostream &f, INT verbose_level);
	void conic_type_randomized(INT nb_times, 
		INT *set, INT set_size, 
		INT **&Pts_on_conic, INT *&nb_pts_on_conic, INT &len, 
		INT verbose_level);
	void conic_intersection_type(INT f_randomized, INT nb_times, 
		INT *set, INT set_size, 
		INT *&intersection_type, INT &highest_intersection_number, 
		INT f_save_largest_sets, set_of_sets *&largest_sets, 
		INT verbose_level);
	void conic_type(
		INT *set, INT set_size, 
		INT **&Pts_on_conic, INT *&nb_pts_on_conic, INT &len, 
		INT verbose_level);
	void find_nucleus(INT *set, INT set_size, INT &nucleus, INT verbose_level);
};


// ####################################################################################
// buekenhout_metz.C:
// ####################################################################################

class buekenhout_metz {
public:
	finite_field *FQ, *Fq;
	INT q;
	INT Q;

	INT f_classical;
	INT f_Uab;
	INT parameter_a;
	INT parameter_b;

	projective_space *P2; // PG(2,q^2), where the unital lives
	projective_space *P3; // PG(3,q), where the ovoid lives

	INT *v; // [3]
	INT *w1; // [6]
	INT *w2; // [6]
	INT *w3; // [6]
	INT *w4; // [6]
	INT *w5; // [6]
	INT *components;
	INT *embedding;
	INT *pair_embedding;
	INT *ovoid;
	INT *U;
	INT sz;
	INT alpha, t0, t1, T0, T1, theta_3, minus_t0, sz_ovoid;
	INT e1, one_1, one_2;


	// compute_the_design:
	INT *secant_lines;
	INT nb_secant_lines;
	INT *tangent_lines;
	INT nb_tangent_lines;
	INT *Intersection_sets;
	INT *Design_blocks;
	INT *block;
	INT block_size;
	INT *idx_in_unital;
	INT *idx_in_secants;
	INT *tangent_line_at_point;
	INT *point_of_tangency;
	INT *f_is_tangent_line;
	INT *f_is_Baer;

#if 0
	// compute_automorphism_group
	action *A;
	longinteger_object ago, ago2;
	sims *S;
	vector_ge *gens;
	INT *tl;
	BYTE fname_stab[1000];

	// compute_orbits:
	INT f_prefered_line_reps;
	INT *prefered_line_reps;
	INT nb_prefered_line_reps;
	schreier *Orb;
	schreier *Orb2;

	
	// investigate_line_orbit:
	sims *Stab;
	choose_points_or_lines *C;
#endif

	// the block that we choose:
	INT nb_good_points;
	INT *good_points; // = q + 1


	buekenhout_metz();
	~buekenhout_metz();
	void null();
	void freeself();
	void init(finite_field *Fq, finite_field *FQ, 
		INT f_Uab, INT a, INT b, 
		INT f_classical, INT verbose_level);
	void init_ovoid(INT verbose_level);
	void init_ovoid_Uab_even(INT a, INT b, INT verbose_level);
	void create_unital(INT verbose_level);
	void create_unital_tex(INT verbose_level);
	void create_unital_Uab_tex(INT verbose_level);
	void compute_the_design(INT verbose_level);
#if 0
	void compute_automorphism_group(INT verbose_level);
	void compute_orbits(INT verbose_level);
	void investigate_line_orbit(INT h, INT verbose_level);
#endif
	void write_unital_to_file();
	void get_name(BYTE *name);

};


INT buekenhout_metz_check_good_points(INT len, INT *S, void *data, INT verbose_level);


// ####################################################################################
// geometric_object.C:
// ####################################################################################

void do_cone_over(INT n, finite_field *F, 
	INT *set_in, INT set_size_in, INT *&set_out, INT &set_size_out, 
	INT verbose_level);
void do_blocking_set_family_3(INT n, finite_field *F, 
	INT *set_in, INT set_size, 
	INT *&the_set_out, INT &set_size_out, 
	INT verbose_level);
void create_hyperoval(finite_field *F, 
	INT f_translation, INT translation_exponent, 
	INT f_Segre, INT f_Payne, INT f_Cherowitzo, INT f_OKeefe_Penttila, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_subiaco_oval(finite_field *F, 
	INT f_short, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_subiaco_hyperoval(finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_adelaide_hyperoval(subfield_structure *S, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_ovoid(finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_Baer_substructure(INT n, finite_field *FQ, finite_field *Fq, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_BLT_from_database(INT f_embedded, finite_field *F, INT BLT_k, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_BLT(INT f_embedded, finite_field *FQ, finite_field *Fq, 
	INT f_Linear,
	INT f_Fisher,
	INT f_Mondello,
	INT f_FTWKB,
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_orthogonal(INT epsilon, INT n, finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_hermitian(INT n, finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_twisted_cubic(finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_ttp_code(finite_field *FQ, finite_field *Fq, 
	INT f_construction_A, INT f_hyperoval, INT f_construction_B, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_unital_XXq_YZq_ZYq(finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_desarguesian_line_spread_in_PG_3_q(finite_field *FQ, finite_field *Fq, 
	INT f_embedded_in_PG_4_q, 
	BYTE *fname, INT &nb_lines, INT *&Lines, 
	INT verbose_level);
void create_whole_space(INT n, finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_hyperplane(INT n, finite_field *F, 
	INT pt, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_segre_variety(finite_field *F, INT a, INT b, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);
void create_Maruta_Hamada_arc(finite_field *F, 
	BYTE *fname, INT &nb_pts, INT *&Pts, 
	INT verbose_level);

// ####################################################################################
// geometric_operations.C:
// ####################################################################################

void do_Klein_correspondence(INT n, finite_field *F, 
	INT *set_in, INT set_size,
	INT *&the_set_out, INT &set_size_out, 
	INT verbose_level);
void do_m_subspace_type(INT n, finite_field *F, INT m, 
	INT *set, INT set_size, 
	INT f_show, INT verbose_level);
void do_m_subspace_type_fast(INT n, finite_field *F, INT m, 
	INT *set, INT set_size, 
	INT f_show, INT verbose_level);
void do_line_type(INT n, finite_field *F, 
	INT *set, INT set_size, 
	INT f_show, INT verbose_level);
void do_plane_type(INT n, finite_field *F, 
	INT *set, INT set_size, 
	INT *&intersection_type, INT &highest_intersection_number, 
	INT verbose_level);
void do_plane_type_failsafe(INT n, finite_field *F, 
	INT *set, INT set_size, 
	INT verbose_level);
void do_conic_type(INT n, finite_field *F, INT f_randomized, INT nb_times, 
	INT *set, INT set_size, 
	INT *&intersection_type, INT &highest_intersection_number, 
	INT verbose_level);
void do_test_diagonal_line(INT n, finite_field *F, 
	INT *set_in, INT set_size, 
	BYTE *fname_orbits_on_quadrangles, 
	INT verbose_level);
void do_andre(finite_field *FQ, finite_field *Fq, 
	INT *the_set_in, INT set_size_in, 
	INT *&the_set_out, INT &set_size_out, 
	INT verbose_level);
void do_print_lines_in_PG(INT n, finite_field *F, 
	INT *set_in, INT set_size);
void do_print_points_in_PG(INT n, finite_field *F, 
	INT *set_in, INT set_size);
void do_print_points_in_orthogonal_space(INT epsilon, INT n, finite_field *F,  
	INT *set_in, INT set_size, INT verbose_level);
void do_print_points_on_grassmannian(INT n, INT k, finite_field *F, 
	INT *set_in, INT set_size);
void do_embed_orthogonal(INT epsilon, INT n, finite_field *F, 
	INT *set_in, INT *&set_out, INT set_size, INT verbose_level);
void do_embed_points(INT n, finite_field *F, 
	INT *set_in, INT *&set_out, INT set_size, INT verbose_level);

#if 0
void do_move_line_in_PG(INT n, finite_field *F, 
	INT from_line, INT to_line, 
	INT *the_set_in, INT set_size_in, 
	INT *&the_set_out, INT &set_size_out, 
	INT verbose_level);
void do_group_in_PG(INT n, finite_field *F, 
	INT *the_set_in, INT set_size_in, INT f_list_group_elements, INT verbose_level);
#endif


// ####################################################################################
// a_domain.C:
// ####################################################################################

enum domain_kind {
	not_applicable, domain_the_integers, domain_integer_fractions
};

class a_domain {
public:
	domain_kind kind;
	INT size_of_instance_in_INT;
	
	a_domain();
	~a_domain();
	void null();
	void freeself();
	
	void init_integers(INT verbose_level);
	void init_integer_fractions(INT verbose_level);
	INT as_INT(INT *elt, INT verbose_level);
	void make_integer(INT *elt, INT n, INT verbose_level);
	void make_zero(INT *elt, INT verbose_level);
	void make_zero_vector(INT *elt, INT len, INT verbose_level);
	INT is_zero_vector(INT *elt, INT len, INT verbose_level);
	INT is_zero(INT *elt, INT verbose_level);
	void make_one(INT *elt, INT verbose_level);
	INT is_one(INT *elt, INT verbose_level);
	void copy(INT *elt_from, INT *elt_to, INT verbose_level);
	void copy_vector(INT *elt_from, INT *elt_to, INT len, INT verbose_level);
	void swap_vector(INT *elt1, INT *elt2, INT n, INT verbose_level);
	void swap(INT *elt1, INT *elt2, INT verbose_level);
	void add(INT *elt_a, INT *elt_b, INT *elt_c, INT verbose_level);
	void add_apply(INT *elt_a, INT *elt_b, INT verbose_level);
	void subtract(INT *elt_a, INT *elt_b, INT *elt_c, INT verbose_level);
	void negate(INT *elt, INT verbose_level);
	void negate_vector(INT *elt, INT len, INT verbose_level);
	void mult(INT *elt_a, INT *elt_b, INT *elt_c, INT verbose_level);
	void mult_apply(INT *elt_a, INT *elt_b, INT verbose_level);
	void power(INT *elt_a, INT *elt_b, INT n, INT verbose_level);
	void divide(INT *elt_a, INT *elt_b, INT *elt_c, INT verbose_level);
	void inverse(INT *elt_a, INT *elt_b, INT verbose_level);
	void print(INT *elt);
	void print_vector(INT *elt, INT n);
	void print_matrix(INT *A, INT m, INT n);
	void make_element_from_integer(INT *elt, INT n, INT verbose_level);
	void mult_by_integer(INT *elt, INT n, INT verbose_level);
	void divide_by_integer(INT *elt, INT n, INT verbose_level);
	INT *offset(INT *A, INT i);
	INT Gauss_echelon_form(INT *A, INT f_special, INT f_complete, INT *base_cols, 
		INT f_P, INT *P, INT m, INT n, INT Pn, INT verbose_level);
	// returns the rank which is the number of entries in base_cols
	// A is a m x n matrix,
	// P is a m x Pn matrix (if f_P is TRUE)
	void Gauss_step(INT *v1, INT *v2, INT len, INT idx, INT verbose_level);
	// afterwards: v2[idx] = 0 and v1,v2 span the same space as before
	// v1 is not changed if v1[idx] is nonzero
	void matrix_get_kernel(INT *M, INT m, INT n, INT *base_cols, INT nb_base_cols, 
		INT &kernel_m, INT &kernel_n, INT *kernel, INT verbose_level);
	// kernel must point to the appropriate amount of memory! (at least n * (n - nb_base_cols) INT's)
	// kernel is stored as column vectors, i.e. kernel_m = n and kernel_n = n - nb_base_cols.
	void matrix_get_kernel_as_row_vectors(INT *M, INT m, INT n, INT *base_cols, INT nb_base_cols, 
		INT &kernel_m, INT &kernel_n, INT *kernel, INT verbose_level);
	// kernel must point to the appropriate amount of memory! (at least n * (n - nb_base_cols) INT's)
	// kernel is stored as row vectors, i.e. kernel_m = n - nb_base_cols and kernel_n = n.
	void get_image_and_kernel(INT *M, INT n, INT &rk, INT verbose_level);
	void complete_basis(INT *M, INT m, INT n, INT verbose_level);
	void mult_matrix(INT *A, INT *B, INT *C, INT ma, INT na, INT nb, INT verbose_level);
	void mult_matrix3(INT *A, INT *B, INT *C, INT *D, INT n, INT verbose_level);
	void add_apply_matrix(INT *A, INT *B, INT m, INT n, INT verbose_level);
	void matrix_mult_apply_scalar(INT *A, INT *s, INT m, INT n, INT verbose_level);
	void make_block_matrix_2x2(INT *Mtx, INT n, INT k, INT *A, INT *B, INT *C, INT *D, INT verbose_level);
	// A is k x k, B is k x (n - k), C is (n - k) x k, D is (n - k) x (n - k), Mtx is n x n
	void make_identity_matrix(INT *A, INT n, INT verbose_level);
	void matrix_inverse(INT *A, INT *Ainv, INT n, INT verbose_level);
	void matrix_invert(INT *A, INT *T, INT *basecols, INT *Ainv, INT n, INT verbose_level);

};


// #################################################################################
// diophant.C:
// #################################################################################



class diophant {
public:
	void *operator new(size_t bytes);
	void *operator new[](size_t bytes);
	void operator delete(void *ptr, size_t bytes);
	void operator delete[](void *ptr, size_t bytes);
	static INT cntr_new;
	static INT cntr_objects;
	static INT f_debug_memory;

	BYTE label[1000];
	INT m; // number of equations or inequalities
	INT n; // number of indeterminates
	INT sum; /* constraint: sum(i=0..(n-1); x[i]) = sum */
	INT sum1;
	INT f_x_max;
	/* with constraints: x[i] <= x_max[i] for i=0..(n-1) */
	
	INT *A; // [m][n] the coefficient matrix
	INT *G; // [m][n] matrix of gcd values
	INT *x_max; // [n] upper bounds for x
	INT *x; // [n]  current value of x
	INT *RHS; // [m] the right hand sides
	INT *RHS1; // [m] he current values of the RHS (=RHS - what is chose on the left
	diophant_equation_type *type;
	//INT *f_le; // [m] TRUE the i-th equation an inequality
	BYTE **eqn_label; // [m] a label for each equation / inequality

	INT *X; // [n]
	INT *Y; // [m]
	
	deque<vector<int> > _results;
	int _maxresults;
	int _resultanz;
	int _cur_result;
	INT nb_steps_betten;
	INT f_max_time;
	INT f_broken_off_because_of_maxtime;
	INT max_time_in_sec;
	INT max_time_in_ticks;
	INT t0;


	diophant();
	~diophant();
	void null();
	void freeself();
	
	void open(INT m, INT n);
	void join_problems(diophant *D1, diophant *D2, INT verbose_level);
	void init_problem_of_Steiner_type_with_RHS(INT nb_rows, INT nb_cols, INT *Inc, INT nb_to_select, 
		INT *Rhs, INT verbose_level);
	void init_problem_of_Steiner_type(INT nb_rows, INT nb_cols, INT *Inc, INT nb_to_select, INT verbose_level);
	void init_clique_finding_problem(INT *Adj, INT nb_pts, INT nb_to_select, INT verbose_level);
	void fill_coefficient_matrix_with(INT a);
	INT &Aij(INT i, INT j);
	INT &Gij(INT i, INT j);
	INT &RHSi(INT i);
	void init_eqn_label(INT i, BYTE *label);
	void print();
	void print_tight();
	void print_dense();
	void print2(INT f_with_gcd);
	void print_compressed();
	void print_eqn(INT i, INT f_with_gcd);
	void print_eqn_compressed(INT i);
	void print_eqn_dense(INT i);
	void print_x_long();
	void print_x(INT header);
	INT RHS_ge_zero();
	INT solve_first(INT verbose_level);
	INT solve_next();
	INT solve_first_wassermann(INT verbose_level);
	void write_solutions(const BYTE *fname, INT verbose_level);
	void read_solutions_from_file(const BYTE *fname_sol, INT verbose_level);
	void get_solutions(INT *&Sol, INT &nb_sol, INT verbose_level);
	INT solve_all_DLX(INT f_write_tree, const BYTE *fname_tree, INT verbose_level);
	INT solve_all_DLX_with_RHS(INT f_write_tree, const BYTE *fname_tree, INT verbose_level);
	INT solve_all_betten(INT verbose_level);
	INT solve_all_betten_with_conditions(INT verbose_level, 
		INT f_max_sol, INT max_sol, 
		INT f_max_time, INT max_time_in_seconds);
	INT solve_first_betten(INT verbose_level);
	INT solve_next_betten(INT verbose_level);
	INT j_fst(INT j, INT verbose_level);
	INT j_nxt(INT j, INT verbose_level);
	void latex_it();
	void latex_it(ostream &ost);
	void trivial_row_reductions(INT &f_no_solution, INT verbose_level);
	INT count_non_zero_coefficients_in_row(INT i);
	void coefficient_values_in_row(INT i, INT &nb_values, 
		INT *&values, INT *&multiplicities, INT verbose_level);
	INT maximum_number_of_non_zero_coefficients_in_row();
	void save_in_compact_format(const BYTE *fname, INT verbose_level);
	void read_compact_format(const BYTE *fname, INT verbose_level);
	void save_in_general_format(const BYTE *fname, INT verbose_level);
	void read_general_format(const BYTE *fname, INT verbose_level);
	void save_in_wassermann_format(const BYTE *fname, INT verbose_level);
	void solve_wassermann(INT verbose_level);
	void eliminate_zero_rows_quick(INT verbose_level);
	void eliminate_zero_rows(INT *&eqn_number, INT verbose_level);
	INT is_zero_outside(INT first, INT len, INT i);
	void project(diophant *D, INT first, INT len, INT *&eqn_number, INT &nb_eqns_replaced, INT *&eqns_replaced, INT verbose_level);
	void multiply_A_x_to_RHS1();
	void write_xml(ostream &ost, const BYTE *label);
	void read_xml(ifstream &f, BYTE *label);
		// label will be set to the label that is in the file
		// therefore, label must point to sufficient memory
	void append_equation();
	void delete_equation(INT I);
	void write_gurobi_binary_variables(const BYTE *fname);
	void draw_it(const BYTE *fname_base, INT xmax_in, INT ymax_in, INT xmax_out, INT ymax_out);
	INT test_solution(INT *sol, INT len, INT verbose_level);
	void get_columns(INT *col, INT nb_col, set_of_sets *&S, INT verbose_level);
	void test_solution_file(const BYTE *solution_file, INT verbose_level);
	void analyze(INT verbose_level);
	INT is_of_Steiner_type();
	void make_clique_graph_adjacency_matrix(UBYTE *&Adj, INT verbose_level);
	void make_clique_graph(colored_graph *&CG, INT verbose_level);
	void make_clique_graph_and_save(const BYTE *clique_graph_fname, INT verbose_level);
};

void diophant_callback_solution_found(INT *sol, INT len, INT nb_sol, void *data);


// ##################################################################################################
// mindist.C:
// ##################################################################################################

int mindist(int n, int k, int q, int *G, 
	int f_v, int f_vv, int idx_zero, int idx_one, 
	INT *add_table, INT *mult_table);
//Main routine for the code minimum distance computation.
//The tables are only needed if $q = p^f$ with $f > 1$. 
//In the GF(p) case, just pass a NULL pointer. 


// ##################################################################################################
// null_polarity_generator.C:
// ##################################################################################################

class null_polarity_generator {
public:

	finite_field *F; // no ownership, do not destroy
	INT n, q;
	INT qn; // = q^n

	INT *nb_candidates; // [n + 1]
	INT *cur_candidate; // [n]
	INT **candidates; // [n + 1][q^n]
	
	INT *Mtx; // [n * n]
	INT *v; // [n]
	INT *w; // [n]
	INT *Points; // [qn * n]

	INT nb_gens;
	INT *Data;
	INT *transversal_length;

	null_polarity_generator();
	~null_polarity_generator();
	void null();
	void freeself();
	void init(finite_field *F, INT n, INT verbose_level);
	INT count_strong_generators(INT &nb, INT *transversal_length, INT &first_moved, INT depth, INT verbose_level);
	INT get_strong_generators(INT *Data, INT &nb, INT &first_moved, INT depth, INT verbose_level);
	void backtrack_search(INT &nb_sol, INT depth, INT verbose_level);
	void create_first_candidate_set(INT verbose_level);
	void create_next_candidate_set(INT level, INT verbose_level);
	INT dot_product(INT *u1, INT *u2);
};

// ##################################################################################################
// klein_correspondence.C:
// ##################################################################################################


class klein_correspondence {
public:

	projective_space *P3;
	projective_space *P5;
	orthogonal *O;
	finite_field *F;
	INT q;
	INT nb_Pts; // number of points on the klein quadric
	INT nb_pts_PG; // number of points in PG(5,q)

	grassmann *Gr63;
	grassmann *Gr62;

	
	INT *Form; // [d * d]
	INT *Line_to_point_on_quadric; // [P3->N_lines]
	INT *Point_on_quadric_to_line; // [P3->N_lines]
	INT *Point_on_quadric_embedded_in_P5; // [P3->N_lines]
	INT *coordinates_of_quadric_points; // [P3->N_lines * d]
	INT *Pt_rk; // [P3->N_lines]
	INT *Pt_idx; // [nb_pts_PG]

	klein_correspondence();
	~klein_correspondence();
	void null();
	void freeself();
	void init(finite_field *F, orthogonal *O, INT verbose_level);
	void plane_intersections(INT *lines_in_PG3, INT nb_lines, 
		longinteger_object *&R,
		INT **&Pts_on_plane, 
		INT *&nb_pts_on_plane, 
		INT &nb_planes, 
		INT verbose_level);
};

// ##################################################################################################
// file_output.C:
// ##################################################################################################


class file_output {
public:
	BYTE fname[1000];
	INT f_file_is_open;
	ofstream *fp;
	void *user_data;
	
	file_output();
	~file_output();
	void null();
	void freeself();
	void open(const BYTE *fname, void *user_data, INT verbose_level);
	void close();
	void write_line(INT nb, INT *data, INT verbose_level);
};



