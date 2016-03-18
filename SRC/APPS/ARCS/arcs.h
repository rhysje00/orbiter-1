// arcs.h
//
// Anton Betten
// December 6, 2004


typedef class arc_generator arc_generator;

extern INT t0; // the system time when the program started


class arc_generator {

public:

	INT q;
	finite_field *F;
	int argc;
	const char **argv;

	exact_cover_arguments *ECA;
	isomorph_arguments *IA;

	INT verbose_level;
	INT f_starter;
	INT f_draw_poset;


	INT nb_points_total;
	INT f_target_size;
	INT target_size;

	BYTE starter_directory_name[1000];
	BYTE prefix[1000];
	BYTE prefix_with_directory[1000];
	INT starter_size;


	INT f_recognize;
	INT *recognize_set;
	INT recognize_set_sz;



	INT f_no_arc_testing;


	INT f_semilinear;

	action *A;
	
	grassmann *Grass;
	action_on_grassmannian *AG;
	action *A_on_lines;
	
	projective_space *P2;
	
	rank_checker rc;
		
	generator *gen;

	INT *Data1; // [P2->N_points * 3]
	INT *Data2; // [P2->N_points * 3]
	INT *Data3; // [9]
	


	arc_generator();
	~arc_generator();
	void null();
	void freeself();
	void read_arguments(int argc, const char **argv);
	void main(INT verbose_level);
	void init(finite_field *F,
		const BYTE *input_prefix, 
		const BYTE *base_fname,
		INT starter_size,  
		int argc, const char **argv, 
		INT verbose_level);
	void prepare_generator(INT verbose_level);
	void compute_starter(
		INT f_recognize, INT *set, INT sz, 
		INT verbose_level);


	void early_test_func(INT *S, INT len, 
		INT *candidates, INT nb_candidates, 
		INT *good_candidates, INT &nb_good_candidates, 
		INT verbose_level);
	INT check_arc(INT *S, INT len, INT verbose_level);
	void print_set_in_affine_plane(INT len, INT *S);
	void point_unrank(INT *v, INT rk);
	INT point_rank(INT *v);
	void lifting_prepare_function_new(exact_cover *E, INT starter_case, 
		INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
		diophant *&Dio, INT *&col_labels, 
		INT &f_ruled_out, 
		INT verbose_level);
	// compute the incidence matrix of tangent lines versus candidate points
	// extended by external lines versus candidate points
	INT arc_test(INT *S, INT len, INT verbose_level);
	void report(isomorph &Iso, INT verbose_level);
	void report_decompositions(isomorph &Iso, ofstream &f, INT orbit, 
		INT *data, INT verbose_level);
	void report_stabilizer(isomorph &Iso, ofstream &f, INT orbit, INT verbose_level);
};


INT callback_arc_test(exact_cover *EC, INT *S, INT len, void *data, INT verbose_level);
INT check_arc(INT len, INT *S, void *data, INT verbose_level);
INT placebo_test_function(INT len, INT *S, void *data, INT verbose_level);
void arc_generator_early_test_function(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level);
void placebo_early_test_function(INT *S, INT len, 
	INT *candidates, INT nb_candidates, 
	INT *good_candidates, INT &nb_good_candidates, 
	void *data, INT verbose_level);
void arc_generator_lifting_prepare_function_new(exact_cover *EC, INT starter_case, 
	INT *candidates, INT nb_candidates, strong_generators *Strong_gens, 
	diophant *&Dio, INT *&col_labels, 
	INT &f_ruled_out, 
	INT verbose_level);
void print_arc(INT len, INT *S, void *data);
void print_point(INT pt, void *data);
void callback_arc_report(isomorph *Iso, void *data, INT verbose_level);



