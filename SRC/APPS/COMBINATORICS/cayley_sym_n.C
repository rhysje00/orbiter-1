// cayley_sym_n.C
//
// Anton Betten
// February 16, 2015
//

#include "orbiter.h"


INT t0;

void do_it(INT n, INT f_special, INT verbose_level);



int main(int argc, char **argv)
{
	INT i;
	INT verbose_level = 0;
	INT f_n = FALSE;
	INT n = 0;
	INT f_star = FALSE;
	
	t0 = os_ticks();

	
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-n") == 0) {
			f_n = TRUE;
			n = atoi(argv[++i]);
			cout << "-n " << n << endl;
			}
		else if (strcmp(argv[i], "-star") == 0) {
			f_star = TRUE;
			cout << "-star" << endl;
			}
		}
	


	if (!f_n) {
		cout << "please specify -n <n>" << endl;
		exit(1);
		}

	do_it(n, f_star, verbose_level);
	
}

void do_it(INT n, INT f_star, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "do_it" << endl;
		}

	action *A;
	longinteger_object go;
	INT goi;
	INT **Elt;
	INT *v;
	INT i, j, h;
	INT nb_gens;

	A = new action;
	A->init_symmetric_group(n, verbose_level);
	A->group_order(go);

	goi = go.as_INT();
	cout << "Created group Sym(" << n << ") of size " << goi << endl;


	nb_gens = n - 1;
	Elt = NEW_PINT(nb_gens);
	for (i = 0; i < nb_gens; i++) {
		Elt[i] = NEW_INT(A->elt_size_in_INT);
		}
	v = NEW_INT(n);


	if (f_star) {
		for (i = 0; i < nb_gens; i++) {
			for (j = 0; j < n; j++) {
				v[j] = j;
				}
			v[0] = i + 1;
			v[i + 1] = 0;
			A->make_element(Elt[i], v, 0 /* verbose_level */);
			}
		}
	else {
		for (i = 0; i < nb_gens; i++) {
			for (j = 0; j < n; j++) {
				v[j] = j;
				}
			v[i] = i + 1;
			v[i + 1] = i;
			A->make_element(Elt[i], v, 0 /* verbose_level */);
			}
		}
	cout << "generators:" << endl;
	for (i = 0; i < nb_gens; i++) {
		cout << "generator " << i << ":" << endl;
		A->element_print(Elt[i], cout);
		}

	

	sims *Sims;

	Sims = A->Sims;	

	INT *Adj;
	INT *Elt1, *Elt2;

	Elt1 = NEW_INT(A->elt_size_in_INT);
	Elt2 = NEW_INT(A->elt_size_in_INT);
	Adj = NEW_INT(goi * goi);

	INT_vec_zero(Adj, goi * goi);

	cout << "Computing the Cayley graph:" << endl;
	for (i = 0; i < goi; i++) {
		Sims->element_unrank_INT(i, Elt1);
		//cout << "i=" << i << endl;
		for (h = 0; h < nb_gens; h++) {
			A->element_mult(Elt1, Elt[h], Elt2, 0);
#if 0
			cout << "i=" << i << " h=" << h << endl;
			cout << "Elt1=" << endl;
			A->element_print_quick(Elt1, cout);
			cout << "g_h=" << endl;
			A->element_print_quick(gens->ith(h), cout);
			cout << "Elt2=" << endl;
			A->element_print_quick(Elt2, cout);
#endif
			j = Sims->element_rank_INT(Elt2);
			Adj[i * goi + j] = Adj[j * goi + i] = 1;
			if (i == 0) {
				cout << "edge " << i << " " << j << endl;
				}
			}
		}

	cout << "The adjacency matrix of a graph with " << goi << " vertices has been computed" << endl;
	//INT_matrix_print(Adj, goi, goi);


	colored_graph *CG;
	BYTE fname[1000];

	CG = new colored_graph;
	CG->init_adjacency_no_colors(goi, Adj, verbose_level);

	if (f_star) {
		sprintf(fname, "Cayley_Sym_%ld_star.colored_graph", n);
		}
	else {
		sprintf(fname, "Cayley_Sym_%ld_coxeter.colored_graph", n);
		}

	CG->save(fname, verbose_level);

	cout << "Written file " << fname << " of size " << file_size(fname) << endl;


	delete CG;
	FREE_INT(Adj);
	FREE_INT(Elt1);
	FREE_INT(Elt2);
	FREE_INT(v);
	for (i = 0; i < n - 1; i++) {
		FREE_INT(Elt[i]);
		}
	FREE_PINT(Elt);
		
}

