// elliptic_curve.C
// 
// Anton Betten
// Oct 27, 2009
//
//
// pulled out of crypto.C: November 19, 2014
//
//

#include "orbiter.h"

static void get_ab(INT q, INT x1, INT x2, INT x3, INT &a, INT &b);

elliptic_curve::elliptic_curve()
{
	null();
}

elliptic_curve::~elliptic_curve()
{
	freeself();
}


void elliptic_curve::null()
{
	T = NULL;
	A = NULL;
}

void elliptic_curve::freeself()
{
	if (T) {
		FREE_INT(T);
		}
	if (A) {
		FREE_INT(A);
		}
}


void elliptic_curve::init(finite_field *F, INT b, INT c, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_v) {
		cout << "elliptic_curve::init q=" << F->q << " b=" << b << " c=" << c << endl;
		}
	elliptic_curve::F = F;
	q = F->q;
	p = F->p;
	e = F->e;
	elliptic_curve::b = b;
	elliptic_curve::c = c;


	if (f_v) {
		cout << "elliptic_curve::init before compute_points" << endl;
		}

	compute_points(verbose_level);

	if (f_v) {
		cout << "elliptic_curve::init after compute_points" << endl;
		}

#if 0
	if (E.nb < 20) {
		print_integer_matrix_width(cout, E.A, E.nb, E.nb, E.nb, 3);
		}
	
#endif


	
#if 0
	cout << "point : order" << endl;
	for (i = 0; i < E.nb; i++) {
		j = order_of_point(E, i);
		cout << setw(4) << i << " : " << setw(4) << j << endl;
		}
	
	cout << "the curve has " << E.nb << " points" << endl;
#endif

#if 0
	{
	INT a, b, c;
	a = 1;
	b = 2;
	c = E.A[a * E.nb + b];
	cout << "P_" << a << " + P_" << b << " = P_" << c << endl;
	}

	j = multiple_of_point(E, 0, 37);
	cout << "37 * P_0 = P_" << j << endl;
#endif
	if (f_v) {
		cout << "elliptic_curve::init done" << endl;
		}
}


void elliptic_curve::compute_points(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT x, y;
	INT r, l;
	INT bound;

	if (f_v) {
		cout << "elliptic_curve::compute_points" << endl;
		}
	bound = q + 1 + 2 * ((INT)(sqrt(q)) + 1);// Hasse Weil bound
	


	T = new INT[bound * 3];
	nb = 0;


	for (x = 0; x < q; x++) {
		r = evaluate_RHS(x);
		if (r == 0) {
			T[nb * 3 + 0] = x;
			T[nb * 3 + 1] = 0;
			T[nb * 3 + 2] = 1;
			nb++;
			if (nb == bound) {
				cout << "The number of points exceeds the bound" << endl;
				exit(1);
				}
			//cout << nb++ << " : (" << x << "," << 0 << ",1)" << endl;
			}
		else {
			if (p != 2) {
				l = Legendre(r, q, 0);
					// GALOIS/number_theory.C

				if (l == 1) {
					y = sqrt_mod_involved(r, q);
						// DISCRETA/global.C

					if (F->mult(y, y) != r) {
						cout << "There is a problem with the square root" << endl;
						exit(1);
						}
					T[nb * 3 + 0] = x;
					T[nb * 3 + 1] = y;
					T[nb * 3 + 2] = 1;
					nb++;
					if (nb == bound) {
						cout << "The number of points exceeds the bound" << endl;
						exit(1);
						}
					T[nb * 3 + 0] = x;
					T[nb * 3 + 1] = F->negate(y);
					T[nb * 3 + 2] = 1;
					nb++;
					if (nb == bound) {
						cout << "The number of points exceeds the bound" << endl;
						exit(1);
						}
					//cout << nb++ << " : (" << x << "," << y << ",1)" << endl;
					//cout << nb++ << " : (" << x << "," << F.negate(y) << ",1)" << endl;
					}
				}
			else {
				y = F->frobenius_power(r, e - 1);
				T[nb * 3 + 0] = x;
				T[nb * 3 + 1] = y;
				T[nb * 3 + 2] = 1;
				nb++;
				if (nb == bound) {
					cout << "The number of points exceeds the bound" << endl;
					exit(1);
					}
				//cout << nb++ << " : (" << x << "," << y << ",1)" << endl;
				}
			}
		}

	// the point at infinity comes last:
	T[nb * 3 + 0] = 0;
	T[nb * 3 + 1] = 1;
	T[nb * 3 + 2] = 0;
	nb++;

	if (nb == bound) {
		cout << "The number of points exceeds the bound" << endl;
		exit(1);
		}


	if (f_v) {
		cout << "elliptic_curve::compute_points done, we found " << nb << " points" << endl;
		}
}

INT elliptic_curve::evaluate_RHS(INT x)
// evaluates x^3 + bx + c
{
	INT x2, x3, t;
	
	x2 = F->mult(x, x);
	x3 = F->mult(x2, x);
	t = F->add(x3, F->mult(b, x));
	t = F->add(t, c);
	return t;
}

void elliptic_curve::print_points()
{
	INT i;
	
	cout << "i : point (x,y,x)" << endl;
	for (i = 0; i < nb; i++) {
		cout << setw(4) << i << " : " << T[i * 3 + 0] << "," << T[i * 3 + 1] << "," << T[i * 3 + 2] << endl;
		}
}


void elliptic_curve::addition(
	INT x1, INT x2, INT x3, 
	INT y1, INT y2, INT y3, 
	INT &z1, INT &z2, INT &z3, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT a, two, three, top, bottom, m;
	
	if (f_v) {
		cout << "elliptic_curve::addition: ";
		cout << "(" << x1 << "," << x2 << "," << x3 << ")";
		cout << " + ";
		cout << "(" << y1 << "," << y2 << "," << y3 << ")";
		cout << endl;
		}
	if (x3 == 0) {
		z1 = y1;
		z2 = y2;
		z3 = y3;
		return;
		}
	if (y3 == 0) {
		z1 = x1;
		z2 = x2;
		z3 = x3;
		return;
		}
	if (x3 != 1) {
		a = F->inverse(x3);
		x1 = F->mult(x1, a);
		x2 = F->mult(x2, a);
		}
	if (y3 != 1) {
		a = F->inverse(y3);
		y1 = F->mult(y1, a);
		y2 = F->mult(y2, a);
		}
	if (x1 == y1 && x2 != y2) {
		if (F->negate(x2) != y2) {
			cout << "x1 == y1 && x2 != y2 && F.negate(x2) != y2" << endl;
			exit(1);
			}
		z1 = 0;
		z2 = 1;
		z3 = 0;
		return;
		}
	if (x1 == y1 && x2 == 0 && y2 == 0) {
		z1 = 0;
		z2 = 1;
		z3 = 0;
		return;
		}
	if (x1 == y1 && x2 == y2) {
		two = F->add(1, 1);
		three = F->add(two, 1);
		top = F->add(F->mult(three, F->mult(x1, x1)), b);
		bottom = F->mult(two, x2);
		a = F->inverse(bottom);  // this does not work in characteristic two !!!
		m = F->mult(top, a);
		}
	else {
		top = F->add(y2, F->negate(x2));
		bottom = F->add(y1, F->negate(x1));
		a = F->inverse(bottom);
		m = F->mult(top, a);
		}
	z1 = F->add(F->add(F->mult(m, m), F->negate(x1)), F->negate(y1));
	z2 = F->add(F->mult(m, F->add(x1, F->negate(z1))), F->negate(x2));
	z3 = 1;
}

void elliptic_curve::draw_grid(char *fname, INT xmax, INT ymax, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT x_min = 0, x_max = 1000;
	INT y_min = 0, y_max = 1000;
	INT factor_1000 = 1000;
	BYTE fname_full[1000];
	INT f_embedded = TRUE;
	INT f_sideways = FALSE;
	
	if (f_v) {
		cout << "draw_grid" << endl;
		}
	sprintf(fname_full, "%s.mp", fname);
	{
	mp_graphics G(fname_full, x_min, y_min, x_max, y_max, f_embedded, f_sideways);
	G.out_xmin() = 0;
	G.out_ymin() = 0;
	G.out_xmax() = xmax;
	G.out_ymax() = ymax;
	cout << "xmax/ymax = " << xmax << " / " << ymax << endl;
	
	G.header();
	G.begin_figure(factor_1000);
	
	draw_grid_(G, verbose_level);


	G.draw_boxes_final();
	G.end_figure();
	G.footer();
	}
	cout << "written file " << fname_full << " of size " << file_size(fname_full) << endl;
	if (f_v) {
		cout << "draw_grid done" << endl;
		}
	
}

#define Xcoord(x) (2 + 2 * (x))
#define Ycoord(y) (2 + 2 * (y))

void elliptic_curve::draw_grid_(mp_graphics &G, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT x, y, Q, a, b, c, d;
	INT *Px, *Py;
	INT u;
	INT x1, x2, x3;
	//INT y1, y2, y3;
	INT rad = 20;
	INT i;
	
	if (f_v) {
		cout << "draw_grid_" << endl;
		}
	u = 500 / q;
	if (q == 4) {
		u = 400 / q;
		}
	
	Q = 2 * (q + 3) + 1;
	if (f_v) {
		cout << "Q=" << Q << endl;
		}
	Px = NEW_INT(Q * Q);
	Py = NEW_INT(Q * Q);
	
	for (x = 0; x < Q; x++) {
		for (y = 0; y < Q; y++) {
			Px[x * Q + y] = x * u;
			Py[x * Q + y] = y * u;
			}
		}
	
		


	if (f_v) {
		cout << "drawing grid" << endl;
		}
	//G.polygon2(Px, Py, qq, n - 1);
	for (x = 0; x < q; x++) {
		a = Xcoord(x);
		b = Ycoord(0);
		c = Xcoord(x);
		d = Ycoord(q - 1);
		//cout << a << "," << b << "," << c << "," << d << endl;
		G.polygon2(Px, Py, a * Q + b, c * Q + d);
		}
	for (y = 0; y < q; y++) {
		a = Xcoord(0);
		b = Ycoord(y);
		c = Xcoord(q - 1);
		d = Ycoord(y);
		//cout << a << "," << b << "," << c << "," << d << endl;
		G.polygon2(Px, Py, a * Q + b, c * Q + d);
		}
	if (f_v) {
		cout << "drawing text" << endl;
		}
	for (x = 0; x < q; x++) {
		BYTE str[1000];
		sprintf(str, "%ld", x);
		a = Xcoord(x);
		b = Ycoord(-1);
		G.aligned_text(Px[a * Q + b], Py[a * Q + b], "t", str);
		}
	for (y = 0; y < q; y++) {
		BYTE str[1000];
		sprintf(str, "%ld", y);
		a = Xcoord(-1);
		b = Ycoord(y);
		G.aligned_text(Px[a * Q + b], Py[a * Q + b], "r", str);
		}

	if (f_v) {
		cout << "drawing points" << endl;
		}
	for (i = 0; i < nb; i++) {
		x1 = T[3 * i + 0];
		x2 = T[3 * i + 1];
		x3 = T[3 * i + 2];
		get_ab(q, x1, x2, x3, a, b);
		G.nice_circle(Px[a * Q + b], Py[a * Q + b], rad);
		}

	if (f_v) {
		cout << "drawing point labels" << endl;
		}
	// drawing point labels:
	for (i = 0; i < nb; i++) {
		BYTE str[1000];
		sprintf(str, "%ld", i);
		x1 = T[3 * i + 0];
		x2 = T[3 * i + 1];
		x3 = T[3 * i + 2];
		get_ab(q, x1, x2, x3, a, b);
		G.aligned_text(Px[a * Q + b], Py[a * Q + b], "", str);
		}



#if 0
	if (start_idx < 0) {
		goto done;
		}
	G.sl_ends(0 /* line_beg_style */, 1 /* line_end_style*/);
	
	INT ord;

	if (f_v) {
		cout << "drawing multiples of point" << endl;
		}
	i = start_idx;
	ord = 1;
	while (TRUE) {
		x1 = E->T[3 * i + 0];
		x2 = E->T[3 * i + 1];
		x3 = E->T[3 * i + 2];
		j = E->A[i * E->nb + start_idx];
		ord++;
		y1 = E->T[3 * j + 0];
		y2 = E->T[3 * j + 1];
		y3 = E->T[3 * j + 2];
		get_ab(q, x1, x2, x3, a, b);
		get_ab(q, y1, y2, y3, c, d);
		G.polygon2(Px, Py, a * Q + b, c * Q + d);
		if (j == E->nb - 1) {
			cout << "point P_" << start_idx << " has order " << ord << endl;
			break;
			}
		i = j;
		}
#endif



#if 0
	G.sl_udsty(100);
	a = Xcoord(-1);
	b = Ycoord(-1);
	c = Xcoord(q + 1);
	d = Ycoord(q + 1);
	cout << a << "," << b << "," << c << "," << d << endl;
	G.polygon2(Px, Py, a * Q + b, c * Q + d);
	a = Xcoord(q + 1);
	b = Ycoord(-1);
	c = Xcoord(-1);
	d = Ycoord(q + 1);
	cout << a << "," << b << "," << c << "," << d << endl;
	G.polygon2(Px, Py, a * Q + b, c * Q + d);


	q2 = q >> 1;
	if (ODD(q))
		r = 1;
	else 
		r = 0;
	
	a = Xcoord(q2) + r;
	b = Ycoord(-1);
	c = Xcoord(q2) + r;
	d = Ycoord(q + 1);
	cout << a << "," << b << "," << c << "," << d << endl;
	G.polygon2(Px, Py, a * Q + b, c * Q + d);

	a = Xcoord(-1);
	b = Ycoord(q2) + r;
	c = Xcoord(q + 1);
	d = Ycoord(q2) + r;
	cout << a << "," << b << "," << c << "," << d << endl;
	G.polygon2(Px, Py, a * Q + b, c * Q + d);
#endif
//done:
	FREE_INT(Px);
	FREE_INT(Py);
}

static void get_ab(INT q, INT x1, INT x2, INT x3, INT &a, INT &b)
{
	if (x3 == 0) {
		a = Xcoord(q >> 1);
		b = Ycoord(q);
		}
	else {
		a = Xcoord(x1);
		b = Ycoord(x2);
		}
}




#if 0
INT order_of_point(elliptic_curve &E, INT i)
{
	INT j;
	INT ord = 1;

	j = i;
	while (j != E.nb - 1) {
		j = E.A[i * E.nb + j];
		ord++;
		}
	return ord;
}

INT multiple_of_point(elliptic_curve &E, INT i, INT n)
{
	INT j, a;

	a = E.nb - 1;
	for (j = 0; j < n; j++) {
		a = E.A[a * E.nb + i];
		}
	return a;
}
#endif

void elliptic_curve::compute_addition_table(INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_v3 = (verbose_level >= 3);
	INT i, j, k;
	INT x1, x2, x3;
	INT y1, y2, y3;
	INT z1, z2, z3;

	if (f_v) {
		cout << "elliptic_curve::compute_addition_table" << endl;
		}
	
	A = new INT[nb * nb];
	for (i = 0; i < nb; i++) {
		x1 = T[3 * i + 0];
		x2 = T[3 * i + 1];
		x3 = T[3 * i + 2];
		for (j = 0; j < nb; j++) {
			y1 = T[3 * j + 0];
			y2 = T[3 * j + 1];
			y3 = T[3 * j + 2];
			if (f_v3) {
				cout << "add " << i << " " << j << endl;
				}
			addition(
				x1, x2, x3, 
				y1, y2, y3,
				z1, z2, z3, verbose_level - 1);
			k = index_of_point(z1, z2, z3);
			A[i * nb + j] = k;
			}
		}
	if (f_v) {
		cout << "elliptic_curve::compute_addition_table done" << endl;
		}
}

INT elliptic_curve::index_of_point(INT x1, INT x2, INT x3)
{
	INT a, i;
	
	if (x3 == 0) {
		return nb - 1;
		}
	if (x3 != 1) {
		a = F->inverse(x3);
		x1 = F->mult(x1, a);
		x2 = F->mult(x2, a);
		x3 = 1;
		}
	for (i = 0; i < nb - 1; i++) {
		if (T[3 * i + 0] == x1 && T[3 * i + 1] == x2) {
			return i;
			}
		}
	cout << "did not find point " << x1 << "," << x2 << "," << x3 << " in table" << endl;
	exit(1);
}


