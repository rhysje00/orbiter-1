// projective.C
//
// Anton Betten
//
// started:  April 2, 2003




#include "galois.h"

INT nb_PG_elements(INT n, INT q)
// $\frac{q^{n+1} - 1}{q-1} = \sum_{i=0}^{n} q^i $
{
	INT qhl, l, deg;
	
	l = 0;
	qhl = 1;
	deg = 0;
	while (l <= n) {
		deg += qhl;
		qhl *= q;
		l++;
		}	
	return deg;
}

INT nb_PG_elements_not_in_subspace(INT n, INT m, INT q)
// |PG_n(q)| - |PG_m(q)|
{
	INT a, b;
	
	a = nb_PG_elements(n, q);
	b = nb_PG_elements(m, q);
	return a - b;
}

INT nb_AG_elements(INT n, INT q)
// $q^n$
{
	return i_power_j(q, n);
}

void all_PG_elements_in_subspace(finite_field *F, INT *genma, INT k, INT n, INT *&point_list, INT &nb_points, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = FALSE; //(verbose_level >= 2);
	INT *message;
	INT *word;
	INT i, j;

	if (f_v) {
		cout << "all_PG_elements_in_subspace" << endl;
		}
	message = NEW_INT(k);
	word = NEW_INT(n);
	nb_points = generalized_binomial(k, 1, F->q);
	point_list = NEW_INT(nb_points);
	
	for (i = 0; i < nb_points; i++) {
		PG_element_unrank_modified(*F, message, 1, k, i);
		if (f_vv) {
			cout << "message " << i << " / " << nb_points << " is ";
			INT_vec_print(cout, message, k);
			cout << endl;
			}
		F->mult_vector_from_the_left(message, genma, word, k, n);
		if (f_vv) {
			cout << "yields word ";
			INT_vec_print(cout, word, n);
			cout << endl;
			}
		PG_element_rank_modified(*F, word, 1, n, j);
		if (f_vv) {
			cout << "which has rank " << j << endl;
			}
		point_list[i] = j;
		}
	if (f_v) {
		cout << "before FREE_INT(message);" << endl;
		}

	FREE_INT(message);
	FREE_INT(word);
	if (f_v) {
		cout << "all_PG_elements_in_subspace done" << endl;
		}
}

void all_PG_elements_in_subspace_array_is_given(finite_field *F, INT *genma, INT k, INT n, INT *point_list, INT &nb_points, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = FALSE; //(verbose_level >= 2);
	INT *message;
	INT *word;
	INT i, j;

	if (f_v) {
		cout << "all_PG_elements_in_subspace_array_is_given" << endl;
		}
	message = NEW_INT(k);
	word = NEW_INT(n);
	nb_points = generalized_binomial(k, 1, F->q);
	//point_list = NEW_INT(nb_points);
	
	for (i = 0; i < nb_points; i++) {
		PG_element_unrank_modified(*F, message, 1, k, i);
		if (f_vv) {
			cout << "message " << i << " / " << nb_points << " is ";
			INT_vec_print(cout, message, k);
			cout << endl;
			}
		F->mult_vector_from_the_left(message, genma, word, k, n);
		if (f_vv) {
			cout << "yields word ";
			INT_vec_print(cout, word, n);
			cout << endl;
			}
		PG_element_rank_modified(*F, word, 1, n, j);
		if (f_vv) {
			cout << "which has rank " << j << endl;
			}
		point_list[i] = j;
		}
	if (f_v) {
		cout << "before FREE_INT(message);" << endl;
		}

	FREE_INT(message);
	FREE_INT(word);
	if (f_v) {
		cout << "all_PG_elements_in_subspace_array_is_given done" << endl;
		}
}

void display_all_PG_elements(INT n, finite_field &GFq)
{
	INT *v = NEW_INT(n + 1);
	INT l = nb_PG_elements(n, GFq.q);
	INT i, j, a;
	
	for (i = 0; i < l; i++) {
		PG_element_unrank_modified(GFq, v, 1, n + 1, i);
		cout << i << " : ";
		for (j = 0; j < n + 1; j++) {
			cout << v[j] << " ";
			}
		PG_element_rank_modified(GFq, v, 1, n + 1, a);
		cout << " : " << a << endl;
		}
	FREE_INT(v);
}

void display_all_PG_elements_not_in_subspace(INT n, INT m, finite_field &GFq)
{
	INT *v = NEW_INT(n + 1);
	INT l = nb_PG_elements_not_in_subspace(n, m, GFq.q);
	INT i, j, a;
	
	for (i = 0; i < l; i++) {
		PG_element_unrank_modified_not_in_subspace(GFq, v, 1, n + 1, m, i);
		cout << i << " : ";
		for (j = 0; j < n + 1; j++) {
			cout << v[j] << " ";
			}
		PG_element_rank_modified_not_in_subspace(GFq, v, 1, n + 1, m, a);
		cout << " : " << a << endl;
		}
	FREE_INT(v);
}

void display_all_AG_elements(INT n, finite_field &GFq)
{
	INT *v = NEW_INT(n);
	INT l = nb_AG_elements(n, GFq.q);
	INT i, j;
	
	for (i = 0; i < l; i++) {
		AG_element_unrank(GFq.q, v, 1, n + 1, i);
		cout << i << " : ";
		for (j = 0; j < n; j++) {
			cout << v[j] << " ";
			}
		cout << endl;
		}
	FREE_INT(v);
}

void PG_element_apply_frobenius(INT n, finite_field &GFq, INT *v, INT f)
{
	INT i;
	
	for (i = 0; i < n; i++) {
		v[i] = GFq.frobenius_power(v[i], f);
		}
}

void PG_element_normalize(finite_field &GFq, INT *v, INT stride, INT len)
// last non-zero element made one
{
	INT i, j, a;
	
	for (i = len - 1; i >= 0; i--) {
		a = v[i * stride];
		if (a) {
			if (a == 1)
				return;
			a = GFq.inverse(a);
			v[i * stride] = 1;
			for (j = i - 1; j >= 0; j--) {
				v[j * stride] = GFq.mult(v[j * stride], a);
				}
			return;
			}
		}
	cout << "PG_element_normalize() zero vector()" << endl;
	exit(1);
}

void PG_element_normalize_from_front(finite_field &GFq, INT *v, INT stride, INT len)
// first non zero element made one
{
	INT i, j, a;
	
	for (i = 0; i < len; i++) {
		a = v[i * stride];
		if (a) {
			if (a == 1)
				return;
			a = GFq.inverse(a);
			v[i * stride] = 1;
			for (j = i + 1; j < len; j++) {
				v[j * stride] = GFq.mult(v[j * stride], a);
				}
			return;
			}
		}
	cout << "PG_element_normalize() zero vector()" << endl;
	exit(1);
}

void PG_element_rank_modified(finite_field &GFq, INT *v, INT stride, INT len, INT &a)
{
	INT i, j, q_power_j, b, sqj;
	INT f_v = FALSE;
	
	if (len <= 0) {
		cout << "PG_element_rank_modified() len <= 0" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "the vector before normalization is ";
		for (i = 0; i < len; i++) {
			cout << v[i * stride] << " ";
			}
		cout << endl;
		}
	PG_element_normalize(GFq, v, stride, len);
	if (f_v) {
		cout << "the vector after normalization is ";
		for (i = 0; i < len; i++) {
			cout << v[i * stride] << " ";
			}
		cout << endl;
		}
	for (i = 0; i < len; i++) {
		if (v[i * stride])
			break;
		}
	if (i == len) {
		cout << "PG_element_rank_modified() zero vector" << endl;
		exit(1);
		}
	for (j = i + 1; j < len; j++) {
		if (v[j * stride])
			break;
		}
	if (j == len) {
		// we have the unit vector vector e_i
		a = i;
		return;
		}
	
	// test for the all one vector: 
	if (i == 0 && v[i * stride] == 1) {
		for (j = i + 1; j < len; j++) {
			if (v[j * stride] != 1)
				break;
			}
		if (j == len) {
			a = len;
			return;
			}
		}
	
	
	for (i = len - 1; i >= 0; i--) {
		if (v[i * stride])
			break;
		}
	if (i < 0) {
		cout << "PG_element_rank_modified() zero vector" << endl;
		exit(1);
		}
	if (v[i * stride] != 1) {
		cout << "PG_element_rank_modified() vector not normalized" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "i=" << i << endl;
		}

	b = 0;
	q_power_j = 1;
	sqj = 0;
	for (j = 0; j < i; j++) {
		b += q_power_j - 1;
		sqj += q_power_j;
		q_power_j *= GFq.q;
		}
	if (f_v) {
		cout << "b=" << b << endl;
		cout << "sqj=" << sqj << endl;
		}


	a = 0;
	for (j = i - 1; j >= 0; j--) {
		a += v[j * stride];
		if (j > 0)
			a *= GFq.q;
		if (f_v) {
			cout << "j=" << j << ", a=" << a << endl;
			}
		}

	if (f_v) {
		cout << "a=" << a << endl;
		}
	
	// take care of 1111 vector being left out
	if (i == len - 1) {
		//cout << "sqj=" << sqj << endl;
		if (a >= sqj)
			a--;
		}
	
	a += b;
	a += len;
}

void PG_element_unrank_modified(finite_field &GFq, INT *v, INT stride, INT len, INT a)
{
	INT n, l, ql, sql, k, j, r, a1 = a;
	
	n = len;
	if (n <= 0) {
		cout << "PG_element_unrank_modified() len <= 0" << endl;
		exit(1);
		}
	if (a < n) {
		for (k = 0; k < n; k++) {
			if (k == a)
				v[k * stride] = 1;
			else
				v[k * stride] = 0;
			}
		return;
		}
	a -= n;
	if (a == 0) {
		for (k = 0; k < n; k++)
			v[k * stride] = 1;
		return;
		}
	a--;
	
	l = 1;
	ql = GFq.q;
	sql = 1;
	// sql = q^0 + q^1 + \cdots + q^{l-1}
	while (l < n) {
		if (a >= ql - 1) {
			a -= (ql - 1);
			sql += ql;
			ql *= GFq.q;
			l++;
			continue;
			}
		v[l * stride] = 1;
		for (k = l + 1; k < n; k++) {
			v[k * stride] = 0;
			}
		a++; // take into account that we do not want 00001
		if (l == n - 1 && a >= sql)
			a++; // take int account that the vector 11111 has already been listed
		j = 0;
		while (a != 0) {
			r = a % GFq.q;
			v[j * stride] = r;
			j++;
			a -= r;
			a /= GFq.q;
			}
		for ( ; j < l; j++)
			v[j * stride] = 0;
		return;
		}
	cout << "PG_element_unrank_modified() a too large" << endl;
	cout << "len = " << len << endl;
	cout << "a = " << a1 << endl;
	exit(1);
}

void PG_element_rank_modified_not_in_subspace(finite_field &GFq, INT *v, INT stride, INT len, INT m, INT &a)
{
	INT s, qq, i;
	
	qq = 1;
	s = qq;
	for (i = 0; i < m; i++) {
		qq *= GFq.q;
		s += qq;
		}
	s -= (m + 1);
	
	PG_element_rank_modified(GFq, v, stride, len, a);
	if (a > len + s)
		a -= s;
	a -= (m + 1);
}

void PG_element_unrank_modified_not_in_subspace(finite_field &GFq, INT *v, INT stride, INT len, INT m, INT a)
{
	INT s, qq, i;
	
	qq = 1;
	s = qq;
	for (i = 0; i < m; i++) {
		qq *= GFq.q;
		s += qq;
		}
	s -= (m + 1);
	
	a += (m + 1);
	if (a > len)
		a += s;
	
	PG_element_unrank_modified(GFq, v, stride, len, a);
}

void AG_element_rank(INT q, INT *v, INT stride, INT len, INT &a)
{
	INT i;
	
	if (len <= 0) {
		cout << "AG_element_rank() len <= 0" << endl;
		exit(1);
		}
	a = 0;
	for (i = len - 1; i >= 0; i--) {
		a += v[i * stride];
		if (i > 0)
			a *= q;
		}
}

void AG_element_unrank(INT q, INT *v, INT stride, INT len, INT a)
{
	INT i, b;
	
#if 1
	if (len <= 0) {
		cout << "AG_element_unrank() len <= 0" << endl;
		exit(1);
		}
#endif
	for (i = 0; i < len; i++) {
		b = a % q;
		v[i * stride] = b;
		a /= q;
		}
}

void AG_element_rank_longinteger(INT q, INT *v, INT stride, INT len, longinteger_object &a)
{
	longinteger_domain D;
	longinteger_object Q, a1;
	INT i;
	
	if (len <= 0) {
		cout << "AG_element_rank_longinteger() len <= 0" << endl;
		exit(1);
		}
	a.create(0);
	Q.create(q);
	for (i = len - 1; i >= 0; i--) {
		a.add_INT(v[i * stride]);
		//cout << "AG_element_rank_longinteger after add_INT " << a << endl;
		if (i > 0) {
			D.mult(a, Q, a1);
			a.swap_with(a1);
			//cout << "AG_element_rank_longinteger after mult " << a << endl;
			}
		}
}

void AG_element_unrank_longinteger(INT q, INT *v, INT stride, INT len, longinteger_object &a)
{
	INT i, r;
	longinteger_domain D;
	longinteger_object Q, a1;
	
	if (len <= 0) {
		cout << "AG_element_unrank_longinteger() len <= 0" << endl;
		exit(1);
		}
	for (i = 0; i < len; i++) {
		D.integral_division_by_INT(a, q, a1, r);
		//r = a % q;
		v[i * stride] = r;
		//a /= q;
		a.swap_with(a1);
		}
}


INT PG_element_modified_is_in_subspace(INT n, INT m, INT *v)
{
	INT j;
	
	for (j = m + 1; j < n + 1; j++) {
		if (v[j])
			return FALSE;
		}
	return TRUE;
}

void PG_element_modified_not_in_subspace_perm(INT n, INT m, 
	finite_field &GFq, INT *orbit, INT *orbit_inv, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *v = NEW_INT(n + 1);
	INT l = nb_PG_elements(n, GFq.q);
	INT ll = nb_PG_elements_not_in_subspace(n, m, GFq.q);
	INT i, j1 = 0, j2 = ll, f_in, j;
	
	for (i = 0; i < l; i++) {
		PG_element_unrank_modified(GFq, v, 1, n + 1, i);
		f_in = PG_element_modified_is_in_subspace(n, m, v);
		if (f_v) {
			cout << i << " : ";
			for (j = 0; j < n + 1; j++) {
				cout << v[j] << " ";
				}
			}
		if (f_in) {
			if (f_v) {
				cout << " is in the subspace" << endl;
				}
			orbit[j2] = i;
			orbit_inv[i] = j2;
			j2++;
			}
		else {
			if (f_v) {
				cout << " is not in the subspace" << endl;
				}
			orbit[j1] = i;
			orbit_inv[i] = j1;
			j1++;
			}
		}
	if (j1 != ll) {
		cout << "j1 != ll" << endl;
		exit(1);
		}
	if (j2 != l) {
		cout << "j2 != l" << endl;
		exit(1);
		}
	FREE_INT(v);
}

INT PG2_line_on_point_unrank(finite_field &GFq, INT *v1, INT rk)
{
	INT v2[3];
	
	PG2_line_on_point_unrank_second_point(GFq, v1, v2, rk);
	return PG2_line_rank(GFq, v1, v2, 1);
}

void PG2_line_on_point_unrank_second_point(finite_field &GFq, INT *v1, INT *v2, INT rk)
{
	INT V[2];
	INT idx0, idx1, idx2;
	
	PG_element_normalize(GFq, v1, 1/* stride */, 3);
	if (v1[2] == 1) {
		idx0 = 2;
		idx1 = 0;
		idx2 = 1;
		}
	else if (v1[1] == 1) {
		idx0 = 1;
		idx1 = 0;
		idx2 = 2;
		}
	else {
		idx0 = 0;
		idx1 = 1;
		idx2 = 2;
		}
	PG_element_unrank_modified(GFq, V, 1/* stride */, 2, rk);
	v2[idx0] = v1[idx0];
	v2[idx1] = GFq.add(v1[idx1], V[0]);
	v2[idx2] = GFq.add(v1[idx2], V[1]);
}

INT PG2_line_rank(finite_field &GFq, INT *v1, INT *v2, INT stride)
{
	INT A[9];
	INT base_cols[3];
	INT kernel_m, kernel_n;
	INT kernel[9];
	INT rk, line_rk;
	
	A[0] = v1[0];
	A[1] = v1[1];
	A[2] = v1[2];
	A[3] = v2[0];
	A[4] = v2[1];
	A[5] = v2[2];
	rk = GFq.Gauss_INT(A, FALSE /* f_special */, TRUE /* f_complete */, base_cols, 
		FALSE /* f_P */, NULL /*P*/, 2, 3, 3, 0 /* verbose_level */);
	if (rk != 2) {
		cout << "PG2_line_rank rk != 2" << endl;
		exit(1);
		}
	GFq.matrix_get_kernel(A, 2, 3, base_cols, rk /*nb_base_cols*/, 
		kernel_m, kernel_n, kernel);
	if (kernel_m != 3) {
		cout << "PG2_line_rank kernel_m != 3" << endl;
		exit(1);
		}
	if (kernel_n != 1) {
		cout << "PG2_line_rank kernel_n != 1" << endl;
		exit(1);
		}
	PG_element_rank_modified(GFq, kernel, 1 /* stride*/, 3, line_rk);
	return line_rk;
}

void PG2_line_unrank(finite_field &GFq, INT *v1, INT *v2, INT stride, INT line_rk)
{
	INT A[9];
	INT base_cols[3];
	INT kernel_m, kernel_n;
	INT kernel[9];
	INT rk;
	
	PG_element_unrank_modified(GFq, A, 1 /* stride*/, 3, line_rk);
	rk = GFq.Gauss_INT(A, FALSE /* f_special */, TRUE /* f_complete */, 
		base_cols, 
		FALSE /* f_P */, NULL /*P*/, 1, 3, 3, 0 /* verbose_level */);
	if (rk != 1) {
		cout << "PG2_line_unrank rk != 1" << endl;
		exit(1);
		}
	GFq.matrix_get_kernel(A, 1, 3, base_cols, rk /*nb_base_cols*/, 
		kernel_m, kernel_n, kernel);
	if (kernel_m != 3) {
		cout << "PG2_line_rank kernel_m != 3" << endl;
		exit(1);
		}
	if (kernel_n != 2) {
		cout << "PG2_line_rank kernel_n != 2" << endl;
		exit(1);
		}
	v1[0] = kernel[0];
	v2[0] = kernel[1];
	v1[1] = kernel[2];
	v2[1] = kernel[3];
	v1[2] = kernel[4];
	v2[2] = kernel[5];
}

void diagonal_orbit_perm(INT n, finite_field &GFq, INT *orbit, INT *orbit_inv, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *v = NEW_INT(n + 1);
	INT l = nb_PG_elements(n - 1, GFq.q);
	INT ll = nb_AG_elements(n - 1, GFq.q - 1);
	INT a, b, c;
	INT i, j;
	
	//cout << "l = " << l << endl;
	for (i = 0; i < l; i++) {
		orbit[i] = i;
		orbit_inv[i] = i;
		}
	for (i = 0; i < ll; i++) {
		v[0] = 1;
		AG_element_unrank(GFq.q - 1, v + 1, 1, n - 1, i);
		for (j = 1; j < n; j++) {
			v[j]++;
			}
		if (f_v) {
			cout << i << " : ";
			for (j = 0; j < n; j++) {
				cout << v[j] << " ";
				}
			}
		PG_element_rank_modified(GFq, v, 1, n, a);
		if (f_v) {
			cout << " : " << a << endl;
			}
		b = orbit_inv[a];
		c = orbit[i];
		orbit[i] = a;
		orbit[b] = c;
		orbit_inv[a] = i;
		orbit_inv[c] = b;
		}
	FREE_INT(v);
}

void frobenius_orbit_perm(INT n, finite_field &GFq, INT *orbit, INT *orbit_inv, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT *v = NEW_INT(n);
	INT l = nb_PG_elements(n - 1, GFq.q);
	INT ll = GFq.e;
	INT a, b, c;
	INT i, j;
	
	if (f_v) {
		cout << "frobenius_orbit_perm n=" << n << " (vector space dimension)" << endl;
		cout << "l=" << l << endl;
		}
	if (GFq.e == 1) {
		cout << "frobenius_orbit_perm GFq.e == 1" << endl;
		exit(1);
		}
	//cout << "l = " << l << endl;
	for (i = 0; i < l; i++) {
		orbit[i] = i;
		orbit_inv[i] = i;
		}
	if (f_v) {
		cout << "before PG_element_unrank_modified(" << n + GFq.p << ")" << endl;
		}
	PG_element_unrank_modified(GFq, v, 1, n, n + GFq.p);
	if (f_v) {
		cout << "after PG_element_unrank_modified(" << n + GFq.p << ")" << endl;
		}
	for (i = 0; i < ll; i++) {
		if (f_v) {
			cout << i << " : ";
			for (j = 0; j < n; j++) {
				cout << v[j] << " ";
				}
			}
		PG_element_rank_modified(GFq, v, 1, n, a);
		if (f_v) {
			cout << " : " << a << endl;
			}
		b = orbit_inv[a];
		c = orbit[i];
		orbit[i] = a;
		orbit[b] = c;
		orbit_inv[a] = i;
		orbit_inv[c] = b;
		PG_element_apply_frobenius(n, GFq, v, 1);
		}
	FREE_INT(v);
}

void test_PG(INT n, INT q)
{
	finite_field GFq;
	INT m;
	INT verbose_level = 1;
	
	GFq.init(q, verbose_level);
	
	cout << "all elements of PG_" << n << "(" << q << ")" << endl;
	display_all_PG_elements(n, GFq);
	
	for (m = 0; m < n; m++) {
		cout << "all elements of PG_" << n << "(" << q << "), not in a subspace of dimension " << m << endl;
		display_all_PG_elements_not_in_subspace(n, m, GFq);
		}
	
}


void line_through_two_points(finite_field &GFq, INT len, INT pt1, INT pt2, INT *line)
{
	INT v1[100], v2[100], v3[100], alpha, a, ii;
	
	if (len > 100) {
		cout << "line_through_two_points() len >= 100" << endl;
		exit(1);
		}
	PG_element_unrank_modified(GFq, v1, 1 /* stride */, len, pt1);
	PG_element_unrank_modified(GFq, v2, 1 /* stride */, len, pt2);
	line[0] = pt1;
	line[1] = pt2;
	for (alpha = 1; alpha < GFq.q; alpha++) {
		for (ii = 0; ii < len; ii++) {
			a = GFq.mult(v1[ii], alpha);
			v3[ii] = GFq.add(a, v2[ii]);
			}
		PG_element_normalize(GFq, v3, 1 /* stride */, len);
		PG_element_rank_modified(GFq, v3, 1 /* stride */, len, line[1 + alpha]);
		}
}

void print_set_in_affine_plane(finite_field &GFq, INT len, INT *S)
{
	INT *A;
	INT i, j, x, y, v[3];
	
	
	A = NEW_INT(GFq.q * GFq.q);
	for (x = 0; x < GFq.q; x++) {
		for (y = 0; y < GFq.q; y++) {
			A[(GFq.q - 1 - y) * GFq.q + x] = 0;
			}
		}
	for (i = 0; i < len; i++) {
		PG_element_unrank_modified(GFq, v, 1 /* stride */, 3 /* len */, S[i]);
		if (v[2] != 1) {
			//cout << "my_generator::print_set_in_affine_plane() not an affine point" << endl;
			cout << "(" << v[0] << "," << v[1] << "," << v[2] << ")" << endl;
			continue;
			}
		x = v[0];
		y = v[1];
		A[(GFq.q - 1 - y) * GFq.q + x] = 1;
		}
	for (i = 0; i < GFq.q; i++) {
		for (j = 0; j < GFq.q; j++) {
			cout << A[i * GFq.q + j];
			}
		cout << endl;
		}
	FREE_INT(A);
}

INT consecutive_ones_property_in_affine_plane(ostream &ost, finite_field &GFq, INT len, INT *S)
{
	INT i, y, v[3];
	
	
	for (i = 0; i < len; i++) {
		PG_element_unrank_modified(GFq, v, 1 /* stride */, 3 /* len */, S[i]);
		if (v[2] != 1) {
			cout << "my_generator::consecutive_ones_property_in_affine_plane() not an affine point" << endl;
			ost << "(" << v[0] << "," << v[1] << "," << v[2] << ")" << endl;
			exit(1);
			}
		y = v[1];
		if (y != i)
			return FALSE;
		}
	return TRUE;
}

void oval_polynomial(finite_field &GFq, INT *S, unipoly_domain &D, unipoly_object &poly, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i, v[3], x, y;
	INT *map;
	
	if (f_v) {
		cout << "oval_polynomial" << endl;
		}
	map = NEW_INT(GFq.q);
	for (i = 0; i < GFq.q; i++) {
		PG_element_unrank_modified(GFq, v, 1 /* stride */, 3 /* len */, S[2 + i]);
		if (v[2] != 1) {
			cout << "oval_polynomial not an affine point" << endl;
			exit(1);
			}
		x = v[0];
		y = v[1];
		//cout << "map[" << i << "] = " << xx << endl;
		map[i] = x;
		}
	if (f_v) {
		cout << "the map" << endl;
		for (i = 0; i < GFq.q; i++) {
			cout << map[i] << " ";
			}
		cout << endl;
		}
	
	D.create_Dickson_polynomial(poly, map);
	
	FREE_INT(map);
	if (f_v) {
		cout << "oval_polynomial done" << endl;
		}
}


#if 0
void regular_oval_point_list(finite_field &GFq, INT &len, INT *&pts)
{
	INT i;
	
	len = GFq.q + 2;
	pts = NEW_INT(len);
	for (i = 0; i < len; i++) {
		pts[i] = regular_oval_point(GFq, i);
		}
}

void segre_oval_point_list(finite_field &GFq, INT &len, INT *&pts)
{
	INT i;
	
	if ((GFq.e % 2) == 0) {
		cout << "segre oval needs odd extension degree" << endl;
		exit(1);
		}
	len = GFq.q + 2;
	pts = NEW_INT(len);
	for (i = 0; i < len; i++) {
		pts[i] = segre_oval_point(GFq, i);
		}
}

void translation_oval_point_list(finite_field &GFq, INT &len, INT *&pts)
{
	INT i;
	
	len = GFq.q + 2;
	pts = NEW_INT(len);
	
	if (GFq.q == 32) {
		INT data[] = { 0, 1, 2, 3, 100, 133, 184, 199, 236, 270, 311, 328, 380, 414, 443, 468, 496, 531, 555, 586, 617, 655, 678, 729, 765, 789, 831, 864, 867, 918, 945, 978, 1018, 1037 };
		for (i = 0; i < len; i++) {
			pts[i] = data[i];
			}
		}
	else {
		cout << "don't have translation oval for q=" << GFq.q << endl;
		exit(1);
		}
}

void penttila_okeefe_oval_point_list(finite_field &GFq, INT &len, INT *&pts)
{
	INT i;
	
	len = GFq.q + 2;
	pts = NEW_INT(len);
	
	if (GFq.q == 32) {
		INT data[] = { 0, 1, 2, 3, 100, 131, 169, 202, 238, 283, 316, 350, 383, 400, 423, 453, 505, 520, 556, 596, 634, 661, 690, 728, 755, 779, 815, 849, 896, 925, 951, 982, 998, 1037 };
		for (i = 0; i < len; i++) {
			pts[i] = data[i];
			}
		}
	else {
		cout << "Penttila O'Keefe oval exists only for q=32, but here q=" << GFq.q << endl;
		exit(1);
		}
}

void oval_32_5_point_list(finite_field &GFq, INT &len, INT *&pts)
{
	INT i;
	
	len = GFq.q + 2;
	pts = NEW_INT(len);
	
	if (GFq.q == 32) {
		INT data[] = { 0, 1, 2, 3, 100, 131, 169, 202, 241, 281, 318, 333, 378, 415, 443, 469, 494, 536, 556, 591, 619, 656, 700, 724, 742, 787, 805, 864, 886, 914, 951, 967, 1000, 1053 };
		for (i = 0; i < len; i++) {
			pts[i] = data[i];
			}
		}
	else {
		cout << "oval 32_5 exists only for q=32, but here q=" << GFq.q << endl;
		exit(1);
		}
}

void oval_32_10_point_list(finite_field &GFq, INT &len, INT *&pts)
{
	INT i;
	
	len = GFq.q + 2;
	pts = NEW_INT(len);
	
	if (GFq.q == 32) {
		INT data[] = { 0, 1, 2, 3, 100, 131, 169, 202, 236, 274, 313, 351, 378, 412, 446, 464, 500, 526, 552, 597, 635, 669, 678, 727, 749, 787, 815, 843, 886, 903, 960, 977, 997, 1048 };
		for (i = 0; i < len; i++) {
			pts[i] = data[i];
			}
		}
	else {
		cout << "oval 32_10 exists only for q=32, but here q=" << GFq.q << endl;
		exit(1);
		}
}

INT regular_oval_point(finite_field &GFq, INT i)
{
	INT t, t2, rk, xyz[3];
	
	xyz[0] = xyz[1] = xyz[2] = 0;
	if (i <= 2) {
		return i;
		}
	t = i - 2;
	xyz[2] = 1;
	xyz[1] = t;
	t2 = GFq.power(t, 2);
	xyz[0] = t2;
	PG_element_rank_modified(GFq, xyz, 1, 3, rk);
	return rk;
}

INT segre_oval_point(finite_field &GFq, INT i)
{
	INT t, t6, rk, xyz[3];
	
	xyz[0] = xyz[1] = xyz[2] = 0;
	if (i <= 2) {
		return i;
		}
	t = i - 2;
	xyz[2] = 1;
	xyz[1] = t;
	t6 = GFq.power(t, 6);
	xyz[0] = t6;
	PG_element_rank_modified(GFq, xyz, 1, 3, rk);
	return rk;
}
#endif

INT line_intersection_with_oval(finite_field &GFq, INT *f_oval_point, INT line_rk, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT line[3], base_cols[3], K[6], points[6], rk, kernel_m, kernel_n;
	INT j, w[2], a, b, pt[3], nb = 0;
	
	if (f_v) {
		cout << "intersecting line " << line_rk << " with the oval" << endl;
		}
	PG_element_unrank_modified(GFq, line, 1, 3, line_rk);
	rk = GFq.Gauss_INT(line, FALSE /* f_special */, TRUE /* f_complete */, base_cols, 
		FALSE /* f_P */, NULL /* P */, 1 /* m */, 3 /* n */, 0 /* Pn */, 
		0 /* verbose_level */);
	if (f_vv) {
		cout << "after Gauss:" << endl;
		print_integer_matrix(cout, line, 1, 3);
		}
	GFq.matrix_get_kernel(line, 1, 3, base_cols, rk, 
		kernel_m, kernel_n, K);
	INT_matrix_transpose(K, kernel_m, kernel_n, points);
	if (f_vv) {
		cout << "points:" << endl;
		print_integer_matrix(cout, points, 2, 3);
		}
	for (j = 0; j < GFq.q + 1; j++) {
		PG_element_unrank_modified(GFq, w, 1, 2, j);
		a = GFq.mult(points[0], w[0]);
		b = GFq.mult(points[3], w[1]);
		pt[0] = GFq.add(a, b);
		a = GFq.mult(points[1], w[0]);
		b = GFq.mult(points[4], w[1]);
		pt[1] = GFq.add(a, b);
		a = GFq.mult(points[2], w[0]);
		b = GFq.mult(points[5], w[1]);
		pt[2] = GFq.add(a, b);
		PG_element_rank_modified(GFq, pt, 1, 3, rk);
		//cout << j << " : " << pt[0] << "," << pt[1] << "," << pt[2] << " : " << rk << " : ";
		if (f_oval_point[rk]) {
			//cout << "yes" << endl;
			if (f_v) {
				cout << "found oval point" << endl;
				}
			nb++;
			}
		else {
			//cout << "no" << endl;
			}
		}
	return nb;
}

INT get_base_line(finite_field &GFq, INT plane1, INT plane2, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);

	if (f_v) {
		cout << "get_base_line()" << endl;
		}
	INT planes[8], base_cols[4], rk;
	INT intersection[16], intersection2[16], lines[6], line[3], kernel_m, kernel_n, line_rk;

	AG_element_unrank(GFq.q, planes, 1, 3, plane1);
	planes[3] = 1;
	AG_element_unrank(GFq.q, planes + 4, 1, 3, plane2);
	planes[7] = 1;
	if (f_v) {
		cout << "planes:" << endl;
		print_integer_matrix(cout, planes, 2, 4);
		}
	rk = GFq.Gauss_INT(planes, FALSE /* f_special */, TRUE /* f_complete */, base_cols, 
		FALSE /* f_P */, NULL /* P */, 2 /* m */, 4 /* n */, 0 /* Pn */, 
		0 /* verbose_level */);
	if (f_v) {
		cout << "after Gauss:" << endl;
		print_integer_matrix(cout, planes, 2, 4);
		}
	GFq.matrix_get_kernel(planes, 2, 4, base_cols, rk, 
		kernel_m, kernel_n, intersection);
	INT_matrix_transpose(intersection, kernel_m, kernel_n, intersection2);
	if (f_v) {
		cout << "kernel:" << endl;
		print_integer_matrix(cout, intersection2, kernel_n, kernel_m);
		}
	lines[0] = intersection2[0];
	lines[1] = intersection2[1];
	lines[2] = intersection2[2];
	lines[3] = intersection2[4];
	lines[4] = intersection2[5];
	lines[5] = intersection2[6];
	if (f_v) {
		cout << "lines:" << endl;
		print_integer_matrix(cout, lines, 2, 3);
		}
	rk = GFq.Gauss_INT(lines, FALSE /* f_special */, TRUE /* f_complete */, base_cols, 
		FALSE /* f_P */, NULL /* P */, 2 /* m */, 3 /* n */, 0 /* Pn */, 
		0 /* verbose_level */);
	if (rk != 2) {
		cout << "rk != 2" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "after Gauss:" << endl;
		print_integer_matrix(cout, lines, 2, 3);
		}
	GFq.matrix_get_kernel(lines, 2, 3, base_cols, rk, kernel_m, kernel_n, line);
	if (f_v) {
		cout << "the line:" << endl;
		print_integer_matrix(cout, line, 1, 3);
		}
	PG_element_rank_modified(GFq, line, 1, 3, line_rk);
	if (f_v) {
		cout << "has rank " << line_rk << endl;
		}
	return line_rk;
}

void translation_in_AG(finite_field &GFq, INT n, INT i, INT a, INT *perm, INT *v, INT verbose_level)
// v[n] needs to be allocated 
// p[q^n] needs to be allocated
{
	INT f_v = (verbose_level >= 1);
	INT ii, j, l, q;
	
	q = GFq.q;
	l = nb_AG_elements(n, q);
	for (ii = 0; ii < l; ii++) {
		AG_element_unrank(q, v, 1 /* stride */, n, ii);
		// cout << "ii=" << ii << " v=" << v;
		v[i] = GFq.add(v[i], a);
		
		AG_element_rank(q, v, 1 /* stride */, n, j);
		perm[ii] = j;
		// cout << " j=" << j << endl;
		}
	if (f_v) {
		cout << "translation_in_AG() i=" << i << " a=" << a << " : ";
		perm_print(cout, perm, l);
		cout << endl;
		}
}

void frobenius_in_AG(finite_field &GFq, INT n, INT *perm, INT *v, INT verbose_level)
// v[n] needs to be allocated 
// p[q^n] needs to be allocated
{
	INT f_v = (verbose_level >= 1);
	INT i, j, l, q, p;
	
	q = GFq.q;
	p = GFq.p;
	l = nb_AG_elements(n, q);
	for (i = 0; i < l; i++) {
		AG_element_unrank(q, v, 1 /* stride */, n, i);
		for (j = 0; j < n; j++) {
			v[j] = GFq.power(v[j], p);
			}
		AG_element_rank(q, v, 1 /* stride */, n, j);
		perm[i] = j;
		}
	if (f_v) {
		cout << "frobenius_in_AG() : ";
		perm_print(cout, perm, l);
		cout << endl;
		}
}

void frobenius_in_PG(finite_field &GFq, INT n, INT *perm, INT *v, INT verbose_level)
// v[n + 1] needs to be allocated 
// p[q^n+...+q+1] needs to be allocated
{
	INT f_v = (verbose_level >= 1);
	INT i, j, l, q, p;
	
	q = GFq.q;
	p = GFq.p;
	l = nb_PG_elements(n, q);
	for (i = 0; i < l; i++) {
		PG_element_unrank_modified(GFq, v, 1 /* stride */, n + 1, i);
		for (j = 0; j <= n; j++) {
			v[j] = GFq.power(v[j], p);
			}
		PG_element_unrank_modified(GFq, v, 1 /* stride */, n + 1, j);
		perm[i] = j;
		}
	if (f_v) {
		cout << "frobenius_in_PG() : ";
		perm_print(cout, perm, l);
		cout << endl;
		}
}

void AG_representation_of_matrix(finite_field &GFq, INT n, INT f_from_the_right, 
	INT *M, INT *v, INT *w, INT *perm, INT verbose_level)
// perm[q^n] needs to be already allocated
{
	INT f_v = (verbose_level >= 1);
	INT i, j, l, q;
	
	q = GFq.q;
	l = nb_AG_elements(n, q);
	for (i = 0; i < l; i++) {
		AG_element_unrank(q, v, 1 /* stride */, n, i);
		if (f_from_the_right) {
			GFq.mult_matrix_matrix(v, M, w, 1, n, n);
			}
		else {
			GFq.mult_matrix_matrix(M, v, w, n, n, 1);
			}
		AG_element_rank(q, w, 1 /* stride */, n, j);
		perm[i] = j;
		}
	if (f_v) {
		cout << "AG_representation_of_matrix() : ";
		perm_print(cout, perm, l);
		cout << endl;
		}
	
}

void AG_representation_one_dimensional(finite_field &GFq, 
	INT a, INT *perm, INT verbose_level)
// perm[q] needs to be already allocated
{
	INT f_v = (verbose_level >= 1);
	INT i, j, l, q, v, w;
	
	q = GFq.q;
	l = q;
	if (f_v) {
		cout << "AG_representation_one_dimensional() : q = " << q << " a=" << a << endl;
		}
	for (i = 0; i < q; i++) {
		AG_element_unrank(q, &v, 1 /* stride */, 1, i);
		w = GFq.mult(a, v);
		AG_element_rank(q, &w, 1 /* stride */, 1, j);
		perm[i] = j;
		}
	if (f_v) {
		cout << "AG_representation_one_dimensional() : ";
		perm_print(cout, perm, l);
		cout << endl;
		}
	
}

INT nb_generators_affine_translations(finite_field &GFq, INT n)
{
	return n * GFq.e;
}

void generators_affine_translations(finite_field &GFq, INT n, INT *perms, INT verbose_level)
// primes[n * d] needs to be allocated, where d = q^n
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *v, i, j, l, k = 0, a = 1;
	
	l = nb_AG_elements(n, GFq.q);
	
	if (f_v) {
		cout << "computing generators for affine translations, q=" << GFq.q << " n = " << n << endl;
		}
	v = NEW_INT(n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < GFq.e; j++) {
			translation_in_AG(GFq, n, i, a, perms + k * l, v, f_vv);
			k++;
			a *= GFq.p;
			}
		}
	FREE_INT(v);
}

void generators_AGL1xAGL1_subdirect1(finite_field &GFq1, finite_field &GFq2, 
	INT u, INT v, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *perms1;
	INT *perms2;
	INT nb1, nb2, q1, q2, q12, i, k = 0;
	
	q1 = GFq1.q;
	q2 = GFq2.q;
	q12 = q1 * q2;
	nb1 = nb_generators_affine_translations(GFq1, 1);
	nb2 = nb_generators_affine_translations(GFq2, 1);
	perms1 = NEW_INT((nb1 + 3) * q1);
	perms2 = NEW_INT((nb2 + 3) * q2);
	nb_perms = nb1 + nb2 + 1;
	perms = NEW_INT(nb_perms * q12);

	perm_identity(perms1, q1);
	perm_identity(perms2, q2);

	generators_affine_translations(GFq1, 1, perms1 + q1, verbose_level - 2);
	generators_affine_translations(GFq2, 1, perms2 + q2, verbose_level - 2);
	if (f_v) {
		cout << "affine translations created" << endl;
		}
	
	AG_representation_one_dimensional(GFq1, GFq1.alpha, perms1 + (nb1 + 1) * q1, 
		verbose_level - 2);
	AG_representation_one_dimensional(GFq2, GFq2.alpha, perms2 + (nb2 + 1) * q2, 
		verbose_level - 2);
	if (f_v) {
		cout << "AG_representation_one_dimensional created" << endl;
		if (f_vv) {
			perm_print(cout, perms1 + (nb1 + 1) * q1, q1); cout << endl;
			perm_print(cout, perms2 + (nb2 + 1) * q2, q2); cout << endl;
			}
		}
	
	perm_raise(perms1 + (nb1 + 1) * q1, perms1 + (nb1 + 2) * q1, u, q1);
	perm_raise(perms2 + (nb2 + 1) * q2, perms2 + (nb2 + 2) * q2, v, q2);
	if (f_v) {
		cout << "raised to the powers u and v" << endl;
		if (f_vv) {
			perm_print(cout, perms1 + (nb1 + 2) * q1, q1); cout << endl;
			perm_print(cout, perms2 + (nb2 + 2) * q2, q2); cout << endl;
			}
		}
	
	for (i = 0; i < nb1; i++) {
		perm_direct_product(q1, q2, perms1 + (i + 1) * q1, perms2, perms + k * q12);
		k++;
		}
	for (i = 0; i < nb2; i++) {
		perm_direct_product(q1, q2, perms1, perms2 + (i + 1) * q2, perms + k * q12);
		k++;
		}
	perm_direct_product(q1, q2, 
		perms1 + (nb1 + 2) * q1, 
		perms2 + (nb2 + 2) * q2, 
		perms + k * q12);
	k++;
	if (f_v) {
		cout << "generators for subdirect product AGL(1," << q1 << ") x AGL(1," << q2 << ") created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * q12, q12);
			cout << endl;
			}
		}
	FREE_INT(perms1);
	FREE_INT(perms2);
}

void generators_AGL1q(finite_field &GFq, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT nb, q, i;
	
	q = GFq.q;
	nb = nb_generators_affine_translations(GFq, 1);
	nb_perms = nb + 1;
	perms = NEW_INT(nb_perms * q);

	generators_affine_translations(GFq, 1, perms, verbose_level - 2);
	if (f_v) {
		cout << "affine translations created" << endl;
		}
	
	AG_representation_one_dimensional(GFq, GFq.alpha, perms + nb * q, verbose_level - 2);
	if (f_v) {
		cout << "AG_representation_one_dimensional created" << endl;
		}
	
	if (f_v) {
		cout << "generators for AGL(1," << q << ") created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * q, q);
			cout << endl;
			}
		}
}

void generators_AGL1q_subgroup(finite_field &GFq, INT index_in_multiplicative_group, 
	INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT nb, q, i, a, b;
	
	q = GFq.q;
	nb = nb_generators_affine_translations(GFq, 1);
	nb_perms = nb + 1;
	perms = NEW_INT(nb_perms * q);

	generators_affine_translations(GFq, 1, perms, verbose_level - 2);
	if (f_v) {
		cout << "affine translations created" << endl;
		}
	
	a = GFq.alpha;
	b = GFq.power(a, index_in_multiplicative_group);
	AG_representation_one_dimensional(GFq, b, perms + nb * q, verbose_level - 2);
	if (f_v) {
		cout << "AG_representation_one_dimensional created" << endl;
		}
	
	if (f_v) {
		cout << "generators for AGL(1," << q << ") created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * q, q);
			cout << endl;
			}
		}
}

void generators_AGL1_x_AGL1(finite_field &GFq1, finite_field &GFq2, INT &deg, 
	INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT deg1, nb_perms1, *perms1;
	INT deg2, nb_perms2, *perms2;
	INT i;
	
	deg1 = GFq1.q;
	deg2 = GFq2.q;
	
	generators_AGL1q(GFq1, nb_perms1, perms1, verbose_level - 1);
	generators_AGL1q(GFq2, nb_perms2, perms2, verbose_level - 1);
	
	generators_direct_product(deg1, nb_perms1, perms1, deg2, 
		nb_perms2, perms2, deg, nb_perms, perms, verbose_level - 1);
	
	FREE_INT(perms1);
	FREE_INT(perms2);
	if (f_v) {
		cout << "generators for AGL(1," << deg1 << ") x AGL(1," << deg2 << ") created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_AGL1_x_AGL1_extension(finite_field &GFq1, finite_field &GFq2, INT u, INT v, 
	INT &deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT *perms1, *perms2;
	INT q1, q2, i;
	
	q1 = GFq1.q;
	q2 = GFq2.q;
	
	perms1 = NEW_INT(2 * q1);
	perms2 = NEW_INT(2 * q2);
	
	deg = q1 * q2;
	nb_perms = 1;
	perms = NEW_INT(nb_perms * deg);

	AG_representation_one_dimensional(GFq1, GFq1.alpha, perms1, f_vv);
	AG_representation_one_dimensional(GFq2, GFq2.alpha, perms2, f_vv);
	if (f_v) {
		cout << "AG_representation_one_dimensional created" << endl;
		}
	
	perm_raise(perms1, perms1 + q1, u, q1);
	perm_raise(perms2, perms2 + q2, v, q2);
	if (f_v) {
		cout << "raised to the powers u and v" << endl;
		}
	
	perm_direct_product(q1, q2, perms1 + q1, perms2 + q2, perms);
	FREE_INT(perms1);
	FREE_INT(perms2);
	
	if (f_v) {
		cout << "generators for a^" << u << "b^" << v << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_AGL1_x_AGL1_extended_once(finite_field &F1, finite_field &F2, INT u, INT v, 
	INT &deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT deg1, nb_perms1, *perms1;
	INT deg2, nb_perms2, *perms2;
	INT q1, q2, i;
	
	q1 = F1.q;
	q2 = F2.q;
	generators_AGL1_x_AGL1(F1, F2, deg1, nb_perms1, perms1, verbose_level - 1);
	generators_AGL1_x_AGL1_extension(F1, F2, u, v, deg2, nb_perms2, perms2, verbose_level - 1);
	
	generators_concatenate(deg1, nb_perms1, perms1, 
		deg2, nb_perms2, perms2, 
		deg, nb_perms, perms, verbose_level - 1);
	
	FREE_INT(perms1);
	FREE_INT(perms2);
	
	if (f_v) {
		cout << "generators for AGL(1," << q1 << ") x AGL(1," << q2 << ") extended by a^" << u << "b^" << v << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_AGL1_x_AGL1_extended_twice(finite_field &F1, finite_field &F2, INT u1, INT v1, INT u2, INT v2, INT &deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT deg1, nb_perms1, *perms1;
	INT deg2, nb_perms2, *perms2;
	INT q1, q2, i;
	
	q1 = F1.q;
	q2 = F2.q;
	generators_AGL1_x_AGL1_extended_once(F1, F2, u1, v1, deg1, nb_perms1, perms1, verbose_level - 1);
	generators_AGL1_x_AGL1_extension(F1, F2, u2, v2, deg2, nb_perms2, perms2, verbose_level - 1);
	
	generators_concatenate(deg1, nb_perms1, perms1, deg2, nb_perms2, perms2, deg, nb_perms, perms, verbose_level - 1);
	
	FREE_INT(perms1);
	FREE_INT(perms2);
	
	if (f_v) {
		cout << "generators for AGL(1," << q1 << ") x AGL(1," << q2 << ") extended by a^" << u1 << "b^" << v1 << " and by a^" << u2 << "b^" << v2 << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_symmetric_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i;
	
	if (deg <= 1) {
		nb_perms = 0;
		perms = NULL;
		return;
		}
	nb_perms = deg - 1;
	perms = NEW_INT(nb_perms * deg);
	for (i = 0; i < nb_perms; i++) {
		perm_identity(perms + i * deg, deg);
		perms[i * deg + i] = i + 1;
		perms[i * deg + i + 1] = i;
		}
	if (f_v) {
		cout << "generators for symmetric group of degree " << deg << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_cyclic_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i = 0, j;
	
	if (deg <= 1) {
		nb_perms = 0;
		perms = NULL;
		return;
		}
	nb_perms = 1;
	perms = NEW_INT(nb_perms * deg);
	for (j = 0; j < deg; j++) {
		perms[i * deg + j] = j + 1;
		}
	perms[i * deg + i + deg - 1] = 0;
	if (f_v) {
		cout << "generators for cyclic group of degree " << deg << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_dihedral_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i = 0, j, d2;
	
	if (deg <= 1) {
		nb_perms = 0;
		perms = NULL;
		return;
		}
	d2 = deg >> 1;
	nb_perms = 2;
	perms = NEW_INT(nb_perms * deg);
	for (j = 0; j < deg; j++) {
		perms[i * deg + j] = j + 1;
		}
	perms[i * deg + i + deg - 1] = 0;
	i++;
	for (j = 0; j <= d2; j++) {
		perms[i * deg + j] = deg - 1 - j;
		perms[i * deg + deg - 1 - j] = j;
		}
	if (f_v) {
		cout << "generators for dihedral group of degree " << deg << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_dihedral_involution(INT deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i = 0, j, d2;
	
	if (deg <= 1) {
		nb_perms = 0;
		perms = NULL;
		return;
		}
	d2 = deg >> 1;
	nb_perms = 1;
	perms = NEW_INT(nb_perms * deg);
	i = 0;
	for (j = 0; j <= d2; j++) {
		perms[i * deg + j] = deg - 1 - j;
		perms[i * deg + deg - 1 - j] = j;
		}
	if (f_v) {
		cout << "generators for dihedral involution of degree " << deg << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_identity_group(INT deg, INT &nb_perms, INT *&perms, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i = 0, j;
	
	if (deg <= 1) {
		nb_perms = 0;
		perms = NULL;
		return;
		}
	nb_perms = 1;
	perms = NEW_INT(nb_perms * deg);
	for (j = 0; j < deg; j++) {
		perms[j] = j;
		}
	if (f_v) {
		cout << "generators for identity group of degree " << deg << " created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms; i++) {
			perm_print(cout, perms + i * deg, deg);
			cout << endl;
			}
		}
}

void generators_direct_product(INT deg1, INT nb_perms1, INT *perms1, 
	INT deg2, INT nb_perms2, INT *perms2, 
	INT &deg3, INT &nb_perms3, INT *&perms3, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, k = 0;
	INT *id1, *id2;
	
	deg3 = deg1 * deg2;
	nb_perms3 = nb_perms1 + nb_perms2;
	perms3 = NEW_INT(nb_perms3 * deg3);
	id1 = NEW_INT(deg1);
	id2 = NEW_INT(deg2);
	perm_identity(id1, deg1);
	perm_identity(id2, deg2);
	
	for (i = 0; i < nb_perms1; i++) {
		perm_direct_product(deg1, deg2, perms1 + i * deg1, id2, perms3 + k * deg3);
		k++;
		}
	for (i = 0; i < nb_perms2; i++) {
		perm_direct_product(deg1, deg2, id1, perms2 + i * deg2, perms3 + k * deg3);
		k++;
		}
	FREE_INT(id1);
	FREE_INT(id2);
	if (f_v) {
		cout << "generators for direct product created" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms3; i++) {
			perm_print(cout, perms3 + i * deg3, deg3);
			cout << endl;
			}
		}
}

void generators_concatenate(INT deg1, INT nb_perms1, INT *perms1, 
	INT deg2, INT nb_perms2, INT *perms2, 
	INT &deg3, INT &nb_perms3, INT *&perms3, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i, k = 0;
	
	if (deg1 != deg2) {
		cout << "generators_concatenate: deg1 != deg2" << endl;
		exit(1);
		}
	deg3 = deg1;
	nb_perms3 = nb_perms1 + nb_perms2;
	perms3 = NEW_INT(nb_perms3 * deg3);
	
	k = 0;
	for (i = 0; i < nb_perms1; i++) {
		perm_move(perms1 + i * deg1, perms3 + k * deg3, deg3);
		k++;
		}
	for (i = 0; i < nb_perms2; i++) {
		perm_move(perms2 + i * deg1, perms3 + k * deg3, deg3);
		k++;
		}
	if (f_v) {
		cout << "generators concatenated" << endl;
		}
	if (f_vv) {
		for (i = 0; i < nb_perms3; i++) {
			perm_print(cout, perms3 + i * deg3, deg3);
			cout << endl;
			}
		}
}

void O4_isomorphism_4to2(finite_field *F, INT *At, INT *As, INT &f_switch, INT *B, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT a, b, c, d, e, f, g, h;
	INT ev, fv;
	INT P[4], Q[4], R[4], S[4];
	INT Rx, Ry, Sx, Sy;
	INT b11, b12, b13, b14;
	INT b21, b22, b23, b24;
	INT b31, b32, b33, b34;
	INT b41, b42, b43, b44;

	if (f_v) {
		cout << "O4_isomorphism_4to2" << endl;
		}
	b11 = B[0 * 4 + 0];
	b12 = B[0 * 4 + 1];
	b13 = B[0 * 4 + 2];
	b14 = B[0 * 4 + 3];
	b21 = B[1 * 4 + 0];
	b22 = B[1 * 4 + 1];
	b23 = B[1 * 4 + 2];
	b24 = B[1 * 4 + 3];
	b31 = B[2 * 4 + 0];
	b32 = B[2 * 4 + 1];
	b33 = B[2 * 4 + 2];
	b34 = B[2 * 4 + 3];
	b41 = B[3 * 4 + 0];
	b42 = B[3 * 4 + 1];
	b43 = B[3 * 4 + 2];
	b44 = B[3 * 4 + 3];
	O4_grid_coordinates_unrank(*F, P[0], P[1], P[2], P[3], 0, 0, verbose_level);
	if (f_vv) {
		cout << "grid point (0,0) = ";
		INT_vec_print(cout, P, 4);
		cout << endl;
		}
	O4_grid_coordinates_unrank(*F, Q[0], Q[1], Q[2], Q[3], 1, 0, verbose_level);
	if (f_vv) {
		cout << "grid point (1,0) = ";
		INT_vec_print(cout, Q, 4);
		cout << endl;
		}
	F->mult_vector_from_the_left(P, B, R, 4, 4);
	F->mult_vector_from_the_left(Q, B, S, 4, 4);
	O4_grid_coordinates_rank(*F, R[0], R[1], R[2], R[3], Rx, Ry, verbose_level);
	O4_grid_coordinates_rank(*F, S[0], S[1], S[2], S[3], Sx, Sy, verbose_level);
	if (f_vv) {
		cout << "Rx=" << Rx << " Ry=" << Ry << " Sx=" << Sx << " Sy=" << Sy << endl;
		}
	if (Ry == Sy) {
		f_switch = FALSE;
		}
	else {
		f_switch = TRUE;
		}
	if (f_vv) {
		cout << "f_switch=" << f_switch << endl;
		}
	if (f_switch) {
		if (b22 == 0 && b24 == 0 && b32 == 0 && b34 == 0) {
			a = 0;
			b = 1;
			f = b12;
			h = b14;
			e = b42;
			g = b44;
			if (e == 0) {
				fv = F->inverse(f);
				c = F->mult(fv, b33);
				d = F->negate(F->mult(fv, b13));
				}
			else {
				ev = F->inverse(e);
				c = F->negate(F->mult(ev, b23));
				d = F->negate(F->mult(ev, b43));
				}
			}
		else {
			a = 1;
			e = b22;
			g = b24;
			f = F->negate(b32);
			h = F->negate(b34);
			if (e == 0) {
				fv = F->inverse(f);
				b = F->mult(fv, b12);
				c = F->mult(fv, b33);
				d = F->negate(F->mult(fv, b13));
				}
			else {
				ev = F->inverse(e);
				b = F->mult(ev, b42);
				c = F->negate(F->mult(ev, b23));
				d = F->negate(F->mult(ev, b43));
				}
			}
		}
	else {
		// no switch
		if (b22 == 0 && b24 == 0 && b42 == 0 && b44 == 0) {
			a = 0;
			b = 1;
			f = b12;
			h = b14;
			e = F->negate(b32);
			g = F->negate(b34);
			if (e == 0) {
				fv = F->inverse(f);
				c = F->negate(F->mult(fv, b43));
				d = F->negate(F->mult(fv, b13));
				}
			else {
				ev = F->inverse(e);
				c = F->negate(F->mult(ev, b23));
				d = F->mult(ev, b33);
				}
			}
		else {
			a = 1;
			e = b22;
			g = b24;
			f = b42;
			h = b44;
			if (e == 0) {
				fv = F->inverse(f);
				b = F->mult(fv, b12);
				c = F->negate(F->mult(fv, b43));
				d = F->negate(F->mult(fv, b13));
				}
			else {
				ev = F->inverse(e);
				b = F->negate(F->mult(ev, b32));
				c = F->negate(F->mult(ev, b23));
				d = F->mult(ev, b33);
				}
			}
		}
	if (f_vv) {
		cout << "a=" << a << " b=" << b << " c=" << c << " d=" << d << endl;
		cout << "e=" << e << " f=" << f << " g=" << g << " h=" << h << endl;
		}
	At[0] = d;
	At[1] = b;
	At[2] = c;
	At[3] = a;
	As[0] = h;
	As[1] = f;
	As[2] = g;
	As[3] = e;
	if (f_v) {
		cout << "At:" << endl;
		print_integer_matrix_width(cout, At, 2, 2, 2, F->log10_of_q);
		cout << "As:" << endl;
		print_integer_matrix_width(cout, As, 2, 2, 2, F->log10_of_q);
		}
	
}

void O4_isomorphism_2to4(finite_field *F, INT *At, INT *As, INT f_switch, INT *B)
{
	INT a, b, c, d, e, f, g, h;

	a = At[3];
	b = At[1];
	c = At[2];
	d = At[0];
	e = As[3];
	f = As[1];
	g = As[2];
	h = As[0];
	if (f_switch) {
		B[0 * 4 + 0] = F->mult(h, d);
		B[0 * 4 + 1] = F->mult(f, b);
		B[0 * 4 + 2] = F->negate(F->mult(f, d));
		B[0 * 4 + 3] = F->mult(h, b);
		B[1 * 4 + 0] = F->mult(g, c);
		B[1 * 4 + 1] = F->mult(e, a);
		B[1 * 4 + 2] = F->negate(F->mult(e, c));
		B[1 * 4 + 3] = F->mult(g, a);
		B[2 * 4 + 0] = F->negate(F->mult(h, c));
		B[2 * 4 + 1] = F->negate(F->mult(f, a));
		B[2 * 4 + 2] = F->mult(f, c);
		B[2 * 4 + 3] = F->negate(F->mult(h, a));
		B[3 * 4 + 0] = F->mult(g, d);
		B[3 * 4 + 1] = F->mult(e, b);
		B[3 * 4 + 2] = F->negate(F->mult(e, d));
		B[3 * 4 + 3] = F->mult(g, b);
		}
	else {
		B[0 * 4 + 0] = F->mult(h, d);
		B[0 * 4 + 1] = F->mult(f, b);
		B[0 * 4 + 2] = F->negate(F->mult(f, d));
		B[0 * 4 + 3] = F->mult(h, b);
		B[1 * 4 + 0] = F->mult(g, c);
		B[1 * 4 + 1] = F->mult(e, a);
		B[1 * 4 + 2] = F->negate(F->mult(e, c));
		B[1 * 4 + 3] = F->mult(g, a);
		B[2 * 4 + 0] = F->negate(F->mult(g, d));
		B[2 * 4 + 1] = F->negate(F->mult(e, b));
		B[2 * 4 + 2] = F->mult(e, d);
		B[2 * 4 + 3] = F->negate(F->mult(g, b));
		B[3 * 4 + 0] = F->mult(h, c);
		B[3 * 4 + 1] = F->mult(f, a);
		B[3 * 4 + 2] = F->negate(F->mult(f, c));
		B[3 * 4 + 3] = F->mult(h, a);
		}
}

void O4_grid_coordinates_rank(finite_field &F, INT x1, INT x2, INT x3, INT x4, INT &grid_x, INT &grid_y, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT a, b, c, d, av, e;
	INT v[2], w[2];
	
	a = x1;
	b = x4;
	c = F.negate(x3);
	d = x2;
	
	if (a) {
		if (a != 1) {
			av = F.inverse(a);
			b = F.mult(b, av);
			c = F.mult(c, av);
			d = F.mult(d, av);
			}
		v[0] = 1;
		w[0] = 1;
		v[1] = c;
		w[1] = b;
		e = F.mult(c, b);
		if (e != d) {
			cout << "O4_grid_coordinates e != d" << endl;
			exit(1);
			}
		}
	else if (b == 0) {
		v[0] = 0;
		v[1] = 1;
		w[0] = c;
		w[1] = d;
		}
	else {
		if (c) {
			cout << "a is zero, b and c are not" << endl;
			exit(1);
			}
		w[0] = 0;
		w[1] = 1;
		v[0] = b;
		v[1] = d;
		}
	PG_element_normalize_from_front(F, v, 1, 2);
	PG_element_normalize_from_front(F, w, 1, 2);
	if (f_v) {
		INT_vec_print(cout, v, 2);
		INT_vec_print(cout, w, 2);
		cout << endl;
		}
	
	PG_element_rank_modified(F, v, 1, 2, grid_x);
	PG_element_rank_modified(F, w, 1, 2, grid_y);
}

void O4_grid_coordinates_unrank(finite_field &F, INT &x1, INT &x2, INT &x3, INT &x4, INT grid_x, INT grid_y, INT verbose_level)
{	
	INT f_v = (verbose_level >= 1);
	INT a, b, c, d;
	INT v[2], w[2];

	PG_element_unrank_modified(F, v, 1, 2, grid_x);
	PG_element_unrank_modified(F, w, 1, 2, grid_y);
	PG_element_normalize_from_front(F, v, 1, 2);
	PG_element_normalize_from_front(F, w, 1, 2);
	if (f_v) {
		INT_vec_print(cout, v, 2);
		INT_vec_print(cout, w, 2);
		cout << endl;
		}

	a = F.mult(v[0], w[0]);
	b = F.mult(v[0], w[1]);
	c = F.mult(v[1], w[0]);
	d = F.mult(v[1], w[1]);
	x1 = a;
	x2 = d;
	x3 = F.negate(c);
	x4 = b;
}

void O4_find_tangent_plane(finite_field &F, INT pt_x1, INT pt_x2, INT pt_x3, INT pt_x4, INT *tangent_plane, INT verbose_level)
{
	INT A[4];
	INT C[3 * 4];
	INT q, size, x, y, z, xx, yy, zz, h, k;
	INT x1, x2, x3, x4;
	INT y1, y2, y3, y4;
	INT f_special = FALSE;
	INT f_complete = FALSE;
	INT base_cols[4];
	INT f_P = FALSE;
	INT rk, det;
	INT vec2[2];

	
	cout << "O4_find_tangent_plane pt_x1=" << pt_x1 
		<< " pt_x2=" << pt_x2 
		<< " pt_x3=" << pt_x3 
		<< " pt_x4=" << pt_x4 << endl; 
	q = F.q;
	size = q + 1;
	A[0] = pt_x1;
	A[3] = pt_x2;
	A[2] = F.negate(pt_x3);
	A[1] = pt_x4;
	
	INT *secants1;
	INT *secants2;
	INT nb_secants = 0;
	INT *complement;
	INT nb_complement = 0;
	
	secants1 = new INT[size * size];
	secants2 = new INT[size * size];
	complement = new INT[size * size];
	for (x = 0; x < size; x++) {
		for (y = 0; y < size; y++) {
			z = x * size + y;
			
			//cout << "trying grid point (" << x << "," << y << ")" << endl;
			//cout << "nb_secants=" << nb_secants << endl;
			O4_grid_coordinates_unrank(F, x1, x2, x3, x4, x, y, 0);

			//cout << "x1=" << x1 << " x2=" << x2 << " x3=" << x3 << " x4=" << x4 << endl;
			
			


			for (k = 0; k < size; k++) {
				PG_element_unrank_modified(F, vec2, 1, 2, k);
				y1 = F.add(F.mult(pt_x1, vec2[0]), F.mult(x1, vec2[1]));
				y2 = F.add(F.mult(pt_x2, vec2[0]), F.mult(x2, vec2[1]));
				y3 = F.add(F.mult(pt_x3, vec2[0]), F.mult(x3, vec2[1]));
				y4 = F.add(F.mult(pt_x4, vec2[0]), F.mult(x4, vec2[1]));
				det = F.add(F.mult(y1, y2), F.mult(y3, y4));
				if (det != 0) {
					continue;
					}
				O4_grid_coordinates_rank(F, y1, y2, y3, y4, xx, yy, 0);
				zz = xx * size + yy;
				if (zz == z) 
					continue;
				C[0] = pt_x1;
				C[1] = pt_x2;
				C[2] = pt_x3;
				C[3] = pt_x4;

				C[4] = x1;
				C[5] = x2;
				C[6] = x3;
				C[7] = x4;

				C[8] = y1;
				C[9] = y2;
				C[10] = y3;
				C[11] = y4;

				rk = F.Gauss_INT(C, f_special, f_complete, base_cols, 
					f_P, NULL, 3, 4, 4, 0);
				if (rk < 3) {
					secants1[nb_secants] = z;
					secants2[nb_secants] = zz;
					nb_secants++;
					}
				
				}

#if 0

			for (xx = 0; xx < size; xx++) {
				for (yy = 0; yy < size; yy++) {
					zz = xx * size + yy;
					if (zz == z) 
						continue;
					O4_grid_coordinates_unrank(F, y1, y2, y3, y4, xx, yy, 0);
					//cout << "y1=" << y1 << " y2=" << y2 << " y3=" << y3 << " y4=" << y4 << endl;
					C[0] = pt_x1;
					C[1] = pt_x2;
					C[2] = pt_x3;
					C[3] = pt_x4;

					C[4] = x1;
					C[5] = x2;
					C[6] = x3;
					C[7] = x4;

					C[8] = y1;
					C[9] = y2;
					C[10] = y3;
					C[11] = y4;

					rk = F.Gauss_INT(C, f_special, f_complete, base_cols, 
						f_P, NULL, 3, 4, 4, 0);
					if (rk < 3) {
						secants1[nb_secants] = z;
						secants2[nb_secants] = zz;
						nb_secants++;
						}
					}
				}
#endif


			}
		}
	cout << "nb_secants=" << nb_secants << endl;
	INT_vec_print(cout, secants1, nb_secants);
	cout << endl;
	INT_vec_print(cout, secants2, nb_secants);
	cout << endl;
	h = 0;
	for (zz = 0; zz < size * size; zz++) {
		if (secants1[h] > zz) {
			complement[nb_complement++] = zz;
			}
		else {
			h++;
			}
		}
	cout << "complement = tangents:" << endl;
	INT_vec_print(cout, complement, nb_complement);
	cout << endl;

	INT *T;
	T = new INT[4 * nb_complement];

	for (h = 0; h < nb_complement; h++) {
		z = complement[h];
		x = z / size;
		y = z % size;
		cout << setw(3) << h << " : " << setw(4) << z << " : " << x << "," << y << " : ";
		O4_grid_coordinates_unrank(F, y1, y2, y3, y4, x, y, verbose_level);
		cout << "y1=" << y1 << " y2=" << y2 << " y3=" << y3 << " y4=" << y4 << endl;
		T[h * 4 + 0] = y1;
		T[h * 4 + 1] = y2;
		T[h * 4 + 2] = y3;
		T[h * 4 + 3] = y4;
		}


	rk = F.Gauss_INT(T, f_special, f_complete, base_cols, 
		f_P, NULL, nb_complement, 4, 4, 0);
	cout << "the rank of the tangent space is " << rk << endl;
	cout << "basis:" << endl;
	print_integer_matrix_width(cout, T, rk, 4, 4, F.log10_of_q);

	if (rk != 3) {
		cout << "rk = " << rk << " not equal to 3" << endl;
		exit(1);
		}
	INT i;
	for (i = 0; i < 12; i++) {
		tangent_plane[i] = T[i];
		}
	delete [] secants1;
	delete [] secants2;
	delete [] complement;
	delete [] T;
	
#if 0
	for (h = 0; h < nb_secants; h++) {
		z = secants1[h];
		zz = secants2[h];
		x = z / size;
		y = z % size;
		xx = zz / size;
		yy = zz % size;
		cout << "(" << x << "," << y << "),(" << xx << "," << yy << ")" << endl;
		O4_grid_coordinates_unrank(F, x1, x2, x3, x4, x, y, verbose_level);
		cout << "x1=" << x1 << " x2=" << x2 << " x3=" << x3 << " x4=" << x4 << endl;
		O4_grid_coordinates_unrank(F, y1, y2, y3, y4, xx, yy, verbose_level);
		cout << "y1=" << y1 << " y2=" << y2 << " y3=" << y3 << " y4=" << y4 << endl;
		}
#endif
}

INT PHG_element_normalize(finite_ring &R, INT *v, INT stride, INT len)
// last unit element made one
{
	INT i, j, a;
	
	for (i = len - 1; i >= 0; i--) {
		a = v[i * stride];
		if (R.is_unit(a)) {
			if (a == 1)
				return i;
			a = R.inverse(a);
			for (j = len - 1; j >= 0; j--) {
				v[j * stride] = R.mult(v[j * stride], a);
				}
			return i;
			}
		}
	cout << "PHG_element_normalize() vector is not free" << endl;
	exit(1);
}


INT PHG_element_normalize_from_front(finite_ring &R, INT *v, INT stride, INT len)
// first non unit element made one
{
	INT i, j, a;
	
	for (i = 0; i < len; i++) {
		a = v[i * stride];
		if (R.is_unit(a)) {
			if (a == 1)
				return i;
			a = R.inverse(a);
			for (j = 0; j < len; j++) {
				v[j * stride] = R.mult(v[j * stride], a);
				}
			return i;
			}
		}
	cout << "PHG_element_normalize_from_front() vector is not free" << endl;
	exit(1);
}

INT PHG_element_rank(finite_ring &R, INT *v, INT stride, INT len)
{
	INT i, j, idx, a, b, r1, r2, rk, N;
	INT f_v = FALSE;
	INT *w;
	INT *embedding;
	
	if (len <= 0) {
		cout << "PHG_element_rank() len <= 0" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "the vector before normalization is ";
		for (i = 0; i < len; i++) {
			cout << v[i * stride] << " ";
			}
		cout << endl;
		}
	idx = PHG_element_normalize(R, v, stride, len);
	if (f_v) {
		cout << "the vector after normalization is ";
		for (i = 0; i < len; i++) {
			cout << v[i * stride] << " ";
			}
		cout << endl;
		}
	w = NEW_INT(len - 1);
	embedding = NEW_INT(len - 1);
	for (i = 0, j = 0; i < len - 1; i++, j++) {
		if (i == idx) {
			j++;
			}
		embedding[i] = j;
		}
	for (i = 0; i < len - 1; i++) {
		w[i] = v[embedding[i] * stride];
		}
	for (i = 0; i < len - 1; i++) {
		a = w[i];
		b = a % R.p;
		v[embedding[i] * stride] = b;
		w[i] = (a - b) / R.p;
		}
	if (f_v) {
		cout << "w=";
		INT_vec_print(cout, w, len - 1);
		cout << endl;
		}
	AG_element_rank(R.e, w, 1, len - 1, r1);
	PG_element_rank_modified(*R.Fp, v, stride, len, r2);

	N = nb_PG_elements(len - 1, R.p);
	rk = r1 * N + r2;

	FREE_INT(w);
	FREE_INT(embedding);

	return rk;
}

void PHG_element_unrank(finite_ring &R, INT *v, INT stride, INT len, INT rk)
{
	INT i, j, idx, r1, r2, N;
	INT f_v = FALSE;
	INT *w;
	INT *embedding;
	
	if (len <= 0) {
		cout << "PHG_element_unrank() len <= 0" << endl;
		exit(1);
		}

	w = NEW_INT(len - 1);
	embedding = NEW_INT(len - 1);

	N = nb_PG_elements(len - 1, R.p);
	r2 = rk % N;
	r1 = (rk - r2) / N;
	
	AG_element_unrank(R.e, w, 1, len - 1, r1);
	PG_element_unrank_modified(*R.Fp, v, stride, len, r2);

	if (f_v) {
		cout << "w=";
		INT_vec_print(cout, w, len - 1);
		cout << endl;
		}

	idx = PHG_element_normalize(R, v, stride, len);
	for (i = 0, j = 0; i < len - 1; i++, j++) {
		if (i == idx) {
			j++;
			}
		embedding[i] = j;
		}
	
	for (i = 0; i < len - 1; i++) {
		v[embedding[i] * stride] += w[i] * R.p;
		}



	FREE_INT(w);
	FREE_INT(embedding);

}

INT nb_PHG_elements(INT n, finite_ring &R)
{
	INT N1, N2;
	
	N1 = nb_PG_elements(n, R.p);
	N2 = nb_AG_elements(n, R.e);
	return N1 * N2;
}

void display_all_PHG_elements(INT n, INT q)
{
	INT *v = NEW_INT(n + 1);
	INT l;
	INT i, j, a;
	finite_ring R;

	R.init(q, 0);
	l = nb_PHG_elements(n, R);
	for (i = 0; i < l; i++) {
		PHG_element_unrank(R, v, 1, n + 1, i);
		cout << i << " : ";
		for (j = 0; j < n + 1; j++) {
			cout << v[j] << " ";
			}
		a = PHG_element_rank(R, v, 1, n + 1);
		cout << " : " << a << endl;
		}
	FREE_INT(v);
}

INT matrix_group_base_len_projective_group(INT n, INT q, INT f_semilinear, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT base_len;

	base_len = n;
	if (q > 2) {
		base_len++;
		}
	if (f_semilinear) {
		base_len++;
		}
	if (f_v) {
		cout << "matrix_group_base_len_projective_group: n=" << n << " q=" << q << " f_semilinear=" << f_semilinear << " base_len = " << base_len << endl;
		}
	return base_len;
}

INT matrix_group_base_len_affine_group(INT n, INT q, INT f_semilinear, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT base_len;

	base_len = 1; // the point 0 takes care of killing the translations
	base_len += n;
	if (f_semilinear) {
		base_len++;
		}
	if (f_v) {
		cout << "matrix_group_base_len_affine_group: n=" << n << " q=" << q << " f_semilinear=" << f_semilinear << " base_len = " << base_len << endl;
		}
	return base_len;
}

INT matrix_group_base_len_general_linear_group(INT n, INT q, INT f_semilinear, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT base_len;

	base_len = 0; // no need to kill translations
	base_len += n;
	if (f_semilinear) {
		base_len++;
		}
	if (f_v) {
		cout << "matrix_group_base_len_general_linear_group: n=" << n << " q=" << q << " f_semilinear=" << f_semilinear << " base_len = " << base_len << endl;
		}
	return base_len;
}



void projective_matrix_group_base_and_orbits(INT n, 
	finite_field *F, INT f_semilinear, 
	INT base_len, INT degree, 
	INT *base, INT *transversal_length, 
	INT **orbit, INT **orbit_inv, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT i;
	INT q;


	if (f_v) {
		cout << "projective_matrix_group_base_and_orbits" << endl;
		}
	q = F->q;
	for (i = 0; i < base_len; i++) {
		base[i] = i;
		}
	if (f_semilinear) {
		base[base_len - 1] = n + F->p;
			// here was an error: the -1 was missing 
			// A.B. 11/11/05
			// no that -1 needs to go 
			// A.B. 3/9/2006
		}
	//transversal_length[0] = nb_PG_elements(n - 1, q);
	for (i = 0; i < n; i++) {
		transversal_length[i] = nb_PG_elements_not_in_subspace(n - 1, i - 1, q);
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits transversal " << i << " of length " << transversal_length[i] << endl;
			}
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits before PG_element_modified_not_in_subspace_perm" << endl;
			}
		PG_element_modified_not_in_subspace_perm(n - 1, i - 1, 
			*F, orbit[i], orbit_inv[i], 0);
			// global function in GALOIS/projective.C

		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits after PG_element_modified_not_in_subspace_perm" << endl;
			}
		
		if (FALSE) {
			print_set(cout, degree, orbit[i]);
			cout << endl;
			print_set(cout, degree, orbit_inv[i]);
			cout << endl;
			}
		}
	if (q > 2) {
		transversal_length[i] = nb_AG_elements(n - 1, q - 1);
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits: diagonal transversal " << i << " of length " << transversal_length[i] << endl;
			}
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits before diagonal_orbit_perm" << endl;
			}
		diagonal_orbit_perm(n, *F, orbit[i], orbit_inv[i], 0);
			// global function in GALOIS/projective.C
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits after diagonal_orbit_perm" << endl;
			}

		if (FALSE) {
			print_set(cout, degree, orbit[i]);
			cout << endl;
			print_set(cout, degree, orbit_inv[i]);
			cout << endl;
			}
		i++;
		}
	if (f_semilinear) {
		transversal_length[i] = F->e;
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits: frobenius transversal " << i << " of length " << transversal_length[i] << endl;
			}
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits before frobenius_orbit_perm" << endl;
			}
		frobenius_orbit_perm(n, *F, orbit[i], orbit_inv[i], verbose_level - 2);
			// global function in GALOIS/projective.C
		if (f_vv) {
			cout << "projective_matrix_group_base_and_orbits after frobenius_orbit_perm" << endl;
			}

		if (FALSE) {
			print_set(cout, degree, orbit[i]);
			cout << endl;
			print_set(cout, degree, orbit_inv[i]);
			cout << endl;
			}
		i++;
		}
	if (f_v) {
		cout << "projective_matrix_group_base_and_orbits base: ";
		INT_vec_print(cout, base, base_len);
		cout << endl;
		cout << "projective_matrix_group_base_and_orbits transversal_length: ";
		INT_vec_print(cout, transversal_length, base_len);
		cout << endl;
		}
	if (f_v) {
		cout << "projective_matrix_group_base_and_orbits done" << endl;
		}
}

void affine_matrix_group_base_and_transversal_length(INT n, 
	finite_field *F, INT f_semilinear, 
	INT base_len, INT degree, 
	INT *base, INT *transversal_length, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT i, c;
	INT q;


	if (f_v) {
		cout << "affine_matrix_group_base_and_transversal_length" << endl;
		}
	q = F->q;
	c = 0;
	base[c] = 0;
	transversal_length[c] = i_power_j(q, n);
	c++;
	for (i = 0; i < n; i++) {
		base[c] = i_power_j(q, i);
		transversal_length[c] = i_power_j(q, n) - i_power_j(q, i);
		c++;
		}
	if (f_semilinear) {
		base[c] = F->q + F->p;
		transversal_length[c] = F->e;
		c++;
		}
	if (c != base_len) {
		cout << "affine_matrix_group_base_and_transversal_length c != base_len" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "affine_matrix_group_base_and_transversal_length base: ";
		INT_vec_print(cout, base, base_len);
		cout << endl;
		cout << "affine_matrix_group_base_and_transversal_length transversal_length: ";
		INT_vec_print(cout, transversal_length, base_len);
		cout << endl;
		}
	if (f_v) {
		cout << "affine_matrix_group_base_and_transversal_length done" << endl;
		}
}


void general_linear_matrix_group_base_and_transversal_length(INT n, 
	finite_field *F, INT f_semilinear, 
	INT base_len, INT degree, 
	INT *base, INT *transversal_length, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	INT i, c;
	INT q;


	if (f_v) {
		cout << "general_linear_matrix_group_base_and_transversal_length" << endl;
		}
	q = F->q;
	c = 0;
	for (i = 0; i < n; i++) {
		base[c] = i_power_j(q, i);
		transversal_length[c] = i_power_j(q, n) - i_power_j(q, i);
		c++;
		}
	if (f_semilinear) {
		base[c] = F->q + F->p;
		transversal_length[c] = F->e;
		c++;
		}
	if (c != base_len) {
		cout << "general_linear_matrix_group_base_and_transversal_length c != base_len" << endl;
		cout << "c=" << c << endl;
		cout << "base_len=" << base_len << endl;
		exit(1);
		}
	if (f_v) {
		cout << "general_linear_matrix_group_base_and_transversal_length base: ";
		INT_vec_print(cout, base, base_len);
		cout << endl;
		cout << "general_linear_matrix_group_base_and_transversal_length transversal_length: ";
		INT_vec_print(cout, transversal_length, base_len);
		cout << endl;
		}
	if (f_v) {
		cout << "general_linear_matrix_group_base_and_transversal_length done" << endl;
		}
}


void strong_generators_for_projective_linear_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT *&data, INT &size, INT &nb_gens, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT h, u, cur;
	
	if (f_v) {
		cout << "strong_generators_for_projective_linear_group" << endl;
		}
	size = n * n;
	if (f_semilinear) {
		size++;
		}
	nb_gens = 0;
	if (f_semilinear) {
		nb_gens++;
		}
	if (F->q > 2) {
		nb_gens += n - 1;
		}
	nb_gens += (n - 1) * F->e;
	nb_gens += n - 1;
	data = NEW_INT(size * nb_gens);

	cur = 0;
	if (f_semilinear) {
		F->identity_matrix(data + cur * size, n);
		data[cur * size + n * n] = 1;
		cur++;
		}
	if (F->q > 2) {
		// the primitive elements on the diagonal:
		for (h = 0; h < n - 1; h++) {
			if (f_vv) {
				cout << "generators for primitive elements on the diagonal:" << endl;
				}
			INT_vec_zero(data + cur * size, size);
			F->identity_matrix(data + cur * size, n);

			data[cur * size + h * n + h] = F->primitive_root();
			if (f_semilinear) {
				data[cur * size + n * n] = 0;
				}
			cur++;
			} // next h
		} // if
	// the entries in the last row:
	for (h = 0; h < n - 1; h++) {
		if (f_vv) {
			cout << "generators for entries in the last row (e=" << F->e << "):" << endl;
			}
		for (u = 0; u < F->e; u++) {
			//INT_vec_zero(data + cur * size, size);
			F->identity_matrix(data + cur * size, n);

			data[cur * size + (n - 1) * n + h] = i_power_j(F->p, u);
			if (f_semilinear) {
				data[cur * size + n * n] = 0;
				}
			cur++;
			} // next u
		} // next h
	// the swaps along the diagonal:
	for (h = n - 2; h >= 0; h--) {
		if (f_vv) {
			cout << "generators for swaps along the diagonal:" << endl;
			}
		//INT_vec_zero(data + cur * size, size);
		F->identity_matrix(data + cur * size, n);
		data[cur * size + h * n + h] = 0;
		data[cur * size + h * n + h + 1] = 1;
		data[cur * size + (h + 1) * n + h] = 1;
		data[cur * size + (h + 1) * n + h + 1] = 0;
		if (f_semilinear) {
			data[cur * size + n * n] = 0;
			}
		cur++;
		} // next h
	if (cur != nb_gens) {
		cout << "strong_generators_for_projective_linear_group cur != nb_gens" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "strong_generators_for_projective_linear_group done" << endl;
		}
}


void strong_generators_for_affine_linear_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT *&data, INT &size, INT &nb_gens, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT h, u, cur;
	
	if (f_v) {
		cout << "strong_generators_for_affine_linear_group" << endl;
		}
	size = n * n + n;
	if (f_semilinear) {
		size++;
		}
	nb_gens = 0;
	if (f_semilinear) {
		nb_gens++; // the field automorphism
		}
	nb_gens += (n - 1) * F->e; // the bottom layer

	if (F->q > 2) {
		nb_gens++;
		}

	nb_gens += n - 1; // the transpositions

	nb_gens += n * F->e; // the translations

	data = NEW_INT(size * nb_gens);

	cur = 0;
	if (f_semilinear) {
		INT_vec_zero(data + cur * size, size);
		F->identity_matrix(data + cur * size, n);
		data[cur * size + n * n + n] = 1;
		cur++;
		}

	// the entries in the last row:
	for (h = 0; h < n - 1; h++) {
		if (f_vv) {
			cout << "generators for entries in the last row (e=" << F->e << "):" << endl;
			}
		for (u = 0; u < F->e; u++) {
			INT_vec_zero(data + cur * size, size);
			F->identity_matrix(data + cur * size, n);

			data[cur * size + (n - 1) * n + h] = i_power_j(F->p, u);
			if (f_semilinear) {
				data[cur * size + n * n + n] = 0;
				}
			cur++;
			} // next u
		} // next h

	if (F->q > 2) {
		// the primitive element on the last diagonal:
		h = n - 1;
		if (f_vv) {
			cout << "generators for primitive element on the last diagonal:" << endl;
			}
		INT_vec_zero(data + cur * size, size);
		F->identity_matrix(data + cur * size, n);

		data[cur * size + h * n + h] = F->primitive_root();
		if (f_semilinear) {
			data[cur * size + n * n + n] = 0;
			}
		cur++;
		} // if


	// the swaps along the diagonal:
	for (h = n - 2; h >= 0; h--) {
		if (f_vv) {
			cout << "generators for swaps along the diagonal:" << endl;
			}
		INT_vec_zero(data + cur * size, size);
		F->identity_matrix(data + cur * size, n);
		data[cur * size + h * n + h] = 0;
		data[cur * size + h * n + h + 1] = 1;
		data[cur * size + (h + 1) * n + h] = 1;
		data[cur * size + (h + 1) * n + h + 1] = 0;
		if (f_semilinear) {
			data[cur * size + n * n + n] = 0;
			}
		cur++;
		} // next h

	// the translations:
	for (h = 0; h < n; h++) {
		for (u = 0; u < F->e; u++) {
			INT_vec_zero(data + cur * size, size);
			F->identity_matrix(data + cur * size, n);

			data[cur * size + n * n + h] = i_power_j(F->p, u);
			if (f_semilinear) {
				data[cur * size + n * n + n] = 0;
				}
			cur++;
			} // next u
		} // next h

	if (cur != nb_gens) {
		cout << "strong_generators_for_affine_linear_group cur != nb_gens" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "strong_generators_for_affine_linear_group done" << endl;
		}
}

void strong_generators_for_general_linear_group(INT n, finite_field *F, 
	INT f_semilinear, 
	INT *&data, INT &size, INT &nb_gens, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT f_vv = (verbose_level >= 2);
	INT h, u, cur;
	
	if (f_v) {
		cout << "strong_generators_for_general_linear_group" << endl;
		}
	size = n * n;
	if (f_semilinear) {
		size++;
		}
	nb_gens = 0;
	if (f_semilinear) {
		nb_gens++; // the field automorphism
		}
	nb_gens += (n - 1) * F->e; // the bottom layer

	if (F->q > 2) {
		nb_gens++;
		}

	nb_gens += n - 1; // the transpositions


	data = NEW_INT(size * nb_gens);

	cur = 0;
	if (f_semilinear) {
		INT_vec_zero(data + cur * size, size);
		F->identity_matrix(data + cur * size, n);
		data[cur * size + n * n] = 1;
		cur++;
		}

	// the entries in the last row:
	for (h = 0; h < n - 1; h++) {
		if (f_vv) {
			cout << "generators for entries in the last row (e=" << F->e << "):" << endl;
			}
		for (u = 0; u < F->e; u++) {
			INT_vec_zero(data + cur * size, size);
			F->identity_matrix(data + cur * size, n);

			data[cur * size + (n - 1) * n + h] = i_power_j(F->p, u);
			if (f_semilinear) {
				data[cur * size + n * n] = 0;
				}
			cur++;
			} // next u
		} // next h

	if (F->q > 2) {
		// the primitive element on the last diagonal:
		h = n - 1;
		if (f_vv) {
			cout << "generators for primitive element on the last diagonal:" << endl;
			}
		INT_vec_zero(data + cur * size, size);
		F->identity_matrix(data + cur * size, n);

		data[cur * size + h * n + h] = F->primitive_root();
		if (f_semilinear) {
			data[cur * size + n * n] = 0;
			}
		cur++;
		} // if


	// the swaps along the diagonal:
	for (h = n - 2; h >= 0; h--) {
		if (f_vv) {
			cout << "generators for swaps along the diagonal:" << endl;
			}
		INT_vec_zero(data + cur * size, size);
		F->identity_matrix(data + cur * size, n);
		data[cur * size + h * n + h] = 0;
		data[cur * size + h * n + h + 1] = 1;
		data[cur * size + (h + 1) * n + h] = 1;
		data[cur * size + (h + 1) * n + h + 1] = 0;
		if (f_semilinear) {
			data[cur * size + n * n] = 0;
			}
		cur++;
		} // next h


	if (cur != nb_gens) {
		cout << "strong_generators_for_general_linear_group cur != nb_gens" << endl;
		exit(1);
		}
	if (f_v) {
		cout << "strong_generators_for_general_linear_group done" << endl;
		}
}



