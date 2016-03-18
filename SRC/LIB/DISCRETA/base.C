// base.C
//
// Anton Betten
// 18.12.1998
// moved from D2 to ORBI Nov 15, 2007


#include "orbiter.h"

#undef BASE_SETTYPE_VERBOSE


base::base()
{
	k = BASE;
	clearself();
}

base::base(const base &x)
	// copy constructor:    this := x
{
	// cout << "base::copy constructor for object: " << x << "\n";
	clearself();
#if 0
	if (x.k != x.s_virtual_kind()) {
		x.c_kind(k);
		}
#endif
	// cout << "base::copy constructor, calling copyobject_to()\n";
	const_cast<base &>(x).copyobject_to(*this);
	// cout << "base::copy constructor finished\n";
}

base& base::operator = (const base &x)
	// copy assignment
{
	// cout << "base::operator = (copy assignment)" << endl;
	// cout << "source=" << x << endl;
	copyobject(const_cast<base &>(x));
	return *this;
}

base::~base()
{
	// cout << "base::~base()\n";
	// printobjectkindln(cout);
	// cout << endl;
	freeself_kind(k); // virtual kind may be different form k ! */
}

void base::freeself_base()
{
	// cout << "base::freeself_base()\n";
	// printobjectkindln(cout);
	// cout << "self=" << self.vector_pointer << endl;
	if (s_kind() != BASE) {
		cout << "freeself() not implemented for class ";
		printobjectkindln(cout);
		exit(1);
		// freeself();
		// return;
		}
	clearself();
}

void base::freeself()
{
	freeself_kind(s_kind());
}

void base::freeself_kind(kind k)
{
	switch (k) {
		case BASE: freeself_base(); break;
		case INTEGER: as_integer().freeself_integer(); break;
		case VECTOR: as_vector().freeself_vector(); break;
		case NUMBER_PARTITION: as_number_partition().freeself_number_partition(); break;
		case PERMUTATION: as_permutation().freeself_permutation(); break;
		case MATRIX: as_matrix().freeself_matrix(); break;
		case LONGINTEGER: as_longinteger().freeself_longinteger(); break;
		case MEMORY: as_memory().freeself_memory(); break;
		//case PERM_GROUP: as_perm_group().freeself_perm_group(); break;
		//case PERM_GROUP_STAB_CHAIN: as_perm_group_stab_chain().freeself_perm_group_stab_chain(); break;
		case UNIPOLY: as_unipoly().freeself_unipoly(); break;
		case SOLID: as_solid().freeself_solid(); break;
		case BITMATRIX: as_bitmatrix().freeself_bitmatrix(); break;
		//case PC_PRESENTATION: as_pc_presentation().freeself_pc_presentation(); break;
		//case PC_SUBGROUP: as_pc_subgroup().freeself_pc_subgroup(); break;
		//case GROUP_WORD: as_group_word().freeself_group_word(); break;
		//case GROUP_TABLE: as_group_table().freeself_group_table(); break;
		// case ACTION: as_action().freeself_action(); break;
		case GEOMETRY: as_geometry().freeself_geometry(); break;
		case HOLLERITH: as_hollerith().freeself_hollerith(); break;
		case GROUP_SELECTION: as_group_selection().freeself_group_selection(); break;
		case BT_KEY: as_bt_key().freeself_bt_key(); break;
		case DATABASE: as_database().freeself_database(); break;
		case BTREE: as_btree().freeself_btree(); break;
		case DESIGN_PARAMETER_SOURCE: as_design_parameter_source().freeself_design_parameter_source(); break;
		case DESIGN_PARAMETER: as_design_parameter().freeself_design_parameter(); break;
		default: cout << "base::freeself_kind(), unknown kind: k= " << kind_ascii(k) << "\n";
		}
}

void base::settype_base()
{
#ifdef BASE_SETTYPE_VERBOSE
	if (s_kind() != BASE) {
		cout << "warning: base::settype_base() converting from " 
			<< kind_ascii(s_kind()) << " to BASE\n";
		}
#endif
	new(this) base;
}

kind base::s_kind()
{
	kind kv;
	
	kv = s_virtual_kind();
	if (k != kv) {
		cout << "base::s_kind(): kind != virtual kind\n";
		cout << "k=" << kind_ascii(k) << ", virtual kind = " << kind_ascii(kv) << endl;
		exit(1);
		}
	return k;
}

kind base::s_virtual_kind()
{
	return BASE;
}

void base::c_kind(kind k)
{
	// cout << "base::c_kind(), k= " << kind_ascii(k) << "\n";
	switch (k) {
		case BASE: settype_base(); break;
		case INTEGER: as_integer().settype_integer(); break;
		case VECTOR: as_vector().settype_vector(); break;
		case NUMBER_PARTITION: as_number_partition().settype_number_partition(); break;
		case PERMUTATION: as_permutation().settype_permutation(); break;
		case MATRIX: as_matrix().settype_matrix(); break;
		case LONGINTEGER: as_longinteger().settype_longinteger(); break;
		case MEMORY: as_memory().settype_memory(); break;
		//case PERM_GROUP: as_perm_group().settype_perm_group(); break;
		//case PERM_GROUP_STAB_CHAIN: as_perm_group_stab_chain().settype_perm_group_stab_chain(); break;
		case UNIPOLY: as_unipoly().settype_unipoly(); break;
		case SOLID: as_solid().settype_solid(); break;
		case BITMATRIX: as_bitmatrix().settype_bitmatrix(); break;
		//case PC_PRESENTATION: as_pc_presentation().settype_pc_presentation(); break;
		//case PC_SUBGROUP: as_pc_subgroup().settype_pc_subgroup(); break;
		//case GROUP_WORD: as_group_word().settype_group_word(); break;
		//case GROUP_TABLE: as_group_table().settype_group_table(); break;
		// case ACTION: as_action().settype_action(); break;
		case GEOMETRY: as_geometry().settype_geometry(); break;
		case HOLLERITH: as_hollerith().settype_hollerith(); break;
		case GROUP_SELECTION: as_group_selection().settype_group_selection(); break;
		case BT_KEY: as_bt_key().settype_bt_key(); break;
		case DATABASE: as_database().settype_database(); break;
		case BTREE: as_btree().settype_btree(); break;
		case DESIGN_PARAMETER_SOURCE: as_design_parameter_source().settype_design_parameter_source(); break;
		case DESIGN_PARAMETER: as_design_parameter().settype_design_parameter(); break;
		default: cout << "base::c_kind(), k= " << kind_ascii(k) << " unknown\n";
		}
	if (s_kind() != k) {
		cout << "base::c_kind() did not work\n";
		}
	// cout << "base::c_kind() finished \n";
}

void base::swap(base &a)
{
	kind k, ka;
	OBJECTSELF s, sa;
	
	k = s_kind();
	ka = a.s_kind();
	s = self;
	sa = a.self;
	c_kind(ka);
	self = sa;
	a.c_kind(k);
	a.self = s;
}

void base::copyobject(base &x)
// this := x
{
	// cout << "base::copyobject\n";
	// cout << "source=" << x << endl;
	x.copyobject_to(*this);
}

void base::copyobject_to(base &x)
{
	kind k = s_kind();
	OBJECTSELF s = self;
	
	if (k != BASE) {
		cout << "error: base::copyobject_to() for object of kind " << kind_ascii(k) << endl;
		exit(1);
		}
	cout << "warning: base::copyobject_to() for object: " << *this << "\n";
	x.freeself();
	x.c_kind(k);
	x.self = s;
}

ostream& base::print(ostream& ost)
{
	ost << "object of kind BASE";
	return ost;
}

ostream& base::println(ostream &ost)
{
	print(ost) << endl;
	return ost;
}

void base::print_to_hollerith(hollerith& h)
{
	ostringstream s;
	int l;
	char *p;
	
	s << *this << ends;
	l = s.str().length();
	p = new char [l + 1];
	s.str().copy(p, l, 0);
	p[l] = 0;
	h.init(p);
	delete [] p;
}

ostream& base::printobjectkind(ostream& ost)
{
	::printobjectkind(ost, s_kind());
	return ost;
}

ostream& base::printobjectkindln(ostream& ost)
{
	printobjectkind(ost) << "\n";
	return ost;
}

INT& base::s_i_i()
{
	if (s_kind() != INTEGER) {
		cout << "base::s_i_i() not an integer, objectkind=";
		printobjectkindln(cout);
		exit(1);
		}
	return as_integer().s_i();
}

void base::m_i_i(INT i)
{
	change_to_integer().m_i(i);
}


INT base::compare_with(base &a)
{
	if (s_kind() != BASE) {
		cout << "compare_with() not implemented for class ";
		printobjectkindln(cout);
		exit(1);
		// return compare_with(a);
		}
	NOT_EXISTING_FUNCTION("base::compare_with");
	exit(1);
	return 0;
}

INT base::eq(base &a)
{
	INT r = compare_with(a);
	if (r == 0)
		return TRUE;
	else
		return FALSE;
}

INT base::neq(base &a)
{
	INT r = compare_with(a);
	if (r != 0)
		return TRUE;
	else
		return FALSE;
}

INT base::le(base &a)
{
	INT r = compare_with(a);
	if (r <= 0)
		return TRUE;
	else
		return FALSE;
}

INT base::lt(base &a)
{
	//cout << "lt(): " << *this << ", " << a;
	INT r = compare_with(a);
	//cout << " r=" << r << endl;
	if (r < 0)
		return TRUE;
	else
		return FALSE;
}

INT base::ge(base &a)
{
	INT r = compare_with(a);
	if (r >= 0)
		return TRUE;
	else
		return FALSE;
}

INT base::gt(base &a)
{
	INT r = compare_with(a);
	if (r > 0)
		return TRUE;
	else
		return FALSE;
}

INT base::is_even()
{
	base a, q, r;
	
	a.m_i_i(2);
	integral_division(a, q, r, 0);
	if (r.is_zero())
		return TRUE;
	else
		return FALSE;
}

INT base::is_odd()
{
	if (is_even())
		return FALSE;
	else
		return TRUE;
}



// mathematical functions:

void base::mult(base &x, base &y)
{
	x.mult_to(y, *this);
}

void base::mult_mod(base &x, base &y, base &p)
{
	base z;
	
	x.mult_to(y, z);
	z.modulo(p);
	swap(z);
}

void base::mult_to(base &x, base &y)
{
	if (s_kind() != BASE) {
		cout << "mult_to() not implemented for class ";
		printobjectkindln(cout);
		exit(1);
		// mult_to(x, y);
		// return;
		}
	NOT_EXISTING_FUNCTION("base::mult_to");
	exit(1);
}

INT base::invert()
{
	base a;
	INT ret;
	
	ret = invert_to(a);
	// cout << "base::invert() a="; a.println();
	// freeself();
	swap(a);
	return ret;
}

INT base::invert_mod(base &p)
{
	base u, v, g;
	
	// cout << "base::invert_mod() this=" << *this << endl;
	extended_gcd(p, u, v, g, 0);
	// cout << "base::invert_mod(): gcd = " << g << " = " << u << " * " << *this << " + " << v << " * " << p << endl;
	if (!g.is_one()) {
		return FALSE;
		}
	swap(u);
	// normalize(p);
	return TRUE;
}

INT base::invert_to(base &x)
{
	if (s_kind() != BASE) {
		// cout << "invert_to() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		return invert_to(x);
		}
	NOT_EXISTING_FUNCTION("base::invert_to");
	exit(1);
}

void base::mult_apply(base &x)
{
	base a;
	
	// cout << "base::mult_apply() calling mult_to()\n";
	mult_to(x, a);
	freeself();
	swap(a);
}

#if 1
base& base::power_int(INT l)
{
	base a, b;
	
	if (l < 0) {
		invert();
		l *= -1;
		}
	a = *this;
	a.one();
	b = *this;
	while (l) {
		if (EVEN(l)) {
			b *= b;
			l >>= 1;
			}
		if (ODD(l)) {
			a *= b;
			l--;
			}
		}
	*this = a;
	return *this;
}
#endif

#if 0
base& base::power_int(INT l)
{
	base *a = callocobject(BASE);
	base *b = callocobject(BASE);
	
	*a = *this;
	a->one();
	*b = *this;
	while (l) {
		if (EVEN(l)) {
			*b *= *b;
			l >>= 1;
			}
		if (ODD(l)) {
			*a *= *b;
			l--;
			}
		}
	*this = *a;
	freeobject(a);
	freeobject(b);
	return *this;
}
#endif

base& base::power_int_mod(INT l, base &p)
{
	base a, b, c;
	
	a = *this;
	a.one();
	b = *this;
	// cout << "base:power_int_mod() x=" << *this << " l=" << l << " p=" << p << endl;
	while (l) {
		// cout << "= " << a << " * " << b << "^" << l << endl;
		if (EVEN(l)) {
			c.mult_mod(b, b, p);
			c.swap(b);
			l >>= 1;
			}
		// cout << "= " << a << " * " << b << "^" << l << endl;
		if (ODD(l)) {
			c.mult_mod(a, b, p);
			c.swap(a);
			l--;
			}
		}
	// cout << "= " << a << " * " << b << "^" << l << endl;
	*this = a;
	return *this;
}

base& base::power_longinteger(longinteger& l)
{
	base a, b, c;
	
	a = *this;
	a.one();
	b = *this;
	while (!l.is_zero()) {
		if (a.s_kind() == LONGINTEGER) {
			longinteger &B = b.as_longinteger();
			INT d;
			
			d = B.s_len();
			cout << "l=" << l << " digits=" << d << endl;
			}
		if (l.is_even()) {
			b *= b;
			// c.mult(b, b);
			// b.swap(c);
			l.divide_out_int(2);
			}
		if (l.is_odd()) {
			a *= b;
			// c.mult(a, b);
			// a.swap(c);
			l.dec();
			}
		}
	*this = a;
	return *this;
}

base& base::power_longinteger_mod(longinteger& l, base &p)
{
	base a, b, c;
	
	a = *this;
	a.one();
	b = *this;
	while (!l.is_zero()) {
		if (a.s_kind() == LONGINTEGER) {
			longinteger &B = a.as_longinteger();
			INT d;
			
			d = B.s_len();
			cout << "l=" << l << " digits=" << d << endl;
			}
		if (l.is_even()) {
			c.mult_mod(b, b, p);
			c.swap(b);
			l.divide_out_int(2);
			}
		if (l.is_odd()) {
			c.mult_mod(a, b, p);
			c.swap(a);
			l.dec();
			}
		}
	*this = a;
	return *this;
}

base& base::commutator(base &x, base &y)
{
	base xv, yv, a;
	
	x.invert_to(xv);
	y.invert_to(yv);
	a.mult(xv, yv);
	a *= x;
	a *= y;
	swap(a);
	xv.freeself();
	yv.freeself();
	return *this;
}

base& base::conjugate(base &x, base &y)
{
	base yv, a;
	
	// cout << "base::conjugate: y.invert_to(yv)\n";
	y.invert_to(yv);
	// cout << "yv= " << yv << endl;
	// cout << "x= " << x << endl;
	// cout << "base::conjugate: a.mult(yv, x)\n";
	a.mult(yv, x);
	// cout << "a=" << a << endl;
	// cout << "base::conjugate: a *= y\n";
	a *= y;
	swap(a);
	return *this;
}

base& base::divide_by(base& x)
{
	base q, r;
	integral_division(x, q, r, 0);
	swap(q);
	return *this;
}

base& base::divide_by_exact(base& x)
{
	base q;
	integral_division_exact(x, q);
	swap(q);
	return *this;
}

#undef DEBUG_ORDER

INT base::order()
{
	base a, b;
	INT i = 1;
	
	copyobject_to(a);
	copyobject_to(b);
	while (!b.is_one()) {
#ifdef DEBUG_ORDER
		cout << "base::order b^" << i << "=" << b << endl;
#endif
		b *= a;
		i++;
		}
#ifdef DEBUG_ORDER
	cout << "base::order b^" << i << "=" << b << " is one " << endl;
#endif
	return i;
}

INT base::order_mod(base &p)
{
	base a, b, c;
	INT i = 1;
	
	copyobject_to(a);
	copyobject_to(b);
	while (!b.is_one()) {
		c.mult_mod(a, b, p);
		b.swap(c);
		i++;
		}
	return i;
}

void base::add(base &x, base &y)
{
	// cout << "base::add() x=" << x << ", y=" << y << endl;
	x.add_to(y, *this);
}

void base::add_mod(base &x, base &y, base &p)
{
	base z;
	
	x.add_to(y, z);
	z.modulo(p);
	swap(z);
}

void base::add_to(base &x, base &y)
{
	if (s_kind() != BASE) {
		// cout << "add_to() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		add_to(x, y);
		return;
		}
	NOT_EXISTING_FUNCTION("base::add_to");
	exit(1);
}

void base::negate()
{
	base a;
	
	negate_to(a);
	swap(a);
}

void base::negate_to(base &x)
{
	if (s_kind() != BASE) {
		cout << "negate_to() not implemented for class ";
		printobjectkindln(cout);
		exit(1);
		// negate_to(x);
		// return;
		}
	NOT_EXISTING_FUNCTION("base::negate_to");
	exit(1);
}

void base::add_apply(base &x)
{
	base a;
	
	add_to(x, a);
	swap(a);
}

void base::normalize(base &p)
{
	if (s_kind() != BASE) {
		cout << "normalize() not implemented for class ";
		printobjectkindln(cout);
		exit(1);
		// normalize(p);
		// return;
		}
	NOT_EXISTING_FUNCTION("base::normalize");
	exit(1);
}

void base::zero()
{
	if (s_kind() != BASE) {
		// cout << "zero() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		zero();
		return;
		}
	NOT_EXISTING_FUNCTION("base::zero");
	exit(1);
}

void base::one()
{
	if (s_kind() != BASE) {
		// cout << "one() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		one();
		return;
		}
	NOT_EXISTING_FUNCTION("base::one");
	exit(1);
}

void base::m_one()
{
	if (s_kind() != BASE) {
		// cout << "m_one() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		m_one();
		return;
		}
	NOT_EXISTING_FUNCTION("base::m_one");
	exit(1);
}

void base::homo_z(INT z)
{
	if (s_kind() != BASE) {
		// cout << "homo_z() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		homo_z(z);
		return;
		}
	NOT_EXISTING_FUNCTION("base::homo_z");
	exit(1);
}

void base::inc()
{
	if (s_kind() != BASE) {
		// cout << "inc() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		inc();
		return ;
		}
	NOT_EXISTING_FUNCTION("base::inc");
	exit(1);
}

void base::dec()
{
	if (s_kind() != BASE) {
		// cout << "dec() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		dec();
		return;
		}
	NOT_EXISTING_FUNCTION("base::dec");
	exit(1);
}

INT base::is_zero()
{
	if (s_kind() != BASE) {
		// cout << "is_zero() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		return is_zero();
		}
	NOT_EXISTING_FUNCTION("base::is_zero");
	exit(1);
}

INT base::is_one()
{
	if (s_kind() != BASE) {
		// cout << "is_one() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		return is_one();
		}
	NOT_EXISTING_FUNCTION("base::is_one");
	exit(1);
}

INT base::is_m_one()
{
	if (s_kind() != BASE) {
		// cout << "is_m_one() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		return is_m_one();
		}
	NOT_EXISTING_FUNCTION("base::is_m_one");
	exit(1);
}

base& base::factorial(INT z)
{
	base a, b;
	
	a.m_i_i(1);
	while (z) {
		b.m_i_i(z);
		a *= b;
		z--;
		}
	*this = a;
	return *this;
}

base& base::i_power_j(INT i, INT j)
{
	m_i_i(i);
	power_int(j);
	return *this;
}

INT base::compare_with_euklidean(base &a)
{
	if (s_kind() != BASE) {
		// cout << "compare_with_euklidean() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		return compare_with_euklidean(a);
		}
	NOT_EXISTING_FUNCTION("base::compare_with_euklidean");
	exit(1);
}

void base::integral_division(base &x, base &q, base &r, INT verbose_level)
{
	if (s_kind() != BASE) {
		// cout << "integral_division() not implemented for class ";
		// printobjectkindln(cout);
		// exit(1);
		integral_division(x, q, r, verbose_level);
		return;
		}
	NOT_EXISTING_FUNCTION("base::integral_division");
	exit(1);
}

void base::integral_division_exact(base &x, base &q)
{
	base r;
	
	if (s_kind() != BASE) {
		integral_division(x, q, r, 0);
		if (r.is_zero())
			return;
		cout << "integral_division_exact() remainder not zero\n";
		cout << "this=" << *this << " divided by " << x << " gives remainder " << r << endl;
		exit(1);
		}
	NOT_EXISTING_FUNCTION("base::integral_division");
	exit(1);
}

void base::integral_division_by_integer(INT x, base &q, base &r)
{
	base a;
	
	a.m_i_i(x);
	integral_division(a, q, r, 0);
}

void base::integral_division_by_integer_exact(INT x, base &q)
{
	base a;
	
	a.m_i_i(x);
	integral_division_exact(a, q);
}

void base::integral_division_by_integer_exact_apply(INT x)
{
	base a, q;
	
	a.m_i_i(x);
	integral_division_exact(a, q);
	swap(q);
}

INT base::is_divisor(base& y)
{
	base q, r;
	
	y.integral_division(*this, q, r, 0);
	if (r.is_zero())
		return TRUE;
	else
		return FALSE;
}

void base::modulo(base &p)
{
	base q, r;
	
	integral_division(p, q, r, 0);
	swap(r);
}

void base::extended_gcd(base &n, base &u, base &v, base &g, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	//INT f_vv = (verbose_level >= 2);
	base sign1, sign2;
	INT c;
	
	if (f_v) {
		cout << "base::extended_gcd() m=" << *this << " n=" << n << endl;
		}
	c = compare_with_euklidean(n);
	if (c < 0) {
		n.extended_gcd(*this, v, u, g, verbose_level);
		return;
		}
	if (f_v) {
		cout << "base::extended_gcd() m=" << *this << "(" << kind_ascii(s_kind()) << ")" 
			<< "n=" << n << "(" << kind_ascii(n.s_kind()) << ")" << endl;
		}
	u = *this;
	v = n;
	if (/* c == 0 ||*/ n.is_zero()) {
		u.one();
		v.zero();
		g = *this;
		return;
		}

	if (s_kind() == INTEGER) {
		INT a;
		a = s_i_i();
		if (a < 0) {
			sign1.m_i_i(-1);
			m_i_i(-a);
			}
		else
			sign1.m_i_i(1);
		
		a = n.s_i_i();
		if (a < 0) {
			sign2.m_i_i(-1);
			n.m_i_i(-a);
			}
		else
			sign2.m_i_i(1);
		
		}


	base M, N, Q, R;
	base u1, u2, u3, v1, v2, v3;
	
	M = *this;
	N = n;
	u1 = *this; u1.one();
	u2 = *this; u2.zero();
	v1 = n; v1.zero();
	v2 = n; v2.one();
	while (TRUE) {
		if (f_v) {
			cout << "loop:" << endl;
			cout << "M=" << M << "(" << kind_ascii(M.s_kind()) << ") N=" 
				<< N << "(" << kind_ascii(N.s_kind()) << ")" << endl;
			cout << "before integral_division" << endl;
			}
		M.integral_division(N, Q, R, verbose_level);
		if (f_v) {
			cout << "after integral_division" << endl;
			cout << "Q=" << Q << " R=" << R << endl;
			}
		if (R.is_zero()) {
			break;
			}
		// u3 := u1 - Q * u2
		u3 = u2;
		u3 *= Q;
		u3.negate();
		u3 += u1;
		
		// v3 := v1 - Q * v2
		v3 = v2;
		v3 *= Q;
		v3.negate();
		v3 += v1;
		
		M = N;
		N = R;
		u1 = u2;
		u2 = u3;
		v1 = v2;
		v2 = v3;
		}
	u = u2;
	v = v2;
	g = N;
	if (s_kind() == INTEGER) {
		// cout << "sign1=" << sign1 << endl;
		// cout << "sign2=" << sign2 << endl;
		INT a;
		
		a = s_i_i();
		a *= sign1.s_i_i();
		m_i_i(a);

		a = u.s_i_i();
		a *= sign1.s_i_i();
		u.m_i_i(a);
		
		a = n.s_i_i();
		a *= sign2.s_i_i();
		n.m_i_i(a);
		
		a = v.s_i_i();
		a *= sign2.s_i_i();
		v.m_i_i(a);
		
		// *this *= sign1;
		// u *= sign1;
		// n *= sign2;
		// v *= sign2;
		}
	if (f_v) {
		cout << "g=" << g << " =" << u << " * " << *this << " + " << v << " * " << n << endl;
		}
}

void base::write_memory(memory &m, INT debug_depth)
{
	enum kind k;
	INT i;
	char c;
	
	k = s_kind();
	i = (INT) k;
	c = (char) k;
	if (!ONE_BYTE_INT(i)) {
		cout << "write_memory(): kind not 1 byte" << endl;
		exit(1);
		}
	m.write_char(c);
	if (debug_depth > 0) {
		cout << "base::write_memory() object of kind = " << kind_ascii(k) << endl;
		}
	switch (k) {
		case BASE:
			break;
		case INTEGER:
			m.write_int(s_i_i());
			break;
		case VECTOR:
			as_vector().write_mem(m, debug_depth);
			break;
		case NUMBER_PARTITION:
			as_number_partition().write_mem(m, debug_depth);
			break;
		case PERMUTATION:
			as_permutation().write_mem(m, debug_depth);
			break;
		case MATRIX:
			as_matrix().write_mem(m, debug_depth);
			break;
		case LONGINTEGER:
			// as_longinteger().write_mem(m, debug_depth);
			cout << "base::write_mem() no write_mem for LONGINTEGER" << endl;
			break;
		case MEMORY:
			as_memory().write_mem(m, debug_depth);
			break;
		case HOLLERITH:
			as_hollerith().write_mem(m, debug_depth);
			break;
		//case PERM_GROUP:
			//as_perm_group().write_mem(m, debug_depth);
			//break;
		//case PERM_GROUP_STAB_CHAIN:
			//as_perm_group_stab_chain().write_mem(m, debug_depth);
			//break;
		case UNIPOLY:
			as_unipoly().write_mem(m, debug_depth);
			break;
		case SOLID:
			as_solid().write_mem(m, debug_depth);
			break;
		case BITMATRIX:
			as_bitmatrix().write_mem(m, debug_depth);
			break;
		//case PC_PRESENTATION:
			//as_pc_presentation().write_mem(m, debug_depth);
			//break;
		//case PC_SUBGROUP:
			//as_pc_subgroup().write_mem(m, debug_depth);
			//break;
		//case GROUP_WORD:
			//as_group_word().write_mem(m, debug_depth);
			//break;
		//case GROUP_TABLE:
			//as_group_table().write_mem(m, debug_depth);
			//break;
#if 0
		case ACTION:
			as_action().write_mem(m, debug_depth);
			break;
#endif
		case GEOMETRY:
			as_geometry().write_mem(m, debug_depth);
			break;
		case GROUP_SELECTION:
			as_group_selection().write_mem(m, debug_depth);
			break;
		case DESIGN_PARAMETER:
			as_design_parameter().write_mem(m, debug_depth);
			break;
		case DESIGN_PARAMETER_SOURCE:
			as_design_parameter_source().write_mem(m, debug_depth);
			break;
		default:
			cout << "base::write_memory() no write_mem for " << kind_ascii(k) << endl;
			exit(1);
		}
}

void base::read_memory(memory &m, INT debug_depth)
{
	enum kind k;
	INT i;
	char c;
	
	m.read_char(&c);
	k = (enum kind) c;
	c_kind(k);
	switch (k) {
		case BASE:
			break;
		case INTEGER:
			m.read_int(&i);
			m_i_i(i);
			break;
		case VECTOR:
			as_vector().read_mem(m, debug_depth);
			break;
		case NUMBER_PARTITION:
			as_number_partition().read_mem(m, debug_depth);
			break;
		case PERMUTATION:
			as_permutation().read_mem(m, debug_depth);
			break;
		case MATRIX:
			as_matrix().read_mem(m, debug_depth);
			break;
		case LONGINTEGER:
			// as_longinteger().read_mem(m, debug_depth);
			cout << "base::read_mem() no read_mem for LONGINTEGER" << endl;
			break;
		case MEMORY:
			as_memory().read_mem(m, debug_depth);
			break;
		case HOLLERITH:
			as_hollerith().read_mem(m, debug_depth);
			break;
		//case PERM_GROUP:
			//as_perm_group().read_mem(m, debug_depth);
			//break;
		//case PERM_GROUP_STAB_CHAIN:
			//as_perm_group_stab_chain().read_mem(m, debug_depth);
			//break;
		case UNIPOLY:
			as_unipoly().read_mem(m, debug_depth);
			break;
		case SOLID:
			as_vector().read_mem(m, debug_depth);
			break;
		case BITMATRIX:
			as_bitmatrix().read_mem(m, debug_depth);
			break;
		//case PC_PRESENTATION:
			//as_pc_presentation().read_mem(m, debug_depth);
			//break;
		//case PC_SUBGROUP:
			//as_pc_subgroup().read_mem(m, debug_depth);
			//break;
		//case GROUP_WORD:
			//as_group_word().read_mem(m, debug_depth);
			//break;
		//case GROUP_TABLE:
			//as_group_table().read_mem(m, debug_depth);
			//break;
#if 0
		case ACTION:
			as_action().read_mem(m, debug_depth);
			break;
#endif
		case GEOMETRY:
			as_geometry().read_mem(m, debug_depth);
			break;
		case GROUP_SELECTION:
			as_group_selection().read_mem(m, debug_depth);
			break;
		case DESIGN_PARAMETER:
			as_design_parameter().read_mem(m, debug_depth);
			break;
		case DESIGN_PARAMETER_SOURCE:
			as_design_parameter_source().read_mem(m, debug_depth);
			break;
		default:
			cout << "base::read_memory() no read_mem for " << kind_ascii(k) << endl;
			exit(1);
		}
}

INT base::calc_size_on_file()
{
	enum kind k;
	INT i, size;
	char c;
	
	k = s_kind();
	i = (INT) k;
	c = (char) k;
	if (!ONE_BYTE_INT(i)) {
		cout << "write_memory(): kind not 1 byte" << endl;
		exit(1);
		}
	size = 1;
	switch (k) {
		case BASE:
			break;
		case INTEGER:
			size += 4;
			break;
		case VECTOR:
			size += as_vector().csf();
			break;
		case NUMBER_PARTITION:
			size += as_number_partition().csf();
			break;
		case PERMUTATION:
			size += as_permutation().csf();
			break;
		case MATRIX:
			size += as_matrix().csf();
			break;
		case LONGINTEGER:
			// size += as_longinteger().csf();
			cout << "base::write_mem() no csf for LONGINTEGER" << endl;
			break;
		case MEMORY:
			size += as_memory().csf();
			break;
		case HOLLERITH:
			size += as_hollerith().csf();
			break;
		//case PERM_GROUP:
			//size += as_perm_group().csf();
			//break;
		//case PERM_GROUP_STAB_CHAIN:
			//size += as_perm_group_stab_chain().csf();
			//break;
		case UNIPOLY:
			size += as_unipoly().csf();
			break;
		case SOLID:
			size += as_vector().csf();
			break;
		case BITMATRIX:
			size += as_bitmatrix().csf();
			break;
		//case PC_PRESENTATION:
			//size += as_pc_presentation().csf();
			//break;
		//case PC_SUBGROUP:
			//size += as_pc_subgroup().csf();
			//break;
		//case GROUP_WORD:
			//size += as_group_word().csf();
			//break;
		//case GROUP_TABLE:
			//size += as_group_table().csf();
			//break;
#if 0
		case ACTION:
			size += as_action().csf();
			break;
#endif
		case GEOMETRY:
			size += as_geometry().csf();
			break;
		case GROUP_SELECTION:
			size += as_group_selection().csf();
			break;
		case DESIGN_PARAMETER:
			size += as_design_parameter().csf();
			break;
		case DESIGN_PARAMETER_SOURCE:
			size += as_design_parameter_source().csf();
			break;
		default:
			cout << "base::calc_size_on_file() no csf() for " << kind_ascii(k) << endl;
			exit(1);
		}
	return size;
}

void base::pack(memory & M, INT f_v, INT debug_depth)
// used to pack (i.e. to serialize) objects into (binary) strings in memory objects.
{
	INT size, size0;
	
	if (f_v) {
		cout << "base::pack(): calculating memory size" << endl;
		}
	size0 = calc_size_on_file();
	// M.init(0, NULL);
	if (f_v) {
		cout << "base::pack(): allocating memory of size " << size0 << endl;
		}
	M.alloc(size0);
	M.used_length() = 0;
	if (f_v) {
		cout << "base::pack(): calling write_memory()" << endl;
		}
	write_memory(M, debug_depth);
	size = M.used_length();
	if (size != size0) {
		cout << "base::pack(): WARNING!!!: size = " << size << " != size0 = " << size0 << endl;
		}
}

void base::unpack(memory & M, INT f_v, INT debug_depth)
// unpacks an object from a binary representation in a memory object
{
	read_memory(M, debug_depth);
}

void base::save_ascii(ostream & f)
// writes in ASCII text format (uuencoded like) into the stream f. 
{
	memory M;
	INT f_v = FALSE, f_vv = FALSE;
	INT size, debug_depth;
	INT i;
	UINT a, a1, a2;
	UBYTE *pc, c1, c2;

	if (f_v) {
		cout << "base::save_ascii(): calculating memory size" << endl;
		}
	if (f_vv)
		debug_depth = 1;
	else
		debug_depth = 0;
	if (f_v) {
		cout << "base::save_ascii(): packing object" << endl;
		}
	pack(M, f_v, debug_depth);
#ifdef SAVE_ASCII_USE_COMPRESS
	if (f_v) {
		cout << "base::save_ascii(): compressing object" << endl;
		}
	M.compress(f_v);
#endif
	if (f_v) {
		cout << "base::save_ascii(): saving data" << endl;
		}
	size = M.used_length();
	pc = (UBYTE *) M.self.char_pointer;
	
	f << "ASCII " << size << endl;
	for (i = 0; i < size; i++) {
		a = (UINT) pc[i];
		a1 = a % (UINT) 16;
		a2 = a >> 4;
		c1 = '0' + a1;
		c2 = '0' + a2;
		f << c1 << c2;
		if ((i + 1) % 40 == 0)
			f << endl;
		}
	f << endl << "ASCIIEND" << endl;
}

//#define BUFSIZE 10000

void base::load_ascii(istream & f)
// reads ASCII style objects written with save-ascii
{
	memory M;
	BYTE buf[BUFSIZE];
	BYTE str[1024], *p;
	INT f_v = TRUE;
	INT f_vv = FALSE;
	INT size, i, debug_depth;
	UBYTE *pc;
	UBYTE c;
	INT a;
	UINT a1, a2;
	char cc;
		
	f.getline(buf, sizeof(buf));
	p = buf;
	s_scan_token(&p, str);
	if (strcmp(str, "ASCII") != 0) {
		cout << "base::load_ascii(): error reading header: ASCII keyword not found" << endl;
		exit(1);
		}
	s_scan_int(&p, &size);
	if (f_v) {
		cout << "base::load_ascii(): reading ASCII file of size " << size << endl;
		}
	M.alloc(size);
	pc = (UBYTE *) M.self.char_pointer;
	for (i = 0; i < size; i++) {
		while (TRUE) {
			if (f.eof()) {
				cout << "base::load_ascii() primature EOF" << endl;
				exit(1);
				}
			f >> cc;
			if (cc == '\n')
				continue;
			break;
			}
		a1 = (UINT) cc;
		if (f.eof()) {
			cout << "base::load_ascii() primature EOF" << endl;
			exit(1);
			}
		f >> cc;
		a2 = (UINT) cc;
		a1 = a1 - '0';
		a2 = a2 - '0';
		a = a2 << 4;
		a += a1;
		c = (UBYTE) a;
		pc[i] = c; 
		}
#if 1
	while (TRUE) {
		f.getline(buf, sizeof(buf));
		if (strlen(buf))
			break;
		// if (buf[0] != '\n') break;
		}
#endif
	// f.getline(buf, sizeof(buf));
	// cout << "base::load_ascii(): buf = " << buf << endl;
	p = buf;
	s_scan_token(&p, str);
	if (strcmp(str, "ASCIIEND") != 0) {
		cout << "base::load_ascii(): error reading footer: ASCIIEND keyword not found" << endl;
		exit(1);
		}

			
	
	if (f_v) {
		cout << "file read." << endl;
		}
	M.used_length() = size;
#ifdef SAVE_ASCII_USE_COMPRESS
	M.decompress(TRUE /* f_verbose */);
#endif
	M.cur_pointer() = 0;
	if (f_vv)
		debug_depth = 1;
	else
		debug_depth = 0;
	unpack(M, f_v, debug_depth);
}

void base::save_file(char *fname)
// writes in ASCII text format (uuencoded like) into the file.
{
	ofstream f(fname);
	save_ascii(f);
}
 
void base::load_file(char *fname)
// read in ASCII text format (uuencoded like) from the file.
{
	ifstream f(fname);
	load_ascii(f);
}
 

