// mp_graphics.C
//
// Anton Betten
// March 6, 2003

#include "galois.h"

INT mp_graphics::cntr_new = 0;
INT mp_graphics::cntr_objects = 0;
INT mp_graphics::f_debug_memory = FALSE;

void *mp_graphics::operator new(size_t bytes)
{
	cntr_new++;
	cntr_objects++;
	if (f_debug_memory) {
		cout << "mp_graphics::operator new bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	return malloc(bytes);
}

void *mp_graphics::operator new[](size_t bytes)
{
	INT n;
	
	n = bytes / sizeof(mp_graphics);
	cntr_new++;
	cntr_objects += n;
	if (f_debug_memory) {
		cout << "mp_graphics::operator new[] n=" << n 
			<< " bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	return malloc(bytes);
}

void mp_graphics::operator delete(void *ptr, size_t bytes)
{
	if (f_debug_memory) {
		cout << "mp_graphics::operator delete bytes=" << bytes 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects--;
	return free(ptr);
}

void mp_graphics::operator delete[](void *ptr, size_t bytes)
{
	INT n;
	
	n = bytes / sizeof(mp_graphics);
	if (f_debug_memory) {
		cout << "mp_graphics::operator delete[] n=" << n 
			<< " cntr_new=" << cntr_new 
			<< " cntr_objects=" << cntr_objects 
			<< endl;
		}
	cntr_new--;
	cntr_objects -= n;
	return free(ptr);
}

mp_graphics::mp_graphics()
{
	default_values();
}

mp_graphics::mp_graphics(char *file_name, INT xmin, INT ymin, INT xmax, INT ymax, INT f_embedded, INT f_sideways)
{
	
	default_values();
	init(file_name, xmin, ymin, xmax, ymax, f_embedded, f_sideways);
}

mp_graphics::~mp_graphics()
{
	// cout << "mp_graphics::~mp_graphics()\n";
	exit(cout, FALSE);
	//label_pos.freeself();
	//label_idx.freeself();
}

void mp_graphics::default_values()
{
	f_file_open = FALSE;
	f_min_max_set = FALSE;
	txt_halign = 0;
	txt_valign = 0;
	txt_boxed = 0;
	txt_overwrite = 0;
	txt_rotate = 0;
	line_beg_style = 0;
	line_end_style = 0; // if 1, draw an arrow
	fill_interior = 0;
	fill_color = 0;
	fill_shape = 1;
	fill_outline = 0;
	fill_nofill = 0;
	line_dashing = 0; // between 0 and 100  sl_udsty
	line_thickness = 100; // 100 is the old 1
	
	cur_path = 1;
	
	f_embedded = FALSE;
	f_sideways = FALSE;

	tikz_global_scale = .45;
	tikz_global_line_width = 1.5;
		
	//label_pos.m_l(0);
	//label_idx.m_l(0);
}

void mp_graphics::init(char *file_name, INT xmin, INT ymin, INT xmax, INT ymax, 
	INT f_embedded, INT f_sideways)
{
	mp_graphics::f_embedded = f_embedded;
	mp_graphics::f_sideways = f_sideways;
	
	strcpy(fname_mp, file_name);
	
	get_fname_base(file_name, fname_base);
	sprintf(fname_log, "%s.commands", fname_base);
	sprintf(fname_tikz, "%s.tex", fname_base);

	fp_mp.open(fname_mp);
	fp_log.open(fname_log);
	fp_tikz.open(fname_tikz);
	f_file_open = TRUE;
	
	user[0] = xmin;
	user[1] = ymin;
	user[2] = xmax;
	user[3] = ymax;
	
	// the identity transformation:
	
	dev[0] = xmin;
	dev[1] = ymin;
	dev[2] = xmax;
	dev[3] = ymax;
}

void mp_graphics::exit(ostream &ost, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	
	if (f_file_open) {
		fp_mp.close();
		fp_log.close();
		fp_tikz.close();
		if (f_v) {
			ost << "written file " << fname_mp << " of size " << file_size(fname_mp) << endl;
			ost << "written file " << fname_log << " of size " << file_size(fname_log) << endl;
			ost << "written file " << fname_tikz << " of size " << file_size(fname_tikz) << endl;
			}
		f_file_open = FALSE;
		}
}

void mp_graphics::setup(const char *fname_base, 
	INT in_xmin, INT in_ymin, INT in_xmax, INT in_ymax, 
	INT xmax, INT ymax, INT f_embedded, INT f_sideways, 
	double scale, double line_width)
{
	//INT x_min = 0, x_max = 1000;
	//INT y_min = 0, y_max = 1000;
	INT factor_1000 = 1000;
	BYTE fname_full[1000];
	
	sprintf(fname_full, "%s.mp", fname_base);
	init(fname_full, in_xmin, in_ymin, in_xmax, in_ymax, f_embedded, f_sideways);
#if 0
	out_xmin() = -(xmax >> 1);
	out_ymin() = -(ymax >> 1);
	out_xmax() = xmax >> 1;
	out_ymax() = ymax >> 1;
#else
	out_xmin() = 0;
	out_ymin() = 0;
	out_xmax() = xmax;
	out_ymax() = ymax;
#endif
	//cout << "xmax/ymax = " << xmax << " / " << ymax << endl;
	
	tikz_global_scale = scale;
	tikz_global_line_width = line_width;
	header();
	begin_figure(factor_1000);
}

void mp_graphics::set_parameters(double scale, double line_width)
{
	tikz_global_scale = scale;
	tikz_global_line_width = line_width;
}

void mp_graphics::set_scale(double scale)
{
	tikz_global_scale = scale;
}


void mp_graphics::frame(double move_out)
{
	INT Px[30];
	INT Py[30];
	INT x, y, dx, dy, ox, oy, Ox, Oy;
	INT xmin, ymin, xmax, ymax;
	
	xmin = out_xmin();
	xmax = out_xmax();
	ymin = out_ymin();
	ymax = out_ymax();
	
	dx = xmax - xmin;
	dy = ymax - ymin;
	
	ox = (INT)(dx * .05);
	oy = (INT)(dy * .05);
	Ox = (INT)(dx * move_out);
	Oy = (INT)(dy * move_out);
	xmin -= Ox;
	xmax += Ox;
	ymin -= Oy;
	ymax += Oy;

	x = xmin;
	y = ymin;
	dev2user(x, y);
	Px[0] = x;
	Py[0] = y;
	//cout << "Px[0]=" << Px[0] << endl;
	//cout << "Py[0]=" << Py[0] << endl;
	x = xmin + ox;
	y = ymin;
	dev2user(x, y);
	Px[1] = x;
	Py[1] = y;
	x = xmin;
	y = ymin + oy;
	dev2user(x, y);
	Px[2] = x;
	Py[2] = y;
	
	x = xmax;
	y = ymin;
	dev2user(x, y);
	Px[3] = x;
	Py[3] = y;
	x = xmax - ox;
	y = ymin;
	dev2user(x, y);
	Px[4] = x;
	Py[4] = y;
	x = xmax;
	y = ymin + oy;
	dev2user(x, y);
	Px[5] = x;
	Py[5] = y;
	
	x = xmax;
	y = ymax;
	dev2user(x, y);
	Px[6] = x;
	Py[6] = y;
	//cout << "Px[6]=" << Px[6] << endl;
	//cout << "Py[6]=" << Py[6] << endl;
	x = xmax - ox;
	y = ymax;
	dev2user(x, y);
	Px[7] = x;
	Py[7] = y;
	x = xmax;
	y = ymax - oy;
	dev2user(x, y);
	Px[8] = x;
	Py[8] = y;
	
	x = xmin;
	y = ymax;
	dev2user(x, y);
	Px[9] = x;
	Py[9] = y;
	x = xmin + ox;
	y = ymax;
	dev2user(x, y);
	Px[10] = x;
	Py[10] = y;
	x = xmin;
	y = ymax - oy;
	dev2user(x, y);
	Px[11] = x;
	Py[11] = y;

	polygon3(Px, Py, 2, 0, 1);
	polygon3(Px, Py, 4, 3, 5);
	polygon3(Px, Py, 7, 6, 8);
	polygon3(Px, Py, 11, 9, 10);
}

void mp_graphics::frame_constant_aspect_ratio(double move_out)
{
	INT Px[30];
	INT Py[30];
	INT x, y, dx, dy, ox, oy, Ox, Oy;
	INT xmin, ymin, xmax, ymax;
	
	xmin = out_xmin();
	xmax = out_xmax();
	ymin = out_ymin();
	ymax = out_ymax();
	
	cout << "mp_graphics::frame_constant_aspect_ratio:" << endl;
	cout << xmin << "," << xmax << "," << ymin << "," << ymax << endl;
	dx = xmax - xmin;
	dy = ymax - ymin;
	
	INT adjust_x = 0;
	INT adjust_y = 0;
	
	if (dx > dy) {
		adjust_y = (dx - dy) >> 1;
		cout << "mp_graphics::frame_constant_aspect_ratio adjust_y=" << adjust_y << endl;
		}
	else {
		adjust_x = (dy - dx) >> 1;
		cout << "mp_graphics::frame_constant_aspect_ratio adjust_x=" << adjust_x << endl;
		}
	xmin -= adjust_x;
	xmax += adjust_x;
	ymin -= adjust_y;
	ymax += adjust_y;
	cout << "mp_graphics::frame_constant_aspect_ratio after adjustment:" << endl;
	cout << xmin << "," << xmax << "," << ymin << "," << ymax << endl;
	dx = xmax - xmin;
	dy = ymax - ymin;
	
	ox = (INT)(dx * .05);
	oy = (INT)(dy * .05);
	Ox = (INT)(dx * move_out);
	Oy = (INT)(dy * move_out);
	xmin -= Ox;
	xmax += Ox;
	ymin -= Oy;
	ymax += Oy;

	x = xmin;
	y = ymin;
	dev2user(x, y);
	Px[0] = x;
	Py[0] = y;
	//cout << "Px[0]=" << Px[0] << endl;
	//cout << "Py[0]=" << Py[0] << endl;
	x = xmin + ox;
	y = ymin;
	dev2user(x, y);
	Px[1] = x;
	Py[1] = y;
	x = xmin;
	y = ymin + oy;
	dev2user(x, y);
	Px[2] = x;
	Py[2] = y;
	
	x = xmax;
	y = ymin;
	dev2user(x, y);
	Px[3] = x;
	Py[3] = y;
	x = xmax - ox;
	y = ymin;
	dev2user(x, y);
	Px[4] = x;
	Py[4] = y;
	x = xmax;
	y = ymin + oy;
	dev2user(x, y);
	Px[5] = x;
	Py[5] = y;
	
	x = xmax;
	y = ymax;
	dev2user(x, y);
	Px[6] = x;
	Py[6] = y;
	//cout << "Px[6]=" << Px[6] << endl;
	//cout << "Py[6]=" << Py[6] << endl;
	x = xmax - ox;
	y = ymax;
	dev2user(x, y);
	Px[7] = x;
	Py[7] = y;
	x = xmax;
	y = ymax - oy;
	dev2user(x, y);
	Px[8] = x;
	Py[8] = y;
	
	x = xmin;
	y = ymax;
	dev2user(x, y);
	Px[9] = x;
	Py[9] = y;
	x = xmin + ox;
	y = ymax;
	dev2user(x, y);
	Px[10] = x;
	Py[10] = y;
	x = xmin;
	y = ymax - oy;
	dev2user(x, y);
	Px[11] = x;
	Py[11] = y;

	polygon3(Px, Py, 2, 0, 1);
	polygon3(Px, Py, 4, 3, 5);
	polygon3(Px, Py, 7, 6, 8);
	polygon3(Px, Py, 11, 9, 10);
}

void mp_graphics::finish(ostream &ost, INT verbose_level)
{
	draw_boxes_final();
	end_figure();
	footer();
	exit(cout, verbose_level - 1);
}

INT& mp_graphics::out_xmin()
{
	return dev[0];
}

INT& mp_graphics::out_ymin()
{
	return dev[1];
}

INT& mp_graphics::out_xmax()
{
	return dev[2];
}

INT& mp_graphics::out_ymax()
{
	return dev[3];
}

void mp_graphics::user2dev(INT &x, INT &y)
{
	transform_llur(user, dev, x, y);
}

void mp_graphics::dev2user(INT &x, INT &y)
{
	transform_llur(dev, user, x, y);
}

void mp_graphics::user2dev_dist(INT &x, INT &y)
{
	transform_dist(user, dev, x, y);
}

void mp_graphics::user2dev_dist_x(INT &x)
{
	transform_dist_x(user, dev, x);
}

void mp_graphics::user2dev_dist_y(INT &y)
{
	transform_dist_y(user, dev, y);
}

void mp_graphics::dev2user_dist(INT &x, INT &y)
{
	transform_dist(dev, user, x, y);
}

void mp_graphics::nice_circle(INT x, INT y, INT rad)
{
	fp_log << "NiceCircle " << x << " " << y << " " << rad << endl;

	sf_interior(100);
	sf_color(0); // 1 = black, 0 = white
	circle(x, y, rad);
	sf_interior(0);
	sf_color(1); // 1 = black, 0 = white
	circle(x, y, rad);
}

void mp_graphics::grid_polygon2(grid_frame *F, 
	INT x0, INT y0, INT x1, INT y1)
{
	INT *Px, *Py, *Idx;
	INT i;

	Px = NEW_INT(2);
	Py = NEW_INT(2);
	Idx = NEW_INT(2);

	for (i = 0; i < 2; i++) {
		Idx[i] = i;
		}
	if (F->f_matrix_notation) {
		Px[0] = (INT)(F->origin_x + y0 * F->dx);
		Py[0] = (INT)(F->origin_y + (F->m - x0) * F->dy);
		Px[1] = (INT)(F->origin_x + y1 * F->dx);
		Py[1] = (INT)(F->origin_y + (F->m - x1) * F->dy);
		}
	else {
		Px[0] = (INT)(F->origin_x + x0 * F->dx);
		Py[0] = (INT)(F->origin_y + y0 * F->dy);
		Px[1] = (INT)(F->origin_x + x1 * F->dx);
		Py[1] = (INT)(F->origin_y + y1 * F->dy);
		}
	polygon_idx(Px, Py, Idx, 2);
	FREE_INT(Px);
	FREE_INT(Py);
	FREE_INT(Idx);
}

void mp_graphics::grid_polygon4(grid_frame *F, 
	INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3)
{
	INT *Px, *Py, *Idx;
	INT i;

	Px = NEW_INT(4);
	Py = NEW_INT(4);
	Idx = NEW_INT(4);

	for (i = 0; i < 4; i++) {
		Idx[i] = i;
		}
	if (F->f_matrix_notation) {
		Px[0] = (INT)(F->origin_x + y0 * F->dx);
		Py[0] = (INT)(F->origin_y + (F->m - x0) * F->dy);
		Px[1] = (INT)(F->origin_x + y1 * F->dx);
		Py[1] = (INT)(F->origin_y + (F->m - x1) * F->dy);
		Px[2] = (INT)(F->origin_x + y2 * F->dx);
		Py[2] = (INT)(F->origin_y + (F->m - x2) * F->dy);
		Px[3] = (INT)(F->origin_x + y3 * F->dx);
		Py[3] = (INT)(F->origin_y + (F->m - x3) * F->dy);
		}
	else {
		Px[0] = (INT)(F->origin_x + x0 * F->dx);
		Py[0] = (INT)(F->origin_y + y0 * F->dy);
		Px[1] = (INT)(F->origin_x + x1 * F->dx);
		Py[1] = (INT)(F->origin_y + y1 * F->dy);
		Px[2] = (INT)(F->origin_x + x2 * F->dx);
		Py[2] = (INT)(F->origin_y + y2 * F->dy);
		Px[3] = (INT)(F->origin_x + x3 * F->dx);
		Py[3] = (INT)(F->origin_y + y3 * F->dy);
		}
	polygon_idx(Px, Py, Idx, 4);
	FREE_INT(Px);
	FREE_INT(Py);
	FREE_INT(Idx);
}

void mp_graphics::grid_polygon5(grid_frame *F, 
	INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3, INT x4, INT y4)
{
	INT *Px, *Py, *Idx;
	INT i;

	Px = NEW_INT(5);
	Py = NEW_INT(5);
	Idx = NEW_INT(5);

	for (i = 0; i < 5; i++) {
		Idx[i] = i;
		}
	if (F->f_matrix_notation) {
		Px[0] = (INT)(F->origin_x + y0 * F->dx);
		Py[0] = (INT)(F->origin_y + (F->m - x0) * F->dy);
		Px[1] = (INT)(F->origin_x + y1 * F->dx);
		Py[1] = (INT)(F->origin_y + (F->m - x1) * F->dy);
		Px[2] = (INT)(F->origin_x + y2 * F->dx);
		Py[2] = (INT)(F->origin_y + (F->m - x2) * F->dy);
		Px[3] = (INT)(F->origin_x + y3 * F->dx);
		Py[3] = (INT)(F->origin_y + (F->m - x3) * F->dy);
		Px[4] = (INT)(F->origin_x + y4 * F->dx);
		Py[4] = (INT)(F->origin_y + (F->m - x4) * F->dy);
		}
	else {
		Px[0] = (INT)(F->origin_x + x0 * F->dx);
		Py[0] = (INT)(F->origin_y + y0 * F->dy);
		Px[1] = (INT)(F->origin_x + x1 * F->dx);
		Py[1] = (INT)(F->origin_y + y1 * F->dy);
		Px[2] = (INT)(F->origin_x + x2 * F->dx);
		Py[2] = (INT)(F->origin_y + y2 * F->dy);
		Px[3] = (INT)(F->origin_x + x3 * F->dx);
		Py[3] = (INT)(F->origin_y + y3 * F->dy);
		Px[4] = (INT)(F->origin_x + x4 * F->dx);
		Py[4] = (INT)(F->origin_y + y4 * F->dy);
		}
	polygon_idx(Px, Py, Idx, 5);
	FREE_INT(Px);
	FREE_INT(Py);
	FREE_INT(Idx);
}

void mp_graphics::polygon(INT *Px, INT *Py, INT n)
{
	INT *Idx = NEW_INT(n);
	INT i;
	
	fp_log << "Polygon " << n;
	for (i = 0; i < n; i++) {
		fp_log << " " << Px[i] << " " << Py[i];
		}
	fp_log << endl;
	for (i = 0; i < n; i++) {
		Idx[i] = i;
		}
	polygon_idx(Px, Py, Idx, n);
	FREE_INT(Idx);
}

void mp_graphics::polygon2(INT *Px, INT *Py, INT i1, INT i2)
{
	INT Idx[2];
	Idx[0] = i1;
	Idx[1] = i2;
	polygon_idx(Px, Py, Idx, 2);
}

void mp_graphics::polygon3(INT *Px, INT *Py, INT i1, INT i2, INT i3)
{
	INT Idx[3];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	polygon_idx(Px, Py, Idx, 3);
}

void mp_graphics::polygon4(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4)
{
	INT Idx[4];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	polygon_idx(Px, Py, Idx, 4);
}

void mp_graphics::polygon5(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5)
{
	INT Idx[5];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;

	polygon_idx(Px, Py, Idx, 5);
}

void mp_graphics::polygon6(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	polygon_idx(Px, Py, Idx, 6);
}

void mp_graphics::polygon7(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	polygon_idx(Px, Py, Idx, 7);
}

void mp_graphics::polygon8(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	polygon_idx(Px, Py, Idx, 8);
}

void mp_graphics::polygon9(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	Idx[8] = i9;
	polygon_idx(Px, Py, Idx, 9);
}

void mp_graphics::polygon10(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10)
{
	INT Idx[20];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	Idx[8] = i9;
	Idx[9] = i10;
	polygon_idx(Px, Py, Idx, 10);
}

void mp_graphics::polygon11(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10, INT i11)
{
	INT Idx[20];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	Idx[8] = i9;
	Idx[9] = i10;
	Idx[10] = i11;
	polygon_idx(Px, Py, Idx, 11);
}

void mp_graphics::polygon_idx(INT *Px, INT *Py, INT *Idx, INT n)
{
	INT i;
	
	fp_log << "Polygon " << n;
	for (i = 0; i < n; i++) {
		fp_log << " " << Px[Idx[i]] << " " << Py[Idx[i]];
		}
	fp_log << endl;
	polygon_or_bezier_idx(Px, Py, Idx, n, "--", FALSE);
}

void mp_graphics::bezier(INT *Px, INT *Py, INT n)
{
	INT *Idx = NEW_INT(n);
	INT i;
	
	for (i = 0; i < n; i++) {
		Idx[i] = i;
		}
	bezier_idx(Px, Py, Idx, n);
	FREE_INT(Idx);
}

void mp_graphics::bezier2(INT *Px, INT *Py, INT i1, INT i2)
{
	INT Idx[2];
	Idx[0] = i1;
	Idx[1] = i2;
	bezier_idx(Px, Py, Idx, 2);
}

void mp_graphics::bezier3(INT *Px, INT *Py, INT i1, INT i2, INT i3)
{
	INT Idx[3];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	bezier_idx(Px, Py, Idx, 3);
}

void mp_graphics::bezier4(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4)
{
	INT Idx[4];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	bezier_idx(Px, Py, Idx, 4);
}

void mp_graphics::bezier5(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5)
{
	INT Idx[5];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	bezier_idx(Px, Py, Idx, 5);
}

void mp_graphics::bezier6(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6)
{
	INT Idx[6];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	bezier_idx(Px, Py, Idx, 6);
}

void mp_graphics::bezier7(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7)
{
	INT Idx[7];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	bezier_idx(Px, Py, Idx, 7);
}

void mp_graphics::bezier_idx(INT *Px, INT *Py, INT *Idx, INT n)
{
	INT i;
	
	fp_log << "Bezier " << n;
	for (i = 0; i < n; i++) {
		fp_log << " " << Px[Idx[i]] << " " << Py[Idx[i]];
		}
	fp_log << endl;
	polygon_or_bezier_idx(Px, Py, Idx, n, "..", FALSE);
}

void mp_graphics::grid_fill_polygon4(grid_frame *F, 
	INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3)
{
	INT *Px, *Py, *Idx;
	INT i;

	Px = NEW_INT(4);
	Py = NEW_INT(4);
	Idx = NEW_INT(4);

	for (i = 0; i < 4; i++) {
		Idx[i] = i;
		}
	if (F->f_matrix_notation) {
		Px[0] = (INT)(F->origin_x + y0 * F->dx);
		Py[0] = (INT)(F->origin_y + (F->m - x0) * F->dy);
		Px[1] = (INT)(F->origin_x + y1 * F->dx);
		Py[1] = (INT)(F->origin_y + (F->m - x1) * F->dy);
		Px[2] = (INT)(F->origin_x + y2 * F->dx);
		Py[2] = (INT)(F->origin_y + (F->m - x2) * F->dy);
		Px[3] = (INT)(F->origin_x + y3 * F->dx);
		Py[3] = (INT)(F->origin_y + (F->m - x3) * F->dy);
		}
	else {
		Px[0] = (INT)(F->origin_x + x0 * F->dx);
		Py[0] = (INT)(F->origin_y + y0 * F->dy);
		Px[1] = (INT)(F->origin_x + x1 * F->dx);
		Py[1] = (INT)(F->origin_y + y1 * F->dy);
		Px[2] = (INT)(F->origin_x + x2 * F->dx);
		Py[2] = (INT)(F->origin_y + y2 * F->dy);
		Px[3] = (INT)(F->origin_x + x3 * F->dx);
		Py[3] = (INT)(F->origin_y + y3 * F->dy);
		}
	fill_idx(Px, Py, Idx, 4, "--", TRUE);
	FREE_INT(Px);
	FREE_INT(Py);
	FREE_INT(Idx);
}

void mp_graphics::grid_fill_polygon5(grid_frame *F, 
	INT x0, INT y0, INT x1, INT y1, INT x2, INT y2, INT x3, INT y3, INT x4, INT y4)
{
	INT *Px, *Py, *Idx;
	INT i;

	Px = NEW_INT(5);
	Py = NEW_INT(5);
	Idx = NEW_INT(5);

	for (i = 0; i < 5; i++) {
		Idx[i] = i;
		}
	if (F->f_matrix_notation) {
		Px[0] = (INT)(F->origin_x + y0 * F->dx);
		Py[0] = (INT)(F->origin_y + (F->m - x0) * F->dy);
		Px[1] = (INT)(F->origin_x + y1 * F->dx);
		Py[1] = (INT)(F->origin_y + (F->m - x1) * F->dy);
		Px[2] = (INT)(F->origin_x + y2 * F->dx);
		Py[2] = (INT)(F->origin_y + (F->m - x2) * F->dy);
		Px[3] = (INT)(F->origin_x + y3 * F->dx);
		Py[3] = (INT)(F->origin_y + (F->m - x3) * F->dy);
		Px[4] = (INT)(F->origin_x + y4 * F->dx);
		Py[4] = (INT)(F->origin_y + (F->m - x4) * F->dy);
		}
	else {
		Px[0] = (INT)(F->origin_x + x0 * F->dx);
		Py[0] = (INT)(F->origin_y + y0 * F->dy);
		Px[1] = (INT)(F->origin_x + x1 * F->dx);
		Py[1] = (INT)(F->origin_y + y1 * F->dy);
		Px[2] = (INT)(F->origin_x + x2 * F->dx);
		Py[2] = (INT)(F->origin_y + y2 * F->dy);
		Px[3] = (INT)(F->origin_x + x3 * F->dx);
		Py[3] = (INT)(F->origin_y + y3 * F->dy);
		Px[4] = (INT)(F->origin_x + x4 * F->dx);
		Py[4] = (INT)(F->origin_y + y4 * F->dy);
		}
	fill_idx(Px, Py, Idx, 5, "--", TRUE);
	FREE_INT(Px);
	FREE_INT(Py);
	FREE_INT(Idx);
}

void mp_graphics::fill_polygon3(INT *Px, INT *Py, INT i1, INT i2, INT i3)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	fill_idx(Px, Py, Idx, 3, "--", FALSE);
}

void mp_graphics::fill_polygon4(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	fill_idx(Px, Py, Idx, 4, "--", FALSE);
}

void mp_graphics::fill_polygon5(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	fill_idx(Px, Py, Idx, 5, "--", FALSE);
}

void mp_graphics::fill_polygon6(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	fill_idx(Px, Py, Idx, 6, "--", FALSE);
}

void mp_graphics::fill_polygon7(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	fill_idx(Px, Py, Idx, 7, "--", FALSE);
}

void mp_graphics::fill_polygon8(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	fill_idx(Px, Py, Idx, 8, "--", FALSE);
}

void mp_graphics::fill_polygon9(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9)
{
	INT Idx[10];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	Idx[8] = i9;
	fill_idx(Px, Py, Idx, 9, "--", FALSE);
}

void mp_graphics::fill_polygon10(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10)
{
	INT Idx[20];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	Idx[8] = i9;
	Idx[9] = i10;
	fill_idx(Px, Py, Idx, 10, "--", FALSE);
}

void mp_graphics::fill_polygon11(INT *Px, INT *Py, INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7, INT i8, INT i9, INT i10, INT i11)
{
	INT Idx[20];
	Idx[0] = i1;
	Idx[1] = i2;
	Idx[2] = i3;
	Idx[3] = i4;
	Idx[4] = i5;
	Idx[5] = i6;
	Idx[6] = i7;
	Idx[7] = i8;
	Idx[8] = i9;
	Idx[9] = i10;
	Idx[10] = i11;
	fill_idx(Px, Py, Idx, 11, "--", FALSE);
}

void mp_graphics::polygon2_arrow_halfway(INT *Px, INT *Py, INT i1, INT i2)
{
	INT Idx[3];
	INT X[3], Y[3];
	
	X[0] = Px[i1];
	X[1] = Px[i2];
	X[2] = (Px[i1] + Px[i2]) >> 1;
	Y[0] = Py[i1];
	Y[1] = Py[i2];
	Y[2] = (Py[i1] + Py[i2]) >> 1;

	sl_ends(0, 0);
	Idx[0] = 0;
	Idx[1] = 1;
	
	polygon_idx(X, Y, Idx, 2);

	sl_ends(0, 1);

	Idx[0] = 0;
	Idx[1] = 2;
	polygon_idx(X, Y, Idx, 2);

	sl_ends(0, 0);
}

void mp_graphics::polygon2_arrow_halfway_and_label(INT *Px, INT *Py, INT i1, INT i2, 
	const BYTE *alignment, const BYTE *txt)
{
	INT Idx[3];
	INT X[3], Y[3];
	
	X[0] = Px[i1];
	X[1] = Px[i2];
	X[2] = (Px[i1] + Px[i2]) >> 1;
	Y[0] = Py[i1];
	Y[1] = Py[i2];
	Y[2] = (Py[i1] + Py[i2]) >> 1;

	sl_ends(0, 0);
	Idx[0] = 0;
	Idx[1] = 1;
	
	polygon_idx(X, Y, Idx, 2);

	sl_ends(0, 1);

	Idx[0] = 0;
	Idx[1] = 2;
	polygon_idx(X, Y, Idx, 2);

	sl_ends(0, 0);
	
	aligned_text(X[2], Y[2], alignment, txt);
}

void mp_graphics::grid_aligned_text(grid_frame *F, INT x, INT y, const char *alignment, const char *p)
{
	INT *Px, *Py;

	Px = NEW_INT(1);
	Py = NEW_INT(1);

	if (F->f_matrix_notation) {
		Px[0] = (INT)(F->origin_x + y * F->dx);
		Py[0] = (INT)(F->origin_y + (F->m - x) * F->dy);
		}
	else {
		Px[0] = (INT)(F->origin_x + x * F->dx);
		Py[0] = (INT)(F->origin_y + y * F->dy);
		}
	aligned_text(Px[0], Py[0], alignment, p);
	FREE_INT(Px);
	FREE_INT(Py);
}

void mp_graphics::aligned_text(INT x, INT y, const char *alignment, const char *p)
{
	fp_log << "AlignedText " << x << " " << y << " " << alignment << " \"" << p << "\"" << endl;
	aligned_text_with_offset(x, y, 0, 0, alignment, p);
}

void mp_graphics::aligned_text_array(INT *Px, INT *Py, INT idx, const char *alignment, const char *p)
{
	aligned_text(Px[idx], Py[idx], alignment, p);
}

void mp_graphics::aligned_text_with_offset(INT x, INT y, INT xoffset, INT yoffset, 
	const char *alignment, const char *p)
{
	INT h_align = 1, v_align = 1;
	INT l, i;
	BYTE c;
	
	fp_log << "AlignedText " << x << " " << y << " " << xoffset << " " << yoffset << " " << alignment << " \"" << p << "\"" << endl;
	l = strlen(alignment);
	for (i = 0; i < l; i++) {
		c = alignment[i];
		if (c == 'r')
			h_align = 2;
		else if (c == 'l')
			h_align = 0;
		else if (c == 'b')
			v_align = 0;
		else if (c == 't')
			v_align = 2;
		else {
			cout << "mp_graphics::aligned_text: unknown alignment character " << c << endl;
			}
		}
	//cout << "xoffset=" << xoffset << endl;
	//cout << "yoffset=" << yoffset << endl;
	//cout << "text=" << p << endl;
	st_alignment(h_align, v_align);
	text(x + xoffset, y + yoffset, p);
}







void mp_graphics::st_alignment(INT txt_halign, INT txt_valign)
{
	mp_graphics::txt_halign = txt_halign;
	mp_graphics::txt_valign = txt_valign;
}

void mp_graphics::get_alignment(BYTE *align)
{
	if (txt_halign == 2) { // right aligned, text to the left of the current position
		if (txt_valign == 2) 
			strcpy(align, ".llft");
		else if (txt_valign == 1) 
			strcpy(align, ".lft");
		else if (txt_valign == 0) 
			strcpy(align, ".ulft");
		}
	else if (txt_halign == 1) { // horizontally centered
		if (txt_valign == 2) 
			strcpy(align, ".bot");
		else if (txt_valign == 1) 
			strcpy(align, "");
		else if (txt_valign == 0) 
			strcpy(align, ".top");
		}
	else if (txt_halign == 0) {
		if (txt_valign == 2) 
			strcpy(align, ".lrt");
		else if (txt_valign == 1) 
			strcpy(align, ".rt");
		else if (txt_valign == 0) 
			strcpy(align, ".urt");
		}
}


void mp_graphics::sl_udsty(INT line_dashing)
{
	mp_graphics::line_dashing = line_dashing;
}

void mp_graphics::sl_ends(INT line_beg_style, INT line_end_style)
{
	mp_graphics::line_beg_style = line_beg_style;
	mp_graphics::line_end_style = line_end_style;
}

void mp_graphics::sl_thickness(INT line_thickness)
{
	double d;
	mp_graphics::line_thickness = line_thickness;

	d = line_thickness * 0.01;
	fp_mp << "pickup pencircle scaled " << d << "pt;" << endl;

	//cout << "mp_graphics::line_thickness = " << mp_graphics::line_thickness << endl;
}

void mp_graphics::sf_interior(INT fill_interior)
{
	mp_graphics::fill_interior = fill_interior;
	//cout << "fill_interior = " << mp_graphics::fill_interior << endl;
}

void mp_graphics::sf_color(INT fill_color)
{
	mp_graphics::fill_color = fill_color;
	//cout << "fill_color = " << mp_graphics::fill_color << endl;
}

void mp_graphics::sf_shape(INT fill_shape)
{
	mp_graphics::fill_shape = fill_shape;
}

void mp_graphics::sf_outline(INT fill_outline)
{
	mp_graphics::fill_outline = fill_outline;
}

void mp_graphics::sf_nofill(INT fill_nofill)
{
	mp_graphics::fill_nofill = fill_nofill;
}

void mp_graphics::st_boxed(INT txt_boxed)
{
	mp_graphics::txt_boxed = txt_boxed;
}

void mp_graphics::st_overwrite(INT txt_overwrite)
{
	mp_graphics::txt_overwrite = txt_overwrite;
}

void mp_graphics::st_rotate(INT txt_rotate)
{
	mp_graphics::txt_rotate = txt_rotate;
}

void mp_graphics::coords_min_max(INT x, INT y)
{
	if (!f_min_max_set) {
		x_min = x_max = x;
		y_min = y_max = y;
		}
	else {
		x_min = MINIMUM(x_min, x);
		y_min = MINIMUM(y_min, y);
		x_max = MAXIMUM(x_max, x);
		y_max = MAXIMUM(y_max, y);
		}
	f_min_max_set = TRUE;
}

INT mp_graphics::get_label(INT x, INT y)
{
	static int i = 0;
	
	return i++;
#if 0
	vector v;
	integer idx_ob;
	INT i, idx;
	
	v.m_l(2);
	v.m_ii(0, x);
	v.m_ii(1, y);
	if (label_pos.search(v, &i)) {
		idx = label_idx.s_ii(i);
		return idx;
		}
	else {
		label_pos.insert_element(i, v);
		idx = label_pos.s_l();
		idx_ob.m_i(idx);
		label_idx.insert_element(i, idx_ob);
		// cout << "mp_graphics::get_label: new label " << idx << " at " << i << endl;
		return idx;
		}
#endif
}

void mp_graphics::draw_boxes_final()
{
#if 0
	INT i, l, idx;
	char str[256];
	char str1[256];
	
	l = label_pos.s_l();
	for (i = 0; i < l; i++) {
		idx = label_idx.s_ii(i);
		sprintf(str1, "%ld", idx);
		f << "unfill bpath " << str << ";" << endl;
		f << "drawboxed(" << str << ");" << endl;
		}
#endif
}


// ####################################################################################
// output commands:
// ####################################################################################

void mp_graphics::header()
{
	BYTE str[1024];
	
	f_min_max_set = FALSE;
	fp_mp << "% file: " << fname_mp << endl;
	fp_mp << "% created by DISCRETA MetaPost interface" << endl;
	fp_log << "% file: " << fname_log << endl;
	fp_log << "% created by DISCRETA MetaPost interface" << endl;
	fp_tikz << "% file: " << fname_tikz << endl;
	fp_tikz << "% created by DISCRETA MetaPost interface" << endl;
	//system("rm a");
	system("date >a");
	{
	ifstream f1("a");
	f1.getline(str, sizeof(str));
	}
	fp_mp << "% creation date: " << str << endl << endl;
	fp_log << "% creation date: " << str << endl << endl;
	fp_tikz << "% creation date: " << str << endl;

	// no extra spaces in tikz mode so that we do not create a line feed.
	// this allows for multiple tikz pictures on the same line


	fp_log << "DeviceCoordinates " << dev[0] << " " << dev[1] << " " << dev[2] << " " << dev[3] << endl;
	fp_log << "UserCoordinates " << user[0] << " " << user[1] << " " << user[2] << " " << user[3] << endl;
	
	fp_tikz << "% DeviceCoordinates " << dev[0] << " " << dev[1] << " " << dev[2] << " " << dev[3] << endl;
	fp_tikz << "% UserCoordinates " << user[0] << " " << user[1] << " " << user[2] << " " << user[3] << endl;
	
	fp_mp << "input boxes" << endl << endl;
	
	fp_mp << "verbatimtex" << endl;
	fp_mp << "\\documentclass[10pt]{article}" << endl;
	fp_mp << "\\usepackage[noBBpl]{mathpazo}" << endl;
	fp_mp << "\\usepackage{amsmath}" << endl;
	fp_mp << "\\usepackage{amssymb}" << endl;
	fp_mp << "\\begin{document}" << endl;
	fp_mp << "\\scriptsize" << endl;
	fp_mp << "etex;" << endl;
#if 0
	fprintf(fp, "vardef  cuta(suffix a,b) expr p =\n");
	fprintf(fp, "  draw p cutbefore bpath.a cutafter bpath.b;\n");
	fprintf(fp, "enddef;\n\n");
#endif


	if (f_embedded) {
		fp_tikz << "\\documentclass[12pt]{article}" << endl;
		fp_tikz << "\\usepackage{amsmath, amssymb, amsthm}" << endl;
		fp_tikz << "\\usepackage{tikz} " << endl;
		if (f_sideways) {
			fp_tikz << "\\usepackage{rotating} " << endl;
			}
		fp_tikz << "%\\usepackage{anysize}" << endl;
		fp_tikz << "\\begin{document}" << endl;
		fp_tikz << "%\\bibliographystyle{plain}" << endl;
		fp_tikz << "\\pagestyle{empty}" << endl;
		}

	if (f_sideways) {
		fp_tikz << "\\begin{sideways}" << endl;
		}
	fp_tikz << "\\begin{tikzpicture}[scale=" << tikz_global_scale << ",line width = " << tikz_global_line_width << "pt]" << endl;
	//fp_tikz << "\\begin{tikzpicture}[scale=.05,line width = 0.5pt]" << endl;



}

void mp_graphics::footer()
{
	fp_mp << "end" << endl << endl;
	
	fp_log << "END_OF_FILE" << endl << endl;

	fp_tikz << "\\end{tikzpicture}" << endl;
	if (f_sideways) {
		fp_tikz << "\\end{sideways}" << endl;
		}
	if (f_embedded) {
		fp_tikz << "\\end{document}" << endl;
		}
}

void mp_graphics::begin_figure(INT factor_1000)
{
	double d;
	BYTE str[1000];
	INT i, l;
	
	d = (double) factor_1000 * 0.001  * 0.1;
	
	//fp_mp << "defaultfont:=\"cmr7\";\n";

	sprintf(str, "%lf", d);
	l = strlen(str);
	for (i = 0; i < l; i++) {
		if (str[i] == ',')
			str[i] = '.';
		}
	fp_mp << "u=" << str << "mm;\n";
	fp_mp << "beginfig(1);" << endl;
	fp_mp << "path p[];" << endl;
}

void mp_graphics::end_figure()
{
	fp_mp << "endfig;" << endl << endl;
}

void mp_graphics::output_xy_metapost(INT x, INT y)
{
	fp_mp << "(";
	output_x_metapost(x);
	fp_mp << ",";
	output_y_metapost(y);
	fp_mp << ")";
}

void mp_graphics::output_x_metapost(INT x)
{
	double d;

	d = (double) x;
	d /= 1000.;
	fp_mp << d << "u ";
}

void mp_graphics::output_y_metapost(INT y)
{
	double d;

	d = (double) y;
	d /= 1000.;
	fp_mp << d << "u ";
}

void mp_graphics::output_xy_tikz(INT x, INT y)
{
	fp_tikz << "(";
	output_x_tikz(x);
	fp_tikz << ",";
	output_y_tikz(y);
	fp_tikz << ")";
}

void mp_graphics::output_x_tikz(INT x)
{
	double d;

	d = (double) x;
	d /= 30000.;
	fp_tikz << d;
}

void mp_graphics::output_y_tikz(INT y)
{
	double d;

	d = (double) y;
	d /= 30000.;
	fp_tikz << d;
}

void mp_graphics::comment(const BYTE *p)
{
	fp_mp << "% " << p << endl;
	fp_log << "Comment " << p << endl;
	fp_tikz << "% " << p << endl;
}


void mp_graphics::text(INT x, INT y, const char *p)
{
	INT x1, y1;
	BYTE align[64];
	INT lab;
	
	fp_log << "Text " << x << " " << y << " \"" << p << "\"" << endl;

	x1 = x;
	y1 = y;
	coords_min_max(x1, y1);
	user2dev(x1, y1);
	
	get_alignment(align);
	if (txt_boxed) {
		lab = get_label(x1, y1);
		fp_mp << "boxit.l" << lab << "(btex " << p << " etex);" << endl;
		fp_mp << lab << ".c=";
		output_xy_metapost(x1, y1);
		fp_mp << endl;
		if (txt_overwrite) {
			fp_mp << "unfill bpath " << lab << ";" << endl;
			}
		fp_mp << "drawboxed(" << lab << ");" << endl;
		
		}
	else {
		fp_mp << "label" << align << "(btex " << p << " etex";
		if (txt_rotate) {
			fp_mp << " rotated " << txt_rotate;
			}
		fp_mp << ", ";
		output_xy_metapost(x1, y1);
		fp_mp << ");" << endl;
		}


	if (txt_overwrite) {
		fp_tikz << "\\draw ";
		output_xy_tikz(x1, y1);
		fp_tikz << " node[fill=white] {";
		fp_tikz << p;
		fp_tikz << "};" << endl;
		}
	else {
		fp_tikz << "\\draw ";
		output_xy_tikz(x1, y1);
		fp_tikz << " node{";
		fp_tikz << p;
		fp_tikz << "};" << endl;
		}
	

}

void mp_graphics::circle(INT x, INT y, INT rad)
{
	INT X[10], Y[10], i;

	fp_log << "Circle " << x << " " << y << " " << rad << endl;

	coords_min_max(x, y);
	user2dev(x, y);
	user2dev_dist_x(rad);
	
	if (rad <= 0) rad = 1;
	
	X[0] = x +  rad;
	Y[0] = y;
	X[1] = x;
	Y[1] = y +  rad;
	X[2] = x -  rad;
	Y[2] = y;
	X[3] = x;
	Y[3] = y -  rad;
	X[4] = x +  rad;
	Y[4] = y;
	fp_mp << "path pp;" << endl;
	//fp_mp << "pickup pencircle scaled " << line_thickness << "pt;" << endl;
	fp_mp << "pp = ";
	for (i = 0; i < 5; i++) {
		if (i) {
			fp_mp << " .. ";
			}
		output_xy_metapost(X[i], Y[i]);
		}
	fp_mp << " .. cycle;" << endl;
	if (fill_interior > 0) {
		fp_mp << "fill pp withcolor ";
		if (fill_interior > 99) {
			fp_mp << "1 ";
			}
		else {
			fp_mp << "." << fill_interior << " ";
			// fprintf(mp_draw_fp, ".%02ld ", fill_interior);
			}
		if (fill_color == 1)
			fp_mp << "black";
		else
			fp_mp << "white";
		fp_mp << ";" << endl;
		}
	else {
		fp_mp << "draw pp";
		if (line_dashing) {
			fp_mp << " dashed evenly";
			if (line_dashing != 100) {
				fp_mp << " scaled " << (double) line_dashing / 100.;
				}
			}
		fp_mp << ";" << endl;
		}

	if (fill_interior > 0) {
		fp_tikz << "\\filldraw[fill=white] "; 
		output_xy_tikz(X[0] - rad, Y[0]);
		fp_tikz << " circle [radius = ";
		output_x_tikz(rad * 1);
		fp_tikz << "cm];" << endl;
		}
	else {
		fp_tikz << "\\draw "; 
		output_xy_tikz(X[0] - rad, Y[0]);
		fp_tikz << " circle [radius = ";
		output_x_tikz(rad * 1);
		fp_tikz << "cm];" << endl;
		}
}


void mp_graphics::circle_text(INT x, INT y, const char *text)
{
	INT idx;
	
	fp_log << "CircleText " << x << " " << y << " \"" << text << "\"" << endl;
	// cout << "mp_graphics::circle_text() text = " << text << " at " << x << ", " << y << endl;
	idx = get_label(x, y);

	
	coords_min_max(x, y);
	user2dev(x, y);

	output_circle_text_mp(x, y, idx, text);
	output_circle_text_tikz(x, y, idx, text);
}


void mp_graphics::output_circle_text_mp(INT x, INT y, INT idx, const char *text)
{
	fp_mp << "circleit.l" << idx << "(btex " << text << " etex);" << endl;
	fp_mp << "l" << idx << ".c = ";
	output_xy_metapost(x, y);
	fp_mp << endl;
	fp_mp << "unfill bpath l" << idx << ";" << endl;
	fp_mp << "drawboxed(l" << idx << ");" << endl;
}

void mp_graphics::output_circle_text_tikz(INT x, INT y, INT idx, const char *text)
{
	BYTE str[1000];

	sprintf(str, "%ld", idx);

	fp_tikz << "\\draw ";
	output_xy_tikz(x, y);
	fp_tikz << " node[fill=white]{" << text;
	fp_tikz << "};" << endl;
}

void mp_graphics::polygon_or_bezier_idx(INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle)
{
	INT x, y, i;

	//f << "pickup pencircle scaled " << line_thickness << "pt;" << endl;
	if (line_end_style == 1)
		fp_mp << "drawarrow ";
	else
		fp_mp << "draw ";
	for (i = 0; i < n; i++) {
		x = Px[Idx[i]];
		y = Py[Idx[i]];
		coords_min_max(x, y);
		user2dev(x, y);
		
		if (i) {
			fp_mp << symbol;
			}
		//fp_mp << "(" << x << "u," << y << "u)";
		output_xy_metapost(x, y);
		if (((i + 1) % 30) == 0)
			fp_mp << endl;
		}
	if (f_cycle) {
		fp_mp << " " << symbol << " cycle ";
		}

	if (line_dashing) {
		fp_mp << " dashed evenly";
		if (line_dashing != 100) {
			fp_mp << " scaled " << (double) line_dashing / 100.;
			}
		}
	fp_mp << ";" << endl;

	INT f_need_comma = FALSE;
	fp_tikz << "\\draw [";
	if (line_end_style == 1 && line_beg_style == 0) {
		fp_tikz << "->";
		f_need_comma = TRUE;
		}
	else if (line_end_style == 0 && line_beg_style == 1) {
		fp_tikz << "<-";
		f_need_comma = TRUE;
		}
	else if (line_end_style == 1 && line_beg_style == 0) {
		fp_tikz << "<->";
		f_need_comma = TRUE;
		}
	if (line_thickness != 100) {
		if (f_need_comma) {
			fp_tikz << ",";
			}
		fp_tikz << "line width=" << ((double)line_thickness * 0.01) << "mm";
		}
	fp_tikz << "] ";
	for (i = 0; i < n; i++) {
		x = Px[Idx[i]];
		y = Py[Idx[i]];
		coords_min_max(x, y);
		user2dev(x, y);

		if (i) {
			fp_tikz << " -- ";
			}
		output_xy_tikz(x, y);
		}
	fp_tikz << ";" << endl;


}

void mp_graphics::fill_idx(INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle)
{
	fill_idx_mp(Px, Py, Idx, n, symbol, f_cycle);
	fill_idx_tikz(fp_tikz, Px, Py, Idx, n, symbol, f_cycle);
}

void mp_graphics::fill_idx_mp(INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle)
{
	INT x, y, i;

	fp_mp << "fill buildcycle(";
	//fp_mp << "p[" << cur_path << "] = buildcycle(";
	for (i = 0; i < n; i++) {
		x = Px[Idx[i]];
		y = Py[Idx[i]];
		coords_min_max(x, y);
		user2dev(x, y);
		
		if (i) {
			fp_mp << symbol;
			}
		output_xy_metapost(x, y);
		if (((i + 1) % 30) == 0)
			fp_mp << endl;
		}
	if (f_cycle) {
		fp_mp << " " << symbol << " cycle ";
		}
	fp_mp << ")";

	//fp_mp << "fill p" << cur_path << " withcolor ";
	fp_mp << " withcolor ";
	fp_mp << ((double)fill_interior / (double) 100);
	if (fill_color == 0) {
		fp_mp << " white;" << endl;
		}
	else if (fill_color == 1) {
		fp_mp << " black;" << endl;
		}
	cur_path++;
	
}

void mp_graphics::fill_idx_tikz(ofstream &fp, INT *Px, INT *Py, INT *Idx, INT n, const char *symbol, INT f_cycle)
{
	INT f_need_comma;
	INT i, x, y;
	
	fp << "\\fill [fill=";

	if (fill_color == 0)
		fp << "white";
	else if (fill_color == 1)
		fp << "black";
	else if (fill_color == 2)
		fp << "red";
	else if (fill_color == 3)
		fp << "green";

	if (fill_interior < 100) {
		fp << "!" << fill_interior;
		}
	f_need_comma = TRUE;
	if (line_end_style == 1 && line_beg_style == 0) {
		if (f_need_comma) {
			fp << ", ";
			}
		fp << "->";
		}
	else if (line_end_style == 0 && line_beg_style == 1) {
		if (f_need_comma) {
			fp << ", ";
			}
		fp << "<-";
		}
	else if (line_end_style == 1 && line_beg_style == 0) {
		if (f_need_comma) {
			fp << ", ";
			}
		fp << "<->";
		}
	if (line_thickness != 100) {
		if (f_need_comma) {
			fp << ", ";
			}
		fp << "line width=" << ((double)line_thickness * 0.01) << "mm";
		}
	fp << "] ";
	for (i = 0; i < n; i++) {
		x = Px[Idx[i]];
		y = Py[Idx[i]];
		coords_min_max(x, y);
		user2dev(x, y);

		if (i) {
			fp << " -- ";
			}
		output_xy_tikz(x, y);
		}
	fp << ";" << endl;



}


