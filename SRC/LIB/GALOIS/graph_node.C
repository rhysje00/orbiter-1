// graph_node.C
// 
// Anton Betten
// December 30, 2013
//
//
// 
//
//

#include "galois.h"

graph_node::graph_node()
{
	null();
}

graph_node::~graph_node()
{
	freeself();
}

void graph_node::null()
{
	label = NULL;
	id = -1;
	layer = -1;
	f_has_data1 = FALSE;
	data1 = -1;
	f_has_data2 = FALSE;
	data2 = -1;
	f_has_data3 = FALSE;
	data3 = -1;
	f_has_vec_data = FALSE;
	vec_data = NULL;
	vec_data_len = 0;
	f_has_distinguished_element = FALSE;
	distinguished_element_index = -1;
	nb_neighbors = 0;
	neighbor_list = NULL;
	neighbor_list_allocated = 0;
}

void graph_node::freeself()
{
	if (label) {
		FREE_BYTE(label);
		}
	if (neighbor_list) {
		FREE_INT(neighbor_list);
		}
	if (f_has_vec_data) {
		FREE_INT(vec_data);
		}
	null();
}

void graph_node::add_neighbor(INT l, INT n, INT id)
{
	INT i;
	
	if (nb_neighbors >= neighbor_list_allocated) {
		INT new_neighbor_list_allocated;
		INT *new_neighbor_list;
		
		if (neighbor_list_allocated) {
			new_neighbor_list_allocated = 2 * neighbor_list_allocated;
			}
		else {
			new_neighbor_list_allocated = 16;
			}
		new_neighbor_list = NEW_INT(new_neighbor_list_allocated);
		for (i = 0; i < nb_neighbors; i++) {
			new_neighbor_list[i] = neighbor_list[i];
			}
		if (neighbor_list) {
			FREE_INT(neighbor_list);
			}
		neighbor_list = new_neighbor_list;
		neighbor_list_allocated = new_neighbor_list_allocated;
		}
	neighbor_list[nb_neighbors] = id;
	nb_neighbors++;
}

void graph_node::add_text(const BYTE *text)
{
	INT l;
	BYTE *p;

	l = strlen(text);
	p = NEW_BYTE(l + 1);
	strcpy(p, text);
	if (label) {
		FREE_BYTE(label);
		}
	label = p;
}

void graph_node::add_vec_data(INT *v, INT len)
{
	vec_data = NEW_INT(len);
	vec_data_len = len;
	INT_vec_copy(v, vec_data, len);
	f_has_vec_data = TRUE;
}

void graph_node::set_distinguished_element(INT idx)
{
	f_has_distinguished_element = TRUE;
	distinguished_element_index = idx;
}


void graph_node::add_data1(INT data)
{
	f_has_data1 = TRUE;
	data1 = data;
}

void graph_node::add_data2(INT data)
{
	f_has_data2 = TRUE;
	data2 = data;
}

void graph_node::add_data3(INT data)
{
	f_has_data3 = TRUE;
	data3 = data;
}

void graph_node::write_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "graph_node::write_memory_object" << endl;
		}
	if (label == NULL) {
		m->write_string("");
		}
	else {
		m->write_string(label);
		}
	m->write_int(id);
	m->write_int(f_has_data1);
	m->write_int(data1);
	m->write_int(f_has_data2);
	m->write_int(data2);
	m->write_int(f_has_data3);
	m->write_int(data3);
	m->write_int(f_has_vec_data);
	if (f_has_vec_data) {
		m->write_int(vec_data_len);
		for (i = 0; i < vec_data_len; i++) {
			m->write_int(vec_data[i]);
			}
		}
	m->write_int(f_has_distinguished_element);
	m->write_int(distinguished_element_index);
	m->write_int(layer);
	m->write_int(nb_neighbors);
	for (i = 0; i < nb_neighbors; i++) {
		m->write_int(neighbor_list[i]);
		}
	m->write_double(x_coordinate);
	if (f_v) {
		cout << "graph_node::write_memory_object finished, data size (in bytes) = " << m->used_length << endl;
		}
}

void graph_node::read_memory_object(memory_object *m, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;
	
	if (f_v) {
		cout << "graph_node::read_memory_object" << endl;
		}
	m->read_string(label);
	m->read_int(&id);
	m->read_int(&f_has_data1);
	m->read_int(&data1);
	m->read_int(&f_has_data2);
	m->read_int(&data2);
	m->read_int(&f_has_data3);
	m->read_int(&data3);
	m->read_int(&f_has_vec_data);
	if (f_has_vec_data) {
		m->read_int(&vec_data_len);
		vec_data = NEW_INT(vec_data_len);
		for (i = 0; i < vec_data_len; i++) {
			m->read_int(&vec_data[i]);
			}
		}
	else {
		vec_data_len = 0;
		vec_data = NULL;
		}
	
	m->read_int(&f_has_distinguished_element);
	m->read_int(&distinguished_element_index);


	m->read_int(&layer);
	m->read_int(&nb_neighbors);
	neighbor_list = NEW_INT(nb_neighbors);
	for (i = 0; i < nb_neighbors; i++) {
		m->read_int(&neighbor_list[i]);
		}
	
	m->read_double(&x_coordinate);
	if (f_v) {
		cout << "graph_node::read_memory_object finished" << endl;
		}
}


