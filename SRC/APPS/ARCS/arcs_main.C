// arcs_main.C
// 
// Anton Betten
//
// previous version Dec 6, 2004
// revised June 19, 2006
// revised Aug 17, 2008
//
// Searches for arcs in desarguesian projective planes
//
//

#include "orbiter.h"
#include "arcs.h"

// global data:

INT t0; // the system time when the program started


int main(int argc, const char **argv)
{
	t0 = os_ticks();
	
	
	{
	arc_generator *Gen;
	finite_field *F;

	
	Gen = new arc_generator;

	Gen->read_arguments(argc, argv);
	

	F = new finite_field;
	F->init(Gen->q, 0);


	Gen->init(F, 
		Gen->ECA->input_prefix, 
		Gen->ECA->base_fname,
		Gen->ECA->starter_size, 
		argc, argv, 
		Gen->verbose_level);
	


	Gen->main(Gen->verbose_level);
	
	delete Gen;
	delete F;
	
	}
	//the_end(t0);
	the_end_quietly(t0);
}



