// codes.C
//
// Anton Betten
// December 30, 2003

#include "orbiter.h"
#include "discreta.h"
#include "codes.h"

// global data:

INT t0; // the system time when the program started
const BYTE *version = "codes.C version 5/18/2009";

void print_usage()
{
	cout << "usage: codes.out [options] -n <n> -k <k> -q <q> -d <d>" << endl;
	cout << "where options can be:" << endl;

	generator gen;
	
	gen.usage();

}


int main(int argc, const char **argv)
{
	cout << version << endl;
	t0 = os_ticks();
	
	
	{
	code_generator cg;
	
	cg.init(argc, argv);

	cg.main();
	
	//cg.gen->print_tree();
	
	}
	the_end(t0);
	//the_end_quietly(t0);
}

