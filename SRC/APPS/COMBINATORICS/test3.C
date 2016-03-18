// test3.C
// 
// Anton Betten
// May 12, 2010
//
//
// 
// Creates all NSEW paths from (0,0) to (I,J) using H hops.
//

#include <iostream>
#include <fstream>
//#include <sstream>
#include <iomanip>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

using namespace std;


// global data:

int nb_sol = 0;

int generate(int H, int I, int J);
void print_path(int *digits, int length);
void generate_recursion(int H, int I, int J, int *digits, int depth, int max_depth);

int main(int argc, char **argv)
{
	//int verbose_level = 2;
#if 0
	if (argc <= 1) {
		usage(argc, argv);
		exit(1);
		}
#endif
	int H, I, J;

	H = atoi(argv[argc - 3]);
	I = atoi(argv[argc - 2]);
	J = atoi(argv[argc - 1]);
	cout << "H=" << H << " I=" << I << " J=" << J << endl;
	generate(H, I, J);
	cout << "we found " << nb_sol << " solutions for H=" << H << " I=" << I << " J=" << J << endl;
	
}

int generate(int H, int I, int J)
{
	int *digits;
	
	digits = new int[H];
	nb_sol = 0;
	generate_recursion(H, I, J, digits, 0, H);
	return nb_sol;
}

void print_path(int *digits, int length)
{
	int i;
	
	for (i = 0; i < length; i++) {
		if (digits[i] == 0) {
			cout << "N";
			}
		if (digits[i] == 1) {
			cout << "S";
			}
		if (digits[i] == 2) {
			cout << "E";
			}
		if (digits[i] == 3) {
			cout << "W";
			}
		}
	cout << "\\\\" << endl;
}

void generate_recursion(int H, int I, int J, int *digits, int depth, int max_depth)
{
	int i, j;
	
	//cout << "H=" << setw(5) << H << " I=" << setw(5) << I << " J=" << setw(5) << J << " : ";
	//print_path(digits, depth);

	if (H == 0) {
		if (I == 0 && J == 0) {
			nb_sol++;
#if 1
			cout << setw(10) << nb_sol << " & ";
			print_path(digits, max_depth);
#endif
			}
		return;
			
		}
	for (i = 0; i < 4; i++) {



		if (depth) {
			j = digits[depth - 1];

			if (i == 0 && j == 1) {
				continue;
				}
			if (i == 1 && j == 0) {
				continue;
				}
			if (i == 2 && j == 3) {
				continue;
				}
			if (i == 3 && j == 2) {
				continue;
				}
			}
		digits[depth] = i;
		//print_path(digits, depth + 1);
		if (i == 0) {
			// N
			generate_recursion(H - 1, I, J - 1, digits, depth + 1, max_depth);
			}
		if (i == 1) {
			// S
			generate_recursion(H - 1, I, J + 1, digits, depth + 1, max_depth);
			}
		if (i == 2) {
			// E
			generate_recursion(H - 1, I - 1, J, digits, depth + 1, max_depth);
			}
		if (i == 3) {
			// W
			generate_recursion(H - 1, I + 1, J, digits, depth + 1, max_depth);
			}
		}

}


