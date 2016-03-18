// test4.C
// 
// Anton Betten
// October 6, 2015


#include "orbiter.h"


int main(int argc, char **argv)
{
	INT i, j, a, b, c;

	a = 10000;
	for (i = 1; i <= a; i++) {
		if ((a % i) != 0) {
			continue;
			}
		b = a / i;
		for (j = i; j <= b; j++) {
			if ((b % j) != 0) {
				continue;
				}
			c = b / j;
			if (c >= j) {
				cout << i << " . " << j << " . " << c << endl;
				}
			}
		}
}


