//
// Created by revans on 16/11/16.
//

#include "orbiter.h"

int main(int argc, char **argv)
{
    INT verbose_level = 0;
    INT i, j;
    INT *Adj;
    INT f_n = FALSE;
    INT n;

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
    }

    if (!f_n) {
        cout << "Please use option -n <n>" << endl;
        exit(1);
    }

    Adj = NEW_INT(n * n);
    INT_vec_zero(Adj, n * n);

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            Adj[i * n + j] = 1;
            Adj[j * n + 1] = 1;
        }
    }

    colored_graph *CG;
    BYTE fname[1000];

    CG = new colored_graph;
    CG->init_adjacency_no_colors(n, Adj, verbose_level);

    sprintf(fname, "Complete_%ld.colored_graph", n);

    CG->save(fname, verbose_level);

    delete CG;
    FREE_INT(Adj);
}
