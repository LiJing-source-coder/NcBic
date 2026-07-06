#ifndef MAIN_H
#define MAIN_H

#include "struct.h"
/*get_option subroutine prototypes */
void get_options (int argc, char* argv[]);

/* read_array subroutine prototypes  */
extern char** alloc2c (int rows, int cols);
extern ColPair** alloc2m(int rr, int cc);
extern void get_matrix_size (FILE* fp);
extern void read_labels (FILE* fp);
extern void read_expression (FILE* fp);

/* construct ranking matrix*/
extern void construct_col_labels ();
extern void construct_ranking_matrix();
extern void free_input_matrix();

/*make graph*/
extern void make_graph (const char* fn);


#endif
