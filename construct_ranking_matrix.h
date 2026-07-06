#ifndef CONSTRUCT_RANKING_MATRIX
#define CONSTRUCT_RANKING_MATRIX

/************************************************************************/

#include "struct.h"

/***********************************************************************/

/* prototypes  */
bool** alloc2b(int rr, int cc);
void construct_col_labels ();
void construct_ranking_matrix();
void free_input_matrix();
#endif
