#ifndef READ_ARRAY_H
#define READ_ARRAY_H

/************************************************************************/

#include "struct.h"

/***********************************************************************/

/* prototypes  */
float **alloc2d (int rr, int cc);
short **alloc2c (int rr, int cc);
ColPair** alloc2m(int rr, int cc);
void get_matrix_size (FILE *fp);
void read_labels (FILE *fp);
void read_expression (FILE *fp);
#endif
