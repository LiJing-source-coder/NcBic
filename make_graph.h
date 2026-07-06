#ifndef MAKE_GRAPH_H
#define MAKE_GRAPH_H

#include "struct.h"

#include<fstream>

//using namespace std;
/* global data */
Edge **edge_list;
Edge *edge_ptr;

/* prototypes */
int common_num_1 (const bool *s1, const bool *s2);
static int edge_cmpr(void *a, void *b);
static void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min);
static void fh_dump(struct fibheap *a, Edge **res);
void make_graph (const char* fn);

/*from enum_gene_pvalue: enumerate all the p-values for genes*/
extern int compute_confidence(const int *conditionL,int m,double repvalue,int &cons);

/* from cluster */
extern int cluster (FILE *fw, Edge **el, int n, Edge **edge_list = nullptr);
#endif
