#ifndef _STRUCT_H
#define _STRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <omp.h>
#include <stdarg.h>
#include <limits.h>
#include <sys/time.h>
#include <time.h>
#include "fib.h"
#include <vector>
#include <algorithm>
#include<iostream>
#include<numeric>
#include<cmath>
#include<boost/multiprecision/cpp_dec_float.hpp>

using namespace boost::multiprecision;
using namespace std;
/* Constants */
#define LABEL_LEN 64

/* Pretend that C has boolean type */
#define TRUE 1
#define FALSE 0

/* Strings */
/* strcmp: a zero value indicates that both strings are equal.
 * a value greater than zero indicates that the first character that does not match has a greater value in str1 than in str2;
 * And a value less than zero indicates the opposite.
 */
#define sameString(a, b) (strcmp((a), (b))==0)
/* Returns TRUE if two strings are same */

/* global data */
extern float **arr;
extern bool **Marr;
extern char **genes;
extern char **conds;
extern int rows, cols;
extern int Mrows, Mcols;
extern short *lcs_length;
extern char **lcs_tags;

/***** Structures *****/
struct ColPair
  {
   int condA;            /* condition A*/
   int condB;		    /* condition B */
  };

extern ColPair **Mcollabels;     /*column labels of matrix M */
extern int **col_pair_to_idx;    /* fast lookup: col_pair_to_idx[a][b] = column index where condA=a, condB=b */

/* edge between two genes */
typedef struct Edge{
	int gene_one;
	int gene_two;
	int score;
} Edge;
/* one block */
typedef struct Block{
  std::vector<int> BTP = std::vector<int> (cols,0);
  std::vector<bool> genes = std::vector<bool> (Mrows,false);
  cpp_dec_float_100 p_value=1.0;
  int size;
}Block;

/* holds running options */
typedef struct Prog_options{
	char FN[LABEL_LEN];
	bool IS_TFname;
	bool IS_RMAXDefault;
	int COL_WIDTH;
	int ROW_MAXWIDTH;
   	int ROW_MINWIDTH;
	int SEED_MINWEIGHT;
	int SEED_MAXWEIGHT;
	int ROW_MAXINIT;
	double R_THRESHOLD;
	double SEED_RePValue;
	double TOLERANCE;
	double pvalue;
	int SCH_BLOCK;
	int RPT_BLOCK;
	double FILTER;
	char TFname[LABEL_LEN];
	FILE* FP;
	int CHECKPOINT_INTERVAL; /* minutes, 0 = disabled */
	bool IS_DUMP_CKPT;       /* dump checkpoint blocks and exit */
} Prog_options;

extern Prog_options* po;
/***** Helper functions *****/

void progress(const char *format, ...)
/* Print progress message */
     __attribute__((format(printf, 1, 2)));

void uglyTime(const char *label, ...);
/* Print label and how long it's been since last call.  Call with
 * a NULL label to initialize. */

void err(const char *format, ...)
/* Print error message but do not exit */
     __attribute__((format(printf, 1, 2)));

void errAbort(const char *format, ...)
/* Print error message to stderr and exit */
     __attribute__((noreturn, format(printf, 1, 2)));

/* File-related operations */
char *addSuffix(const void *head, const void *suffix);

FILE *mustOpen(const char *fileName, const char *mode);
/* Open a file or die */
#endif
