/*
 * Author:Qin Ma <maqin@csbl.bmb.uga.edu>, Jan. 24, 2010
 * Usage: This is part of the bicluster package. Use, redistribute, modify without limitations.
 *
 * The procedure for file input:
 * ----------------------------------------
 * 	o 	cond1	 cond2	 cond3
 * 	gene1    3.14	  -1.2     0.0
 * 	gene2    0.0       2.8     4.5
 * ----------------------------------------
 */

#include "read_array.h"

/************************************************************************/
/* Helper variables for tokenizer function */

static char *atom = NULL;
static char delims[] = "\t\r\n";

/***********************************************************************/

/* Matrix allocations (2d array) */

float** alloc2d(int rr, int cc)
{
	float** result  = (float **) malloc(sizeof(float *) * rr);
	int i;
	for (i = 0; i < rr; i++)
		result[i] = (float *)malloc(sizeof(float) * cc);
	return result;
}

short** alloc2c(int rr, int cc)
{
	short** result = (short **)malloc(sizeof(short *) * rr);
	int i;
	for (i = 0; i < rr; i++)
          result[i] = (short *)malloc(sizeof(short) * cc);
        return result;
}

ColPair** alloc2m(int rr, int cc)
{
	ColPair** result = (ColPair **)malloc(sizeof(ColPair *) * rr);
	int i;
	for (i = 0; i < rr; i++)
          result[i] = (ColPair *)malloc(sizeof(ColPair) * cc);
        return result;
}
/***********************************************************************/

/* Pre-read the datafile, retrieve gene labels and condition labels
 * as well as determine the matrix size
 */
void get_matrix_size (FILE* fp)
{
	/*size_t is the best type to use if you want to represent sizes of objects.
	 * Using int to represent object sizes is likely to work on most modern systems, but it isn't guaranteed.
	 */
	size_t n = 0;
	char *line = NULL;
	/*getline() reads an entire line, storing the address of the buffer containing the text into *line.
	 *the buffer is null-terminated and includes the newline character, if a newline delimiter was found.
	 */
	if (getline(&line, &n, fp)>=0)
	{
		/*strtok function returns a pointer to the next token in str1, where str2 contains the delimiters that determine the token*/
		atom = strtok(line, delims);
		/*delete the first element in atom because the first element corresponding to description column*/
		atom = strtok(NULL, delims);
		while (atom != NULL)
		{
			/*if we do not set atom = strtok(NULL, delims), here is a infinite loop*/
			atom = strtok(NULL, delims);
			cols++;
		}
	}
	while (getline(&line, &n, fp)>=0)
	{
		atom = strtok(line, delims);
		rows++;
	}
	/*fseed sets the position indicator associated with the stream to a new position defined by adding offset to a reference position specified by origin*/
	fseek(fp, 0, 0);
}
/* Read in the labels on x and y, in microarray/RNA-seq terms, genes(rows) and conditions(cols)*/
void read_labels (FILE* fp)
{
	int row = 0;
	int col;
	size_t n = 0;
	char *line = NULL;
	while (getline(&line, &n, fp)>=0)
	{
		atom = strtok(line, delims);
		/*currently the first element in atom is the gene name of each row when row>=1, the 0 row corresponding to the line of condition names*/
		if (row >= 1)
		{
			strcpy(genes[row-1], atom);
		}
		/*delete the first element in atom because the first element corresponding to description column*/
		atom = strtok(NULL, delims);
		col = 0;
		while (atom != NULL)
		{
			if (row == 0)
				strcpy(conds[col], atom);
			atom = strtok(NULL, delims);
			if (++col == cols) break;
		}
		if (++row == rows+1) break;
	}
	fseek(fp, 0, 0);
}
/* Read in the expression data, and save each expression value as a float.*/
void read_expression (FILE* fp)
{
	int row, col;
	arr = alloc2d(rows,cols);
	/* import data */
	size_t n = 0;
	char *line = NULL;
	row = 1;
	/* ignore header line */
	getline(&line, &n, fp);
	while (getline(&line, &n, fp)>=0)
	{
		atom = strtok(line, delims);
		/*skip the first column*/
		atom = strtok(NULL, delims);
		col = 0;
		while (atom != NULL)
		{
			/*Checks if parameter atom is either an uppercase or a lowercase alphabetic letter*/
			if (isalpha(*atom))
				arr[row-1][col] = 0.0;
			else
				arr[row-1][col] = atof(atom);
			atom = strtok(NULL, delims);
			if (++col == cols) break;
		}
		if (++row == rows+1) break;
	}
	fseek(fp, 0, 0);
}
/***********************************************************************/
