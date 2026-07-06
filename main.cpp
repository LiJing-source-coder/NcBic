/* Author:  Jing Li May 22, 2025
 * Usage: This is part of bicluster package. Use, redistribution, modify without limitations
 * show how does the whole program work
 */

/***********************************************************************/

#include "main.h"
#include "cluster.h"
#include "checkpoint.h"
#include "write_block.h"
/***********************************************************************/
int main(int argc, char* argv[])
{
    /* Start the timer */
        uglyTime(NULL);
        printf("\nNcBic %.1f: greedy biclustering (compiled %s %s)\n\n", (double)VER,__DATE__, __TIME__);
        rows = cols = 0;

	/* get the program options defined in get_options.c */
        get_options(argc, argv);

    /*get the size of input expression matrix*/
        get_matrix_size(po->FP);
        progress("File %s contains %d genes by %d conditions", po -> FN, rows, cols);
        if (rows < 4 || cols < 4)
        {
		  /*neither rows number nor cols number can be too small*/
            errAbort("Not enough genes or conditions to make inference");
        }
        if(po->IS_RMAXDefault)
        {
	//	po->ROW_MAXWIDTH = (int)(0.5*rows);
		po->ROW_MAXWIDTH = rows;
        }
        if (po->ROW_MINWIDTH > po->ROW_MAXWIDTH)
        {
            errAbort("-z or -x minimum row width should be no more than the maximum row width");
        }
        if(po->ROW_MINWIDTH>=rows)
        {
            errAbort("-z minimum row width should be no more than all rows");
        }
        if(po->COL_WIDTH==-1)
        {
            po->COL_WIDTH = (int)(0.05*cols);
        }
        genes = alloc2c(rows, LABEL_LEN);
        conds = alloc2c(cols, LABEL_LEN);

	/* Read in the gene names and condition names */
        read_labels(po -> FP);
    /* Read in the expression data */
     	read_expression(po -> FP);
        fclose(po->FP);
	/* Construct the ranking matrix M */

    /*to find the corresponding columns between matrix M and the input matrix A easily,
      we construct column labels by use of the order of columns of A.*/
        Mrows=2*rows;
        Mcols=cols*(cols-1);
        Mcollabels = alloc2m(Mcols,sizeof(ColPair));
        construct_col_labels();
        construct_ranking_matrix();
        free_input_matrix();
    /* Dump checkpoint blocks if requested */
        if (po->IS_DUMP_CKPT) {
            char* out_fn = addSuffix(po->FN, ".txt");
            FILE* fw = fopen(out_fn, "w");
            if (!fw) errAbort("Cannot open %s for writing", out_fn);
            int n = checkpoint_dump_blocks(fw);
            fclose(fw);
            if (n == 0) {
                remove(out_fn);
            }
            free(out_fn);
            free(po);
            return 0;
        }

    /* Generate the seed list and operate biclustering*/
    /* po->FN: the file that stores all blocks */
        make_graph(addSuffix(po->FN, ".blocks"));
        
    	free(po);
    	return 0;
}
