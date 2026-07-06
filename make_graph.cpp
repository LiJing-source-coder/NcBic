/* Author:  Zhenjia Wang<zhenjia.sdu@gmail.com> Dec. 22, 2014
 * Usage: This is part of bicluster package. Use, redistribution, modify without limitations
 * show how does the whole program work
 *
 * Produces two graphs sequentially, derived from microarray data.
 *
 * The first graph is generated from the raw data where an edge is defined
 * as two genes having common components from same condition, the score on
 * the edge being the number of the same components. The edges in the first/////
 * graph, are used as vertices in the second graph, where edges are defined
 * as the common columns between two co-vertex edges in the first graph,
 * with scores defined as the number of same columns for all three genes.
 *
 */

#include "make_graph.h"
#include "checkpoint.h"
/*we can reduce the HEAP_SIZE when the data contain so many genes so that memory is not enough*/
static const int HEAP_SIZE = 20000000;

/**************************************************************************/

/* String intersection function without string copying, only numbers */
/*caculate the weight of the edge in the first graph*/
int common_num_1 (const bool *s1, const bool *s2)
{
	int common_cnt = 0;
	/* s1 and s2 of equal length, so we check s1 only */
	int i;
    for (i=0;i<Mcols;i++)
    {
        if ((*s1) * (*s2) > 0)
        {
            common_cnt++;
        }

	s1++;
	s2++;
    }
	return common_cnt;
}

/* Fibonacci heap related subroutines */
static int edge_cmpr(void *a, void *b)
{
	float score_a, score_b;
	score_a = ((Edge *)a)->score;
	score_b = ((Edge *)b)->score;
	if (score_a < score_b) return -1;
	if (score_a == score_b) return 0;
	return 1;
}

/* Maintain a fixed size heap */
static void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min)
{
	if (a->fh_n < HEAP_SIZE)
	{
		fh_insert(a, (void *)i);
	}
	else
	{
		if (edge_cmpr(*cur_min, i) < 0)
		{
			/* Remove least value and renew */
			fh_extractmin(a);
			fh_insert(a, (void *)i);
			/* Keep a memory of the current min */
			*cur_min = (Edge *)fh_min(a);

		}
	}
}

/*sort the edges in decrease order so that e1 is the largest edge*/
static void fh_dump(struct fibheap *a, Edge **res)
{
	int i;
	int n = a->fh_n;
	for (i=n-1; i>=0; i--)
		res[i] = (Edge *) fh_extractmin(a);
}

/*********************************************************************/
void make_graph (const char* fn)
{
    FILE *fw = mustOpen(fn, "w");
	int i, j,step;
	int rec_num = 0;
	int PART= 4;
        int conditionL[cols]={0};
	int cons;

	int ckpt_phase = checkpoint_enabled() ? checkpoint_probe_phase() : 0;
	if (ckpt_phase >= 1) {
		checkpoint_load_seeds(rec_num, edge_list);
		progress("Resuming from checkpoint: %d seeds loaded", rec_num);
	} else {
		/* Set a threshold about the weight of a seed (i.e.,the minimum number of common 1's) */
		for(i=0;i<cols;i++)
    		{
        		conditionL[i]=1;
    		}
    		i=compute_confidence(conditionL,cols,po->SEED_RePValue,cons);
    		po->SEED_MAXWEIGHT = i;
    		po->SEED_MINWEIGHT = cons-i;

    		/* edge_ptr describe edges */
    		edge_list = (Edge **)malloc(sizeof(Edge *) * HEAP_SIZE);

    		/* Allocating heap structure */
    		struct fibheap *heap;
    		heap = fh_makeheap();
    		fh_setcmp(heap, edge_cmpr);

		/* Generating seed list and push into heap */
		progress("Generating seed list (higher than %d 1's)", po->SEED_MINWEIGHT);

		Edge __cur_min = {0, 0, po->SEED_MINWEIGHT};
		Edge *_cur_min = &__cur_min;
		Edge **cur_min = & _cur_min;
		/* iterate over all genes to retrieve all edges */
		step = rows/PART;
		#pragma omp parallel num_threads(4) private(i,j)
		{
        		int k,endi,cnt;
			#pragma omp for nowait
			for(k=0;k<PART;k++)
			{
				endi=(k+1)*step;
				if(k==PART-1)
				endi=rows;
				for (i=k*step; i < endi; i++)
				{
					for (j = i+1; j < endi; j++)
					{
						cnt = common_num_1(Marr[i], Marr[j]);
						if(cnt < po->SEED_MAXWEIGHT)
               						{
                   						 cnt = common_num_1(Marr[i], Marr[j+rows]);
                   					 		j=j+rows;
               					 	}
						#pragma omp critical
                                               {
               				       /*select seed with its weight lower than po->SEED_MINWEIGHT and push it into the seed list "heap"*/
						if(cnt > po->SEED_MINWEIGHT) 
               						{
						edge_ptr = (Edge *)malloc(sizeof(*edge_ptr));
                   						edge_ptr -> gene_one = i;
                   					 	edge_ptr -> gene_two = j;
                   					 	edge_ptr -> score = cnt;
						fh_insert_fixed(heap, edge_ptr, cur_min);
				 				}
						}
						j=j%rows;
					}
				}
			}
		}
		rec_num = heap->fh_n;
		if (rec_num == 0)
			errAbort("Not enough overlap between genes");

		/* sort the seeds */
		uglyTime("%d seeds generated", rec_num);
		edge_list = (Edge **)realloc(edge_list, sizeof(*edge_list) * rec_num);
		fh_dump(heap, edge_list);

		if (checkpoint_enabled()) {
			checkpoint_save(1, edge_list, rec_num, 0, {}, {});
		}
	}

	/* bi-clustering */
	int n_blocks = 0;
	progress("Clustering started");
	n_blocks = cluster(fw, edge_list, rec_num, edge_list);
	uglyTime("%d clusters are written to %s", n_blocks, fn);

	if (checkpoint_enabled()) {
		checkpoint_remove();
	}

	/* clean up */
	for (i=0; i<rec_num; i++)
		free(edge_list[i]);
	free(edge_list);
}

/***************************************************************************/
