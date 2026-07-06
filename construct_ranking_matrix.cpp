#include "construct_ranking_matrix.h"

/***********************************************************************/
/* Matrix allocations (bool 2d array) */
bool** alloc2b(int rr, int cc)
{
	bool** result = (bool **)malloc(sizeof(bool *) * rr);
	int i;
	for (i = 0; i < rr; i++)
          result[i] = (bool *)malloc(sizeof(bool) * cc);
        return result;
}
/* Construct the column labels of ranking matrix M */
void construct_col_labels ()
{
    int i,j,snum=0;
    col_pair_to_idx = (int **)malloc(sizeof(int *) * cols);
    for(i=0;i<cols;i++) {
        col_pair_to_idx[i] = (int *)malloc(sizeof(int) * cols);
        for(j=0;j<cols;j++) {
            col_pair_to_idx[i][j] = -1;
        }
    }
	for(i=0;i<cols-1;i++)
        	for(j=i+1;j<cols;j++)
       	 	{
       	 	    Mcollabels[snum]->condA=i;
                Mcollabels[snum]->condB=j;
                Mcollabels[snum+1]->condA=j;
                Mcollabels[snum+1]->condB=i;
                col_pair_to_idx[i][j] = snum;
                col_pair_to_idx[j][i] = snum + 1;
            	snum=snum+2;
        	}

}
/* Construct the ranking matrix M */
void construct_ranking_matrix ()
{
    int i,j,k,snum;
    Marr = alloc2b(Mrows,Mcols);
    for(i=0;i<rows;i++)
    {
        snum=0;
        for(j=0;j<cols-1;j++)
            for(k=j+1;k<cols;k++)
            {
                if(arr[i][j]==arr[i][k])
                {
                    Marr[i][snum]=0;
                    Marr[i][snum+1]=0;
		    Marr[rows+i][snum]=0;
		    Marr[rows+i][snum+1]=0;
                }
                else if(arr[i][j]<arr[i][k])
                {
                    Marr[i][snum]=1;
                    Marr[i][snum+1]=0;
		    Marr[rows+i][snum]=0;
		    Marr[rows+i][snum+1]=1;

                }
                else
                {
                    Marr[i][snum]=0;
                    Marr[i][snum+1]=1;
		    Marr[rows+i][snum]=1;
		    Marr[rows+i][snum+1]=0;
                }
                snum=snum+2;
            }
    }
}
/*Free the input matrix A*/
void free_input_matrix()
{
	int i;
	for(i=0;i<rows;i++)
        {
            free(arr[i]);
            arr[i]=NULL;
        }

        free(arr);
        arr=NULL;
}

