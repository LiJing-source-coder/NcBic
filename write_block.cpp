/******************************************************************/
/* Author: Jing Li <jlimath@mail.sdu.edu.cn>, Jul. 24, 2025
 * Output the identified bicluster block.
 */

#include "write_block.h"

/*************************************************************************/
static bool condcount(int a)
{
    return a>0;
}
/*************************************************************************/
void print_bc(FILE *fw, const std::unique_ptr<Block> &b,int num){
    int geneindex,layerindex,negeneloc;
  /* block height (genes) */
  const int block_rows = std::count(b->genes.begin(),b->genes.end(),true);
  /* block_width (conditions) */
  const int block_cols = std::count_if(b->BTP.begin(),b->BTP.end(),condcount);
  /* BTP(condition layer) */
  const int block_BTPlayer = *std::max_element(b->BTP.begin(),b->BTP.end());
  
  fprintf(fw, "\n");
  fprintf(fw, " Block: %d", num+1);
  fprintf(fw, "\n");
  fprintf(fw, " Genes [%d]: ", block_rows);
  fprintf(fw, "\n");
  fprintf(fw, " Positive Regulatory Genes: ");
  fprintf(fw, "\n");
  for (geneindex=0;geneindex<rows;geneindex++)
  {
            if(b->genes[geneindex])
	    {
            	fprintf(fw, "%s ", genes[geneindex]);
	    }
  }
  fprintf(fw, "\n");
  fprintf(fw, " Negative Regulatory Genes: ");
  fprintf(fw, "\n");
  for (geneindex=rows;geneindex<Mrows;geneindex++)
  {
            if(b->genes[geneindex])
            {
                negeneloc=geneindex-rows;
                fprintf(fw, "%s ", genes[negeneloc]);
            }

  }
  fprintf(fw, "\n");

  fprintf(fw, " Conds [%d]: ", block_cols);
  fprintf(fw, "\n");
  for(layerindex=1;layerindex<=block_BTPlayer;layerindex++)
  {
        fprintf(fw, " Bucket %d: ", layerindex);
        for (std::size_t BTPindex=0;BTPindex<b->BTP.size();BTPindex++)
        {
            if(b->BTP[BTPindex]==layerindex)
            {
                fprintf(fw, "%s ", conds[BTPindex]);
            }
        }
        fprintf(fw, "\n");
  }
  fputc('\n', fw);
}
