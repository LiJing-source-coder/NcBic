#include "cluster.h"
#include <chrono>
#include <iostream>

#include "checkpoint.h"
#include <ctime>

/*******************************************************/
int** alloc2e(int rr, int cc)
{
	int** result = (int **)malloc(sizeof(int *) * rr);
	int i;
	for (i = 0; i < rr; i++)
          result[i] = (int *)malloc(sizeof(int) * cc);
        return result;
}	
/********************************************************/
static int max_value(std::vector<int> vec)
{
	int maxi=-1;
        for(std::size_t k=0;k<vec.size();k++)
    if(maxi<vec[k])
    {maxi=vec[k];}
            return maxi;
}
/*********************************************************/
/* compute the p_value of a BTP-pattern */

static void get_BTP_genes_and_pvalue(std::vector<int> BTP, std::vector<bool> &genes, cpp_dec_float_100 &p_value)
{

	/* according to the BTP, detect all the genes with the BTP-pattern */
	int maxlay=0;
	int *conditionL;
	cpp_dec_float_100 p=cpp_dec_float_100(0);
	int i, j, k, cons,rloc, cloc, colsig,rowsumm,colsumm;
	std::vector<int> sigcol(Mcols,0);
	std::fill(genes.begin(), genes.end(), false);
	for(i=0;i<cols;i++)
	{
		if(BTP[i]>maxlay)
		{
			maxlay=BTP[i];
		}
	}
	conditionL=(int *)malloc(sizeof(int)*(maxlay+1));
        for(i=0;i<maxlay+1;i++)
        {
                conditionL[i]=0;
        }
        for(i=0;i<cols;i++)
        {
            for(j=0;j<=maxlay;j++)
            {
                if(BTP[i]==j)
                {
                    conditionL[j]=conditionL[j]+1;
                }
            }
        }
	p=compute_p(conditionL,maxlay+1,po->TOLERANCE,cons);
        /* first, determine the corresponding columns in the ranking matrix with the BTP-pattern*/
        for(i=0;i<cols;i++)
        {
                if(BTP[i]<0)
                        continue;
            for(j=0;j<cols;j++)
            {
                    if(BTP[j]<0)
                            continue;
                if(BTP[i]<BTP[j])
                {
                     k = col_pair_to_idx[i][j];
                     if(k >= 0)
                         sigcol[k] = 1;
                }
            }
        }
        /* collect active column indices for sparse scanning */
        std::vector<int> active_cols;
        active_cols.reserve(cols*cols);
        for(cloc=0; cloc<Mcols; cloc++) {
            if(sigcol[cloc])
                active_cols.push_back(cloc);
        }
        /* determine the rows adhering to the BTP-pattern*/
        colsig=(int)std::floor(cons*po->TOLERANCE);
        /* extract the corresponding genes*/
        for(rloc=0;rloc<Mrows;rloc++)
        {
            colsumm=0;
            for(int ac : active_cols)
            {
                if(Marr[rloc][ac])
                    colsumm++;
                if(colsumm >= colsig)
                    break;
            }
            if(colsumm>=colsig)
                genes[rloc]=true;
            else
                genes[rloc]=false;
        }
            /* compute the p_value of the block */
            rowsumm=std::count(genes.begin(), genes.end(), true);
	    if(rowsumm>rows)
	    {
		    p_value=cpp_dec_float_100(1);
	    }else{
	    p_value = pbinom_at_least(rowsumm, rows, p); 
	    }
}

/********************************************************/
static bool detect_bicluster(std::unique_ptr<Block> &bc,std::vector<std::vector<unsigned char>> &allcolcand,const std::vector<int> &active_counts,const int min_width)
{
    if(allcolcand.size()==0)
       {return 0;}
	/* variables for adjacency matrix */
	bool **adj_matrix;
	adj_matrix=alloc2b(cols,cols);
	/* variables for counting */
	int rloc,cloc,j,k;
        /* variables for degree, in-degree and out-degree */
	std::vector<int> degree(cols,0);
	std::vector<int> indegree(cols,0);
        std::vector<int> outdegree(cols,0);
	/* variables for the BTP-pattern */
	std::vector<int> BTP(cols,0);
	int max_deg;
	int whichmax_deg;
	std::vector<int> cand_vet;
	std::vector<float> buck_value;
	std::vector<int> forcom_buck;
	std::vector<int> latcom_buck;
	std::vector<int> incom_buck;
	std::vector<int> BTP_vet;
	std::vector<int> cand_res;
	int innode;
	int max_index;
	float numer,denom;
	/* variables for the optimal BTP-pattern */
	std::vector<int> optBTP(cols,0);
	std::vector<bool> genes(Mrows,false);
	std::vector<bool> optgenes(Mrows,false);
	cpp_dec_float_100 p_value,minp_value;
	int genesum;
	int optgenesum;
	int BTPcount;
	int optBTPcount;

	minp_value=cpp_dec_float_100(1.0);

	for(std::size_t i=0;i<allcolcand.size();i++)
	{
	if (active_counts[i] == 0) continue;
	std::fill(BTP.begin(),BTP.end(),0);
	/* construct an adjacent matrix to record the identified ordered pairs */
	for(rloc=0;rloc<cols;rloc++)
	{
	      for(cloc=0;cloc<cols;cloc++)
	      {
		     adj_matrix[rloc][cloc]=false;
	       }
	 } 
          for(j=0; j<Mcols; j++)
         {
              if(allcolcand[i][j])
              {
                  rloc=Mcollabels[j]->condA;
                  cloc=Mcollabels[j]->condB;
                  adj_matrix[rloc][cloc]=true;
              }
	}
	/* The greedy algorithm for finding the largest CMAT subdigraph */
	/* First, compute the in-degree, out-degree, degree for the graph: flow_matrix*/
	std::fill(indegree.begin(),indegree.end(),0);
	std::fill(outdegree.begin(),outdegree.end(),0);
	std::fill(degree.begin(),degree.end(),0);
	for(rloc=0;rloc<cols;rloc++)
		for(cloc=0;cloc<cols;cloc++)
		{
			if(adj_matrix[rloc][cloc])
			{
			
				indegree[cloc]=indegree[cloc]+1;
				outdegree[rloc]=outdegree[rloc]+1;
			}
		}
	for(rloc=0;rloc<cols;rloc++)
        {
                degree[rloc]=indegree[rloc]+outdegree[rloc];
        }
	cand_vet.clear();
	for(j=0;j<cols;j++)
	{
		cand_vet.push_back(j);
	}
	
	for(j=0;j<cols;j++)
	{
		max_deg=0;
		whichmax_deg=0;
		for(int val:cand_vet){
			if(max_deg < degree[val])
			{
				max_deg=degree[val];
				whichmax_deg=val;
			}
		}
		buck_value.clear();
		if(max_value(BTP)==0)
		{
			BTP[whichmax_deg]=1;
		}
		else
		{
			for(k=1;k<=(2*max_value(BTP)+1);k++)
			{
				forcom_buck.clear();
				latcom_buck.clear();
				incom_buck.clear();
				innode=0;
				if(k%2==0)
				{
					for(rloc=0;rloc<cols;rloc++)
					{
						if(BTP[rloc]>0 && BTP[rloc]<(k/2))
						{forcom_buck.push_back(rloc);}
						if(BTP[rloc]>(k/2))
						{latcom_buck.push_back(rloc);}
						if(BTP[rloc]==k/2)
						{incom_buck.push_back(rloc);}
					}
					for(int val:incom_buck)
					{
						if(adj_matrix[whichmax_deg][val])
						{
							innode++;
						}
						if(adj_matrix[val][whichmax_deg])
						{
							innode++;
						}
					}
					if(forcom_buck.empty() && latcom_buck.empty())
					{
						buck_value.push_back(0.0);
					}
					else
					{
						numer=0.0;
						denom=0.0;
						for(int val:forcom_buck)
						{
							if(adj_matrix[val][whichmax_deg])
							{
								numer=numer+1.0;
							}
						}
						for(int val:latcom_buck)
						{
							if(adj_matrix[whichmax_deg][val])
							{
								numer=numer+1.0;
							}
						}
						numer=numer+incom_buck.size()-innode;
						denom=forcom_buck.size()+incom_buck.size()+latcom_buck.size();
						buck_value.push_back(numer/denom);
					}
				} else if(k % 2==1) {
				forcom_buck.clear();
                                latcom_buck.clear();
				for(rloc=0;rloc<cols;rloc++)
                                {
                                       if(BTP[rloc]>0 && BTP[rloc]<=(k-1)/2)
                                         {forcom_buck.push_back(rloc);}
                                        if(BTP[rloc]>=(k+1)/2)
                                          {latcom_buck.push_back(rloc);}
                                          
                                 }
				if(forcom_buck.empty() && latcom_buck.empty())
				{
					buck_value.push_back(0.0);
				}
				else
				{
					 numer=0.0;
                                         denom=0.0;
                                         for(int val:forcom_buck)
                                         {
                                                 if(adj_matrix[val][whichmax_deg])
                                                 {
                                                           numer=numer+1.0;
                                                 }
                                         }
                                         for(int val:latcom_buck)
                                         {
                                                  if(adj_matrix[whichmax_deg][val])
                                                  {
                                                           numer=numer+1.0;
                                                  }
                                         }
					 denom=forcom_buck.size()+latcom_buck.size();
					 buck_value.push_back(numer/denom);
				}
				}
			}
			if (*std::max_element(buck_value.begin(), buck_value.end()) > 0.8)
			{
				max_index=std::distance(buck_value.begin(), std::max_element(buck_value.begin(), buck_value.end()));
				max_index=max_index+1;
				if(max_index % 2 == 0) {
					BTP[whichmax_deg]=max_index/2;
				}
				else{
					for(cloc=0;cloc<cols;cloc++)
					{
						if(BTP[cloc]>=(max_index+1)/2)
						{
							BTP[cloc]=BTP[cloc]+1;
						}
					}
					BTP[whichmax_deg]=(max_index+1)/2;
				}
			}
		}
		cand_vet.clear();
		BTP_vet.clear();
		for(cloc=0;cloc<cols;cloc++)
		{
			if(BTP[cloc]>0)
			{
				BTP_vet.push_back(cloc);
			}
		}
		for(int val:BTP_vet)
		{
			for(cloc=0;cloc<cols;cloc++)
			{
				if(adj_matrix[cloc][val])
				{
					cand_vet.push_back(cloc);
				}
				if(adj_matrix[val][cloc])
				{
					cand_vet.push_back(cloc);
				}
			}
		}
		if(cand_vet.empty())
			break;
		std::sort(cand_vet.begin(),cand_vet.end());
		auto it=std::unique(cand_vet.begin(),cand_vet.end());
		cand_vet.erase(it,cand_vet.end());
		cand_res.clear();
		for(int val:cand_vet)
		{
			if(std::find(BTP_vet.begin(),BTP_vet.end(),val)==BTP_vet.end())
			{
				cand_res.push_back(val);
			}
		}
		cand_vet=cand_res;
	}
	std::transform(BTP.begin(), BTP.end(), BTP.begin(), [](int val) {return val-1; });
		if(max_value(BTP)<=0)
		{continue;}
		get_BTP_genes_and_pvalue(BTP,genes,p_value);
                std::transform(BTP.begin(), BTP.end(), BTP.begin(), [](int val) {return val+1;});
		genesum=std::count(genes.begin(), genes.end(), true);
                BTPcount=count_if(BTP.begin(),BTP.end(),[](int x) {
                                    return x>0;
                         });
		if(genesum>= po->ROW_MINWIDTH && genesum <=po->ROW_MAXWIDTH && BTPcount>=po->COL_WIDTH)
		{
			if(!isnan(p_value) && !isnan(minp_value))
			{
	                 	if(p_value < minp_value)
	                 	{
	                   		optBTP=BTP;
	                   		minp_value=p_value;
			   		optgenes=genes;
			 	}
			 	if(p_value == minp_value){
			 	optgenesum=std::count(optgenes.begin(), optgenes.end(), true);
			 	optBTPcount=count_if(optBTP.begin(),optBTP.end(),[](int x) {
					 return x > 0;
					 });
					if(genesum*BTPcount>optgenesum*optBTPcount){
					 optBTP=BTP;
					 minp_value=p_value;
					 optgenes=genes;
				 	}
			 	}	
			} else {
                         	optgenesum=std::count(optgenes.begin(), optgenes.end(), true);
                         	optBTPcount=count_if(optBTP.begin(),optBTP.end(),[](int x) {
                                         return x > 0;
                                         });
                                if(genesum*BTPcount>optgenesum*optBTPcount){
                                         optBTP=BTP;
                                         minp_value=p_value;
                                         optgenes=genes;
                                 }
			}
		}
	}
	bc->BTP=optBTP;
	bc->genes=optgenes;
	bc->p_value=minp_value;                       
        optgenesum=std::count(optgenes.begin(), optgenes.end(), true);
        optBTPcount=count_if(optBTP.begin(),optBTP.end(),[](int x) {
                             return x > 0;
                             });
	bc->size=optgenesum*optBTPcount;
        



	{if((isnan(bc->p_value) || bc->p_value < cpp_dec_float_100(po->pvalue)) && ((optgenesum >= po->ROW_MINWIDTH) && (optgenesum <=po->ROW_MAXWIDTH) && (optBTPcount>=po->COL_WIDTH)))
	{
	return 1;}
	else
	{
				return 0;}}


}
/********************************************************/
/* for a gene, update its genecans and genecanNums */
static void update_candidate_rows(std::unique_ptr<Blockbase> &b,const std::vector<unsigned char> &colcand,const std::vector<unsigned char> &candidates){
    std::vector<int> max_cnts(b->genes.size(), -1);
    std::vector<int> max_rows(b->genes.size(), -1);
    if (Mrows > 200) {
        #pragma omp parallel
        {
            std::vector<int> local_max_cnts(b->genes.size(), -1);
            std::vector<int> local_max_rows(b->genes.size(), -1);
            #pragma omp for schedule(dynamic)
            for (int row = 0; row < Mrows; row++) {
                if (!candidates[row])
                    continue;
                for (std::size_t gloc = 0; gloc < b->genes.size(); gloc++) {
                    const bool *genevec = Marr[b->genes[gloc]];
                    int cnt = intersect_row(colcand, genevec, Marr[row]);
                    if (cnt > local_max_cnts[gloc]) {
                        local_max_cnts[gloc] = cnt;
                        local_max_rows[gloc] = row;
                    }
                }
            }
            #pragma omp critical
            {
                for (std::size_t gloc = 0; gloc < b->genes.size(); gloc++) {
                    if (local_max_cnts[gloc] > max_cnts[gloc]) {
                        max_cnts[gloc] = local_max_cnts[gloc];
                        max_rows[gloc] = local_max_rows[gloc];
                    }
                }
            }
        }
    } else {
        for (int row = 0; row < Mrows; row++) {
            if (!candidates[row])
                continue;
            for (std::size_t gloc = 0; gloc < b->genes.size(); gloc++) {
                const bool *genevec = Marr[b->genes[gloc]];
                int cnt = intersect_row(colcand, genevec, Marr[row]);
                if (cnt > max_cnts[gloc]) {
                    max_cnts[gloc] = cnt;
                    max_rows[gloc] = row;
                }
            }
        }
    }
    for (std::size_t gloc = 0; gloc < b->genes.size(); gloc++) {
        b->genecans[gloc] = max_rows[gloc];
        b->genecanNums[gloc] = max_cnts[gloc];
    }
}

/*******************************************************/
static void update_blockcol(std::unique_ptr<Blockbase> &b,std::vector<unsigned char> &colcand,std::vector<unsigned short> &colcount,const int cand_threshold,
                          std::vector<std::vector<unsigned char>> &allcolcand, std::vector<int> &active_counts) {

    int i,j,col,max_row,max_cnt;
    int b_size=2;
  /* maintain a candidate list to avoid looping through all rows */
   std::vector<unsigned char> candidates = std::vector<unsigned char>(Mrows, 1);
    for (auto gene : b->genes) {
        candidates[gene] = false;
    }
    b->genecans.push_back(0);
    b->genecans.push_back(0);
    b->genecanNums.push_back(0);
    b->genecanNums.push_back(0);
  while (b_size < Mrows) {
        update_candidate_rows(b,colcand,candidates);
        max_row=-1;
        max_cnt=-1;
        for (i=0;i<b_size;i++)
        {
            if(b->genecanNums[i]>max_cnt)
            {
                max_row=b->genecans[i];
                max_cnt=b->genecanNums[i];
            }
        }
        /*if the common ordered pairs are too few, stop searching for the next row to add*/
        if(max_cnt<cand_threshold)
	{
            break;
	}
        /*we just get the t (t is from 5 to po->ROW_MAXINIT) largest rows to initial the common trend for a bicluster
         * because we believe that these rows can characterize the expression trend of the bicluster.
         * btw, it can reduce the time complexity*/
         if(b_size > po->ROW_MAXINIT)
            break;

        /* update colcount,colcand and candidate */
        int threshold = 8 * (b_size + 1);
        for (col = 0; col < Mcols; col++)
        {
            colcount[col]=colcount[col]+Marr[max_row][col];
            if(colcount[col] * 10 >= threshold)
            {
                colcand[col]=1;
            }
            else
            {
                colcand[col]=0;
            }
        }
	/* preserve colcand*/
	 if(b_size >= 5)
         {
             int cnt = 0;
             for(j=0;j<Mcols;j++)
             {
                 allcolcand[b_size-5][j]=colcand[j];
                 if(colcand[j]) cnt++;
             }
             active_counts[b_size-5] = cnt;
         }
        b->genes.push_back(max_row);
        b->genecans.push_back(0);
        b->genecanNums.push_back(0);
        candidates[max_row] = false;
	b_size++;
   }
 
   if(b_size>=5)
   allcolcand.resize(b_size-5);
   else
   allcolcand.resize(0);



}

/*************************************************************************/
static int intersect_row(const std::vector<unsigned char> &colcand, const bool *genevec, bool *g)
/*calculate the weight of two rows under all columns and certain columns respectively*/
{
    int cnt = 0;
    int col = 0;
    for (; col + 7 < Mcols; col += 8) {
        uint64_t gv, gg, cc;
        memcpy(&gv, genevec + col, 8);
        memcpy(&gg, g + col, 8);
        memcpy(&cc, colcand.data() + col, 8);
        cnt += __builtin_popcountll(gv & gg & cc);
    }
    for (; col < Mcols; col++) {
        if (colcand[col] && (genevec[col] & g[col]))
            cnt++;
    }
    return cnt;
}

/*************************************************************************/
static void init_colcand_and_colcount(const std::vector<int> &genes,std::vector<unsigned char> &colcand,std::vector<unsigned short> &colcount) {
  bool *g1 = Marr[genes[0]];
  bool *g2 = Marr[genes[1]];

  /*update intial colcand and colcount*/
  for (int col = 0; col < Mcols; col++) {
    if (g1[col] != 0 && g1[col] == g2[col]) {
      colcand[col] = 1;
    }
    colcount[col]=g1[col]+g2[col];
  }
}

/***********************************************************************/
static bool bic_detect(Edge *e, std::unique_ptr<Block> &bc, const int min_width,
                       const int max_rowinit, const int cand_threshold) {

    /* to indicate whether one block has been excavated */
    bool flag;

  /*you must allocate a struct if you want to use the pointers related to it*/
  std::unique_ptr<Blockbase> b(new Blockbase());

  /*initial the block b*/
  b->genes.push_back(e->gene_one);
  b->genes.push_back(e->gene_two);
  std::vector<unsigned char> colcand = std::vector<unsigned char>(Mcols, 0);
  std::vector<std::vector<unsigned char>> allcolcand(max_rowinit-4,std::vector<unsigned char>(Mcols,0));
  std::vector<int> active_counts(max_rowinit-4, 0);
  std::vector<unsigned short> colcount = std::vector<unsigned short>(Mcols, 0);
  init_colcand_and_colcount(b->genes,colcand,colcount);
  
  /* maintain a candidate list to avoid looping through all rows */
  update_blockcol(b, colcand,colcount, cand_threshold, allcolcand, active_counts);
  flag = detect_bicluster(bc,allcolcand,active_counts,min_width);
  return flag;
}

/***********************************************************************/
static bool isInBlocks(const std::vector<bool> &allincluster,int element) {
  return allincluster[element];
}

/************************************************************************/
/* Core algorithm */
int cluster(FILE *fw, Edge **el, int n, Edge **edge_list) {
  Edge *e;
  bool flag;
  std::vector<std::unique_ptr<Block>> bb;
  std::unique_ptr<Block> bc;
  std::size_t allocated = po->SCH_BLOCK;
  bb.reserve(allocated);
  std::vector<bool> allincluster(Mrows, false);

  int cand_threshold = std::max(2,po->COL_WIDTH);
    int i = 0;
    int kcc = 0;

  time_t last_ckpt_time = time(NULL);

  /* resume from checkpoint if available */
  if (checkpoint_enabled() && checkpoint_probe_phase() >= 2) {
      checkpoint_load_cluster_state(i, allincluster, bb);
      el += i;
      last_ckpt_time = time(NULL);
  }

	while (i++ < n)
	{
		/* periodic checkpoint: check at loop top so we don't miss it
		 * when iterations end with continue/break */
		if (checkpoint_enabled() && difftime(time(NULL), last_ckpt_time) >= po->CHECKPOINT_INTERVAL * 60) {
			checkpoint_save(2, edge_list, n, i - 1, allincluster, bb);
			last_ckpt_time = time(NULL);
		}

		e = *el++;

        		/* check if both genes already enumerated in previous blocks */
        		if ( isInBlocks(allincluster,e->gene_one) && isInBlocks(allincluster,e->gene_two) )
				continue;

			//test
          		// std::cout<<"e: "<<e->gene_one<<" "<<e->gene_two<<"\n";
        		// end!

		/*you must allocate a struct if you want to use the pointers related to it*/
			bc=std::unique_ptr<Block>(new Block());
		
        		/* detect the bicluster by figuring out its BTP-pattern */
        		flag=bic_detect(e, bc,po->COL_WIDTH, po->ROW_MAXINIT,cand_threshold);

			for(int genum=0;genum<rows;genum++)
			{
				if(bc->genes[genum]){
					allincluster[genum]=true;
					allincluster[(genum+rows)]=true;
				}
			}
			for(int genum=rows;genum<2*rows;genum++)
        		{
                		if(bc->genes[genum]){
                        		allincluster[genum]=true;
                        		allincluster[(genum-rows)]=true;
                		}
    			}
        		if(!flag)
        		{
            			bc.reset();
            			continue;
        		}

        		/*save the current block b to the block list bb so that we can sort the
        		* blocks by their score*/
        		bb.push_back(std::move(bc));

        		/* reaching the results number limit */
        		if (bb.size() == static_cast<size_t>(po->SCH_BLOCK)) 
			{
				/* save checkpoint before breaking so progress isn't lost */
				if (checkpoint_enabled()) {
					checkpoint_save(2, edge_list, n, i, allincluster, bb);
				}
				break;
			}
		}

    /* writes character to the current position in the standard output (stdout)
   * and advances the internal file position indicator to the next position. It
   * is equivalent to putc(character,stdout).*/

     putchar('\n');

     sort_block_list(bb);
     const int blocks = report_blocks(fw, bb, bb.size());
     return blocks;
}
/************************************************************************/
void print_params(FILE *fw) {
//  fprintf(fw, "# TransBic version %.1f output\n", VER);
  fprintf(fw, "# Parameters: -k %d -x %d -z %d -f %.2f -p %.5f -c %.2f -o %u", po->COL_WIDTH,
		  po->ROW_MAXWIDTH,po->ROW_MINWIDTH,po->FILTER,po->pvalue,po->TOLERANCE, po->RPT_BLOCK);
  fprintf(fw, "\n");
}
/************************************************************************/
static bool condcount(int a)
{
    return a>0;
}

static int interrowcount(std::vector<bool> a,std::vector<bool> b)
{
    int rcount=0;
    for(int m=0;m<rows;m++)
    {
            if(a[m]){
                    if(b[m] || b[(m+rows)])
                    {rcount=rcount+1;}
            }
    }
    for(int m=rows;m<2*rows;m++)
    {
            if(a[m]){
                    if(b[m] || b[(m-rows)])
                    {rcount=rcount+1;}
            }
    }
    return rcount;
}

static int intercondcount(std::vector<int> a,std::vector<int> b)
{
    int ccount=0;
    for(std::size_t m=0;m<a.size();m++)
    {
        if(a[m]*b[m]>0)
            ccount=ccount+1;
    }
    return ccount;
}
/************************************************************************/
static int report_blocks(FILE *fw, const std::vector<std::unique_ptr<Block>> &bb,
                  const std::size_t num) {
  print_params(fw);
  const int n = std::min(static_cast<int>(num),po->RPT_BLOCK);

  std::size_t *output = new std::size_t[n];

  std::size_t *bb_ptr = output;

  /* the major post-processing here, filter overlapping blocks*/
  std::size_t i = 0;
  int j = 0;
  while (i < num && j < n) {
    int index = i;
    const double cur_rows = std::count(bb[index]->genes.begin(),bb[index]->genes.end(),true);
    const double cur_cols = std::count_if(bb[index]->BTP.begin(),bb[index]->BTP.end(),condcount);

    bool flag = TRUE;
    int k = 0;
    while (k < j) {
      const double inter_rows =interrowcount(bb[output[k]]->genes,bb[index]->genes);
      const double inter_cols =intercondcount(bb[output[k]]->BTP,bb[index]->BTP);

      if (inter_rows * inter_cols > po->FILTER * cur_rows * cur_cols) {
        flag = FALSE;
        break;
      }
      k++;
    }
    i++;
    if (flag) {
      print_bc(fw, bb[index], j++);
      *bb_ptr++ = index;
    }
  }
  delete[] output;
  return j;
}
/************************************************************************/
 static bool scoreGreater (const std::unique_ptr<Block> &a,const std::unique_ptr<Block> &b)
{
	if (!isnan(a->p_value) && !isnan(b->p_value) && a->p_value != b->p_value) {
    		return (a->p_value < b->p_value);
	}
	return (a->size > b->size);
}
/************************************************************************/
static void sort_block_list(std::vector<std::unique_ptr<Block>> &el)
{
  std::stable_sort(el.begin(), el.end(), scoreGreater);
}
/************************************************************************/
