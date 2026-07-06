#ifndef CLUSTER_H
#define CLUSTER_H

#include "struct.h"
#include "write_block.h"
#include <memory>
#include<boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace boost;
using namespace boost::multiprecision;
using namespace boost::math;

/* biclustering block */
typedef struct Blockbase{
  std::vector<int> genes;
  std::vector<int> conds;
  std::vector<int> genecans;
  std::vector<int> genecanNums;
}Blockbase;

/*construct ranking matrix*/
extern bool** alloc2b(int rr, int cc);

/*from enum_gene_pvalue: enumerate all the p-values for genes*/
extern cpp_dec_float_100 compute_p(const int *conditionL,int m,double tolerance,int &cons);
extern cpp_dec_float_100 pbinom_at_least(int k, int n, cpp_dec_float_100 p);

/* prototypes */
int** alloc2e(int rr, int cc);
static int max_value(std::vector<int> vec);
static void get_BTP_genes_and_pvalue(std::vector<int> BTP, std::vector<bool> &genes, cpp_dec_float_100 &p_value);
static bool detect_bicluster(std::unique_ptr<Block> &bc,std::vector<std::vector<unsigned char>> &allcolcand,const int min_width);
static void update_candidate_rows(std::unique_ptr<Blockbase> &b,const std::vector<unsigned char> &colcand,const std::vector<unsigned char> &candidates);
static void update_blockcol(std::unique_ptr<Blockbase> &b,std::vector<unsigned char> &colcand,std::vector<unsigned short> &colcount,const int cand_threshold,std::vector<std::vector<unsigned char>> &allcolcand);
static int intersect_row(const std::vector<unsigned char> &colcand, const bool *genevec, bool *g);
static void init_colcand_and_colcount(const std::vector<int> &genes,std::vector<unsigned char> &colcand,std::vector<unsigned short> &colcount);
static bool bic_detect(Edge *e, std::unique_ptr<Block> &bc, const int min_width,const int max_rowinit, const int cand_threshold);
static bool isInBlocks(const std::vector<bool> &allincluster,int element);
int cluster(FILE *fw, Edge **el, int n, Edge **edge_list = nullptr);
void print_params(FILE *fw);
static bool condcount(int a);
static int interrowcount(std::vector<bool> a,std::vector<bool> b);
static int intercondcount(std::vector<int> a,std::vector<int> b);
static int report_blocks(FILE *fw, const std::vector<std::unique_ptr<Block>> &bb,const std::size_t num);
static bool scoreGreater (const std::unique_ptr<Block> &a,const std::unique_ptr<Block> &b);
static void sort_block_list(std::vector<std::unique_ptr<Block>> &el);
#endif
