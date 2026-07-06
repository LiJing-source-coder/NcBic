#include "enum_gene_pvalue.h"
#include <unordered_map>

struct VectorHash {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t hash = v.size();
        for (int i : v) {
            hash = hash * 31 + i;
        }
        return hash;
    }
};

/***********************************************************************/


/*Enumerate all the cases of implanting conditions from buck(one condition bucket) into layer (the union of its front buckets).
Compute all the probabilities of different numbers of 1's (denoted by layer_prob[n]).Each value of layer_prob[n] represents the probability
of the corresponding number of 1's (its location).
In the matrix num_1, its each element records the total number of 1's of these added conditions in sequence, simply called a condition sequence.
To avoid counting the repetitive sequences, it is required that any element in the sequence has the number of 1's not more than its last element.
Different rows correspond to different number of 1's of the last element of the condition sequence.
And count_1 records their occurrence number for these condition sequences at the same locations.
nnum_1 and ncount_1 are used to update num_1 and count_1 respectively.
We iteratively implant conditions and compute num_1, count_1*/
void enum_two_layers(int layer,int buck,cpp_dec_float_100 *layer_prob,int n)
{
    /* Fast path: buck == 1 -> uniform over {0,...,layer} */
    if (buck == 1) {
        std::fill(layer_prob, layer_prob + n, cpp_dec_float_100(0));
        cpp_dec_float_100 p = cpp_dec_float_100(1) / cpp_dec_float_100(layer + 1);
        for (int i = 0; i <= layer && i < n; i++) {
            layer_prob[i] = p;
        }
        return;
    }
    /* Fast path: layer == 1 -> uniform over {0,...,buck} */
    if (layer == 1) {
        std::fill(layer_prob, layer_prob + n, cpp_dec_float_100(0));
        cpp_dec_float_100 p = cpp_dec_float_100(1) / cpp_dec_float_100(buck + 1);
        for (int i = 0; i <= buck && i < n; i++) {
            layer_prob[i] = p;
        }
        return;
    }

    int i,j,x,y,s_1;
    std::vector<cpp_dec_float_100> count_1((layer+1) * n, cpp_dec_float_100(0));
    std::vector<cpp_dec_float_100> ncount_1((layer+1) * n, cpp_dec_float_100(0));
    std::vector<std::vector<int>> active(layer+1);
    
    for(i=0;i<layer+1;i++) {
        count_1[i*n+i] = cpp_dec_float_100(1);
        active[i].push_back(i);
    }
    
    if(buck != 1) {
        for(i=2;i<=buck;i++) {
            std::vector<std::vector<int>> nactive(layer+1);
            for(j=0;j<=layer;j++) {
                for(x=j;x<layer+1;x++) {
                    int base = x * n;
                    int nbase = j * n;
                    for(int y_idx : active[x]) {
                        s_1 = y_idx + j;
                        if (ncount_1[nbase+s_1] == 0) {
                            nactive[j].push_back(s_1);
                        }
                        ncount_1[nbase+s_1] += count_1[base+y_idx];
                    }
                }
            }
            std::copy(ncount_1.begin(), ncount_1.end(), count_1.begin());
            std::fill(ncount_1.begin(), ncount_1.end(), cpp_dec_float_100(0));
            active.swap(nactive);
        }
    }
    
    std::fill(layer_prob, layer_prob + n, cpp_dec_float_100(0));
    for(i=0;i<layer+1;i++) {
        int base = i * n;
        for(int y_idx : active[i]) {
            layer_prob[y_idx] += count_1[base+y_idx];
        }
    }
    cpp_dec_float_100 total = 0;
    for (i = 0; i < n; ++i) {
        total += layer_prob[i];
    }
    for(i=0;i<n;i++) {
        layer_prob[i] = layer_prob[i] / total;
    }
}

/*Enumerate all the cases of different numbers of 1's and compute the probabilities mlayer_prob[n].
Each value of mlayer_prob[n] represents the probability of the corresponding number of 1's (its location).*/
/* Fast helper for compute_confidence: computes all-ones inversion distribution using long double */
static void enum_all_layers_allones_ld(int m, long double *mlayer_prob, int n)
{
    std::vector<long double> old(n);
    mlayer_prob[0] = 1.0L;
    int current_max = 0;
    for (int i = 1; i < m; i++) {
        long double p = 1.0L / (long double)(i + 1);
        std::copy(mlayer_prob, mlayer_prob + n, old.begin());
        std::fill(mlayer_prob, mlayer_prob + n, 0.0L);
        long double sum = 0;
        for (int k = 0; k <= current_max + i; k++) {
            if (k <= current_max) sum += old[k];
            if (k > i) sum -= old[k - i - 1];
            mlayer_prob[k] = sum * p;
        }
        current_max += i;
    }
}

void enum_all_layers(const int *layer_num,int m, cpp_dec_float_100 *mlayer_prob,int n)
{
    /* Fast path: all layers have size 1 -> distribution of inversions */
    bool all_ones = true;
    for (int idx = 0; idx < m; idx++) {
        if (layer_num[idx] != 1) {
            all_ones = false;
            break;
        }
    }
    if (all_ones) {
        std::vector<cpp_dec_float_100> old(n);
        std::fill(mlayer_prob, mlayer_prob + n, cpp_dec_float_100(0));
        mlayer_prob[0] = cpp_dec_float_100(1);
        int current_max = 0;
        for (int i = 1; i < m; i++) {
            cpp_dec_float_100 p = cpp_dec_float_100(1) / cpp_dec_float_100(i + 1);
            std::copy(mlayer_prob, mlayer_prob + n, old.begin());
            std::fill(mlayer_prob, mlayer_prob + n, cpp_dec_float_100(0));
            cpp_dec_float_100 sum = 0;
            for (int k = 0; k <= current_max + i; k++) {
                if (k <= current_max) sum += old[k];
                if (k > i) sum -= old[k - i - 1];
                mlayer_prob[k] = sum * p;
            }
            current_max += i;
        }
        return;
    }

    int i,j,k,buck,prenum_1,num_1,layer;
    std::vector<cpp_dec_float_100> layer_probA(n, cpp_dec_float_100(0));
    std::vector<cpp_dec_float_100> layer_probB(n, cpp_dec_float_100(0));
    std::vector<cpp_dec_float_100> layer_probC(n, cpp_dec_float_100(0));
    /*initialization(i.e.m=2)*/
    enum_two_layers(layer_num[0],layer_num[1],layer_probA.data(),layer_num[0]*layer_num[1]+1);
    prenum_1=layer_num[0]*layer_num[1]+1;
    if(m>2)
    {
        for(i=2;i<m;i++)
        {
            layer=0;
            for(j=0;j<i;j++)
            {
                layer=layer+layer_num[j];
            }
            buck=layer_num[i];
            num_1=layer*buck+1;
            enum_two_layers(layer,buck,layer_probB.data(),num_1);
            for(j=0;j<prenum_1;j++)
            {
                for(k=0;k<num_1;k++)
                {
                    layer_probC[j+k]=layer_probC[j+k]+layer_probA[j]*layer_probB[k];
                }
            }
            prenum_1=prenum_1+num_1-1;
	    std::copy_n(layer_probC.begin(), n, layer_probA.begin());
	    std::fill(layer_probB.begin(), layer_probB.begin() + n, cpp_dec_float_100(0));
	    std::fill(layer_probC.begin(), layer_probC.begin() + n, cpp_dec_float_100(0));
        }
    }
    std::copy_n(layer_probA.begin(), n, mlayer_prob);
}

static std::unordered_map<std::vector<int>, std::pair<cpp_dec_float_100, int>, VectorHash> compute_p_cache;

cpp_dec_float_100 compute_p(const int *conditionL,int m,double tolerance,int &cons)
{
	std::vector<int> key(conditionL, conditionL + m);
	auto it = compute_p_cache.find(key);
	if (it != compute_p_cache.end()) {
		cons = it->second.second;
		return it->second.first;
	}
	cpp_dec_float_100 p_value;
	int i,j,sesum,colsig;
	cons=0;
	for(i=0;i<(m-1);i++)
        {
            sesum=0;
            for(j=i+1;j<=(m-1);j++)
            {
                sesum=sesum+conditionL[j];
            }
            cons=cons+conditionL[i]*sesum;
        }
	std::vector<cpp_dec_float_100> row_prob(cons + 1, cpp_dec_float_100(0));
	enum_all_layers(conditionL,m,row_prob.data(),cons+1);
	colsig=(int)std::floor(cons*tolerance);
	for(j=0;j<cons+1;j++)
        {
                if((cons-j)<colsig)
                {
                    break;
                }
                p_value=p_value+row_prob[j]+row_prob[cons-j];
        }
	if(!isnan(p_value)){
		        if(p_value > cpp_dec_float_100(1.0)){
				p_value = cpp_dec_float_100(1.0);}}
	compute_p_cache[key] = {p_value, cons};
	return p_value;
}

int compute_confidence(const int *conditionL,int m,double repvalue,int &cons)
{
	int i,j,sesum;
	cons=0;
        for(i=0;i<(m-1);i++)
        {
            sesum=0;
            for(j=i+1;j<=(m-1);j++)
            {
                sesum=sesum+conditionL[j];
            }
            cons=cons+conditionL[i]*sesum;
        }
	/* Fast path for all-ones: use long double for 30x speedup, verified to give same integer result */
	bool all_ones = true;
	for (i = 0; i < m; i++) {
		if (conditionL[i] != 1) {
			all_ones = false;
			break;
		}
	}
	if (all_ones) {
		const long double threshold = (long double)1.0 - (long double)repvalue;
		long double pcons = 0;
		std::vector<long double> row_prob(cons + 1, 0.0L);
		enum_all_layers_allones_ld(m, row_prob.data(), cons + 1);
		for (i = 0; i < cons + 1; i++) {
			pcons = pcons + row_prob[i] + row_prob[cons - i];
			if (pcons >= threshold) {
				break;
			}
		}
		return i;
	}
	const cpp_dec_float_100 threshold = cpp_dec_float_100(1) - cpp_dec_float_100(repvalue);
	cpp_dec_float_100 pcons;
	std::vector<cpp_dec_float_100> row_prob(cons + 1, cpp_dec_float_100(0));
	enum_all_layers(conditionL,m,row_prob.data(),cons+1);
	for(i=0;i<cons+1;i++)
        {
                pcons=pcons+row_prob[i]+row_prob[cons-i];
                if(pcons >= threshold)
                {
                        break;
                }
        }
	return i;
}

cpp_dec_float_100 pbinom_at_least(int k, int n, cpp_dec_float_100 p) {
    binomial_distribution<cpp_dec_float_100> dist(n, p);
    /* compute P(X >= k) = P(X > k-1) = 1 - P(X <= k-1) */
    return cdf(complement(dist, k-1));
}
      





