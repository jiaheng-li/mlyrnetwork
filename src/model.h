#pragma once
#include <algorithm>
#include <vector>
#include <string> 
#include <omp.h>
#include <iostream>
#include <stdlib.h>

#ifndef _helper_functions_
#define _helper_functions_
#include "helper_functions.h"
#endif

#ifndef _change_stats_
#define _change_stats_
#include "change_stats.h"
#endif


using namespace std;


class mod {

public:
    vector<string> model_terms;
    vector<double (*)(int i, int j, mlnet& network, int k, int l)> change_stat_funs_ml;


public:
    mod(vector<string> mterms) {
        model_terms.resize(mterms.size());
        model_terms = mterms;
        change_stat_funs_ml.resize(get_num_terms());
        int iter = 0;
        int block_counter = 0;
        int node_counter = 0;
        for (string term : model_terms) {
            
            if (term == "ml_order2") { /// calculate change statistics for ml network with 2nd order cross-layer inetractions
                change_stat_funs_ml[iter] = cs_order2;
                iter += 1;
            }
            if (term == "ml_order3") { /// calculate change statistics for ml network with 3rd order cross-layer inetractions
                change_stat_funs_ml[iter] = cs_order3;
                iter += 1;
            }
        }
    }

    int get_num_terms() {
        int val = model_terms.size();
        return val;
    }

};
