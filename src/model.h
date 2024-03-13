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

#ifndef _basis_net_
#define _basis_net_
#include "basis_net.h"
#endif

using namespace std;


class mod {

public:
    string model_terms;
    void (*populate_basisnet_funs)(mlnet& network, vector<double>& arguments) ;


public:
    mod(string mterms) {
        model_terms = mterms;
        if (model_terms == "SBM") populate_basisnet_funs = generate_SBM;
        if (model_terms == "LSM") populate_basisnet_funs = generate_LSM;
        if (model_terms == "BER") populate_basisnet_funs = generate_BER;
        if (model_terms == "other") populate_basisnet_funs = generate_other;
        
    }

    

    

};
