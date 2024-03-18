#pragma once
#include <algorithm>
#include <ctime>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include "graph_class.h"
using namespace std;

void generate_BER(mlnet& basis_net, vector<double>& arguments) {
    double p = arguments[0];
    double random_seed = arguments[1];

    std::random_device rd;
    std::mt19937 g(rd());

    std::random_device rand_dev;
    std::mt19937 rand(random_seed + 1); /// make this seed different from the seed for 
    //Rcpp::Rcout << "the random seed in initial sampling is: " << rand << "\n";

    
    std::uniform_real_distribution<> runif(0.0, 1.0);
    for (int i = 0; i < basis_net.node_count(); ++i) {
        for (int j = i + 1; j < basis_net.node_count(); ++j) {
            
            if (runif(rand) < p) basis_net.add_edge(i, j, 1);
            
        }
    }
}

void generate_LSM( mlnet& basis_net, vector<double>& arguments) {

    double fe = arguments[0]; 
    double nmu = arguments[1];
    double nsd = arguments[2];
    double random_seed = arguments[3];
    std::random_device rd;
    std::mt19937 g(rd());

    std::random_device rand_dev;
    std::mt19937 rand(random_seed + 2); /// make this seed different from the seed for 
    //Rcpp::Rcout << "the random seed in initial sampling is: " << rand << "\n";
    std::normal_distribution<> rnorm(nmu, nsd);
    std::uniform_real_distribution<> runif(0.0, 1.0);
    int N = basis_net.node_count();

    vector<vector<double> > node_pos;
    node_pos.resize(N);
    for (int i = 0; i < N; ++i) {
        node_pos[i].resize(2);
    }


    for (int i = 0; i < N; ++i) {
        node_pos[i][0] = rnorm(rand);
        node_pos[i][1] = rnorm(rand);
    }

    double odds;
    double p;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            odds = fe - sqrt(pow(node_pos[i][0] - node_pos[j][0], 2) + pow(node_pos[i][1] - node_pos[j][1], 2));
            p = exp(odds) / (1 + exp(odds));
            if (runif(rand) < p) basis_net.add_edge(i, j, 1);

        }
    }

}

void generate_SBM( mlnet& basis_net, vector<double>& arguments) {

    int B = arguments[0];
    double p_within = arguments[1];
    double p_between = arguments[2];
    double random_seed = arguments[3];
    std::random_device rd;
    std::mt19937 g(rd());

    std::random_device rand_dev;
    std::mt19937 rand(random_seed + 2);
    std::uniform_real_distribution<> runif(0.0, 1.0);


    vector<int> memb;
    int N = basis_net.node_count();
    memb.resize(N);
    for (int i = 0; i < N; ++i) {
        memb[i] = 1;
    }
    for (int b = 0; b < B; ++b) {
        for (int i = b * N / B; i < N / B * (b + 1); ++i) {
            memb[i] = b + 1;
        }
    }



    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (memb[i] == memb[j] && runif(rand) < p_within) {
                basis_net.add_edge(i, j, 1);

            }

            if (memb[i] != memb[j] && runif(rand) < p_between) {
                basis_net.add_edge(i, j, 1);

            }

        }
    }

}

void generate_other(mlnet& basis_net, vector<double>& arguments) {
    
    for (int i = 0; i <= arguments.size() - 2; i += 2) {
        basis_net.add_edge(int(arguments[i])-1, int(arguments[i + 1])-1, 1);
    }

}