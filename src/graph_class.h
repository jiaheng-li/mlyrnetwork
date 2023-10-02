#pragma once
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <Rcpp.h>
#include "helper_functions.h"
using namespace std;


class mlnet {

public:
    int num_nodes;
    int num_layers;
    vector<vector<vector<int> > > adj; ///need to store edges in different layers

public:
    mlnet(int N, int K)
    {
        num_nodes = N;
        num_layers = K;
        adj.resize(num_nodes);
        for (int i = 0;i < num_nodes;++i) {
            adj[i].resize(num_nodes);
        }
    }

    int node_count() {
        return num_nodes;
    }

    int layer_count() {
        return num_layers;
    }

    void add_edge(int i, int j, int k) {
        adj[i][j].push_back(k);
        adj[j][i].push_back(k);
    }

    bool is_edge(int i, int j, int k) {
        bool val = false;
        val = is_in(adj[i][j], k);
        return val;
    }

    void delete_edge(int i, int j, int k) {
        vector<int>::iterator ind;
        ind = find(adj[i][j].begin(), adj[i][j].end(), k);
        swap(*ind, adj[i][j].back());
        adj[i][j].pop_back();
        ind = find(adj[j][i].begin(), adj[j][i].end(), k);
        swap(*ind, adj[j][i].back());
        adj[j][i].pop_back();
    }
};


