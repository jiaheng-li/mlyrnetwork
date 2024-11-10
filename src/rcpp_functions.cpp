#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <random>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sstream>
#include <vector>
#include <set>  


#ifndef _estimation_class_
#define _estimation_class_ 
#include "estimation_class.h"
#endif 


#ifndef _simulation_class_
#define _simulation_class_
#include "simulation_class.h"
#endif

using namespace std;
using namespace Rcpp;


//' @title rcpp_compute_dyad_suffstats
//' @description
//' Compute dyadic sufficient statistics for multilayer networks with up to H-way interactions
//' @name rcpp_compute_dyad_suffstats
//' @param
//' @examples
//' 
//'
//' @export
// [[Rcpp::export]]
List rcpp_compute_dyad_suffstats(IntegerMatrix RNETWORK, IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers,
    int rhighest_order, int rand_seed, NumericVector arguments) {


    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int highest_order = rhighest_order;
    int random_seed = rand_seed;
    vector<double> basis_arguments;
    string mterms = Rcpp::as<string>(model_terms);
    int num_of_dyads = num_nodes * (num_nodes - 1) / 2;
    NumericMatrix dyad_suffstats(num_of_dyads, 2+model_dim);
    int node_i = RNETWORK(0, 0) - 1;
    int node_j = RNETWORK(0, 1) - 1;
    set<int> edges;
    string lexi_ind;
    string char_ele;
    set<int> edge_set;
    bool flag;
    // Store the mapping between lexicographical orders (under same order of interactions) and parameter indeces.
    //unordered_map<string, int> indexmap;
    unordered_map<string, set<int>> dyadmap;
    sim_ml_Hway sim_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, highest_order, random_seed, basis_arguments);
    
    
    sim_obj.compute_param_index();

    for (int i = 0; i < RNETWORK.nrow(); ++i) {
        if (node_i == RNETWORK(i, 0) - 1 && node_j == RNETWORK(i, 1) - 1) edges.insert(RNETWORK(i, 2) - 1);
        else {
            lexi_ind = "";
            lexi_ind = to_string(node_i) + "+" + to_string(node_j);
            dyadmap.insert(make_pair(lexi_ind, edges));
            edges.clear();

            node_i = RNETWORK(i, 0) - 1;
            node_j = RNETWORK(i, 1) - 1;
            i--;
        }
        //for (auto set_ele : edges) {
          //  cout << set_ele << ',' << endl;
        //}
    }
    

    // Add the last dyad.
    lexi_ind = to_string(node_i) + "+" + to_string(node_j);
    dyadmap.insert(make_pair(lexi_ind, edges));

    int row = 0;
    string node_ind = "";
    int param_dim = 0;
    for (int i = 0; i < num_nodes - 1; ++i) {
        for (int j = i+1; j < num_nodes; ++j) {
            
            dyad_suffstats(row, 0) = i+1;
            dyad_suffstats(row, 1) = j+1;
            node_ind = to_string(i) + "+" + to_string(j);
            edge_set = dyadmap[node_ind];
            
            for (auto ele = sim_obj.indexmap.begin() ; ele != sim_obj.indexmap.end(); ele++) {
                lexi_ind = ele->first;
                
                for (int iter = 0; iter < lexi_ind.size(); ++iter) {
                    flag = false;
                    while (lexi_ind[iter] != '+') { // use '+' here instead of "" becasue lexi_ind[iter] is a char 
                        char_ele += lexi_ind[iter];
                        iter++;
                    }
                    
                    if (edge_set.find(stoi(char_ele)) == edge_set.end()) {
                        flag = true;
                        char_ele = "";
                        break;
                    }
                    char_ele = "";

                }
                param_dim = sim_obj.indexmap[lexi_ind];
                if (flag) dyad_suffstats(row, 2+param_dim) = 0;
                else dyad_suffstats(row, 2+param_dim) = 1;
            }
            row++;
        }
    }

    

    List return_list = List::create(Named("dyad_suffstats") = dyad_suffstats);

    return return_list;

}




//' @title rcpp_estimate_model_ml_Hway
//' @description
//' Estimate (MPLE) parameters for network-separable multilayer network models with up to H-way interactions
//' @name rcpp_estimate_model_ml_Hway
//' @param
//' @examples
//' 
//'
//' @export
// [[Rcpp::export]]
List rcpp_estimate_model_ml_Hway(IntegerMatrix RNETWORK,
    IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers, IntegerVector rhighest_order, IntegerVector random_seeds, NumericVector arguments , IntegerVector itermax) {

    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int highest_order = rhighest_order[0];
    int random_seed = random_seeds[0];
    int iter_max = itermax[0];
    string mterms = Rcpp::as<string>(model_terms);
    vector<double> basis_arguments;
    for (int i = 0; i < arguments.size(); ++i) {
        basis_arguments.push_back(arguments[i]);
    }

    

    est_ml_Hway est_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, highest_order, random_seed, basis_arguments, iter_max);

    int node_i, node_j, layer_k;
    for (int i = 0; i < RNETWORK.nrow(); ++i) {
        node_i = RNETWORK(i, 0) - 1;
        node_j = RNETWORK(i, 1) - 1;
        layer_k = RNETWORK(i, 2) - 1;
        est_obj.obs_net_ml.add_edge(node_i, node_j, layer_k);
    }
    est_obj.compute_initial_estimate();
    NumericVector Rtheta_est(est_obj.theta_est.size());
    for (int p = 0; p < Rtheta_est.length(); ++p) {
        Rtheta_est(p) = est_obj.theta_est.at(p);
    }
    List return_list = List::create(Named("theta_est") = Rtheta_est,
        Named("model_terms") = model_terms);

    Rcpp::Rcout << "\n" << flush;
    //cout << "\nmade it to estimation check point 6\n";
    return return_list;

}





//' @title rcpp_simulate_ml_Hway
//' @description
//' Simulate a network-separable multilayer network model with arbitrary H-way interactions
//' @name rcpp_simulate_ml_Hway
//' @param
//' @examples
//' 
//' 
//' @export
// [[Rcpp::export]]
List rcpp_simulate_ml_Hway(IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers,
    NumericVector rtheta, int rhighest_order ,int rand_seed, NumericVector arguments) {


    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int highest_order = rhighest_order;
    int random_seed = rand_seed;
    vector<double> basis_arguments;
    string mterms = Rcpp::as<string>(model_terms);


    vector<double> theta;
    theta.resize(model_dim);
    
    for (int p = 0; p < model_dim; ++p) {
        theta.at(p) = rtheta[p];
    }


    for (int i = 0; i < arguments.size(); ++i) {
        basis_arguments.push_back(arguments[i]);
    }

    sim_ml_Hway sim_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, highest_order ,random_seed, basis_arguments);
    sim_obj.set_theta(theta);
    sim_obj.simulate_ml();
    //sim_obj.compute_sample_mean();



    // Make edge list

    // print network output for estimation
    int edge_count = 0;
    for (int i = 0; i < sim_obj.mlnetwork.node_count(); ++i) {
        for (int j = i + 1; j < sim_obj.mlnetwork.node_count(); ++j) {
            edge_count += sim_obj.mlnetwork.adj[i][j].size();
        }
    }


    NumericMatrix elist(edge_count, 3);
    int iter = 0;
    for (int node_i = 0; node_i < sim_obj.mlnetwork.node_count(); ++node_i) {
        for (int node_j = node_i + 1; node_j < sim_obj.mlnetwork.node_count(); ++node_j) {
            for (int loc = 0; loc < sim_obj.mlnetwork.adj[node_i][node_j].size(); ++loc) {

                elist(iter, 0) = node_i + 1;
                elist(iter, 1) = node_j + 1;
                elist(iter, 2) = sim_obj.mlnetwork.adj[node_i][node_j][loc] + 1;
                iter += 1;
            }
        }

    }



    // Print some statistics
    int m;
    int i = 0;
    m = sim_obj.netCount_samp.size() + sim_obj.nCr(sim_obj.mlnetwork.layer_count(), 2) + 1;
    NumericMatrix suff_list(m, 3);
    double count_avg;
    for (int iter = 0; iter < sim_obj.mlnetwork.layer_count(); ++iter) {
        suff_list(i, 0) = iter;
        suff_list(i, 1) = iter;
        count_avg = (double)accumulate(sim_obj.netCount_samp[iter].begin(), sim_obj.netCount_samp[iter].end(), 0) / (double)sim_obj.num_samples;
        suff_list(i, 2) = count_avg;
        i += 1;
    }
    for (int iter = 0; iter < sim_obj.mlnetwork.layer_count(); ++iter) {
        for (int iter2 = iter+1; iter2 < sim_obj.mlnetwork.layer_count(); ++iter2) {
            count_avg = (double)accumulate(sim_obj.crossLayerCount_samp[iter][iter2].begin(), sim_obj.crossLayerCount_samp[iter][iter2].end(), 0) / (double)sim_obj.num_samples;
            suff_list(i, 0) = iter;
            suff_list(i, 1) = iter2;
            suff_list(i, 2) = count_avg;
            i += 1;
        }
    }

    // when K > 3, the threewayCount is counting for the number of all dyads that have 3 or more edges, ignoring the specific dimension.
    count_avg = (double)accumulate(sim_obj.threewayCount_samp.begin(), sim_obj.threewayCount_samp.end(), 0) / (double)sim_obj.num_samples;
    suff_list(i, 0) = 3;
    suff_list(i, 1) = 3;
    suff_list(i, 2) = count_avg;


    // print basis network
    edge_count = 0;
    for (int i = 0; i < sim_obj.basis_net.node_count(); ++i) {
        for (int j = i + 1; j < sim_obj.basis_net.node_count(); ++j) {
            edge_count += sim_obj.basis_net.adj[i][j].size();
        }
    }
    NumericMatrix basis_list(edge_count, 3);
    iter = 0;
    for (int node_i = 0; node_i < sim_obj.basis_net.node_count(); ++node_i) { 
        for (int node_j = node_i + 1; node_j < sim_obj.basis_net.node_count(); ++node_j) {
            for (int loc = 0; loc < sim_obj.basis_net.adj[node_i][node_j].size(); ++loc) {

                basis_list(iter, 0) = node_i + 1;
                basis_list(iter, 1) = node_j + 1;
                basis_list(iter, 2) = sim_obj.basis_net.adj[node_i][node_j][loc];  // in basis network, 1 means a present edge.
                iter += 1;
            }
        }

    }



    // Create return list
    List return_list = List::create(Named("elist") = elist, Named("suff_list") = suff_list, Named("basis_list") = basis_list);

    return return_list;
}
