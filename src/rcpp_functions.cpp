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


vector<int> all_interaction_layer;
vector<int> combination;
vector<vector <int> > selected_layer;
// Select all combinations of k layers, where k is the order of interaction
void select_layer(int offset, int k) {
    if (k == 0) {
        selected_layer.push_back(combination);
        return;
    }
    for (int i = offset; i <= all_interaction_layer.size() - k; ++i) {
        combination.push_back(all_interaction_layer[i]);
        select_layer(i + 1, k - 1);
        combination.pop_back();
    }
}

// Compute the parameter index from lexicographical orders
void compute_param_index(int K, int H, unordered_map<string, int>& indexmap) {
    int param_ind = 0;
    for (int h = 0; h < K; ++h) {
        indexmap.insert(make_pair(to_string(h), param_ind++));
    }
    for (int h = 2; h <= H; ++h) {
        select_layer(0, h);
        for (auto ele : selected_layer)
        {
            string lexi_ind = "";
            for (auto ele2 : ele) {
                lexi_ind += to_string(ele2);
            }
            indexmap.insert(make_pair(lexi_ind, param_ind++));
        }


        // Re-initialize for next use
        int s = selected_layer.size();
        for (int num = 0; num < s; ++num) {
            selected_layer.pop_back();
        }


        s = combination.size();
        for (int num = 0; num < s; ++num) {
            combination.pop_back();
        }
    }
}



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
List rcpp_compute_dyad_suffstats(IntegerMatrix RNETWORK, IntegerVector rmodel_dim,
    IntegerVector rnum_nodes, IntegerVector rnum_layers, IntegerVector rhighest_order) {

    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int highest_order = rhighest_order[0];
    int num_of_dyads = num_nodes * (num_nodes - 1) / 2;
    NumericMatrix dyad_suffstats(num_of_dyads, 2+model_dim);
    int node_i = 0;
    int node_j = 1;
    set<int> edges;
    string lexi_ind;
    char char_ele;
    set<int> edge_set;
    bool flag;
    // Store the mapping between lexicographical orders (under same order of interactions) and parameter indeces.
    unordered_map<string, int> indexmap;
    unordered_map<string, set<int>> dyadmap;
    
    for (int e = 0; e < num_layers; ++e) {
        all_interaction_layer.push_back(e);
    }
    compute_param_index(num_layers, highest_order, indexmap);

    for (int i = 0; i < RNETWORK.nrow(); ++i) {
        if (node_i == RNETWORK(i, 0) - 1 && node_j == RNETWORK(i, 1) - 1) edges.insert(RNETWORK(i, 2) - 1);
        else {
            lexi_ind = "";
            lexi_ind = to_string(node_i) + to_string(node_j);
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
    lexi_ind = to_string(node_i) + to_string(node_j);
    dyadmap.insert(make_pair(lexi_ind, edges));

    int row = 0;
    string node_ind = "";
    int param_dim = 0;
    for (int i = 0; i < num_nodes - 1; ++i) {
        for (int j = i+1; j < num_nodes; ++j) {
            
            dyad_suffstats(row, 0) = i+1;
            dyad_suffstats(row, 1) = j+1;
            node_ind = to_string(i) + to_string(j);
            for (auto ele = indexmap.begin() ; ele != indexmap.end(); ele++) {
                lexi_ind = ele->first;
                //cout << lexi_ind << endl;
                edge_set = dyadmap[node_ind];
                //for (auto set_ele : edge_set) {
                //    cout << set_ele << ',';
                //}
                //cout << endl;
                
                for (int iter = 0; iter < lexi_ind.size(); ++iter) {
                    flag = false;
                    char_ele = lexi_ind[iter];
                    if (edge_set.find(int(char_ele - '0')) == edge_set.end()) {
                        flag = true;
                        break;
                    }

                }
                param_dim = indexmap[lexi_ind];
                if (flag) dyad_suffstats(row, 2+param_dim) = 0;
                else dyad_suffstats(row, 2+param_dim) = 1;
            }
            row++;
        }
    }

    for (int e = 0; e < all_interaction_layer.size(); ++e) {
        all_interaction_layer.pop_back();
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


/*
//' @title rcpp_estimate_model_ml_3way
//' @description
//' Estimate (MPLE) parameters for network-separable multilayer network models with 3-way interactions
//' @name rcpp_estimate_model_ml_3way
//' @param
//' @examples
//' 
//'
//' @export
// [[Rcpp::export]]
List rcpp_estimate_model_ml_3way(IntegerMatrix RNETWORK, NumericVector rNR_tol, IntegerVector rNR_max, IntegerVector rMCMLE_max,
    IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers, LogicalVector rcheck_chull, IntegerVector random_seeds, NumericVector g) {

    int NR_max = rNR_max[0];
    int MCMLE_max = rMCMLE_max[0];
    double NR_tol = rNR_tol[0];
    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    bool check_chull = rcheck_chull[0];
    int random_seed = random_seeds[0];
    double gy = g[0];
    
    vector<vector<int> > MEMB;
    MEMB.resize(num_nodes);
    
    vector<string> mterms;
    mterms.resize(model_dim);
    for (int p = 0; p < model_dim; ++p) {
        mterms.at(p) = model_terms[p];
    }

    est_ml_3way est_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, random_seed, NR_max, NR_tol, MCMLE_max, check_chull, gy);
    int node_i, node_j, layer_k;
    for (int i = 0; i < RNETWORK.nrow(); ++i) {
        node_i = RNETWORK(i, 0) - 1;
        node_j = RNETWORK(i, 1) - 1;
        layer_k = RNETWORK(i, 2) - 1;
        est_obj.obs_net_ml.add_edge(node_i, node_j, layer_k);
    }

   
    est_obj.compute_initial_estimate();
    // Create return list
    NumericVector Rtheta_est(est_obj.theta_est.size());
    for (int p = 0; p < Rtheta_est.length(); ++p) {
        Rtheta_est(p) = est_obj.theta_est.at(p);
    }
    List return_list = List::create(Named("theta_est") = Rtheta_est,
        Named("model_terms") = model_terms);

    Rcpp::Rcout << "\n" << flush;
    return return_list;
}

*/

/*
//' @title rcpp_estimate_model_ml
//' @description
//' Estimate (MPLE) parameters for network-separable multilayer network models with 2-way interactions
//' @name rcpp_estimate_model_ml
//' @param
//' @examples
//' 
//'
//' @export
// [[Rcpp::export]]
List rcpp_estimate_model_ml(IntegerMatrix RNETWORK, NumericVector rIter_max,
    IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers, LogicalVector rcheck_chull, IntegerVector random_seeds, NumericVector g) {

    int Iter_max = rIter_max[0];
    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    bool check_chull = rcheck_chull[0];
    int random_seed = random_seeds[0];
    double gy = g[0];
    vector<vector<int> > MEMB;
    MEMB.resize(num_nodes);
    vector<string> mterms;
    mterms.resize(model_dim);
    for (int p = 0; p < model_dim; ++p) {
        mterms.at(p) = model_terms[p];
    }

    est_ml est_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, random_seed, Iter_max, check_chull, gy);

    // Copy observed network to estimation object
    int node_i, node_j, layer_k;
    for (int i = 0; i < RNETWORK.nrow(); ++i) {
        node_i = RNETWORK(i, 0) - 1;
        node_j = RNETWORK(i, 1) - 1;
        layer_k = RNETWORK(i, 2) - 1;
        est_obj.obs_net_ml.add_edge(node_i, node_j, layer_k);
    }

    est_obj.compute_initial_estimate();
   
    // Create return list
    NumericVector Rtheta_est(est_obj.theta_est.size());
    for (int p = 0; p < Rtheta_est.length(); ++p) {
        Rtheta_est(p) = est_obj.theta_est.at(p);
    }
    List return_list = List::create(Named("theta_est") = Rtheta_est,
        Named("model_terms") = model_terms);

    Rcpp::Rcout << "\n" << flush;
    return return_list;
}

*/



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
    sim_obj.compute_sample_mean();



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
        suff_list(i++, 2) = count_avg;
    }
    for (int iter = 0; iter < sim_obj.mlnetwork.layer_count(); ++iter) {
        for (int iter2 = iter+1; iter2 < sim_obj.mlnetwork.layer_count(); ++iter2) {
            count_avg = (double)accumulate(sim_obj.crossLayerCount_samp[iter][iter2].begin(), sim_obj.crossLayerCount_samp[iter][iter2].end(), 0) / (double)sim_obj.num_samples;
            suff_list(i, 0) = iter;
            suff_list(i, 1) = iter2;
            suff_list(i++, 2) = count_avg;

        }

    }
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

/*
//' @title rcpp_simulate_ml_LSM
//' @description
//' Simulate a network-separable multilayer network model with basis network generated by latent space model 
//' @name rcpp_simulate_ml_LSM
//' @param
//' @examples
//' 
//' 
//' @export
// [[Rcpp::export]]
List rcpp_simulate_ml_LSM(IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers,
    NumericVector rtheta, int rand_seed, NumericVector g, double fe, double nmu, double nsd) {
    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int random_seed = rand_seed;
    double gy = g[0];
    double fixed_effect = fe;
    double norm_mu = nmu;
    double norm_sd = nsd;


    vector<double> theta;
    theta.resize(model_dim);
    vector<string> mterms;
    mterms.resize(model_dim);
    for (int p = 0; p < model_dim; ++p) {
        mterms.at(p) = model_terms[p];
        theta.at(p) = rtheta[p];
    }

    sim_ml_LSM sim_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, random_seed, gy, fixed_effect, norm_mu, norm_sd);
    sim_obj.set_theta(theta);
    sim_obj.simulate_ml();
    sim_obj.compute_sample_mean();
   

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
        for (int iter2 = iter; iter2 < sim_obj.mlnetwork.layer_count(); ++iter2) {

            suff_list(i, 0) = iter;
            suff_list(i, 1) = iter2;
            if (iter == iter2) {
                count_avg = (double)accumulate(sim_obj.netCount_samp[iter].begin(), sim_obj.netCount_samp[iter].end(), 0) / (double)sim_obj.num_samples;

            }
            else {
                count_avg = (double)accumulate(sim_obj.crossLayerCount_samp[iter][iter2].begin(), sim_obj.crossLayerCount_samp[iter][iter2].end(), 0) / (double)sim_obj.num_samples;
            }
            suff_list(i, 2) = count_avg;
            ++i;

        }

    }
    count_avg = (double)accumulate(sim_obj.threewayCount_samp.begin(), sim_obj.threewayCount_samp.end(), 0) / (double)sim_obj.num_samples;
    suff_list(6, 0) = 3;
    suff_list(6, 1) = 3;
    suff_list(6, 2) = count_avg;


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

*/



/*
//' @title rcpp_simulate_ml_SBM
//' @description
//' Simulate a network-separable multilayer network model with basis network generated by stochastic block model
//' @name rcpp_simulate_ml_SBM
//' @param
//' @examples
//' 
//' 
//' @export
// [[Rcpp::export]]
List rcpp_simulate_ml_SBM(IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers,
    NumericVector rtheta, int rand_seed, NumericVector g, int block, double p1, double p2) {
    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int random_seed = rand_seed;
    double gy = g[0];
    int block_num = block;
    double p_in = p1;
    double p_between = p2;


    vector<double> theta;
    theta.resize(model_dim);
    vector<string> mterms;
    mterms.resize(model_dim);
    for (int p = 0; p < model_dim; ++p) {
        mterms.at(p) = model_terms[p];
        theta.at(p) = rtheta[p];
    }

    sim_ml_SBM sim_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, random_seed, gy, block_num, p_in, p_between);
    sim_obj.set_theta(theta);
    sim_obj.simulate_ml();
    sim_obj.compute_sample_mean();
   

    // Make edge list
    
    // print network output for estimation
    int edge_count = 0;
    for (int i = 0; i < sim_obj.mlnetwork.node_count(); ++i) {
        for (int j = i+1; j < sim_obj.mlnetwork.node_count(); ++j) {
            edge_count += sim_obj.mlnetwork.adj[i][j].size();
        }
    }

    
    NumericMatrix elist(edge_count, 3);
    int iter = 0;
    for (int node_i = 0; node_i < sim_obj.mlnetwork.node_count(); ++node_i) {
        for (int node_j = node_i+1; node_j < sim_obj.mlnetwork.node_count(); ++node_j) {
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
    for (int iter = 0;iter < sim_obj.mlnetwork.layer_count();++iter) {
        for (int iter2 = iter; iter2 < sim_obj.mlnetwork.layer_count(); ++iter2) {
            
            suff_list(i, 0) = iter;
            suff_list(i, 1) = iter2;
            if (iter == iter2) {
                count_avg = (double)accumulate(sim_obj.netCount_samp[iter].begin(), sim_obj.netCount_samp[iter].end(), 0) / (double)sim_obj.num_samples;

            } 
            else {
                count_avg = (double)accumulate(sim_obj.crossLayerCount_samp[iter][iter2].begin(), sim_obj.crossLayerCount_samp[iter][iter2].end(), 0) / (double)sim_obj.num_samples;
            }
            suff_list(i, 2) = count_avg;
            ++i;
            
        }
        
    }
    count_avg = (double)accumulate(sim_obj.threewayCount_samp.begin(), sim_obj.threewayCount_samp.end(), 0) / (double)sim_obj.num_samples;
    suff_list(6, 0) = 3;
    suff_list(6, 1) = 3;
    suff_list(6, 2) = count_avg;


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
                basis_list(iter, 2) = sim_obj.basis_net.adj[node_i][node_j][loc] ;  // in basis network, 1 means a present edge.
                iter += 1;
            }
        }

    }
    
    
    
    // Create return list
    List return_list = List::create(Named("elist") = elist, Named("suff_list") = suff_list, Named("basis_list") = basis_list);

    return return_list;
}

*/


/*
//' @title rcpp_simulate_ml
//' @description
//' Simulate a network-separable multilayer network model with 2-way interactions
//' @name rcpp_simulate_ml
//' @param
//' @examples
//' 
//' 
//' @export
// [[Rcpp::export]]
List rcpp_simulate_ml(IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers,
    NumericVector rtheta, int rand_seed, NumericVector g) {

    
    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int random_seed = rand_seed;
    double gy = g[0];


    vector<double> theta;
    theta.resize(model_dim);
    vector<string> mterms;
    mterms.resize(model_dim);
    for (int p = 0; p < model_dim; ++p) {
        mterms.at(p) = model_terms[p];
        theta.at(p) = rtheta[p];
    }

    sim_ml sim_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, random_seed, gy);
    sim_obj.set_theta(theta);
    sim_obj.simulate_ml();
    sim_obj.compute_sample_mean();
   

    
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
        for (int iter2 = iter; iter2 < sim_obj.mlnetwork.layer_count(); ++iter2) {

            suff_list(i, 0) = iter;
            suff_list(i, 1) = iter2;
            if (iter == iter2) {
                count_avg = (double)accumulate(sim_obj.netCount_samp[iter].begin(), sim_obj.netCount_samp[iter].end(), 0) / (double)sim_obj.num_samples;

            }
            else {
                count_avg = (double)accumulate(sim_obj.crossLayerCount_samp[iter][iter2].begin(), sim_obj.crossLayerCount_samp[iter][iter2].end(), 0) / (double)sim_obj.num_samples;
            }
            suff_list(i, 2) = count_avg;
            ++i;

        }

    }
    count_avg = (double)accumulate(sim_obj.threewayCount_samp.begin(), sim_obj.threewayCount_samp.end(), 0) / (double)sim_obj.num_samples;
    suff_list(6, 0) = 3;
    suff_list(6, 1) = 3;
    suff_list(6, 2) = count_avg;


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
*/


/*
//' @title rcpp_simulate_ml_suffstats
//' @description
//' Simulate a network-separable multilayer network model with 2-way interactions and print out sufficient statistics
//' @name rcpp_simulate_ml_suffstats
//' @param
//' @examples
//' rcpp_simulate_ml_suffstats(1,100,1,2,c("withinlayer","crosslayer"),3,2,c(0.4,0.7))
//' 
//' @export
// [[Rcpp::export]]
List rcpp_simulate_ml_suffstats(IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers,
    NumericVector rtheta, int rand_seed, NumericVector g) {


    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int random_seed = rand_seed;
    double gy = g[0];

    vector<double> theta;
    theta.resize(model_dim);
    vector<string> mterms;
    mterms.resize(model_dim);
    for (int p = 0; p < model_dim; ++p) {
        mterms.at(p) = model_terms[p];
        theta.at(p) = rtheta[p];
    }

    sim_ml sim_obj(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, random_seed, gy);
    sim_obj.set_theta(theta);
    sim_obj.simulate_ml();
    sim_obj.compute_sample_mean();
    
    // Make edge list

    

    
    // Print some statistics
    int m;
    int i = 0;
    m = sim_obj.netCount_samp.size() + sim_obj.nCr(sim_obj.mlnetwork.layer_count(), 2); //+ 1;
    NumericMatrix elist(m, 3);
    double count_avg;
    for (int iter = 0;iter < sim_obj.mlnetwork.layer_count();++iter) {
        for (int iter2 = iter; iter2 < sim_obj.mlnetwork.layer_count(); ++iter2) {

            elist(i, 0) = iter;
            elist(i, 1) = iter2;
            if (iter == iter2) {
                count_avg = (double)accumulate(sim_obj.netCount_samp[iter].begin(), sim_obj.netCount_samp[iter].end(), 0) / (double)sim_obj.num_samples;

            }
            else {
                count_avg = (double)accumulate(sim_obj.crossLayerCount_samp[iter][iter2].begin(), sim_obj.crossLayerCount_samp[iter][iter2].end(), 0) / (double)sim_obj.num_samples;
            }
            elist(i, 2) = count_avg;
            ++i;

        }

    }
    
    


    // Create return list
    List return_list = List::create(Named("elist") = elist);

    return return_list;
}
*/


/*
//' @title rcpp_exact_simulate_ml
//' @description
//' Simulate a network-separable multilayer network model by its exact joint density
//' @name rcpp_exact_simulate_ml
//' @param
//' @examples
//' rcpp_exact_simulate_ml(1,1,3,7,rep("ml_order2",7),10,3,c(0,0.693,0.693,0,0.693,0,0),23333)
//'
//' @export
// [[Rcpp::export]]
List rcpp_exact_simulate_ml(IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval,
    IntegerVector rmodel_dim, StringVector model_terms,
    IntegerVector rnum_nodes, IntegerVector rnum_layers,
    NumericVector rtheta, int rand_seed) {


    int samp_num = rsamp_num[0];
    int burnin = rburnin[0];
    int interval = rinterval[0];
    int model_dim = rmodel_dim[0];
    int num_nodes = rnum_nodes[0];
    int num_layers = rnum_layers[0];
    int random_seed = rand_seed;


    vector<double> theta;
    theta.resize(model_dim);
    vector<string> mterms;
    mterms.resize(model_dim);
    for (int p = 0; p < model_dim; ++p) {
        mterms.at(p) = model_terms[p];
        theta.at(p) = rtheta[p];
    }

    sim_ml_exact sim_obj_exact(samp_num, burnin, interval, model_dim, mterms, num_nodes, num_layers, random_seed);
    sim_obj_exact.set_theta(theta);
    sim_obj_exact.simulate_ml_exact();
    

    // Make edge list

    int edge_count = 0;
    for (int i = 0; i < sim_obj_exact.mlnetwork.node_count(); ++i) {
        for (int j = i + 1; j < sim_obj_exact.mlnetwork.node_count(); ++j) {
            edge_count += sim_obj_exact.mlnetwork.adj[i][j].size();
        }
    }

    

    
    // Print some statistics
    int m;
    int i = 0;
    m = sim_obj_exact.netCount_samp_exact.size() + sim_obj_exact.nCr(sim_obj_exact.mlnetwork.layer_count(), 2) + 1;
    NumericMatrix elist(m, 3);
    double count_avg;
    for (int iter = 0;iter < sim_obj_exact.mlnetwork.layer_count();++iter) {
        for (int iter2 = iter; iter2 < sim_obj_exact.mlnetwork.layer_count(); ++iter2) {

            elist(i, 0) = iter;
            elist(i, 1) = iter2;
            if (iter == iter2) {
                count_avg = (double)accumulate(sim_obj_exact.netCount_samp_exact[iter].begin(), sim_obj_exact.netCount_samp_exact[iter].end(), 0) / (double)sim_obj_exact.num_samples;

            }
            else {
                count_avg = (double)accumulate(sim_obj_exact.crossLayerCount_samp_exact[iter][iter2].begin(), sim_obj_exact.crossLayerCount_samp_exact[iter][iter2].end(), 0) / (double)sim_obj_exact.num_samples;
            }
            elist(i, 2) = count_avg;
            ++i;

        }

    }
    count_avg = (double)accumulate(sim_obj_exact.threewayCount_samp.begin(), sim_obj_exact.threewayCount_samp.end(), 0) / (double)sim_obj_exact.num_samples;
    elist(6, 0) = 3;
    elist(6, 1) = 3;
    elist(6, 2) = count_avg;


    cout << "\nmade it check point 5\n";
    
    // Create return list
    List return_list = List::create(Named("elist") = elist);

    return return_list;
}

*/
