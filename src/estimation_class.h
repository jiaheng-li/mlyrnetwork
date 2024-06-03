#pragma once
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include <unordered_map>
#include "helper_functions.h" 
#include "simulation_class.h" 
#include "model.h" 
#include "graph_class.h" 

using namespace std;
using namespace arma;
typedef std::vector<double> stdvec;




class est_ml_Hway {

public:
    sim_ml_Hway simulation_ml;
    mlnet obs_net_ml;
    vec theta_est;
    int model_dim;
    int H;
    int iter_max;



public:
    est_ml_Hway(int nsamp, int burn, int intv, int mdim, string mterms,/// Use corsslayer_est as the cross-layer parameter name for estimation
        int N, int K, int highest_order, double random_seed, vector<double> arguments, int itermax)
        : simulation_ml(nsamp, burn, intv, mdim, mterms, N, K, highest_order, random_seed, arguments),
        obs_net_ml(N, K),
        theta_est(mdim)
    {
        model_dim = mdim;
        H = highest_order;
        iter_max = itermax;
        
    }

    int get_model_dim() {
        return simulation_ml.get_model_dim();
    }

    vector<double> get_theta() {
        return conv_to< stdvec >::from(theta_est);
    }

    vector<double> vec_to_stdvec(vec& vec_) {
        return conv_to< stdvec >::from(vec_);
    }

    vec stdvec_to_vec(vector<double>& vec_) {
        return conv_to< vec >::from(vec_);
    }

 


    void compute_initial_estimate() {
        double obs_val;
        vector<double>  change_stat;   /// Change statistics for each edge
        change_stat.resize(get_model_dim());   //need to be changed to the number of dimensions);
        vector<vector<double> > theta;
        theta.resize(2);
        theta[0].resize(get_model_dim());
        theta[1].resize(get_model_dim());
        vector<double> gradient;
        gradient.resize(get_model_dim());

        simulation_ml.compute_param_index();
        simulation_ml.compute_layer_index();

        // simulation_ml.mlnetwork need to be populated. 
        // in rcpp_estimate_model_ml_Hway, only the est_obj.obs_net_ml was populated.
        // this is because we used simulation_ml.compute_change_stats(i, j, k, layers)
        // instead of defining compute_change_stats function for estimation class.
        for (int i = 0; i < obs_net_ml.node_count() - 1; ++i) {
            for (int j = i + 1; j < obs_net_ml.node_count(); ++j) {
                for (int k = 0; k < obs_net_ml.layer_count(); ++k) {
                    if (obs_net_ml.is_edge(i, j, k)) simulation_ml.mlnetwork.add_edge(i, j, k);
                }
                
            }
        }

        // Estimate initial theta using initial guess of zero vector 
        for (int p = 0; p < get_model_dim(); ++p) {
            theta[0][p] = 0;
            theta[1][p] = 0;
        }

        int s = 0;
        bool conv_flag = false;
        double exp_val, scale_val, proposed_step_sum, gamma;
        double inner_prod = 0.0;
        int max_dim = 0;
        int para_dim = 0;
        set<set<int>> layer_ind;
        vector<int> layers;
        string layer_ele;
        //Rcpp::Rcout << "\n  start of the iteration " << "\n";
        for (int iter_num = 0; iter_num < iter_max; ++iter_num) { /// Maximum number of iterations for gradient descent algm.
            // Compute gradient for edge
            for (int p = 0; p < get_model_dim(); ++p) {
                gradient[p] = 0;
            }
            for (int i = 0; i < (obs_net_ml.node_count() - 1); ++i) {
                for (int j = (i + 1); j < obs_net_ml.node_count(); ++j) {
                    for (int k = 0; k < obs_net_ml.layer_count(); ++k) {
                        for (int p = 0; p < get_model_dim(); ++p) {
                            change_stat[p] = 0;
                        }
                        bool flag_for1 = false;
                        for (int l = 0; l < obs_net_ml.layer_count(); ++l) {
                            if (l != k && obs_net_ml.is_edge(i, j, l)) {
                                flag_for1 = true;
                                break;
                            }
                        }
                        if (!flag_for1) continue;

                        inner_prod = 0.0;

                        layer_ind = simulation_ml.layermap[k];
                        

                        for (auto ele : layer_ind) {
                            for (auto ele2 : ele) {
                                layer_ele += to_string(ele2) + "+";
                            }
                            para_dim = simulation_ml.indexmap[layer_ele];
                            layer_ele = "";

                            for (auto ele2 : ele) {
                                layers.push_back(ele2);
                            }
                            change_stat[para_dim] = simulation_ml.compute_change_stats(i, j, k, layers);
                            inner_prod += theta[0][para_dim] * change_stat[para_dim];
                            layers.clear();
                        }

                        exp_val = exp(inner_prod);

                        if (obs_net_ml.is_edge(i, j, k)) { // obs_net.get_edge_type(i, j))) { 
                            obs_val = 1;
                        }
                        else {
                            obs_val = 0;
                        }

                        scale_val = (exp_val / (1 + exp_val)) - obs_val;
                        //Rcpp::Rcout << "\n  scale is " << scale_val << "\n";

                        for (int d = 0; d < get_model_dim(); ++d) {
                            gradient[d] += change_stat[d] * scale_val;
                            //Rcpp::Rcout << "\n  gradient is " << gradient[d] << ", " << "\n";
                        }

                    }

                }
                //Rcpp::Rcout << "\n  end of the iteration of one dyad" << "\n";
            }

            // Compute gamma scaling parameter and step increment
            proposed_step_sum = 0;
            for (int p = 0; p < get_model_dim(); ++p) {
                proposed_step_sum += pow(gradient[p], 2);
            }
            gamma = 1 / (pow(obs_net_ml.node_count(), 2));


            for (int p = 0; p < get_model_dim(); ++p) {
                theta[1][p] = theta[0][p];  // theta[1] keeps values of last step.
                theta[0][p] -= gamma * gradient[p];
            }


            if (sqrt(proposed_step_sum) < 1e-4 * get_model_dim()) {
                Rcpp::Rcout << "\n    Took " << iter_num << " iterations to converge to initial theta = (";
                for (int p = 0; p < get_model_dim(); ++p) {
                    Rcpp::Rcout << theta[0][p];
                    if (p == (get_model_dim() - 1)) {
                        Rcpp::Rcout << ").\n" << flush;
                    }
                    else {
                        Rcpp::Rcout << ", ";
                    }
                }
                conv_flag = true;
                break;
            }
        }
        if (!conv_flag) {
            Rcpp::Rcout << "\n    Initial estimation procedure failed to converge. Starting estimate may be unstable.";
        }
        for (int p = 0; p < get_model_dim(); ++p) {
            theta_est[p] = theta[0][p];
        }
    }
};

