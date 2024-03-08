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
    unordered_map<string, int> indexmap;
    unordered_map<int, vector<string> > layermap;
    vector<int> all_interaction_layer;
    vector<int> combination;
    vector<vector <int> > selected_layer;
   


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


    int compute_change_stats(int i, int j, int k, vector <int>& layers) {
        bool f = true;
        for (int ii = 0; ii < layers.size(); ii++) {
            if (!obs_net_ml.is_edge(i, j, layers[ii]) && layers[ii] != k) {
                f = false;
                break;
            }
        }
        if (f) return 1;
        else return 0;
    }

    int fact(int n)
    {
        if (n == 0) return 1;
        if (n == 1) return 1;
        return n * fact(n - 1);
    }


    int nCr(int n, int r)
    {
        return fact(n) / (fact(r) * fact(n - r));
    }

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
    void compute_param_index() {
        int param_ind = 0;
        for (int h = 0; h < obs_net_ml.layer_count(); ++h) {
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

    // Compute the index of layers associated with each edge
    void compute_layer_index() {
        int param_ind = 0;
        vector<string> layer_vec;
        for (int k = 0; k < obs_net_ml.layer_count(); k++) {
            layer_vec.clear();
            layer_vec.push_back(to_string(k));
            for (int h = 1; h < H; ++h) {
                select_layer(0, h);
                for (auto ele : selected_layer) {
                    if (binary_search(ele.begin(), ele.end(), k)) continue;
                    string layer_ind = "";
                    vector<int> layer_ele = ele;
                    layer_ele.push_back(k);
                    sort(layer_ele.begin(), layer_ele.end());
                    for (auto ele2 : layer_ele) {
                        layer_ind += to_string(ele2);
                    }
                    layer_vec.push_back(layer_ind);
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


            layermap.insert(make_pair(k, layer_vec));
        }
    }

    void compute_initial_estimate() {
        double obs_val;
        
        for (int e = 0; e < obs_net_ml.layer_count(); ++e) {
            all_interaction_layer.push_back(e);
        }

        vector<double>  change_stat;   /// Change statistics for each edge
        change_stat.resize(get_model_dim());   //need to be changed to the number of dimensions);


        vector<vector<double> > theta;
        theta.resize(2);
        theta[0].resize(get_model_dim());
        theta[1].resize(get_model_dim());

        vector<double> gradient;
        gradient.resize(get_model_dim());

        compute_param_index();
        compute_layer_index();

        
        


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

                        vector<string> layer_ind = layermap[k];
                        vector<int> layers;

                        for (auto ele : layer_ind) {
                            para_dim = indexmap[ele];
                            layers.clear();
                            for (auto char_ele : ele) {
                                layers.push_back(int(char_ele - '0')); /// use char - '0' to convert char to integer.
                            }
                            
                  

                            change_stat[para_dim] = compute_change_stats(i, j, k, layers);
                            inner_prod += theta[0][para_dim] * change_stat[para_dim];
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








/*
class est_ml {

public:
    sim_ml simulation_ml;
    mlnet obs_net_ml;
    vec theta_est;
    vec obs_vec;
    vec obs_vec_hold;
    vector<double> weights;
    mat information_matrix;
    int iteration_max;  
    int model_dim;
    bool check_chull;

public:
    est_ml(int nsamp, int burn, int intv, int mdim, vector<string> mterms,  /// Use corsslayer_est as the cross-layer parameter name for estimation
        int N, int K, double random_seed,
        int iter_max, bool check_ch, double g)
        : simulation_ml(nsamp, burn, intv, mdim, mterms, N, K, random_seed, g),
        obs_net_ml(N, K),
        theta_est(mdim),
        obs_vec(mdim),
        obs_vec_hold(mdim),
        information_matrix(mdim, mdim)
    {
        iteration_max = iter_max;
        model_dim = mdim;
        check_chull = check_ch;
        weights.resize(nsamp);
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

    

    void compute_change_stats(int i, int j, vector<vector<double> >& ch_stat, int k, int l) { 
        ch_stat[k][l] = simulation_ml.m.change_stat_funs_ml[0](i, j, obs_net_ml, k, l);  /// Use cs_crosslayer_est for change stats of 2-layer co-ocurrences
    }

    void compute_initial_estimate() {
        double obs_val;

        vector<vector<double> >  change_stats;   /// Change statistics for each edge
        change_stats.resize(obs_net_ml.layer_count());
        for (int i = 0; i < obs_net_ml.layer_count(); ++i) {
            change_stats[i].resize(obs_net_ml.layer_count());
        }

        vector<vector<double> > theta;
        theta.resize(2);
        theta[0].resize(get_model_dim());
        theta[1].resize(get_model_dim());

        vector<double> gradient;
        gradient.resize(get_model_dim());

        // Compute change statistics for each edge
        
        
 
        // Estimate initial theta using initial guess of zero vector 
        for (int p = 0; p < get_model_dim(); ++p) {
            theta[0][p] = 0;
            theta[1][p] = 0;
        }

        vector<vector<double> >  theta_mat;
        theta_mat.resize(obs_net_ml.layer_count());
        for (int i = 0; i < obs_net_ml.layer_count(); ++i) {
            theta_mat[i].resize(obs_net_ml.layer_count());
        }
        

        bool conv_flag = false;
        double exp_val, scale_val, proposed_step_sum, gamma;
        double inner_prod = 0.0;
        vector<double> delta;
        delta.resize(get_model_dim());
        vector<vector<double> > delta_mat;
        delta_mat.resize(obs_net_ml.layer_count());
        for (int i = 0; i < obs_net_ml.layer_count(); ++i) {
            delta_mat[i].resize(obs_net_ml.layer_count());
        }
        for (int iter_num = 0; iter_num < iteration_max; ++iter_num) { /// Maximum number of iterations for gradient descent algm.
            // Compute gradient for edge
            for (int p = 0; p < get_model_dim(); ++p) {
                gradient[p] = 0;
            }
            for (int i = 0; i < (obs_net_ml.node_count() - 1); ++i) {
                for (int j = (i + 1); j < obs_net_ml.node_count(); ++j) {
                    for (int k = 0;k < obs_net_ml.layer_count();++k) {

                        bool flag_for1 = false;
                        for (int l = 0;l < obs_net_ml.layer_count();++l) {
                            if (l != k && obs_net_ml.is_edge(i, j, l)) {
                                flag_for1 = true;
                                break;
                            }
                        }
                        if (!flag_for1) continue;

                        for (int d = 0; d < get_model_dim(); ++d) {
                            delta[d] = 0;
                        }
                        for (int g = 0; g < obs_net_ml.layer_count(); ++g) {
                            for (int q = 0; q < obs_net_ml.layer_count(); ++q) {
                                delta_mat[g][q] = 0;
                            }
                        }
                        inner_prod = 0;
                        for (int l = 0;l < obs_net_ml.layer_count();++l) {
                            compute_change_stats(i, j, change_stats, k, l);  /// This is change statistics for edge x_{i,j}^(k)                            
                            inner_prod += theta_mat[k][l] * change_stats[k][l];
                            delta_mat[k][l] = change_stats[k][l];
                            //Rcpp::Rcout << "\n  theta, cs and ip are " << theta_mat[k][l] << ", " << change_stats[k][l]  << "\n";
                        }
                        exp_val = exp(inner_prod);
                        
                        if (obs_net_ml.is_edge(i, j, k)) { // obs_net.get_edge_type(i, j))) { 
                            obs_val = 1;
                        }
                        else {
                            obs_val = 0;
                        }

                        scale_val = (exp_val / (1 + exp_val)) - obs_val;
                        int pp;
                        for (int q = 0; q < obs_net_ml.layer_count();++q) {
                            if (q < k) {
                                pp = 0;
                                for (int g = 0; g < obs_net_ml.layer_count(); ++g) {
                                    for (int f = g; f < obs_net_ml.layer_count(); ++f) {
                                        if (g == q && f == k) delta[pp] = delta_mat[k][q];
                                        pp++;
                                    }
                                }
                            }
                            if (k <= q && k >= 1) {
                                delta[(2*obs_net_ml.layer_count() - k + 1)*k/2 + q - k] = delta_mat[k][q];
                                
                            }
                            if (k <= q && k == 0) {
                                delta[q] = delta_mat[k][q];
                            }
                            
                        }

                        for (int d = 0; d < get_model_dim(); ++d) {
                            gradient[d] += delta[d] * scale_val;
                            //Rcpp::Rcout << "\n  gradient is " << gradient[d] << ", " << "\n";
                        }
                        
                        
                        
                    }
                    
                }
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

            ///  Assign 1-d vector theta to 2-d vector theta_mat for next updates
            int pp = 0;
            for (int k = 0; k < obs_net_ml.layer_count(); ++k) {
                for (int l = k; l < obs_net_ml.layer_count(); ++l) {
                    theta_mat[k][l] = theta[0][pp];
                    theta_mat[l][k] = theta_mat[k][l];
                    ++pp;
                }
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


*/


