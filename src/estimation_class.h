#pragma once
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include "helper_functions.h" 
#include "simulation_class.h" 
#include "model.h" 
#include "graph_class.h" 

using namespace std;
using namespace arma;
typedef std::vector<double> stdvec;



class est_ml_Hway {
public:
    sim_ml simulation_ml;
    mlnet obs_net_ml;
    vec theta_est;
    vec obs_vec;
    vec obs_vec_hold;
    vector<double> weights;
    mat information_matrix;
    int model_dim;

public:
    est_ml_Hway(int nsamp, int burn, int intv, int mdim, vector<string> mterms, 
        int N, int K, double random_seed, double gy)
        : simulation_ml(nsamp, burn, intv, mdim, mterms, N, K, random_seed, gy),
        obs_net_ml(N, K),
        theta_est(mdim),
        obs_vec(mdim),
        obs_vec_hold(mdim),
        information_matrix(mdim, mdim)
    {
        model_dim = mdim;
        weights.resize(nsamp);
    }


    int get_model_dim() {
        return simulation_ml.get_model_dim();
    }

    void compute_change_stats(int i, int j, vector<vector<double> >&ch_stat, int k, int l) {
        ch_stat[k][l] = simulation_ml.m.change_stat_funs_ml[0](i, j, obs_net_ml, k, l);  
    }

    void compute_initial_estimate() {

    
    
    
    }


};


class est_ml_3way {

public:
    sim_ml simulation_ml;
    mlnet obs_net_ml;
    vec theta_est;
    vec obs_vec;
    vec obs_vec_hold;
    vector<double> weights;
    mat information_matrix;
    int NR_max_iter;  /// For Newton Raphson method
    int MCMLE_max_iter;
    int model_dim;
    double NR_tol;
    bool check_chull;
    int cs3; // change statistic for 3-way interactions.

public:
    est_ml_3way(int nsamp, int burn, int intv, int mdim, vector<string> mterms,  /// Use corsslayer_est as the cross-layer parameter name for estimation
        int N, int K, double random_seed,
        int NR_max, double NRtol, int MCMLE_max, bool check_ch, double gy)
        : simulation_ml(nsamp, burn, intv, mdim, mterms, N, K, random_seed, gy),
        obs_net_ml(N, K),
        theta_est(mdim),
        obs_vec(mdim),
        obs_vec_hold(mdim),
        information_matrix(mdim, mdim)
    {
        NR_max_iter = NR_max;
        NR_tol = NRtol;
        MCMLE_max_iter = MCMLE_max;
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
        for (int iter_num = 0; iter_num < 30000; ++iter_num) { /// Maximum number of iterations for gradient descent algm.
            // Compute gradient for edge
            for (int p = 0; p < get_model_dim(); ++p) {
                gradient[p] = 0;
            }
            bool flag;
            for (int i = 0; i < (obs_net_ml.node_count() - 1); ++i) {
                for (int j = (i + 1); j < obs_net_ml.node_count(); ++j) {

                    /*
                    /// <summary>
                    /// Added for network separable model, 05/31/2022
                    /// skip for dyad where y_ij = 0
                    /// </summary>
                    bool normcheck = true;
                    for (int k = 0; k < obs_net_ml.layer_count(); ++k) {
                        if (obs_net_ml.is_edge(i, j, k)) {
                            normcheck = false;
                            break;
                        }
                        
                    }

                    if (normcheck) continue;
                    
                    /// <summary>
                    /// Modification for network separable model ended
                    /// </summary>
                    */

                    for (int k = 0;k < obs_net_ml.layer_count();++k) {
                        bool flag_for1 = false;
                        for (int l = 0;l < obs_net_ml.layer_count();++l) {
                            if (l != k && obs_net_ml.is_edge(i, j, l)) {
                                flag_for1 = true;
                                break;
                            }
                        }
                        if (!flag_for1) continue;

                        flag = true;
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
                            if (l != k && (!obs_net_ml.is_edge(i, j, l))) flag = false;
                            //Rcpp::Rcout << "\n  theta, cs and ip are " << theta_mat[k][l] << ", " << change_stats[k][l]  << "\n";
                        }
                        if (flag) cs3 = 1;
                        else cs3 = 0;
                        inner_prod += cs3 * theta[0][6];
                        exp_val = exp(inner_prod);

                        if (obs_net_ml.is_edge(i, j, k)) { // obs_net.get_edge_type(i, j))) { 
                            obs_val = 1;
                        }
                        else {
                            obs_val = 0;
                        }


                        /*
                        /// Correct for y_ij = 1
                        ///
                        bool flag_for1 = false;
                        for (int k2 = 0; k2 < obs_net_ml.layer_count(); ++k2) {
                            if (k2 != k && obs_net_ml.is_edge(i, j, k2)) {
                                flag_for1 = true;
                                break;
                            }
                        }

                        if (flag_for1) scale_val = (exp_val / (1 + exp_val)) - obs_val;
                        else scale_val = 0;


                        /// End of changes on 06/02/2022
                        ///
                        */

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
                                delta[(2 * obs_net_ml.layer_count() - k + 1) * k / 2 + q - k] = delta_mat[k][q];

                            }
                            if (k <= q && k == 0) {
                                delta[q] = delta_mat[k][q];
                            }

                        }
                        delta[6] = cs3;
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





