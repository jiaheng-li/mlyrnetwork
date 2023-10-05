#pragma once
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "graph_class.h"
#include "model.h"

using namespace std;
using namespace arma;




class sim_ml_exact {
public:
    vector<double>  theta; // Allow heterogeneity within and across layers. Should specify the highest order of interactions.
    vector<vector<double> > theta_mat;
    vector<double> chain_state;
    vector<double> mean;
    int num_samples;
    int burnin;
    int interval;
    int model_dim;
    double random_seed;
    mat samples;
    mlnet mlnetwork;
    mod m;
    vector<int> count_exact;
    vector<int> count_net_exact;
    vector<vector<int> >  count_cross_exact;
    vector<vector<int> > netCount_samp_exact;
    vector<vector<vector<int> > > crossLayerCount_samp_exact;
    vector<double> ut;
    double Z1; // normalizing constant for pdf
    int count_threeway;
    vector<int> threewayCount_samp;


public:
    sim_ml_exact(int nsamp, int burn, int intv, int mdim, vector<string> mterms,
        int N, int K, double rand_seed)
        : samples(nsamp, mdim),
        mlnetwork(N, K),
        m(mterms)
    {
        num_samples = nsamp;
        burnin = burn;
        interval = intv;
        model_dim = mdim;
        theta.resize(get_model_dim());
        random_seed = rand_seed;
        ut.resize(2^mlnetwork.layer_count());
    }

    void simulate_ml_exact() {

        double u;
        int samp_num = 0;

        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed);
        std::uniform_real_distribution<> runif(0.0, 1.0);

        //Run burnin
        double p1 = ut[0] / Z1;
        double p2 = (ut[0] + ut[1]) / Z1;
        double p3 = (ut[0] + ut[1] + ut[2]) / Z1;
        double p4 = (ut[0] + ut[1] + ut[2] + ut[3]) / Z1;
        double p5 = (ut[0] + ut[1] + ut[2] + ut[3] + ut[4]) / Z1;
        double p6 = (ut[0] + ut[1] + ut[2] + ut[3] + ut[4] + ut[5]) / Z1;
        double p7 = (ut[0] + ut[1] + ut[2] + ut[3] + ut[4] + ut[5] + ut[6]) / Z1;

        
        count_net_exact.resize(mlnetwork.layer_count());
        count_cross_exact.resize(mlnetwork.layer_count());
        for (int i = 0;i < mlnetwork.layer_count();++i) {
            count_cross_exact[i].resize(mlnetwork.layer_count());
        }
        netCount_samp_exact.resize(mlnetwork.layer_count());
        Rcpp::Rcout << netCount_samp_exact.size() << "\n";
        crossLayerCount_samp_exact.resize(mlnetwork.layer_count());
        threewayCount_samp.resize(num_samples);


        Rcpp::Rcout << "sim_ml check point 1.4\n";


        for (int i = 0;i < mlnetwork.layer_count();++i) {
            netCount_samp_exact[i].resize(num_samples);
            Rcpp::Rcout << "one step before assignment\n";
            netCount_samp_exact[i][num_samples - 1] = 100;
            Rcpp::Rcout << netCount_samp_exact[i].size() << "\n";
            Rcpp::Rcout << i << ", " << netCount_samp_exact[i][num_samples-1] << "\n";
        }
        Rcpp::Rcout << "sim_ml check point 1.5\n";
        for (int i = 0;i < mlnetwork.layer_count();++i) {
            Rcpp::Rcout << "sim_ml check point 1.6\n";
            crossLayerCount_samp_exact[i].resize(mlnetwork.layer_count());
            Rcpp::Rcout << crossLayerCount_samp_exact[i].size() << "\n";
            Rcpp::Rcout << "sim_ml check point 1.7\n";
            for (int j = 0;j < mlnetwork.layer_count();++j) {
                Rcpp::Rcout << "sim_ml check point 1.8\n";
                crossLayerCount_samp_exact[i][j].resize(num_samples);
                Rcpp::Rcout << crossLayerCount_samp_exact[i][j].size() << "\n";
                Rcpp::Rcout << "sim_ml check point 1.9\n";
                cout << i << ", " << j << ", " << crossLayerCount_samp_exact[i][j][num_samples-1] << "\n";
            }
        }
        Rcpp::Rcout << "sim_ml check point 2\n";

        int s = 0;
        while (samp_num < num_samples) {
            
            count_net_exact.assign(mlnetwork.layer_count(), 0);
            for (int i = 0;i < mlnetwork.layer_count();++i) {
                count_cross_exact[i].assign(mlnetwork.layer_count(),0);
                
            }
            cout << "sim_ml check point 3\n";


            for (int i = 0;i < mlnetwork.node_count();++i) {
                for (int j = i + 1;j < mlnetwork.node_count();++j) {
                    
                    u = runif(rand);
                    if (u >= p1 && u < p2) {
                        mlnetwork.add_edge(i, j, 2);
                    }
                    if (u >= p2 && u < p3) {
                        mlnetwork.add_edge(i, j, 1);
                    }
                    if (u >= p3 && u < p4 ) {
                        mlnetwork.add_edge(i, j, 1);
                        mlnetwork.add_edge(i, j, 2);
                    }
                    if (u >= p4 && u < p5) {
                        mlnetwork.add_edge(i, j, 0);
                    }
                    if (u >= p5 && u < p6) {
                        mlnetwork.add_edge(i, j, 0);
                        mlnetwork.add_edge(i, j, 2);
                    }
                    if (u >= p6 && u < p7) {
                        mlnetwork.add_edge(i, j, 0);
                        mlnetwork.add_edge(i, j, 1);
                    }
                    if (u >= p7 && u < 1) {
                        mlnetwork.add_edge(i, j, 0);
                        mlnetwork.add_edge(i, j, 1);
                        mlnetwork.add_edge(i, j, 2);
                    }

                    //count the number of edges in each layer                    
                    for (auto ele : mlnetwork.adj[i][j]) {
                        count_net_exact[ele] += 1;
                    }




                    //count the number of co-occurences between 2 layers

                    for (auto ele : mlnetwork.adj[i][j]) {
                        for (auto ele2 : mlnetwork.adj[i][j]) {
                            if (ele >= ele2) continue;
                            count_cross_exact[ele][ele2] += 1;
                        }

                    }

                    if (mlnetwork.adj[i][j].size() == 3) count_threeway += 1;
                }
            }
            cout << "sim_ml check point 4\n";

            for (int i = 0;i < mlnetwork.layer_count();++i) {
                netCount_samp_exact[i][samp_num] = count_net_exact[i];

            }


            for (int i = 0;i < mlnetwork.layer_count() - 1;++i) {
                for (int j = i + 1;j < mlnetwork.layer_count();++j) {
                    crossLayerCount_samp_exact[i][j][samp_num] = count_cross_exact[i][j];
                }
            }

            threewayCount_samp[samp_num] = count_threeway;
            

            samp_num += 1;

            for (int e1 = 0; e1 < mlnetwork.node_count();++e1) {
                for (int e2 = e1+1; e2 < mlnetwork.node_count();++e2) {
                    mlnetwork.adj[e1][e2].clear();
                    mlnetwork.adj[e2][e1].clear();
                }
            }

        }



    }




    int get_model_dim() {
        return model_dim;
    }


    void set_theta(vector<double> vec) {
        theta = vec;
        vector<double> p;
        p.resize(get_model_dim());

        p[0] = exp(theta[0]) / (1 + exp(theta[0]));
        p[3] = exp(theta[3]) / (1 + exp(theta[3]));
        p[5] = exp(theta[5]) / (1 + exp(theta[5]));

        p[1] = exp(theta[1]);
        p[2] = exp(theta[2]);
        p[4] = exp(theta[4]);
        p[6] = exp(theta[6]);

        ut[0] = (1 - p[0]) * (1 - p[3]) * (1 - p[5]);
        ut[1] = (1 - p[0]) * (1 - p[3]) * p[5];
        ut[2] = (1 - p[0]) * p[3] * (1 - p[5]);
        ut[3] = (1 - p[0]) * p[3] * p[5] * p[4];
        ut[4] = p[0] * (1 - p[3]) * (1 - p[5]);
        ut[5] = p[0] * (1 - p[3]) * p[5] * p[2];
        ut[6] = p[0] * p[3] * (1 - p[5]) * p[1];
        ut[7] = p[0] * p[1] * p[2] * p[3] * p[4] * p[5] * p[6];
        Z1 = ut[0] + ut[1] + ut[2] + ut[3] + ut[4] + ut[5] + ut[6] + ut[7];
    }



    // Returns factorial of n
    int fact(int n)
    {
        int res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i;
        return res;
    }


    int nCr(int n, int r)
    {
        return fact(n) / (fact(r) * fact(n - r));
    }


};


class sim_ml_LSM {

public:
    vector<double>  theta; // Allow heterogeneity within and across layers. Should specify the highest order of interactions.
    vector<vector<double> > theta_mat;
    vector<double> chain_state;
    vector<double> mean;
    int num_samples;
    int burnin;
    int interval;
    int model_dim;
    int cs3; // change statistic for 3-way interactions.
    double random_seed;
    mat samples;
    mlnet mlnetwork;
    mlnet basis_net;
    mod m;
    vector<int> count;
    vector<int> count_net;
    int count_threeway;
    vector<vector<vector<int> > > count_cross;
    vector<vector<int> > netCount_samp;
    vector<int> threewayCount_samp;
    vector<vector<vector<int> > > crossLayerCount_samp;
    double gy;  //specify density g(y) of y
    double fixed_effect;
    double norm_mu;
    double norm_sd;


public:
    sim_ml_LSM(int nsamp, int burn, int intv, int mdim, vector<string> mterms,
        int N, int K, double rand_seed, double g, double fe,double nmu, double nsd)
        : samples(nsamp, mdim),
        mlnetwork(N, K),
        basis_net(N, 1),
        m(mterms)
    {
        num_samples = nsamp;
        burnin = burn;
        interval = intv;
        model_dim = mdim;
        chain_state.resize(get_model_dim());
        theta.resize(get_model_dim());
        mean.resize(get_model_dim());
        random_seed = rand_seed;
        theta_mat.resize(mlnetwork.layer_count());
        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            theta_mat[i].resize(mlnetwork.layer_count());
        }
        gy = g;
        fixed_effect = fe;
        norm_mu = nmu;
        norm_sd = nsd;
        


    }

    void simulate_ml() {
        vector<vector<double> > change_stat;
        change_stat.resize(mlnetwork.layer_count());
        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            change_stat[i].resize(mlnetwork.layer_count());
        }



        double prop_prob, u;
        int samp_num = 0;

        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed);
        std::uniform_real_distribution<> runif(0.0, 1.0);

        initial_samp(); /// Sample the network with prob = 0.5 for each edge independently as an initial network
        generate_LSM(fixed_effect, norm_mu, norm_sd);
        //Run burnin
        bool flag;
        for (int i = 0; i < mlnetwork.node_count(); ++i) {
            for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                if (basis_net.is_edge(i, j, 1)) {
                    for (int upd_num = 0; upd_num < burnin; ++upd_num) {
                        for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                            flag = true;
                            for (int l = 0; l < mlnetwork.layer_count(); ++l) {
                                compute_change_stat(i, j, change_stat, k, l);  /// This is the change statistic for x_{i,j}^(k) given the rest layers, i.e., (k,l) and (l,k) are different.
                                if (l != k && (!mlnetwork.is_edge(i, j, l))) flag = false;
                                //cout << "sim_ml check point 1.1\n";
                            }
                            if (flag) cs3 = 1;
                            else cs3 = 0;


                            prop_prob = compute_cond_prob(i, j, k, change_stat);
                            

                            u = runif(rand);
                            if (u <= prop_prob) {
                                if (mlnetwork.is_edge(i, j, k)) {
                                    mlnetwork.delete_edge(i, j, k);
                                    update_chain(-1, change_stat);
                                }
                            }
                            else {
                                if (!mlnetwork.is_edge(i, j, k)) {
                                    mlnetwork.add_edge(i, j, k);
                                    update_chain(1, change_stat);
                                }
                            }


                            // Enforce h(x,y) = 1, use g(y = 1) = 0.333
                            bool h = 0;
                            for (int kk = 0; kk < mlnetwork.layer_count(); ++kk) {
                                if (mlnetwork.is_edge(i, j, kk)) {
                                    h = 1;
                                    break;
                                }
                            }

                            if (h == 0) {
                                mlnetwork.add_edge(i, j, k);
                                update_chain(1, change_stat);
                            }



                        }
                    }
                }
                else {
                    for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                        if (mlnetwork.is_edge(i, j, k)) {
                            mlnetwork.delete_edge(i, j, k);
                            update_chain(-1, change_stat);
                        }
                    }

                }
            }
        }

        //select dyad for statistics tracking
        int loc1 = 0;
        int loc2 = 0;
        double r1;
        double r2;
        int M;
        int m1;
        int m2;

        count_cross.resize(mlnetwork.layer_count());
        netCount_samp.resize(mlnetwork.layer_count());
        crossLayerCount_samp.resize(mlnetwork.layer_count());
        threewayCount_samp.resize(num_samples);


        r1 = runif(rand);
        r2 = runif(rand);
        


        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            netCount_samp[i].resize(num_samples);
        }

        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            crossLayerCount_samp[i].resize(mlnetwork.layer_count());
            for (int j = 0; j < mlnetwork.layer_count(); ++j) {
                crossLayerCount_samp[i][j].resize(num_samples);
            }
        }

        int s = 0;
        while (samp_num < num_samples) {

            count_net.assign(mlnetwork.layer_count(), 0);
            count_threeway = 0;
            for (int i = 0; i < mlnetwork.layer_count(); ++i) {
                count_cross[i].resize(mlnetwork.layer_count());
                for (int j = 0; j < mlnetwork.layer_count(); ++j) {
                    count_cross[i][j].assign(1, 0);
                }
            }

            for (int i = 0; i < mlnetwork.node_count(); ++i) {
                for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                    if (basis_net.is_edge(i, j, 1)) {
                        for (int upd_num = 0; upd_num < interval; ++upd_num) {  
                            for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                                flag = true;
                                for (int l = 0; l < mlnetwork.layer_count(); ++l) {
                                    compute_change_stat(i, j, change_stat, k, l);  /// This is the change statistic for x_{i,j}^(k) given the rest layers, i.e., (k,l) and (l,k) are different.
                                    if (l != k && (!mlnetwork.is_edge(i, j, l))) flag = false;
                                    //cout << "sim_ml check point 1.1\n";
                                }
                                if (flag) cs3 = 1;
                                else cs3 = 0;


                                prop_prob = compute_cond_prob(i, j, k, change_stat);
                                prop_prob = 1 / (1 + exp(prop_prob)); // + cs3 * theta[6]));  // theta[6] would be 0 anyway bi initialization if interactions = 2 and theta[6] is not defined as an input.

                                //Rcpp::Rcout << "second prob: " << prop_prob << "\n";


                                u = runif(rand);
                                if (u <= prop_prob) {
                                    if (mlnetwork.is_edge(i, j, k)) {
                                        mlnetwork.delete_edge(i, j, k);
                                        update_chain(-1, change_stat);
                                    }
                                }
                                else {
                                    if (!mlnetwork.is_edge(i, j, k)) {
                                        mlnetwork.add_edge(i, j, k);
                                        update_chain(1, change_stat);
                                    }
                                }


                                // Enforce h(x,y) = 1, use g(y = 1) = 0.333
                                bool h = 0;
                                for (int kk = 0; kk < mlnetwork.layer_count(); ++kk) {
                                    if (mlnetwork.is_edge(i, j, kk)) {
                                        h = 1;
                                        break;
                                    }
                                }

                                if (h == 0) {
                                    mlnetwork.add_edge(i, j, k);
                                    update_chain(1, change_stat);
                                }



                            }
                        }
                    }
                    else {
                        for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                            if (mlnetwork.is_edge(i, j, k)) {
                                mlnetwork.delete_edge(i, j, k);
                                update_chain(-1, change_stat);
                            }
                        }

                    }
                }
            }

            for (int i = 0; i < mlnetwork.node_count(); ++i) {
                for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                    //count the number of edges in each layer                    
                    for (auto ele : mlnetwork.adj[i][j]) {
                        count_net[ele] += 1;
                    }




                    //count the number of co-occurences between 2 layers

                    for (auto ele : mlnetwork.adj[i][j]) {
                        for (auto ele2 : mlnetwork.adj[i][j]) {
                            if (ele >= ele2) continue;
                            m1 = min(ele, ele2);
                            m2 = max(ele, ele2);
                            count_cross[m1][m2][0] += 1;
                        }

                    }

                    // count the 3-way interactions
                    if (mlnetwork.adj[i][j].size() == 3) count_threeway += 1;




                }
            }

            for (int i = 0; i < mlnetwork.layer_count(); ++i) {
                netCount_samp[i][samp_num] = count_net[i];

            }

            for (int i = 0; i < mlnetwork.layer_count() - 1; ++i) {
                for (int j = i + 1; j < mlnetwork.layer_count(); ++j) {
                    crossLayerCount_samp[i][j][samp_num] = count_cross[i][j][0];
                }
            }

            threewayCount_samp[samp_num] = count_threeway;






            

            record_sample(samp_num);
            samp_num += 1;

        }


        

    }

    void generate_LSM(double fe, double nmu, double nsd) {
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
                if(runif(rand) < p) basis_net.add_edge(i, j, 1);
             
            }
        }

    }

    void record_sample(int ind) {
        for (int p = 0; p < get_model_dim(); ++p) {
            samples.at(ind, p) = chain_state[p];
        }
    }
    ///uniformly initialize the sampler
    void initial_samp() {
        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed + 1); /// make this seed different from the seed for 
        //Rcpp::Rcout << "the random seed in initial sampling is: " << rand << "\n";
        std::uniform_real_distribution<> runif(0.0, 1.0);
        for (int i = 0; i < mlnetwork.node_count(); ++i) {
            for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                    if (runif(rand) < 0.5) mlnetwork.add_edge(i, j, k);
                }
            }
        }
    }

    /// compute conditional probability for x_{i,j}^(k) = 0 given the rest layers of the same dyad
    double compute_cond_prob(int i, int j, int k, vector<vector<double> >& change_stat) {

        double val = 0.0;
        //double prob;
        for (int l = 0; l < mlnetwork.layer_count(); ++l) {
            val += theta_mat[k][l] * change_stat[k][l];
        }

        //prob = 1 / (1 + exp(val));  /// Prob that x_{i,j}^(k) = 0.
        return val;
    }



    int get_model_dim() {
        return model_dim;
    }

    void compute_change_stat(int i, int j, vector<vector<double> >& change_stat, int k, int l) {
        change_stat[k][l] = m.change_stat_funs_ml[0](i, j, mlnetwork, k, l);
    }

    void set_theta(vector<double> vec) {
        theta = vec;
        int p = 0;
        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            for (int j = i; j < mlnetwork.layer_count(); ++j) {
                theta_mat[i][j] = theta[p];
                theta_mat[j][i] = theta_mat[i][j];
                ++p;
            }
        }
    }

    void compute_sample_mean() {
        rowvec col_mean(get_model_dim());
        col_mean = arma::mean(samples, 0);
        for (int p = 0; p < mean.size(); ++p) {
            mean[p] = col_mean.at(p);
        }
    }

    void update_chain(int sign, vector<vector<double>>& change_stat) {
        for (int p = 0; p < get_model_dim(); ++p) {
            chain_state[p] += sign * change_stat[int(p / mlnetwork.layer_count())][p % mlnetwork.layer_count()];
        }
    }


    // Returns factorial of n
    int fact(int n)
    {
        int res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i;
        return res;
    }


    int nCr(int n, int r)
    {
        return fact(n) / (fact(r) * fact(n - r));
    }


};


class sim_ml_SBM {

public:
    vector<double>  theta; // Allow heterogeneity within and across layers. Should specify the highest order of interactions.
    vector<vector<double> > theta_mat;
    vector<double> chain_state;
    vector<double> mean;
    int num_samples;
    int burnin;
    int interval;
    int model_dim;
    int cs3; // change statistic for 3-way interactions.
    double random_seed;
    mat samples;
    mlnet mlnetwork;
    mlnet basis_net;
    mod m;
    vector<int> count;
    vector<int> count_net;
    int count_threeway;
    vector<vector<vector<int> > > count_cross;
    vector<vector<int> > netCount_samp;
    vector<int> threewayCount_samp;
    vector<vector<vector<int> > > crossLayerCount_samp;
    double gy;  //specify density g(y) of y
    int block_num;
    double p_in;
    double p_between;


public:
    sim_ml_SBM(int nsamp, int burn, int intv, int mdim, vector<string> mterms,
        int N, int K, double rand_seed, double g, int block, double p1, double p2)
        : samples(nsamp, mdim),
        mlnetwork(N, K),
        basis_net(N,1),
        m(mterms)
    {
        num_samples = nsamp;
        burnin = burn;
        interval = intv;
        model_dim = mdim;
        chain_state.resize(get_model_dim());
        theta.resize(get_model_dim());
        mean.resize(get_model_dim());
        random_seed = rand_seed;
        theta_mat.resize(mlnetwork.layer_count());
        for (int i = 0;i < mlnetwork.layer_count();++i) {
            theta_mat[i].resize(mlnetwork.layer_count());
        }
        gy = g;
        block_num = block;
        p_in = p1;
        p_between = p2;

        
    }

    void simulate_ml() {
        vector<vector<double> > change_stat;
        change_stat.resize(mlnetwork.layer_count());
        for (int i = 0;i < mlnetwork.layer_count();++i) {
            change_stat[i].resize(mlnetwork.layer_count());
        }

        
        
        double prop_prob, u;
        int samp_num = 0;

        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed);
        std::uniform_real_distribution<> runif(0.0, 1.0);

        initial_samp(); /// Sample the network with prob = 0.5 for each edge independently as an initial network
        generate_SBM(block_num,p_in,p_between);
        //Run burnin
        bool flag;
        for (int i = 0;i < mlnetwork.node_count();++i) {
            for (int j = i+1;j < mlnetwork.node_count();++j) {
                if (basis_net.is_edge(i,j,1)) {
                    for (int upd_num = 0; upd_num < burnin; ++upd_num) {
                        for (int k = 0;k < mlnetwork.layer_count();++k) {
                            flag = true;
                            for (int l = 0; l < mlnetwork.layer_count();++l) {
                                compute_change_stat(i, j, change_stat, k, l);  /// This is the change statistic for x_{i,j}^(k) given the rest layers, i.e., (k,l) and (l,k) are different.
                                if (l != k && (!mlnetwork.is_edge(i, j, l))) flag = false;
                                //cout << "sim_ml check point 1.1\n";
                            }
                            if (flag) cs3 = 1;
                            else cs3 = 0;


                            prop_prob = compute_cond_prob(i, j, k, change_stat);
                            //prop_prob = 1 / (1 + exp(prop_prob + cs3 * theta[6]));  Only use for 3-way interactions



                            u = runif(rand);
                            if (u <= prop_prob) {
                                if (mlnetwork.is_edge(i, j, k)) {
                                    mlnetwork.delete_edge(i, j, k);
                                    update_chain(-1, change_stat);
                                }
                            }
                            else {
                                if (!mlnetwork.is_edge(i, j, k)) {
                                    mlnetwork.add_edge(i, j, k);
                                    update_chain(1, change_stat);
                                }
                            }


                            bool h = 0;
                            for (int kk = 0; kk < mlnetwork.layer_count(); ++kk) {
                                if (mlnetwork.is_edge(i, j, kk)) {
                                    h = 1;
                                    break;
                                }
                            }

                            if (h == 0) {
                                mlnetwork.add_edge(i, j, k);
                                update_chain(1, change_stat);
                            }



                        }
                    }
                }
                else {
                    for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                        if (mlnetwork.is_edge(i, j, k)) {
                            mlnetwork.delete_edge(i, j, k);
                            update_chain(-1, change_stat);
                        }
                    }

                }
            }
        }

        //select dyad for statistics tracking
        int loc1 = 0;
        int loc2 = 0;
        double r1;
        double r2;
        int M;
        int m1;
        int m2;
        count_cross.resize(mlnetwork.layer_count());
        netCount_samp.resize(mlnetwork.layer_count());
        crossLayerCount_samp.resize(mlnetwork.layer_count());
        threewayCount_samp.resize(num_samples);


        r1 = runif(rand);
        r2 = runif(rand);
        
        

        for (int i = 0;i < mlnetwork.layer_count();++i) {
            netCount_samp[i].resize(num_samples);
        }

        for (int i = 0;i < mlnetwork.layer_count();++i) {
            crossLayerCount_samp[i].resize(mlnetwork.layer_count());
            for (int j = 0;j < mlnetwork.layer_count();++j) {
                crossLayerCount_samp[i][j].resize(num_samples);
            }
        }

        int s = 0;
        while (samp_num < num_samples) {

            count_net.assign(mlnetwork.layer_count(),0);
            count_threeway = 0;
            for (int i = 0;i < mlnetwork.layer_count();++i) {
                count_cross[i].resize(mlnetwork.layer_count());
                for (int j = 0;j < mlnetwork.layer_count();++j) {
                    count_cross[i][j].assign(1,0);
                }
            }

            for (int i = 0;i < mlnetwork.node_count();++i) {
                for (int j = i + 1;j < mlnetwork.node_count();++j) {
                    if (basis_net.is_edge(i,j,1)) {
                        for (int upd_num = 0; upd_num < interval; ++upd_num) {  // !!! In the second update loop, it should use interval instead of burnin here !!!
                            for (int k = 0;k < mlnetwork.layer_count();++k) {
                                flag = true;
                                for (int l = 0; l < mlnetwork.layer_count();++l) {
                                    compute_change_stat(i, j, change_stat, k, l);  /// This is the change statistic for x_{i,j}^(k) given the rest layers, i.e., (k,l) and (l,k) are different.
                                    if (l != k && (!mlnetwork.is_edge(i, j, l))) flag = false;
                                    //cout << "sim_ml check point 1.1\n";
                                }
                                if (flag) cs3 = 1;
                                else cs3 = 0;


                                prop_prob = compute_cond_prob(i, j, k, change_stat);
                                prop_prob = 1 / (1 + exp(prop_prob) ); // + cs3 * theta[6]));  


                                u = runif(rand);
                                if (u <= prop_prob) {
                                    if (mlnetwork.is_edge(i, j, k)) {
                                        mlnetwork.delete_edge(i, j, k);
                                        update_chain(-1, change_stat);
                                    }
                                }
                                else {
                                    if (!mlnetwork.is_edge(i, j, k)) {
                                        mlnetwork.add_edge(i, j, k);
                                        update_chain(1, change_stat);
                                    }
                                }


                                bool h = 0;
                                for (int kk = 0; kk < mlnetwork.layer_count(); ++kk) {
                                    if (mlnetwork.is_edge(i, j, kk)) {
                                        h = 1;
                                        break;
                                    }
                                }

                                if (h == 0) {
                                    mlnetwork.add_edge(i, j, k);
                                    update_chain(1, change_stat);
                                }



                            }
                        }
                    }
                    else {
                        for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                            if (mlnetwork.is_edge(i, j, k)) {
                                mlnetwork.delete_edge(i, j, k);
                                update_chain(-1, change_stat);
                            }
                        }

                    }
                }
            }

            for (int i = 0;i < mlnetwork.node_count();++i) {
                for (int j = i + 1;j < mlnetwork.node_count();++j) {
                    //count the number of edges in each layer                    
                    for (auto ele : mlnetwork.adj[i][j]) {
                        count_net[ele] += 1;
                    }
            
                    
                    

                    //count the number of co-occurences between 2 layers
                    
                    for (auto ele : mlnetwork.adj[i][j]) {
                        for (auto ele2 : mlnetwork.adj[i][j]) {
                            if (ele >= ele2) continue;
                            m1 = min(ele, ele2);
                            m2 = max(ele, ele2);                            
                            count_cross[m1][m2][0] += 1;
                        }
                        
                    }

                    // count the 3-way interactions
                    if (mlnetwork.adj[i][j].size() == 3) count_threeway += 1;
                    

                    

                }
            }

            for (int i = 0;i < mlnetwork.layer_count();++i) {
                netCount_samp[i][samp_num] = count_net[i];
                
            }

            for (int i = 0;i < mlnetwork.layer_count() - 1;++i) {
                for (int j = i + 1;j < mlnetwork.layer_count();++j) {
                    crossLayerCount_samp[i][j][samp_num] = count_cross[i][j][0];
                }
            }
            
            threewayCount_samp[samp_num] = count_threeway;


            



            


            record_sample(samp_num);
            samp_num += 1;

        }

        
    }

    void generate_SBM(int B, double p_within, double p_between){
        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed+2); 
        std::uniform_real_distribution<> runif(0.0, 1.0);

        
        vector<int> memb;  
        int N = mlnetwork.node_count();
        memb.resize(N);
        for (int i = 0; i < N; ++i){
            memb[i] = 1;
        }
        for (int b = 0; b < B; ++b){
            for (int i = b * N/B; i < N/B * (b+1); ++i){
                memb[i] = b+1;
            }
        }


        
        for (int i = 0; i < N - 1; ++i){
            for (int j = i+1 ; j < N ; ++j){
                if(memb[i] == memb[j] && runif(rand) < p_within){
                    basis_net.add_edge(i,j,1);

                }

                if(memb[i] != memb[j] && runif(rand) < p_between){
                    basis_net.add_edge(i,j,1);

                }
                
            }
        }
        
    }

    void record_sample(int ind) {
        for (int p = 0; p < get_model_dim(); ++p) {
            samples.at(ind, p) = chain_state[p];
        }
    }
    ///uniformly initialize the sampler
    void initial_samp() {
        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed+1); 
        std::uniform_real_distribution<> runif(0.0, 1.0);
        for (int i = 0;i < mlnetwork.node_count();++i) {
            for (int j = i + 1;j < mlnetwork.node_count();++j) {
                for (int k = 0;k < mlnetwork.layer_count();++k) {
                    if (runif(rand) < 0.5) mlnetwork.add_edge(i, j, k);
                }
            }
        }
    }

    /// compute conditional probability for x_{i,j}^(k) = 0 given the rest layers of the same dyad
    double compute_cond_prob(int i, int j, int k, vector<vector<double> >& change_stat) {
        
        double val = 0.0;
        for (int l = 0;l < mlnetwork.layer_count();++l) {
            val += theta_mat[k][l] * change_stat[k][l];
        }

        return val;
    }

    

    int get_model_dim() {
        return model_dim;
    }

    void compute_change_stat(int i, int j, vector<vector<double> >& change_stat, int k, int l) {
        change_stat[k][l] = m.change_stat_funs_ml[0](i, j, mlnetwork, k,l);
    }

    void set_theta(vector<double> vec) {
        theta = vec;
        int p = 0;
        for (int i = 0;i < mlnetwork.layer_count();++i) {
            for (int j = i;j < mlnetwork.layer_count();++j) {
                theta_mat[i][j] = theta[p];
                theta_mat[j][i] = theta_mat[i][j];
                ++p;
            }
        }
    }

    void compute_sample_mean() {
        rowvec col_mean(get_model_dim());
        col_mean = arma::mean(samples, 0);
        for (int p = 0; p < mean.size(); ++p) {
            mean[p] = col_mean.at(p);
        }
    }

    void update_chain(int sign, vector<vector<double>>& change_stat) {
        for (int p = 0; p < get_model_dim(); ++p) {
            chain_state[p] += sign * change_stat[int(p/mlnetwork.layer_count())][p%mlnetwork.layer_count()];
        }
    }


    // Returns factorial of n
    int fact(int n)
    {
        int res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i;
        return res;
    }


    int nCr(int n, int r)
    {
        return fact(n) / (fact(r) * fact(n - r));
    }


};



class sim_ml {

public:
    vector<double>  theta; // Allow heterogeneity within and across layers. Should specify the highest order of interactions.
    vector<vector<double> > theta_mat;
    vector<double> chain_state;
    vector<double> mean;
    int num_samples;
    int burnin;
    int interval;
    int model_dim;
    int cs3; // change statistic for 3-way interactions.
    double random_seed;
    mat samples;
    mlnet mlnetwork;
    mlnet basis_net;
    mod m;
    vector<int> count;
    vector<int> count_net;
    int count_threeway;
    vector<vector<vector<int> > > count_cross;
    vector<vector<int> > netCount_samp;
    vector<int> threewayCount_samp;
    vector<vector<vector<int> > > crossLayerCount_samp;
    double gy;  //specify density g(y) of y


public:
    sim_ml(int nsamp, int burn, int intv, int mdim, vector<string> mterms,
        int N, int K, double rand_seed, double g)
        : samples(nsamp, mdim),
        mlnetwork(N, K),
        basis_net(N, 1),
        m(mterms)
    {
        num_samples = nsamp;
        burnin = burn;
        interval = intv;
        model_dim = mdim;
        chain_state.resize(get_model_dim());
        theta.resize(get_model_dim());
        mean.resize(get_model_dim());
        random_seed = rand_seed;
        theta_mat.resize(mlnetwork.layer_count());
        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            theta_mat[i].resize(mlnetwork.layer_count());
        }
        gy = g;

    }

    void simulate_ml() {
        vector<vector<double> > change_stat;
        change_stat.resize(mlnetwork.layer_count());
        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            change_stat[i].resize(mlnetwork.layer_count());
        }



        double prop_prob, u;
        int samp_num = 0;

        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed);
        std::uniform_real_distribution<> runif(0.0, 1.0);

        initial_samp(); /// Sample the network with prob = 0.5 for each edge independently as an initial network
        //Run burnin
        bool flag;
        generate_Bernoulli(gy);
        for (int i = 0; i < mlnetwork.node_count(); ++i) {
            for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                if (basis_net.is_edge(i, j, 1)) {
                    for (int upd_num = 0; upd_num < burnin; ++upd_num) {
                        for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                            flag = true;
                            for (int l = 0; l < mlnetwork.layer_count(); ++l) {
                                compute_change_stat(i, j, change_stat, k, l);  /// This is the change statistic for x_{i,j}^(k) given the rest layers, i.e., (k,l) and (l,k) are different.
                                if (l != k && (!mlnetwork.is_edge(i, j, l))) flag = false;
                            }
                            if (flag) cs3 = 1;
                            else cs3 = 0;


                            prop_prob = compute_cond_prob(i, j, k, change_stat);
                            //prop_prob = 1 / (1 + exp(prop_prob + cs3 * theta[6]));  Only use for 3-way interactions



                            u = runif(rand);
                            if (u <= prop_prob) {
                                if (mlnetwork.is_edge(i, j, k)) {
                                    mlnetwork.delete_edge(i, j, k);
                                    update_chain(-1, change_stat);
                                }
                            }
                            else {
                                if (!mlnetwork.is_edge(i, j, k)) {
                                    mlnetwork.add_edge(i, j, k);
                                    update_chain(1, change_stat);
                                }
                            }


                            // Enforce h(x,y) = 1
                            bool h = 0;
                            for (int kk = 0; kk < mlnetwork.layer_count(); ++kk) {
                                if (mlnetwork.is_edge(i, j, kk)) {
                                    h = 1;
                                    break;
                                }
                            }

                            if (h == 0) {
                                mlnetwork.add_edge(i, j, k);
                                update_chain(1, change_stat);
                            }



                        }
                    }
                }
                else {
                    for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                        if (mlnetwork.is_edge(i, j, k)) {
                            mlnetwork.delete_edge(i, j, k);
                            update_chain(-1, change_stat);
                        }
                    }

                }
            }
        }

        //select dyad for statistics tracking
        int loc1 = 0;
        int loc2 = 0;
        double r1;
        double r2;
        int M;
        int m1;
        int m2;
        count_cross.resize(mlnetwork.layer_count());
        netCount_samp.resize(mlnetwork.layer_count());
        crossLayerCount_samp.resize(mlnetwork.layer_count());
        threewayCount_samp.resize(num_samples);


        r1 = runif(rand);
        r2 = runif(rand);
        


        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            netCount_samp[i].resize(num_samples);
        }

        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            crossLayerCount_samp[i].resize(mlnetwork.layer_count());
            for (int j = 0; j < mlnetwork.layer_count(); ++j) {
                crossLayerCount_samp[i][j].resize(num_samples);
            }
        }

        int s = 0;
        while (samp_num < num_samples) {

            count_net.assign(mlnetwork.layer_count(), 0);
            count_threeway = 0;
            for (int i = 0; i < mlnetwork.layer_count(); ++i) {
                count_cross[i].resize(mlnetwork.layer_count());
                for (int j = 0; j < mlnetwork.layer_count(); ++j) {
                    count_cross[i][j].assign(1, 0);
                }
            }

            for (int i = 0; i < mlnetwork.node_count(); ++i) {
                for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                    if (basis_net.is_edge(i, j, 1)) {
                        for (int upd_num = 0; upd_num < interval; ++upd_num) {  
                            for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                                flag = true;
                                for (int l = 0; l < mlnetwork.layer_count(); ++l) {
                                    compute_change_stat(i, j, change_stat, k, l);  /// This is the change statistic for x_{i,j}^(k) given the rest layers, i.e., (k,l) and (l,k) are different.
                                    if (l != k && (!mlnetwork.is_edge(i, j, l))) flag = false;
                                }
                                if (flag) cs3 = 1;
                                else cs3 = 0;


                                prop_prob = compute_cond_prob(i, j, k, change_stat);
                                prop_prob = 1 / (1 + exp(prop_prob + cs3 * theta[6]));



                                u = runif(rand);
                                if (u <= prop_prob) {
                                    if (mlnetwork.is_edge(i, j, k)) {
                                        mlnetwork.delete_edge(i, j, k);
                                        update_chain(-1, change_stat);
                                    }
                                }
                                else {
                                    if (!mlnetwork.is_edge(i, j, k)) {
                                        mlnetwork.add_edge(i, j, k);
                                        update_chain(1, change_stat);
                                    }
                                }


                                // Enforce h(x,y) = 1
                                bool h = 0;
                                for (int kk = 0; kk < mlnetwork.layer_count(); ++kk) {
                                    if (mlnetwork.is_edge(i, j, kk)) {
                                        h = 1;
                                        break;
                                    }
                                }

                                if (h == 0) {
                                    mlnetwork.add_edge(i, j, k);
                                    update_chain(1, change_stat);
                                }



                            }
                        }
                    }
                    else {
                        for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                            if (mlnetwork.is_edge(i, j, k)) {
                                mlnetwork.delete_edge(i, j, k);
                                update_chain(-1, change_stat);
                            }
                        }

                    }
                }
            }

            for (int i = 0; i < mlnetwork.node_count(); ++i) {
                for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                    //count the number of edges in each layer                    
                    for (auto ele : mlnetwork.adj[i][j]) {
                        count_net[ele] += 1;
                    }




                    //count the number of co-occurences between 2 layers

                    for (auto ele : mlnetwork.adj[i][j]) {
                        for (auto ele2 : mlnetwork.adj[i][j]) {
                            if (ele >= ele2) continue;
                            m1 = min(ele, ele2);
                            m2 = max(ele, ele2);
                            count_cross[m1][m2][0] += 1;
                        }

                    }

                    // count the 3-way interactions
                    if (mlnetwork.adj[i][j].size() == 3) count_threeway += 1;




                }
            }
            //cout << "sim_ml check point 4\n";

            for (int i = 0; i < mlnetwork.layer_count(); ++i) {
                netCount_samp[i][samp_num] = count_net[i];

            }

            for (int i = 0; i < mlnetwork.layer_count() - 1; ++i) {
                for (int j = i + 1; j < mlnetwork.layer_count(); ++j) {
                    crossLayerCount_samp[i][j][samp_num] = count_cross[i][j][0];
                }
            }

            threewayCount_samp[samp_num] = count_threeway;






            

            record_sample(samp_num);
            samp_num += 1;

        }


        

    }

    void generate_Bernoulli(double gy) {
        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed + 2); /// make this seed different from the seed for 
        std::uniform_real_distribution<> runif(0.0, 1.0);

        int N = basis_net.node_count();
        for (int i = 0; i < N - 1; ++i) {
            for (int j = i + 1; j < N; ++j) {
                if (runif(rand) < gy) {
                    basis_net.add_edge(i, j, 1);

                }
            }
        }

    }

    void record_sample(int ind) {
        for (int p = 0; p < get_model_dim(); ++p) {
            samples.at(ind, p) = chain_state[p];
        }
    }
    ///uniformly initialize the sampler
    void initial_samp() {
        std::random_device rd;
        std::mt19937 g(rd());

        std::random_device rand_dev;
        std::mt19937 rand(random_seed + 1); /// make this seed different from the seed for 
        std::uniform_real_distribution<> runif(0.0, 1.0);
        for (int i = 0; i < mlnetwork.node_count(); ++i) {
            for (int j = i + 1; j < mlnetwork.node_count(); ++j) {
                for (int k = 0; k < mlnetwork.layer_count(); ++k) {
                    if (runif(rand) < 0.5) mlnetwork.add_edge(i, j, k);
                }
            }
        }
    }

    /// compute conditional probability for x_{i,j}^(k) = 0 given the rest layers of the same dyad
    double compute_cond_prob(int i, int j, int k, vector<vector<double> >& change_stat) {

        double val = 0.0;
        //double prob;
        for (int l = 0; l < mlnetwork.layer_count(); ++l) {
            val += theta_mat[k][l] * change_stat[k][l];
        }

        //prob = 1 / (1 + exp(val));  /// Prob that x_{i,j}^(k) = 0.
        return val;
    }



    int get_model_dim() {
        return model_dim;
    }

    void compute_change_stat(int i, int j, vector<vector<double> >& change_stat, int k, int l) {
        change_stat[k][l] = m.change_stat_funs_ml[0](i, j, mlnetwork, k, l);
    }

    void set_theta(vector<double> vec) {
        theta = vec;
        int p = 0;
        for (int i = 0; i < mlnetwork.layer_count(); ++i) {
            for (int j = i; j < mlnetwork.layer_count(); ++j) {
                theta_mat[i][j] = theta[p];
                theta_mat[j][i] = theta_mat[i][j];
                ++p;
            }
        }
    }

    void compute_sample_mean() {
        rowvec col_mean(get_model_dim());
        col_mean = arma::mean(samples, 0);
        for (int p = 0; p < mean.size(); ++p) {
            mean[p] = col_mean.at(p);
        }
    }

    void update_chain(int sign, vector<vector<double>>& change_stat) {
        for (int p = 0; p < get_model_dim(); ++p) {
            chain_state[p] += sign * change_stat[int(p / mlnetwork.layer_count())][p % mlnetwork.layer_count()];
        }
    }


    // Returns factorial of n
    int fact(int n)
    {
        int res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i;
        return res;
    }


    int nCr(int n, int r)
    {
        return fact(n) / (fact(r) * fact(n - r));
    }


};
