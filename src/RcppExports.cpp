// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_compute_dyad_suffstats
List rcpp_compute_dyad_suffstats(IntegerMatrix RNETWORK, IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval, IntegerVector rmodel_dim, StringVector model_terms, IntegerVector rnum_nodes, IntegerVector rnum_layers, int rhighest_order, int rand_seed, NumericVector arguments);
RcppExport SEXP _mlyrnetwork_rcpp_compute_dyad_suffstats(SEXP RNETWORKSEXP, SEXP rsamp_numSEXP, SEXP rburninSEXP, SEXP rintervalSEXP, SEXP rmodel_dimSEXP, SEXP model_termsSEXP, SEXP rnum_nodesSEXP, SEXP rnum_layersSEXP, SEXP rhighest_orderSEXP, SEXP rand_seedSEXP, SEXP argumentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type RNETWORK(RNETWORKSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rsamp_num(rsamp_numSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rburnin(rburninSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rinterval(rintervalSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rmodel_dim(rmodel_dimSEXP);
    Rcpp::traits::input_parameter< StringVector >::type model_terms(model_termsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rnum_nodes(rnum_nodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rnum_layers(rnum_layersSEXP);
    Rcpp::traits::input_parameter< int >::type rhighest_order(rhighest_orderSEXP);
    Rcpp::traits::input_parameter< int >::type rand_seed(rand_seedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type arguments(argumentsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_compute_dyad_suffstats(RNETWORK, rsamp_num, rburnin, rinterval, rmodel_dim, model_terms, rnum_nodes, rnum_layers, rhighest_order, rand_seed, arguments));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_estimate_model_ml_Hway
List rcpp_estimate_model_ml_Hway(IntegerMatrix RNETWORK, IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval, IntegerVector rmodel_dim, StringVector model_terms, IntegerVector rnum_nodes, IntegerVector rnum_layers, IntegerVector rhighest_order, IntegerVector random_seeds, NumericVector arguments, IntegerVector itermax);
RcppExport SEXP _mlyrnetwork_rcpp_estimate_model_ml_Hway(SEXP RNETWORKSEXP, SEXP rsamp_numSEXP, SEXP rburninSEXP, SEXP rintervalSEXP, SEXP rmodel_dimSEXP, SEXP model_termsSEXP, SEXP rnum_nodesSEXP, SEXP rnum_layersSEXP, SEXP rhighest_orderSEXP, SEXP random_seedsSEXP, SEXP argumentsSEXP, SEXP itermaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type RNETWORK(RNETWORKSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rsamp_num(rsamp_numSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rburnin(rburninSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rinterval(rintervalSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rmodel_dim(rmodel_dimSEXP);
    Rcpp::traits::input_parameter< StringVector >::type model_terms(model_termsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rnum_nodes(rnum_nodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rnum_layers(rnum_layersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rhighest_order(rhighest_orderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type random_seeds(random_seedsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type arguments(argumentsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type itermax(itermaxSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_estimate_model_ml_Hway(RNETWORK, rsamp_num, rburnin, rinterval, rmodel_dim, model_terms, rnum_nodes, rnum_layers, rhighest_order, random_seeds, arguments, itermax));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_simulate_ml_Hway
List rcpp_simulate_ml_Hway(IntegerVector rsamp_num, IntegerVector rburnin, IntegerVector rinterval, IntegerVector rmodel_dim, StringVector model_terms, IntegerVector rnum_nodes, IntegerVector rnum_layers, NumericVector rtheta, int rhighest_order, int rand_seed, NumericVector arguments);
RcppExport SEXP _mlyrnetwork_rcpp_simulate_ml_Hway(SEXP rsamp_numSEXP, SEXP rburninSEXP, SEXP rintervalSEXP, SEXP rmodel_dimSEXP, SEXP model_termsSEXP, SEXP rnum_nodesSEXP, SEXP rnum_layersSEXP, SEXP rthetaSEXP, SEXP rhighest_orderSEXP, SEXP rand_seedSEXP, SEXP argumentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type rsamp_num(rsamp_numSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rburnin(rburninSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rinterval(rintervalSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rmodel_dim(rmodel_dimSEXP);
    Rcpp::traits::input_parameter< StringVector >::type model_terms(model_termsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rnum_nodes(rnum_nodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rnum_layers(rnum_layersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rtheta(rthetaSEXP);
    Rcpp::traits::input_parameter< int >::type rhighest_order(rhighest_orderSEXP);
    Rcpp::traits::input_parameter< int >::type rand_seed(rand_seedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type arguments(argumentsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_simulate_ml_Hway(rsamp_num, rburnin, rinterval, rmodel_dim, model_terms, rnum_nodes, rnum_layers, rtheta, rhighest_order, rand_seed, arguments));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _mlyrnetwork_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _mlyrnetwork_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _mlyrnetwork_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _mlyrnetwork_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mlyrnetwork_rcpp_compute_dyad_suffstats", (DL_FUNC) &_mlyrnetwork_rcpp_compute_dyad_suffstats, 11},
    {"_mlyrnetwork_rcpp_estimate_model_ml_Hway", (DL_FUNC) &_mlyrnetwork_rcpp_estimate_model_ml_Hway, 12},
    {"_mlyrnetwork_rcpp_simulate_ml_Hway", (DL_FUNC) &_mlyrnetwork_rcpp_simulate_ml_Hway, 11},
    {"_mlyrnetwork_rcpparma_hello_world", (DL_FUNC) &_mlyrnetwork_rcpparma_hello_world, 0},
    {"_mlyrnetwork_rcpparma_outerproduct", (DL_FUNC) &_mlyrnetwork_rcpparma_outerproduct, 1},
    {"_mlyrnetwork_rcpparma_innerproduct", (DL_FUNC) &_mlyrnetwork_rcpparma_innerproduct, 1},
    {"_mlyrnetwork_rcpparma_bothproducts", (DL_FUNC) &_mlyrnetwork_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mlyrnetwork(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
