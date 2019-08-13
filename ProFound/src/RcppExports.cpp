// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_median
double get_median(Rcpp::NumericVector clip, Rcpp::Function median);
RcppExport SEXP _ProFound_get_median(SEXP clipSEXP, SEXP medianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type clip(clipSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type median(medianSEXP);
    rcpp_result_gen = Rcpp::wrap(get_median(clip, median));
    return rcpp_result_gen;
END_RCPP
}
// subset_cpp_inplace
void subset_cpp_inplace(Rcpp::NumericMatrix image, const int scol, const int ecol, const int srow, const int erow, const int coffset, const int roffset, Rcpp::NumericMatrix oimage);
RcppExport SEXP _ProFound_subset_cpp_inplace(SEXP imageSEXP, SEXP scolSEXP, SEXP ecolSEXP, SEXP srowSEXP, SEXP erowSEXP, SEXP coffsetSEXP, SEXP roffsetSEXP, SEXP oimageSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type image(imageSEXP);
    Rcpp::traits::input_parameter< const int >::type scol(scolSEXP);
    Rcpp::traits::input_parameter< const int >::type ecol(ecolSEXP);
    Rcpp::traits::input_parameter< const int >::type srow(srowSEXP);
    Rcpp::traits::input_parameter< const int >::type erow(erowSEXP);
    Rcpp::traits::input_parameter< const int >::type coffset(coffsetSEXP);
    Rcpp::traits::input_parameter< const int >::type roffset(roffsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type oimage(oimageSEXP);
    subset_cpp_inplace(image, scol, ecol, srow, erow, coffset, roffset, oimage);
    return R_NilValue;
END_RCPP
}
// water_cpp
Rcpp::IntegerMatrix water_cpp(Rcpp::NumericVector image, const int nx, const int ny, const double abstol, const double reltol, const double cliptol, const int ext, const double skycut, const int pixcut, const bool verbose, const int Ncheck);
RcppExport SEXP _ProFound_water_cpp(SEXP imageSEXP, SEXP nxSEXP, SEXP nySEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP cliptolSEXP, SEXP extSEXP, SEXP skycutSEXP, SEXP pixcutSEXP, SEXP verboseSEXP, SEXP NcheckSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type image(imageSEXP);
    Rcpp::traits::input_parameter< const int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< const int >::type ny(nySEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< const double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< const double >::type cliptol(cliptolSEXP);
    Rcpp::traits::input_parameter< const int >::type ext(extSEXP);
    Rcpp::traits::input_parameter< const double >::type skycut(skycutSEXP);
    Rcpp::traits::input_parameter< const int >::type pixcut(pixcutSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type Ncheck(NcheckSEXP);
    rcpp_result_gen = Rcpp::wrap(water_cpp(image, nx, ny, abstol, reltol, cliptol, ext, skycut, pixcut, verbose, Ncheck));
    return rcpp_result_gen;
END_RCPP
}
// order_cpp
IntegerVector order_cpp(NumericVector x);
RcppExport SEXP _ProFound_order_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// tabulate_cpp
IntegerVector tabulate_cpp(const IntegerVector& x, const int max);
RcppExport SEXP _ProFound_tabulate_cpp(SEXP xSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(tabulate_cpp(x, max));
    return rcpp_result_gen;
END_RCPP
}
// water_cpp_old
IntegerVector water_cpp_old(const NumericVector image, const int nx, const int ny, const double abstol, const double reltol, const double cliptol, const int ext, const double skycut, const int pixcut, const bool verbose, const int Ncheck);
RcppExport SEXP _ProFound_water_cpp_old(SEXP imageSEXP, SEXP nxSEXP, SEXP nySEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP cliptolSEXP, SEXP extSEXP, SEXP skycutSEXP, SEXP pixcutSEXP, SEXP verboseSEXP, SEXP NcheckSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type image(imageSEXP);
    Rcpp::traits::input_parameter< const int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< const int >::type ny(nySEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< const double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< const double >::type cliptol(cliptolSEXP);
    Rcpp::traits::input_parameter< const int >::type ext(extSEXP);
    Rcpp::traits::input_parameter< const double >::type skycut(skycutSEXP);
    Rcpp::traits::input_parameter< const int >::type pixcut(pixcutSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type Ncheck(NcheckSEXP);
    rcpp_result_gen = Rcpp::wrap(water_cpp_old(image, nx, ny, abstol, reltol, cliptol, ext, skycut, pixcut, verbose, Ncheck));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ProFound_get_median", (DL_FUNC) &_ProFound_get_median, 2},
    {"_ProFound_subset_cpp_inplace", (DL_FUNC) &_ProFound_subset_cpp_inplace, 8},
    {"_ProFound_water_cpp", (DL_FUNC) &_ProFound_water_cpp, 11},
    {"_ProFound_order_cpp", (DL_FUNC) &_ProFound_order_cpp, 1},
    {"_ProFound_tabulate_cpp", (DL_FUNC) &_ProFound_tabulate_cpp, 2},
    {"_ProFound_water_cpp_old", (DL_FUNC) &_ProFound_water_cpp_old, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_ProFound(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
