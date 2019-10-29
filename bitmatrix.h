#ifndef rcpp_bitmatrix_module_h
#define rcpp_bitmatrix_module_h
#include <Rcpp.h>
//#include "bitmatrix_RcppExports.h"
using namespace Rcpp;

typedef struct {
  int yOffset;
  int xOffset1;
  int xOffset2;
  int n;
} chord;

typedef struct {
  chord *C;
  int CLength;
  int minYoffset;
  int maxYoffset;
  int minXoffset;
  int maxXoffset;
  int maxN;
} chordSet;

class BitMatrix {
public:
    BitMatrix();
    BitMatrix(int nrows, int ncols);
    BitMatrix(Rcpp::IntegerMatrix mask);
    void fill(bool yesno);
    bool settrue(uint32_t row, uint32_t col);       // 0 relative
    bool _settrue(uint32_t row, uint32_t col);      // 1 relative
    bool setfalse(uint32_t row, uint32_t col);
    bool _setfalse(uint32_t row, uint32_t col);
    bool istrue(uint32_t row, uint32_t col) const;
    bool _istrue(uint32_t row, uint32_t col) const;
    bool isfalse(uint32_t row, uint32_t col) const;
    bool _isfalse(uint32_t row, uint32_t col) const;
    int nrow() const;
    int ncol() const;
    void setnull(bool yesno);
    bool isnull() const;

    std::vector<int> which(NumericVector x, Function f, List args);

    void maskNaN(NumericMatrix x);
    void maskValue(NumericMatrix x, double value);
    void clearValue(IntegerMatrix x, int value);
    void copyTo(IntegerMatrix x);
    void dilate(SEXP kernel);

    std::vector<int> trues(int32_t offset=0) const; // 0 relative
    std::vector<int> _trues() const;                // 1 relative

private:
  void dilate_line(int ***, BitMatrix & destination, chordSet *, int, int);
  void compute_lookup_table_for_line_dilate(int ***T, int yOff, int line, chordSet *set, int nx, int ny);
  bool _isnull;
  uint32_t _nrows;
  uint32_t _ncols;
  uint32_t _npts;
  uint32_t _n32bitwords;
  std::vector<uint32_t> _data;
};

class AdacsHistogram {
public:
  AdacsHistogram();
  void accumulate(Rcpp::NumericVector x,int nbins=16384, double minv=R_NaN, double maxv=R_NaN);
  void accumulateLO(Rcpp::NumericVector x,double offset=0, int nbins=16384, double minv=R_NaN, double maxv=R_NaN);
  void accumulateHI(Rcpp::NumericVector x,double offset=0, int nbins=16384, double minv=R_NaN, double maxv=R_NaN);
  double quantile(double quantile, double offset=0) const;
private:
  int _nbins;
  double _min;
  double _max;
  int _non_null_sample_count;
  int _null_sample_count;
  std::vector<int> _histogram;
  int _toolow;
  int _toohigh;
};
/**
 *  Represents an Akima spline
 *  Given six input (X,Z) pairs this akima spline fits a smooth curve between points 3 and 4.
 *  This akima spline that is valid for interpolations between the third and fourth points (inclusive) of the six (X,Z) pairs
 *  used in its construction
 *  Reference :  A new method of interpolation and smooth curve fitting based on local procedures, Hiroshi Akima,
 *               Journal Of The Association Of Computing Machinery, Volume 17, number 4, October, 1970, pp. 589 - 602.
 *               https://dl.acm.org/citation.cfm?id=321609
 *
 */
class Coeff
{
public:
  Coeff ()
  {
    w0 = 0.;
    w1 = 0.;
    w2 = 0.;
    w3 = 0.;
    x3 = 0.;
    x4 = 0.;
    y3 = 0.;
  }
  static double_t calcSlopeAtMiddle(const double_t *x, const double_t *z);
  double_t interpValue(double_t x) const;
  
  // valid interpolation range
  double_t x3;
  double_t x4;
  // extra coeffs for akima
  double_t w0;
  double_t w1;
  double_t w2;
  double_t w3;
  double_t y3;
};
/**
 *  Represents an array of Akima splines covering a set of (at least 5) input (X,Z) pairs.
 *  Given six input (X,Z) pairs this akima spline fits a smooth curve between points 3 and 4.
 *  The array of akima splines are valid for interpolations between the first and last points (inclusive)
 *  used in its construction.
 *  The evaluation of the Akima splines for the first and last two intervals are adapted to make use of less than 5 (X,Z) pairs.
 *  Reference :  A new method of interpolation and smooth curve fitting based on local procedures, Hiroshi Akima,
 *               Journal Of The Association Of Computing Machinery, Volume 17, number 4, October, 1970, pp. 589 - 602.
 *               https://dl.acm.org/citation.cfm?id=321609
 *
 */
class adacsakima {
public:
  adacsakima();
  adacsakima(int npts, const double_t *xorig, const double_t *yorig);
  bool Initialise(int npts,const double_t *xorig, const double_t *yorig);
  
  bool isValid() const;
  double_t InterpValue(double_t x) const;
private:
  bool isvalid;
  int ncoeffs;
  std::vector<Coeff> coeffs;
};
class Adacs {
public:
    Adacs() {}
    void callwitharg(int nrows,int ncols, BitMatrix & mask);
  void subset_cpp_inplace(
      Rcpp::NumericMatrix image = 0, const int scol=1, const int ecol=1, const int srow=1, const int erow=1, const int coffset=0, const int roffset=0, Rcpp::NumericMatrix oimage = 0);
  void subset_cpp_inplaceI(
      Rcpp::IntegerMatrix image = 0, const int scol=1, const int ecol=1, const int srow=1, const int erow=1, const int coffset=0, const int roffset=0, Rcpp::LogicalMatrix oimage = 0);
#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)
#define ABS(a) (a)<0?(-a):(a)
  Rcpp::NumericVector Cadacs_FindSkyCellValues(Rcpp::NumericMatrix image,
                                               BitMatrix & bobjects, BitMatrix & bmask,
                                               const double loc1, const double loc2,
                                               const double box1, const double box2,
                                               const double boxadd1, const double boxadd2,
                                               const int skypixmin, const int boxiters);
  Rcpp::IntegerVector Cadacs_FindSkyCellValuesBoxC(Rcpp::NumericMatrix image, Rcpp::Nullable<Rcpp::IntegerMatrix> objects, Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
                                                   const double loc1, const double loc2,
                                                   const double box1, const double box2,
                                                   const double boxadd1, const double boxadd2,
                                                   const int skypixmin, const int boxiters);
  Rcpp::NumericVector adacsmagclip(Rcpp::NumericMatrix x, const int sigma, const int clipiters, const double sigmasel, const int estimate);
  Rcpp::NumericVector Cadacs_magclip(Rcpp::NumericVector x, const int sigma, const int clipiters, const double sigmasel, const int estimate);
  
  //==================================================================================
  void interpolateAkimaGrid(Rcpp::NumericVector xseq,Rcpp::NumericVector yseq,Rcpp::NumericMatrix tempmat_sky,Rcpp::NumericMatrix output);
  
  void interpolateLinearGrid(Rcpp::NumericVector xseq,Rcpp::NumericVector yseq,Rcpp::NumericMatrix tempmat_sky,Rcpp::NumericMatrix output);
  
  //==================================================================================
  double_t Cadacs_quantile(Rcpp::NumericVector x, double quantile, int bins=16384, double minv=R_NaN, double maxv=R_NaN);
  double_t Cadacs_quantileLO(Rcpp::NumericVector x, double quantile, const double offset, int bins=16384, double minv=R_NaN, double maxv=R_NaN);
  double_t Cadacs_quantileHI(Rcpp::NumericVector x, double quantile, const double offset, int bins=16384, double minv=R_NaN, double maxv=R_NaN);
  double_t Cadacs_mean(Rcpp::NumericVector x);
  double_t Cadacs_population_variance(Rcpp::NumericVector x, const double offset);
  double_t Cadacs_sample_variance(Rcpp::NumericVector x, const double offset);
  double_t Cadacs_median(Rcpp::NumericVector x);
  double_t Cadacs_mode(Rcpp::NumericVector x);
  void adacsBothFromHistogram(Rcpp::NumericVector x, double quantile,Rcpp::NumericVector results);
  void adacsBothFromHistogramV2(Rcpp::NumericVector x, double quantile,Rcpp::NumericVector results);
  Rcpp::NumericVector Cadacs_SkyEstLoc(Rcpp::NumericMatrix image,
                                       BitMatrix & bobjects, BitMatrix & bmask,
                                       const double loc1, const double loc2,
                                       const double box1, const double box2,
                                       const double boxadd1, const double boxadd2,
                                       const int skypixmin, const int boxiters,
                                       const int doclip, const int skytype, const int skyRMStype, const double sigmasel
  );
  void Cadacs_MakeSkyGrid(Rcpp::NumericMatrix image,
                          BitMatrix & bobjects, BitMatrix & bmask,
                          const int box1, const int box2,
                          const int grid1, const int grid2,
                          const int boxadd1, const int boxadd2,
                          const int type, const int skypixmin, const int boxiters,
                          const int doclip, const int skytype, const int skyRMStype, const double sigmasel,
                          Rcpp::NumericMatrix sky, Rcpp::NumericMatrix skyRMS
  );
};
RCPP_MODULE(adacs){
    using namespace Rcpp ;

    class_<BitMatrix>("BitMatrix")
    // expose the default constructor
    .constructor()
    .constructor<int, int>()
    .constructor<IntegerMatrix>()

    .method("fill", &BitMatrix::fill     , "set or clear all bits")
    .method("settrue", &BitMatrix::_settrue     , "set the bit.  1 relative")
    .method("setfalse", &BitMatrix::_setfalse     , "clear the bit.  1 relative")
    .const_method("istrue", &BitMatrix::_istrue     , "test if bit set.  1 relative")
    .const_method("isfalse", &BitMatrix::_isfalse     , "test if bit clear.  1 relative")
    .method("which", &BitMatrix::which , "C implementation of which")
    .method("maskNaN", &BitMatrix::maskNaN , "Mask any NaN's as true")
    .method("maskValue", &BitMatrix::maskValue, "Mask any cell with the given value")
    .method("clearValue", &BitMatrix::clearValue, "Clears any cell with the given value")
    .method("copyTo", &BitMatrix::copyTo, "exports to IntegerMatrix")
    .const_method("trues", &BitMatrix::_trues, "which are true.  1 relative")
    .const_method("nrow", &BitMatrix::nrow, "return matrix nrows")
    .const_method("ncol", &BitMatrix::ncol, "return matrix ncols")
    .const_method("isnull", &BitMatrix::isnull, "is this object null or not null")
    .method("setnull", &BitMatrix::setnull, "set this object as null or not null")
    .method("dilate", &BitMatrix::dilate, "apply the morphological dilate operation")
    ;
    
    class_<AdacsHistogram>("AdacsHistogram")
    .constructor()
    .method("accumulate", &AdacsHistogram::accumulate, "Acquire the histogram")
    .const_method("quantile", &AdacsHistogram::quantile, "Return the value of the given quantile")
    ;

    class_<Adacs>("Adacs")
    .constructor()
    .method("Cadacs_MakeSkyGrid", &Adacs::Cadacs_MakeSkyGrid, "entry point for adacs_MakeSkyGrid")
    ;
}

#endif

