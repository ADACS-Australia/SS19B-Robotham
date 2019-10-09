//
// Optimisations contributed by ADACS
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#define adacs_BOTH 1
#define adacs_LO 2
#define adacs_HI 3
#define adacs_SD 4
#define adacs_RBOTH 5
#define adacs_RLO 6
#define adacs_RHI 7
#define adacs_RSD 8

#define adacs_AUTO 1
#define adacs_SET 2

#define adacs_MEDIAN 1
#define adacs_MEAN 2
#define adacs_MODE 3
#define adacs_RMEDIAN 4
#define adacs_RMEAN 5
#define adacs_RMODE 6

#define adacs_CLASSIC_BILINEAR 1
#define adacs_AKIMA_BICUBIC 2
// An example of how to call an R function from C/C++
// [[Rcpp::export]]
double get_median(
    Rcpp::NumericVector clip, Rcpp::Function median
) {
  Rcpp::NumericVector value=median(clip);
  return value[0];
}
// An image subsetting method that uses less memory
// [[Rcpp::export]]
void subset_cpp_inplace(
    Rcpp::NumericMatrix image = 0, const int scol=1, const int ecol=1, const int srow=1, const int erow=1, const int coffset=0, const int roffset=0, Rcpp::NumericMatrix oimage = 0)
{
  // 0 relative
  int sscol = scol-1;
  int eecol = ecol-1;
  int ssrow = srow-1;
  int eerow = erow-1;
  
  std::fill( oimage.begin(), oimage.end(), NA_REAL ) ;
  for (int i = ssrow; i <= eerow; i++) {
    for (int j = sscol; j <= eecol; j++) {
      //Rcpp::Rcout << "   in " << i-roffset-ssrow << " " << j-roffset-sscol << "\n";
      oimage(i-roffset,j-coffset) = image(i, j);
    }
  }
}
// [[Rcpp::export]]
void subset_cpp_inplaceI(
    Rcpp::IntegerMatrix image = 0, const int scol=1, const int ecol=1, const int srow=1, const int erow=1, const int coffset=0, const int roffset=0, Rcpp::LogicalMatrix oimage = 0)
{
  // 0 relative
  int sscol = scol-1;
  int eecol = ecol-1;
  int ssrow = srow-1;
  int eerow = erow-1;
  
  std::fill( oimage.begin(), oimage.end(), NA_LOGICAL ) ;
  for (int i = ssrow; i <= eerow; i++) {
    for (int j = sscol; j <= eecol; j++) {
      //Rcpp::Rcout << "   in " << i-roffset-ssrow << " " << j-roffset-sscol << "\n";
      oimage(i-roffset,j-coffset) = image(i, j)==0;
    }
  }
}
#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)
#define ABS(a) (a)<0?(-a):(a)
#define flt64Null -999
// [[Rcpp::export]]
Rcpp::NumericVector Cadacs_FindSkyCellValues(Rcpp::NumericMatrix image, Rcpp::Nullable<Rcpp::IntegerMatrix> objects, Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
    const double loc1, const double loc2,
    const double box1, const double box2,
    const double boxadd1, const double boxadd2,
    const int skypixmin, const int boxiters)
{
  // R is 1 relative
  int iloc1 = (int)(loc1+0.5);
  int iloc2 = (int)(loc2+0.5);
  int ibox1 = (int)(box1/2);
  int ibox2 = (int)(box2/2);
  int nrow = image.nrow();
  int ncol = image.ncol();
  
  //Rcpp::Rcout << "\nCbox "<<ssrow<<" "<<eerow<<" "<<sscol<<" "<<eecol<<"\n";
  
  const double_t* iiimage=REAL(image);
  Rcpp::IntegerMatrix iobjects;
  const int32_t* iiobjects=NULL;
  if (objects.isNotNull()) {
    iobjects = Rcpp::as<Rcpp::IntegerMatrix>(objects);
    iiobjects=INTEGER(objects.get());
  }
  Rcpp::IntegerMatrix imask;
  const int32_t* iimask=NULL;
  if (mask.isNotNull()) {
    imask = Rcpp::as<Rcpp::IntegerMatrix>(mask);
    iimask=INTEGER(mask.get());
  }

  int iboxadd1=0;
  int iboxadd2=0;
  int skyN=0;
  int iterN=0;
  int ssrow = 1;
  int eerow = 0;
  int sscol = 1;
  int eecol = 0;
  
  while(skyN<skypixmin & iterN<=boxiters){
    skyN = 0;
    ibox1 += iboxadd1;
    ibox2 += iboxadd2;
    ssrow = MAX(1,iloc1-ibox1);
    eerow = MIN(nrow,iloc1+ibox1);
    sscol = MAX(1,iloc2-ibox2);
    eecol = MIN(ncol,iloc2+ibox2);
    
    for (int j = sscol; j <= eecol; j++) {
      int ii=(j-1)*ncol+(ssrow-1);
      for (int i = ssrow; i <= eerow; i++,ii++) {
        // Count sky cells (sky cells are those NOT masked out and NOT objects)
        if ((iiobjects!=NULL)) {
          if (iiobjects[ii]==0 && (iimask==NULL || iimask[ii]==0)) {
            skyN++;
          }
        } else if (iimask!=NULL) {
          if (iimask[ii]==0) {
            skyN++;
          }
        } else {
          skyN++;
        }
      }
    }
    iterN++;
    iboxadd1 = (int)(boxadd1/2);
    iboxadd2 = (int)(boxadd2/2);
    
  }
  // copy sky cell values to vec and return
  Rcpp::NumericVector vec(skyN);
  int k=0;
  for (int j = sscol; j <= eecol; j++) {
    int ii=(j-1)*ncol+(ssrow-1);
    for (int i = ssrow; i <= eerow; i++,ii++) {
      if ((iiobjects!=NULL)) {
        if (iiobjects[ii]==0 && (iimask==NULL || iimask[ii]==0)) {
          vec[k++] = iiimage[ii];
        }
      } else if (iimask!=NULL) {
        if (iimask[ii]==0) {
          vec[k++] = iiimage[ii];
        }
      } else {
        vec[k++] = iiimage[ii];
      }
    }
  }
  return vec;
}
// [[Rcpp::export]]
Rcpp::IntegerVector Cadacs_FindSkyCellValuesBoxC(Rcpp::NumericMatrix image, Rcpp::Nullable<Rcpp::IntegerMatrix> objects, Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
                                            const double loc1, const double loc2,
                                            const double box1, const double box2,
                                            const double boxadd1, const double boxadd2,
                                            const int skypixmin, const int boxiters)
{
  // This method is only needed for plotting (and therefore only for testing)
  // R is 1 relative
  int iloc1 = (int)(loc1+0.5);
  int iloc2 = (int)(loc2+0.5);
  int ibox1 = (int)(box1/2);
  int ibox2 = (int)(box2/2);
  int nrow = image.nrow();
  int ncol = image.ncol();
  
  Rcpp::IntegerMatrix iobjects;
  const int32_t* iiobjects=NULL;
  if (objects.isNotNull()) {
    iobjects = Rcpp::as<Rcpp::IntegerMatrix>(objects);
    iiobjects=INTEGER(objects.get());
  }
  Rcpp::IntegerMatrix imask;
  const int32_t* iimask=NULL;
  if (mask.isNotNull()) {
    imask = Rcpp::as<Rcpp::IntegerMatrix>(mask);
    iimask=INTEGER(objects.get());
  }
  
  int iboxadd1=0;
  int iboxadd2=0;
  int skyN=0;
  int iterN=0;
  int ssrow = 1;
  int eerow = 0;
  int sscol = 1;
  int eecol = 0;
  
  while(skyN<skypixmin & iterN<=boxiters){
    skyN = 0;
    ibox1 += iboxadd1;
    ibox2 += iboxadd2;
    ssrow = MAX(1,iloc1-ibox1);
    eerow = MIN(nrow,iloc1+ibox1);
    sscol = MAX(1,iloc2-ibox2);
    eecol = MIN(ncol,iloc2+ibox2);
    
    for (int j = sscol; j <= eecol; j++) {
      int ii=(j-1)*ncol+(ssrow-1);
      for (int i = ssrow; i <= eerow; i++,ii++) {
        // Count sky cells (sky cells are those NOT masked out and NOT objects)
        if ((iiobjects!=NULL && iiobjects[ii]==0) && (iimask==NULL || iimask[ii]==0)) {
          skyN++;
        }
      }
    }
    iterN++;
    iboxadd1 = (int)(boxadd1/2);
    iboxadd2 = (int)(boxadd2/2);
    
  }
  // copy sky cell values to vec and return
  Rcpp::IntegerVector vec(2);
  vec[0] = ibox1;
  vec[1] = ibox2;
  return vec;
}
// [[Rcpp::export]]
Rcpp::NumericVector adacsmagclip(Rcpp::NumericMatrix x, const int sigma, const int clipiters, const double sigmasel, const int estimate){
  const double_t* iiix=REAL(x);
  int nrow = x.nrow();
  int ncol = x.ncol();
  std::vector<double_t> myx (iiix, iiix+nrow*ncol);
  int length=0;
  for (int i=0;i<nrow*ncol;i++)
  {
    if (!std::isnan(myx[i])) {
      myx[length++] = myx[i];
    }
  }
  std::sort (myx.begin(), myx.begin()+length, std::less<double_t>()); // ascending
  
  int newlen = length;
  if(clipiters>0 & length>0){
    double sigcut=R::pnorm(sigmasel, 0.0, 1.0, 1, 0);
    //double sigcut=0.8;
    
    for(int iteration=0; iteration<clipiters; iteration++){
      if(newlen<=1)
        break;
      int oldlen=newlen;
      double_t roughmed=myx[newlen/2-1];
      double_t clipsigma=sigma;
      if (sigma==1) {
        double l1=MAX(newlen,2);
        double_t y=1.0-2.0/l1;
        clipsigma = R::qnorm(y, 0.0, 1.0, 1, 0);
      }
      
      double_t vallims = 0;
      switch(estimate) {
      case 1:
        vallims = clipsigma*(myx[sigcut*newlen-1]-myx[(1-sigcut)*newlen-1])/2/sigmasel;
        break;
      case 2:
        vallims = clipsigma*(roughmed-myx[(1-sigcut)*newlen-1])/sigmasel;
        break;
      case 3:
        vallims = clipsigma*(myx[sigcut*newlen-1]-roughmed)/sigmasel;
        break;
      }
      newlen = 0;
      for (int i=0;i<oldlen;i++)
      {
        if(myx[i]>=(roughmed-vallims) && myx[i]<=(roughmed+vallims))
        {
          myx[newlen++] = myx[i];
        }
      }
      if(oldlen==newlen)
        break;
    }
  }
  // copy sky cell values to vec and return
  Rcpp::NumericVector vec(newlen);
  for (int i=0;i<newlen;i++)
  {
      vec[i] = myx[i];
  }
  return vec;
}
// [[Rcpp::export]]
Rcpp::NumericVector Cadacs_magclip(Rcpp::NumericVector x, const int sigma, const int clipiters, const double sigmasel, const int estimate){
  const double_t* iiix=REAL(x);
  int nb = x.length();
  std::vector<double_t> myx (iiix, iiix+nb);
  int length=0;
  for (int i=0;i<nb;i++)
  {
    if (!std::isnan(myx[i])) {
      myx[length++] = myx[i];
    }
  }
  std::sort (myx.begin(), myx.begin()+length, std::less<double_t>()); // ascending
  
  int newlen = length;
  if(clipiters>0 & length>0){
    double sigcut=R::pnorm(sigmasel, 0.0, 1.0, 1, 0);
    //double sigcut=0.8;
    
    for(int iteration=0; iteration<clipiters; iteration++){
      if(newlen<=1)
        break;
      int oldlen=newlen;
      double_t roughmed=myx[newlen/2-1];
      double_t clipsigma=sigma;
      if (sigma==1) {
        double l1=MAX(newlen,2);
        double_t y=1.0-2.0/l1;
        clipsigma = R::qnorm(y, 0.0, 1.0, 1, 0);
      }
      
      double_t vallims = 0;
      switch(estimate) {
      case 1:
        vallims = clipsigma*(myx[sigcut*newlen-1]-myx[(1-sigcut)*newlen-1])/2/sigmasel;
        break;
      case 2:
        vallims = clipsigma*(roughmed-myx[(1-sigcut)*newlen-1])/sigmasel;
        break;
      case 3:
        vallims = clipsigma*(myx[sigcut*newlen-1]-roughmed)/sigmasel;
        break;
      }
      newlen = 0;
      for (int i=0;i<oldlen;i++)
      {
        if(myx[i]>=(roughmed-vallims) && myx[i]<=(roughmed+vallims))
        {
          myx[newlen++] = myx[i];
        }
      }
      if(oldlen==newlen)
        break;
    }
  }
  // copy sky cell values to vec and return
  Rcpp::NumericVector vec(newlen);
  for (int i=0;i<newlen;i++)
  {
    vec[i] = myx[i];
  }
  return vec;
}

// Akima spline used to interpolate regular 2D grid
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
/*
 *   Finds the slope of the middle point in a sequence of five (x,z) pairs.
 */
double_t Coeff::calcSlopeAtMiddle(const double_t *x, const double_t *z)
{
  int i;
  double_t S13, S34, S12, S24, W2, W3, Z, DEN;
  double_t A[4], B[4];
  
  // Calculate differences between input points
  for (i = 0; i < 4; i++)
  {
    A[i] = x[i + 1] - x[i];
    B[i] = z[i + 1] - z[i];
  }
  
  //
  S13 = A[0] * B[2] - A[2] * B[0];
  S34 = A[2] * B[3] - A[3] * B[2];
  S12 = A[0] * B[1] - A[1] * B[0];
  S24 = A[1] * B[3] - A[3] * B[1];
  W2  = sqrt(ABS(S13 * S34));
  W3  = sqrt(ABS(S12 * S24));
  //
  Z   = W2 * B[1] + W3 * B[2];
  DEN = W2 * A[1] + W3 * A[2];
  
  if (DEN == 0.0)
  {
    return 0.0;
  }
  return Z / DEN;
}
/*
 * Evaluate this akima spline interpolating polynomial at x
 */
double_t Coeff::interpValue(double_t x) const{
  double_t xx = x - x3;
  return y3 + xx * (w0 + xx * (w2 + xx * w3));
}

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
adacsakima::adacsakima() {
  ncoeffs = 0;
  isvalid = false;
}
bool adacsakima::isValid() const {
  return isvalid;
}

/**---------------------------------------------------------------
*   PURPOSE:
*       To determine the array of akima spline that fit a set of points (Xi,Yi) i = 1,n.  The table of coefficients stored
* in coeffs can be used with InterpValue to perform spline interpolations.
*
*/
adacsakima::adacsakima(int npts, const double_t *xorig, const double_t *yorig)
{
  Initialise(npts,xorig,yorig);
}
/**---------------------------------------------------------------
 *   PURPOSE:
 *       To determine the akima spline that fits a set of points (Xi,Yi) i = 1,n.  The table of coefficients stored
 * in wkarea can be used withInterp to perform spline interpolations.
 *
 *  s3 slope at start of interval s4 slope at end
 */
/**
 * Computes the array of akima splines - notes the special treatment of intervals within 2 of the start and end.
 */
bool adacsakima::Initialise(int npts, const double_t *xorig, const double_t *yorig)
{
  double_t r2 = 2, r3 = 3;
  double_t s3 = 0, s4 = 0, dx = 0., dy = 0., p2, p3, x3, x4, y3, y4;
  int i;
  
  if (npts < 5)	// Error status if validity test fails
  {
    return false;
  }
  
  // One spline for each interval (xorig[i], xorig[i+1])
  ncoeffs = npts-1;
  coeffs.reserve(ncoeffs);
  
  for (i = 0; i < ncoeffs; i++)
  {
    x3 = xorig[i];
    y3 = yorig[i];
    
    if (i != npts - 1)
    {
      x4 = xorig[i + 1];
      y4 = yorig[i + 1];
      dx = x4 - x3;
      dy = y4 - y3;
    }
    
    // check for boundary conditions
    if (i == 0)
    {
      // do first interval
      s3 = dy / dx;
      s4 = Coeff::calcSlopeAtMiddle(&(xorig[0]), &(yorig[0]));
      s4 = (s4 + s3) / 2;
    }
    else if (i == 1)
    {
      // do second interval
      s3 = dy / dx;
      s4 = Coeff::calcSlopeAtMiddle(&(xorig[0]), &(yorig[0]));
      s3 = (s4 + s3) / 2;
    }
    else if (i == ncoeffs - 2)
    {
      // to second last interval
      s3 = Coeff::calcSlopeAtMiddle(&xorig[i - 2], &yorig[i - 2]);
      s4 = dy / dx;
      s4 = (s4 + s3) / 2;
    }
    else if (i == ncoeffs - 1)
    {
      // do last interval
      s3 = Coeff::calcSlopeAtMiddle(&xorig[ncoeffs - 4], &yorig[ncoeffs - 4]);
      x3 = xorig[npts - 2];
      y3 = yorig[npts - 2];
      x4 = xorig[npts - 1];
      y4 = yorig[npts - 1];
      dx = x4 - x3;
      dy = y4 - y3;
      s4 = dy / dx;
      s3 = (s4 + s3) / 2;
    }
    else
    {
      // do the "pure" akima intervals
      // Determine slope at beginning and end of current interval.
      s3 = Coeff::calcSlopeAtMiddle(&xorig[i - 2], &yorig[i - 2]);
      s4 = Coeff::calcSlopeAtMiddle(&xorig[i - 1], &yorig[i - 1]);
    }
    
    //
    // Compute coefficients of cubic equation for this akima
    //
    p2 = (r3 * dy / dx - r2 * s3 - s4) / dx;
    p3 = (s3 + s4 - r2 * dy / dx) / (dx * dx);
    //
    // store slopes and interploting coefficients
    Coeff coeff;
    coeff.w0 = s3;
    coeff.w1 = s4;
    coeff.w2 = p2;
    coeff.w3 = p3;
    coeff.x3 = x3;
    coeff.x4 = x4;
    coeff.y3 = y3;
    coeffs.push_back(coeff);
  }
  return true;
}
double_t adacsakima::InterpValue(double_t x) const
{
  for (int spline = 0; spline < ncoeffs; spline++)
  {
    // Determine if this spline's interval covers x
    if (x >= coeffs[spline].x3 && x <= coeffs[spline].x4)
    {
      // Evaluate this akima spline interpolating polynomial at x
      return coeffs[spline].interpValue(x);
    }
  }
  return 0.0;
}

//==================================================================================
// [[Rcpp::export]]
void interpolateAkimaGrid(Rcpp::NumericVector xseq,Rcpp::NumericVector yseq,Rcpp::NumericMatrix tempmat_sky,Rcpp::NumericMatrix output) {
/*
 * An Matrix element is at (row,col)
 * The elements of a row stack vertically
 * Any row I is to the right of row I-1
 */
  int myxnpts = output.nrow();
  int myynpts = output.ncol();
  const double_t* myx=REAL(xseq);
  const double_t* myy=REAL(yseq);
  int ncol=tempmat_sky.ncol();
  int nrow=tempmat_sky.nrow();
  
  std::vector<double_t> xin,zin;
  std::vector<adacsakima> akimaCOL;
  xin.reserve(ncol);
  zin.reserve(nrow);
  akimaCOL.reserve(ncol);
  for (int i = 1; i <= nrow; i++)
  {
    xin.push_back(0);
    zin.push_back(0);
  }
  for (int j = 1; j <= ncol; j++) {
    int ii=(j-1)*ncol+(1-1);
    for (int i = 1; i <= nrow; i++,ii++) {
      xin[i-1] = myx[i-1];
      zin[i-1] = tempmat_sky(i-1,j-1);
    }
    adacsakima thisspline;
    thisspline.Initialise(nrow,xin.data(),zin.data());
    akimaCOL.push_back(thisspline);
  }
  
  // For each vertical row 
  for (int i = 1; i <= myxnpts; i++) {
    // For a spline to interpolate vertically along the elements of the row
    double_t x = -0.5+i;
    for (int j = 1; j <= ncol; j++) {
      xin[j-1] = myy[j-1];
      zin[j-1] = akimaCOL[j-1].InterpValue(x);
    }
    adacsakima thisspline;
    thisspline.Initialise(ncol,xin.data(),zin.data());
    
    // Interpolate vertically for each element (j) in the current (i) output row
    for (int j = 1; j <= myynpts; j++) {
      double_t y = -0.5+j;
      output(i-1,j-1) = thisspline.InterpValue(y);
    }
  }
}

// [[Rcpp::export]]
void interpolateLinearGrid(Rcpp::NumericVector xseq,Rcpp::NumericVector yseq,Rcpp::NumericMatrix tempmat_sky,Rcpp::NumericMatrix output) {
  /*
  * An Matrix element is at (row,col)
  * The elements of a row stack vertically
  * Any row I is to the right of row I-1
  */
  int myxnpts = output.nrow();
  int myynpts = output.ncol();
  const double_t* myx=REAL(xseq);
  const double_t* myy=REAL(yseq);
  int ncol=tempmat_sky.ncol();
  int nrow=tempmat_sky.nrow();
  
  // For each vertical row 
  for (int i = 1; i <= myxnpts; i++) {
    // For a spline to interpolate vertically along the elements of the row
    double_t x = -0.5+i;
    // find the left and right index ibnto xseq
    int left_index = -1;
    int right_index = -1;
    for (int ii = 1; ii < nrow; ii++) {
      if (myx[ii-1] <= x && myx[ii] >= x) {
        left_index = ii-1;
        right_index = ii;
        break;
      }
    }
    
    //Rcpp::Rcout << "x="<<x<<" xindex="<<left_index<<" "<<right_index<<"\n";
    //Rcpp::Rcout << "x="<<x<<" xleft="<<myx[left_index]<<" "<<myx[right_index]<<"\n";
    int top_index = -1;
    int bottom_index = -1;
    for (int j = 1; j < myynpts; j++) {
      double y = -0.5+j;
      for (int jj = 1; jj < ncol; jj++) {
        if (myy[jj-1] <= y && myy[jj] >= y) {
          top_index = jj-1;
          bottom_index = jj;
          // p1...p2
          // .     .
          // .     .
          // p3...p4
          double p1 = tempmat_sky(left_index,top_index);
          double p2 = tempmat_sky(right_index,top_index);
          double p3 = tempmat_sky(left_index,bottom_index);
          double p4 = tempmat_sky(right_index,bottom_index);
        
          double xlambda = (x-myx[left_index])/(myx[right_index]-myx[left_index]);
          double ylambda = (y-myy[top_index])/(myy[bottom_index]-myy[top_index]);
          double ztop = p1 * (1.0-xlambda) + p2 * xlambda;
          double zbottom = p3 * (1.0-xlambda) + p4 * xlambda;
          output(i-1,j-1) = ztop * (1.0-ylambda) +zbottom * ylambda;
          //Rcpp::Rcout << "y="<<y<<" yindex="<<top_index<<" "<<bottom_index<<" "<<p1<<" "<<p2<<" "<<p3<<" "<<p4<<" result="<<output(i-1,j-1)<<"\n";
          break;
        }
      }
    }
  }
}

//==================================================================================
// [[Rcpp::export]]
double_t Cadacs_quantile(Rcpp::NumericVector x, double quantile) {
  const double_t* iiix=REAL(x);
  int size = x.size();
  std::vector<double_t> myx (iiix, iiix+size);
  double min=std::numeric_limits<double>::max();
  double max=-min;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      min = MIN(min,myx[i]);
      max = MAX(max,myx[i]);
    }
  }
  
  if (non_null_sample_count<1)
    return R_NaN;
  if (min==max) {
    return min;
  }
  
  // histogram
  int levels = 16384;
  std::vector<int> histogram;
  histogram.resize(levels);
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  double value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      int index = (myx[i] - min)*value_to_bin_index;
      histogram[index]++;
    }
  }
  int count=0;
  double quantileValue = min;
  double binwidth = (max - min)/levels;
  int target_count = non_null_sample_count*quantile;
  for (int i=0;i<levels;i++) {
    if (count>=target_count)
      return quantileValue;
    quantileValue += binwidth;
    count += histogram[i];
  }
  return quantileValue;
}
// [[Rcpp::export]]
double_t Cadacs_quantileLO(Rcpp::NumericVector x, double quantile, const double offset) {
  // The population we want the quantile for is x-offset where x<offset
  const double_t* iiix=REAL(x);
  int size = x.size();
  std::vector<double_t> myx (iiix, iiix+size);
  double min=std::numeric_limits<double>::max();
  double max=-min;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i]) && myx[i]<offset) {
      non_null_sample_count++;
      min = MIN(min,myx[i]);
      max = MAX(max,myx[i]);
    }
  }
  
  if (non_null_sample_count<1)
    return R_NaN;
  if (min==max) {
    return min;
  }
  
  // histogram
  int levels = 16384;
  std::vector<int> histogram;
  histogram.resize(levels);
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  double value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i]) && myx[i]<offset) {
      int index = (myx[i] - min)*value_to_bin_index;
      histogram[index]++;
    }
  }
  int count=0;
  double quantileValue = min-offset;
  double binwidth = (max - min)/levels;
  int target_count = non_null_sample_count*quantile;
  for (int i=0;i<levels;i++) {
    if (count>=target_count)
      return quantileValue;
    quantileValue += binwidth;
    count += histogram[i];
  }
  return quantileValue;
}
double_t Cadacs_quantileHI(Rcpp::NumericVector x, double quantile, const double offset) {
  // The population we want the quantile for is x-offset where x>offset
  const double_t* iiix=REAL(x);
  int size = x.size();
  std::vector<double_t> myx (iiix, iiix+size);
  double min=std::numeric_limits<double>::max();
  double max=-min;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i]) && myx[i]>offset) {
      non_null_sample_count++;
      min = MIN(min,myx[i]);
      max = MAX(max,myx[i]);
    }
  }
  
  if (non_null_sample_count<1)
    return R_NaN;
  if (min==max) {
    return min;
  }
  
  // histogram
  int levels = 16384;
  std::vector<int> histogram;
  histogram.resize(levels);
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  double value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i]) && myx[i]>offset) {
      int index = (myx[i] - min)*value_to_bin_index;
      histogram[index]++;
    }
  }
  int count=0;
  double quantileValue = min-offset;
  double binwidth = (max - min)/levels;
  int target_count = non_null_sample_count*quantile;
  for (int i=0;i<levels;i++) {
    if (count>=target_count)
      return quantileValue;
    quantileValue += binwidth;
    count += histogram[i];
  }
  return quantileValue;
}
// [[Rcpp::export]]
double_t Cadacs_mean(Rcpp::NumericVector x) {
  const double_t* myx=REAL(x);
  int size = x.size();
  double mean=0;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      mean += myx[i];
    }
  }
  if (non_null_sample_count==0)
    return R_NaN;
  return mean/non_null_sample_count;
}
// [[Rcpp::export]]
double_t Cadacs_population_variance(Rcpp::NumericVector x, const double offset) {
  const double_t* myx=REAL(x);
  int size = x.size();
  double v=0;
  double sum=0;
  double sum_sq=0;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      v = myx[i]-offset;
      sum += v;
      v *= v;
      sum_sq += v;
    }
  }
  if (non_null_sample_count<=1)
    return R_NaN;
  double N=non_null_sample_count;
  return sum_sq/N;
}
// [[Rcpp::export]]
double_t Cadacs_sample_variance(Rcpp::NumericVector x, const double offset) {
  const double_t* myx=REAL(x);
  int size = x.size();
  double v=0;
  double sum=0;
  double sum_sq=0;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      v = myx[i]-offset;
      sum += v;
      v *= v;
      sum_sq += v;
    }
  }
  if (non_null_sample_count<=1)
    return R_NaN;
  double N=non_null_sample_count;
  return (N*sum_sq - sum*sum)/(N * (N - 1));
  //return sqrt(sum_sq);
}
// [[Rcpp::export]]
double_t Cadacs_median(Rcpp::NumericVector x) {
  return Cadacs_quantile(x, 0.5);
}
double_t Cadacs_mode(Rcpp::NumericVector x) {
  const double_t* iiix=REAL(x);
  int size = x.size();
  std::vector<double_t> myx (iiix, iiix+size);
  double min=std::numeric_limits<double>::max();
  double max=std::numeric_limits<double>::min();
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      min = MIN(min,myx[i]);
      max = MAX(max,myx[i]);
    }
  }
  
  // histogram
  int levels = 16384*2;
  std::vector<int> histogram;
  histogram.resize(levels);
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  double value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      int index = (myx[i] - min)*value_to_bin_index;
      histogram[index]++;
    }
  }
  
  double current_bin_lower=min;
  double mode=current_bin_lower;
  
  double binwidth = (max - min)/levels;
  int max_count = 0;
  for (int i=0;i<levels;i++) {
    if (histogram[i]>max_count)
    {
      max_count = histogram[i];
      mode = current_bin_lower;
    }
    current_bin_lower += binwidth;
  }
  return mode;
}
// [[Rcpp::export]]
void adacsBothFromHistogram(Rcpp::NumericVector x, double quantile,Rcpp::NumericVector results) {
  const double_t* iiix=REAL(x);
  double_t* iresults=REAL(results);
  int size = x.size();
  std::vector<double_t> myx (iiix, iiix+size);
  double min=std::numeric_limits<double>::max();
  double max=std::numeric_limits<double>::min();
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      min = MIN(min,myx[i]);
      max = MAX(max,myx[i]);
    }
  }
  
  // histogram
  int levels = 16384*2;
  std::vector<int> histogram;
  histogram.resize(levels);
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  double value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      int index = (myx[i] - min)*value_to_bin_index;
      histogram[index]++;
    }
  }
  int count=0;
  double median = min;
  double binwidth = (max - min)/levels;
  int target_count = non_null_sample_count*0.5;
  for (int i=0;i<levels;i++) {
    if (count>=target_count)
    {
      levels = i;
      iresults[0] = median;
      iresults[2] = binwidth;
      break;
    }
    median += binwidth;
    count += histogram[i];
  }
  
  non_null_sample_count = count - histogram[levels];
  count=0;
  median = min;
  target_count = non_null_sample_count*quantile;
  for (int i=0;i<levels;i++) {
    if (count>=target_count)
    {
      iresults[1] = median;
      iresults[3] = binwidth;
      break;
    }
    median += binwidth;
    count += histogram[i];
  }
}
// [[Rcpp::export]]
void adacsBothFromHistogramV2(Rcpp::NumericVector x, double quantile,Rcpp::NumericVector results) {
  const double_t* iiix=REAL(x);
  double_t* iresults=REAL(results);
  int size = x.size();
  std::vector<double_t> myx (iiix, iiix+size);
  double min=std::numeric_limits<double>::max();
  double max=std::numeric_limits<double>::min();
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      min = MIN(min,myx[i]);
      max = MAX(max,myx[i]);
    }
  }
  
  // histogram
  int levels = 16384;
  std::vector<int> histogram;
  histogram.resize(levels);
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  double value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      int index = (myx[i] - min)*value_to_bin_index;
      histogram[index]++;
    }
  }
  int count=0;
  double median = min;
  double binwidth = (max - min)/levels;
  int target_count = non_null_sample_count*0.5;
  for (int i=0;i<levels;i++) {
    if (count>=target_count)
    {
      iresults[0] = median;
      iresults[2] = binwidth;
      break;
    }
    median += binwidth;
    count += histogram[i];
  }
  
  
  non_null_sample_count = 0;
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  max = median;
  value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      if (myx[i]<median) {
        non_null_sample_count++;
        int index = (myx[i] - min)*value_to_bin_index;
        histogram[index]++;
      }
    }
  }
  count=0;
  median = min;
  binwidth = (max - min)/levels;
  target_count = non_null_sample_count*quantile;
  for (int i=0;i<levels;i++) {
    if (count>=target_count)
    {
      iresults[1] = median-iresults[0];
      iresults[3] = binwidth;
      break;
    }
    median += binwidth;
    count += histogram[i];
  }
}
// [[Rcpp::export]]
Rcpp::NumericVector Cadacs_SkyEstLoc(Rcpp::NumericMatrix image,Rcpp::Nullable<Rcpp::IntegerMatrix> objects,Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
                      const double loc1, const double loc2,
                      const double box1, const double box2,
                      const double boxadd1, const double boxadd2,
                      const int skypixmin, const int boxiters,
                      const int doclip, const int skytype, const int skyRMStype, const double sigmasel
                      ) {
  Rcpp::NumericVector select = Cadacs_FindSkyCellValues(image, objects, mask, loc1, loc2, box1, box2, boxadd1, boxadd2, skypixmin, boxiters);
  Rcpp::NumericVector clip;
  if(doclip) {
    clip = Cadacs_magclip(select,adacs_AUTO,5,sigmasel,adacs_LO);
  } else {
    clip = select;
  }
  double skyloc=0.0;
  switch (skytype) {
  case adacs_MEDIAN:
    skyloc = Cadacs_median(clip);
    //Rcpp::Rcout<<"Cadacs_SkyEstLoc clip.size="<<clip.size()<<" select.size="<<select.size()<<" skyloc="<<skyloc<<"\n";
    break;
  case adacs_MEAN:
    skyloc = Cadacs_mean(clip);
    break;
  case adacs_MODE:
    skyloc = Cadacs_mode(clip);
    break;
  }
  
  double skyRMSloc=0.0;
  switch (skyRMStype) {
  case adacs_LO:
    skyRMSloc = fabs(Cadacs_quantileLO(clip,R::pnorm(-sigmasel, 0.0, 1.0, 1, 0)*2, skyloc))/sigmasel;
    //Rcpp::Rcout<<"Cadacs_SkyEstLoc clip.size="<<clip.size()<<" select.size="<<select.size()<<" skyRMSloc="<<skyRMSloc<<"\n";
    break;
  case adacs_HI:
    skyRMSloc = fabs(Cadacs_quantileHI(clip,R::pnorm(-sigmasel, 0.0, 1.0, 1, 0)*2, skyloc))/sigmasel;
    break;
  case adacs_BOTH:
  {
    double lo=fabs(Cadacs_quantileLO(clip,R::pnorm(-sigmasel, 0.0, 1.0, 1, 0)*2, skyloc))/sigmasel;
    double hi=fabs(Cadacs_quantileHI(clip,R::pnorm(-sigmasel, 0.0, 1.0, 1, 0)*2, skyloc))/sigmasel;
    skyRMSloc = (lo+hi)/2;
  }
    break;
  case adacs_SD:
    skyRMSloc = sqrt(Cadacs_population_variance(clip, skyloc));
    break;
  }
  Rcpp::NumericVector result(2);
  result[0] = skyloc;
  result[1] = skyRMSloc;
  return result;
}
// [[Rcpp::export]]
void Cadacs_MakeSkyGrid(Rcpp::NumericMatrix image,Rcpp::Nullable<Rcpp::IntegerMatrix> objects,Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
                                     const int box1, const int box2,
                                     const int grid1, const int grid2,
                                     const int boxadd1, const int boxadd2,
                                     const int type, const int skypixmin, const int boxiters,
                                     const int doclip, const int skytype, const int skyRMStype, const double sigmasel,
                                     Rcpp::NumericMatrix sky, Rcpp::NumericMatrix skyRMS
) {
  // box MUST NOT be larger than the input image
  double box[2] = {(double)box1, (double)box2};
  if(box[0]>image.nrow())
    box[0]=image.nrow();
  if(box[1]>image.ncol())
    box[1]=image.ncol();
  
  double grid[2] = {(double)grid1, (double)grid2};
  if(grid[0]>image.nrow())
    grid[0]=image.nrow();
  if(grid[1]>image.ncol())
    grid[1]=image.ncol();
  
  // tile over input image with tile size (grid) and no overlap
  // xseq,yseq give the centres of each tile
  int tile_nrows=0;
  double x_tile_centre=grid[0]/2;
  while (x_tile_centre<image.nrow()) {
    tile_nrows++;
    x_tile_centre += grid[0];
  }
  int tile_ncols=0;
  double y_tile_centre=grid[1]/2;
  while (y_tile_centre<image.ncol()) {
    tile_ncols++;
    y_tile_centre += grid[1];
  }
  
  // add room for linearly extrapolated padding
  tile_nrows += 2;
  tile_ncols += 2;
  
  // Construct the vector of tile centroids
  Rcpp::NumericVector xseq(tile_nrows);
  Rcpp::NumericVector yseq(tile_ncols);
  x_tile_centre=grid[0]/2 - grid[0];
  for (int i=0; i<tile_nrows; i++) {
    xseq[i] = x_tile_centre;
    x_tile_centre += grid[0];
  }
  y_tile_centre=grid[1]/2 - grid[1];
  for (int i=0; i<tile_ncols; i++) {
    yseq[i] = y_tile_centre;
    y_tile_centre += grid[1];
  }
  
  Rcpp::NumericMatrix z_sky_centre(tile_nrows, tile_ncols);
  Rcpp::NumericMatrix z_skyRMS_centre(tile_nrows, tile_ncols);
  
  bool hasNaNs=false;
  x_tile_centre=grid[0]/2;
  for (int i=1; i<tile_nrows-1; i++) {
    x_tile_centre = xseq[i];
    for (int j=1; j<tile_ncols-1; j++) {
      y_tile_centre = yseq[j];
      Rcpp::NumericVector z_tile_centre = Cadacs_SkyEstLoc(image,objects,mask,
                                                           x_tile_centre, y_tile_centre,
                                                           box1, box2,
                                                           boxadd1, boxadd2,
                                                           skypixmin, boxiters,
                                                           doclip, skytype, skyRMStype, sigmasel);
      if (std::isnan(z_tile_centre[0]) || std::isnan(z_tile_centre[1])) {
        hasNaNs = true;
      }
        
      z_sky_centre(i, j) = z_tile_centre[0];
      z_skyRMS_centre(i, j) = z_tile_centre[1];
    }
  }
  
  if (hasNaNs) {
    // Replace any NaN's with reasonable substitute
    // initialise the pad area before getting the medians
    for (int i=0; i<tile_nrows; i++) {
      z_sky_centre(i,0) = R_NaN;
      z_sky_centre(i,tile_ncols-1) = R_NaN;
    }
    for (int i=0; i<tile_ncols; i++) {
      z_sky_centre(0, i) = R_NaN;
      z_sky_centre(tile_nrows-1, i) = R_NaN;
    }
    double medianSkyCentre=Cadacs_median(z_sky_centre);
    double medianSkyRMSCentre=Cadacs_median(z_skyRMS_centre);
    // replace NaN's now
    for (int i=1; i<tile_nrows-1; i++) {
      for (int j=1; j<tile_ncols-1; j++) {
        if (std::isnan(z_sky_centre(i, j)))
          z_sky_centre(i, j) = medianSkyCentre;
        if (std::isnan(z_skyRMS_centre(i, j)))
          z_skyRMS_centre(i, j) = medianSkyRMSCentre;
      }
    }
  }
  
  // Padding
  //work out the second point for linear extrapolation (the first one is at 1+1 and length(seq)-1)
  int xstart=MIN(2,tile_nrows-2);
  int ystart=MIN(2,tile_ncols-2);
  int xend=MAX(tile_nrows-3,1);
  int yend=MAX(tile_ncols-3,1);
  for (int i=0; i<tile_nrows; i++) {
    z_sky_centre(i,0) = z_sky_centre(i, 1)*2 - z_sky_centre(i, ystart);
    z_sky_centre(i,tile_ncols-1) = z_sky_centre(i, tile_ncols-2)*2 - z_sky_centre(i, yend);
  }
  for (int i=0; i<tile_ncols; i++) {
    z_sky_centre(0, i) = z_sky_centre(1, i)*2 - z_sky_centre(xstart, i);
    z_sky_centre(tile_nrows-1, i) = z_sky_centre(tile_nrows-2, i)*2 - z_sky_centre(xend, i);
  }
  
  // Now interpolate for each image cell
  
  switch (type) {
  case adacs_CLASSIC_BILINEAR:
    interpolateLinearGrid(xseq, yseq, z_sky_centre, sky);
    interpolateLinearGrid(xseq, yseq, z_skyRMS_centre, skyRMS);
    break;
  case adacs_AKIMA_BICUBIC:
    interpolateAkimaGrid(xseq, yseq, z_sky_centre, sky);
    interpolateAkimaGrid(xseq, yseq, z_skyRMS_centre, skyRMS);
    break;
  }
  
  // Apply mask
  if (mask.isNotNull()) {
    Rcpp::IntegerMatrix imask = Rcpp::as<Rcpp::IntegerMatrix>(mask);
    int nrows=image.nrow();
    int ncols=image.ncol();
    for (int i=0; i<ncols; i++) {
      for (int j=0; j<nrows; j++) {
        if (imask(j, i)==1) {
          sky(j, i) = R_NaN;
          skyRMS(j, i) = R_NaN;
        }
      }
    }
  }
}
