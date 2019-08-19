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
// [[Rcpp::export]]
Rcpp::NumericVector adacsFindSkyCellValuesC(Rcpp::NumericMatrix image, Rcpp::Nullable<Rcpp::IntegerMatrix> objects, Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
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
  Rcpp::NumericVector vec(skyN);
  int k=0;
  for (int j = sscol; j <= eecol; j++) {
    int ii=(j-1)*ncol+(ssrow-1);
    for (int i = ssrow; i <= eerow; i++,ii++) {
      if ((iiobjects!=NULL && iiobjects[ii]==0) && (iimask==NULL || iimask[ii]==0)) {
        vec[k++] = iiimage[ii];
      }
    }
  }
  return vec;
}
// [[Rcpp::export]]
Rcpp::IntegerVector adacsFindSkyCellValuesBoxC(Rcpp::NumericMatrix image, Rcpp::Nullable<Rcpp::IntegerMatrix> objects, Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
                                            const double loc1, const double loc2,
                                            const double box1, const double box2,
                                            const double boxadd1, const double boxadd2,
                                            const int skypixmin, const int boxiters)
{
  //Rcpp::Rcout << "loc(" << loc1 << "," << loc2 << ")\n";
  //Rcpp::Rcout << "box(" << box1 << "," << box2 << ")\n";
  //Rcpp::Rcout << "boxadd(" << boxadd1 << "," << boxadd2 << ")\n";
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