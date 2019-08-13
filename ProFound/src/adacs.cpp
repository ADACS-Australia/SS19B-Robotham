//
// Optimisation's contributed by ADACS
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
  //Rcpp::Rcout << " cols   " << scol << " " << ecol << "\n";
  //Rcpp::Rcout << " rows   " << srow << " " << erow << "\n";
  //Rcpp::Rcout << " offset " << coffset << " " << roffset << "\n";
  int onrow = oimage.nrow();
  int oncol = oimage.ncol();
  // 0 relative
  int sscol = scol-1;
  int eecol = ecol-1;
  int ssrow = srow-1;
  int eerow = erow-1;
  
  for (int i = 0; i < onrow; i++) {
    for (int j = 0; j < oncol; j++) {
      oimage(i,j) = NA_REAL;
    }
  }
  for (int i = ssrow; i <= eerow; i++) {
    for (int j = sscol; j <= eecol; j++) {
      //Rcpp::Rcout << "   in " << i-roffset-ssrow << " " << j-roffset-sscol << "\n";
      oimage(i-roffset,j-coffset) = image(i, j);
    }
  }
}
