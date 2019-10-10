// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// rcpp_module.cpp: Rcpp R/C++ interface class library -- Rcpp Module examples
//
// Copyright (C) 2010 - 2012  Dirk Eddelbuettel and Romain Francois
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>

class BitMatrix {
public:
    BitMatrix() {}
    BitMatrix(int nrows, int ncols) {
        // rows are "vertical"
        _nrows = nrows;
        _ncols = ncols;
        _npts = nrows*ncols;
        _n32bitwords = 1 + (_npts >> 5);
        _data.resize(_n32bitwords);
    }
    void fill(bool yesno) {
        for (int i=0;i<_nrows;i++) {
            for (int j=0;j<_ncols;j++) {
                if (yesno) {
                    settrue(i,j);
                } else {
                    setfalse(i,j);
                }
            }
        }
    }
    bool settrue(uint32_t row, uint32_t col) {
        uint32_t index = col*_nrows+row;
        uint32_t word, bit;
        bool current;

        word    = index >> 5;
        bit     = index & 31;   // % 32
        current = _data[word] & (1 << bit);
        _data[word] |= (1 << bit);
        return current;
    }
    bool setfalse(uint32_t row, uint32_t col) {
        uint32_t index = col*_nrows+row;
        uint32_t word, bit;
        bool current;

        word    = index >> 5;
        bit     = index & 31;   // % 32
        current = _data[word] & (1 << bit);
        _data[word] &= (~(1 << bit));
        return current;
    }
    bool istrue(uint32_t row, uint32_t col) const {
        uint32_t index = col*_nrows+row;
        uint32_t word, bit;
        bool current;

        word        = index >> 5;
        bit         = index & 31;       // % 32
        current     = _data[word] & (1 << bit);
        return current;
    }
    bool isfalse(uint32_t row, uint32_t col) const {
        uint32_t index = col*_nrows+row;
        uint32_t word, bit;
        bool current;

        word        = index >> 5;
        bit         = index & 31;       // % 32
        current     = _data[word] & (1 << bit);
        return ~current;
    }

private:
    uint32_t _nrows;
    uint32_t _ncols;
    uint32_t _npts;
    uint32_t _n32bitwords;
    std::vector<uint32_t> _data;
};



RCPP_MODULE(yada){
    using namespace Rcpp ;

    class_<BitMatrix>("BitMatrix")
    // expose the default constructor
    .constructor()
    .constructor<int, int>()

    .method("fill", &BitMatrix::fill     , "set or clear all bits")
    .method("settrue", &BitMatrix::settrue     , "set the bit")
    .method("setfalse", &BitMatrix::setfalse     , "clear the bit")
    .const_method("istrue", &BitMatrix::istrue     , "test if bit set")
    .const_method("isfalse", &BitMatrix::istrue     , "test if bit clear")
    ;
}


