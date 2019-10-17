// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// rcpp_module.cpp: Rcpp R/C++ interface class library -- Rcpp Module implementation for Bit Matrix
//
// Copyright (C) 2019 Ray Seikel
//
// This file is part of BitMatrix.
//
// BitMatrix is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// BitMatrix is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.
#ifndef rcpp_module_hpp
#define rcpp_module_hpp
#include <Rcpp.h>
using namespace Rcpp;
RCPP_EXPOSED_CLASS(BitMatrix)
#include "bitmatrix.h"
#endif

