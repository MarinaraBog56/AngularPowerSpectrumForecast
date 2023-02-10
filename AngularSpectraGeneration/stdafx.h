#pragma once

#ifndef STDAFX_H
#define STDAFX_H

// Pre-compiled header file

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <thread>
#include <algorithm>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf.h>

#define SPEED_OF_LIGHT (1.0)
#define MEGA_PARSEC (1.0)
#define SPEED_OF_LIGHT_H0 (299792.458)
#define PI M_PI
#define H_0 ((100) / (MEGA_PARSEC*SPEED_OF_LIGHT_H0))
#define CUTOFF 10000 // 141000 is max value
#define WORKMAX 60000
//const std::string FILEPATH = "/Users/cb607/Documents/Power_Spectrum_Calc";
const std::string FILEPATH = "/cosma8/data/dp203/dc-brow6";
const double K_NON_LIN = 0.2;
const double K_LIMIT_SHEAR = 9.5;
const double K_LIMIT = 0.47362;
const double Z_LOW = 0.0001;

#endif // !STDAFX_H
