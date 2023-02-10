#pragma once

#ifndef ISWINTEGRATOR_HPP
#define ISWINTEGRATOR_HPP


// hpp file for integrating Power spectra for C_ells
// Elizabeth Brown

#include "stdafx.h"
#include "ISWODESolver.hpp"


// structures for use in integral parameters
//struct spline_params {gsl_spline * spline_ptr; gsl_interp_accel * accel_ptr;};
struct A_Bessel_params {gsl_spline * spline_ptr1; gsl_interp_accel * accel_ptr1; 
                        gsl_spline * spline_ptr2; gsl_interp_accel * accel_ptr2; int l; double k;};
struct spline_spline_params {gsl_spline * bessel_ptr; gsl_interp_accel * bessel_accel_ptr; 
                        gsl_spline * spline_ptr; gsl_interp_accel * spline_accel_ptr; 
                        gsl_spline * chi_ptr; gsl_interp_accel * chi_accel_ptr; double k;};
struct normalisation_params {double z0;};
struct parameterFile {std::string inputDirectory; int numParams; std::string* paramNames;
                        bool* calcParam; double h; double hDev; int numGal; double* galFSky; double* galNumDen;
                        double* galNumDenDev; double* galZLow; double* galZHigh; double* galZMed; std::string* galSuffix;
                        int numLens; double* lensFSky; double* lensNumden; double* lensNumDenDev; double* lensZLow;
                        double* lensZHigh; double* lensZMed; std::string* lensSuffix; bool calcCrossCov;
                        int* kernInst; int numObs; bool calcSpeedUp; std::string* obsName; bool* calcObs;
                        int** calcKernel; std::string outputDirectory;};



///////////////////////////
// FUNCTION DECLARATIONS //
///////////////////////////

std::pair<gsl_spline*, gsl_interp_accel*> spline1dData(double * x, double * y, int arrSize);
std::pair<gsl_spline*, gsl_interp_accel*> spline1dSteffenData(double * x, double * y, int arrSize);
double splineToIntegrate(double x, void * p);
double BesselToIntegrate(double x, void * p);
double splineSplineToIntegrate(double x, void * p);
std::pair<double, double> integrateForChi(std::pair<gsl_spline*, gsl_interp_accel*> a2Hinv, double a_low, double a_high=1);
std::pair<gsl_spline *, gsl_interp_accel *> chiSpline(double* a, double* loga, double * a2Hinv, int a_length);
double integrateODE(double k, std::pair<gsl_spline *, gsl_interp_accel *> Bessel, 
                        std::pair<gsl_spline *, gsl_interp_accel *> function, double chiLow, double chiHigh, 
                        int a_length);
void calcKernelIntegrals(parameterFile paramFile, double**& angularPowerKernels, double*** densityKernels, int* kernels,
                    double* k, double* kChi, double* chi, double chiLow, double chiHigh, double* lensingSources, 
                    int k_length, int a_length, int lmax, int l, int BesselArray_length);
void angularPowerSpectrumCurvedSky(double** angularPowerSpectra, double*** densityKernels, parameterFile paramFile,
                                    double* P_primordial, int* kernels, int* kernelInstruments, int** obsIndices,
                                    double* k, double* kChi, double* chi, double chiLow, double chiHigh,
                                    double* lensingSources, double** galaxyBins, int k_length, int a_length, 
                                    int l_max, int multiplier, int offset, int BesselArray_length);

void threadFunctionLimber(double* C_ell, std::pair<gsl_spline *, gsl_interp_accel *>* Power, 
                    double* scale, double* loga, double* HubbleRate, int a_length, 
                    std::pair<gsl_spline *, gsl_interp_accel *> chi, int l_lim, int l_start, int l_max, int multiplier, 
                    int offset, int k_length, double aLow, double aHigh, double kmax);
std::pair<double, double> BesselLimitChecker(double* a, int l, double k,  double aLow,
                                                std::pair<gsl_spline *, gsl_interp_accel *> chi,
                                                std::pair<gsl_spline *, gsl_interp_accel *> growth);                    
double timeFunction(std::pair<gsl_spline *, gsl_interp_accel *> chi, std::pair<gsl_spline *, gsl_interp_accel *>  A, 
                    double k, int l, double a);
double** timeIntegralFunction(std::pair<gsl_spline *, gsl_interp_accel *> chi, double* A, double k, double* a,
                                int a_length, int l);
int kChiLength();
double* kChiArray(int kChi_length);
double* BesselArray(double* kChi, int ell, int BesselArray_length);
bool BesselFileChecker(int ell);
void BesselFileMaker(double* Bessel, int ell, int BesselArray_length);
std::pair<gsl_spline *, gsl_interp_accel *> BesselSpline(double* kChi, int ell, int BesselArray_length);


///////////////
// FUNCTIONS //
///////////////


///////////////////////////////////////////
// SPLINES, CHI, & INTEGRATION FUNCTIONS //
///////////////////////////////////////////

// Creates 1D spline using cubic of data
std::pair<gsl_spline*, gsl_interp_accel*> spline1dData(double * x, double * y, int arrSize) {
    gsl_interp_accel * my_accel_ptr;
    gsl_spline * spline_ptr;
    my_accel_ptr = gsl_interp_accel_alloc();
    spline_ptr = gsl_spline_alloc(gsl_interp_cspline, arrSize);
    gsl_spline_init (spline_ptr, x, y, arrSize);
    std::pair<gsl_spline*, gsl_interp_accel*> Spline;
    Spline.first = spline_ptr;
    Spline.second = my_accel_ptr;
    return Spline;
}

// Creates 1D spline using steffen of data
std::pair<gsl_spline*, gsl_interp_accel*> spline1dSteffenData(double * x, double * y, int arrSize) {
    gsl_interp_accel * my_accel_ptr;
    gsl_spline * spline_ptr;
    my_accel_ptr = gsl_interp_accel_alloc();
    spline_ptr = gsl_spline_alloc(gsl_interp_steffen, arrSize);
    gsl_spline_init(spline_ptr, x, y, arrSize);
    std::pair<gsl_spline*, gsl_interp_accel*> Spline;
    Spline.first = spline_ptr;
    Spline.second = my_accel_ptr;
    return Spline;
}

// Provides a spline function to integrate
double splineToIntegrate(double x, void * p) {
    // sort parameters here
    struct spline_params * params = (struct spline_params *)p;
    gsl_spline * spline = (params->spline_ptr);
    gsl_interp_accel * accelerator = (params->accel_ptr);

    // calculate numerical value of function
    return exp(x) * gsl_spline_eval(spline, exp(x), accelerator);
}

double limberToIntegrate(double x, void* p) {
    // sort parameters here
    struct spline_params * params = (struct spline_params *)p;
    gsl_spline * spline = (params->spline_ptr);
    gsl_interp_accel * accelerator = (params->accel_ptr);

    // calculate numerical value of function
    return gsl_spline_eval(spline, x, accelerator);
}

// Provides a (Bessel x spline) function to integrate
double BesselToIntegrate(double x, void * p) {
    // sort parameters here
    struct A_Bessel_params * params = (struct A_Bessel_params *)p;
    gsl_spline * spline1 = (params->spline_ptr1);
    gsl_interp_accel * accelerator1 = (params->accel_ptr1);
    gsl_spline * spline2 = (params->spline_ptr2);
    gsl_interp_accel * accelerator2 = (params->accel_ptr2);
    int l = (params->l);
    double k = (params->k);

    // calculate numerical value of function
    double Chi = gsl_spline_eval(spline1, x, accelerator1);
    double Bessel = gsl_sf_bessel_jl(l, Chi*k);
    double A = gsl_spline_eval(spline2, exp(x), accelerator2);
    return SPEED_OF_LIGHT * Bessel * A;
}

// Provides a (Bessel x spline) function to integrate
double splineSplineToIntegrate(double x, void * p) {
    // sort parameters here
    struct spline_spline_params * params = (struct spline_spline_params *)p;
    gsl_spline * bessel = (params->bessel_ptr);
    gsl_interp_accel * bessel_acc = (params->bessel_accel_ptr);
    gsl_spline * spline = (params->spline_ptr);
    gsl_interp_accel * spline_acc = (params->spline_accel_ptr);
    gsl_spline * chi = (params->chi_ptr);
    gsl_interp_accel * chi_acc = (params->chi_accel_ptr);
    double k = (params->k);

    // calculate numerical value of function
    double Chi = gsl_spline_eval(chi, x, chi_acc);
    double Bessel;
    if (k*Chi<CUTOFF) {
        Bessel = gsl_spline_eval(bessel, Chi*k, bessel_acc);
    }
    else {
        Bessel=0;
    }
    double A = gsl_spline_eval(spline, exp(x), spline_acc);
    return SPEED_OF_LIGHT * Bessel * A * exp(x);
}

// Provides a normalisation function to integrate
double normalisationToIntegrate(double x, void * p) {
    // Sort parameters
    struct normalisation_params * params = (struct normalisation_params *)p;
    double z0 = (params->z0);

    // calculate selection normalisation
    double exponent = -pow(x/z0, 1.5);
    double y = (1.5*x*x/(z0*z0*z0))*exp(exponent);
    return y;
}

// Integrates radial geodesic distance from a->1
std::pair<double, double> integrateForChi(std::pair<gsl_spline*, gsl_interp_accel*> a2Hinv, double a_low, double a_high) {
    // allocate max number of steps (100,000)
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000000);
    // initialise answer variables
    double result, error;
    std::pair<double, double> Result;
    // initialise function
    gsl_function F;
    struct spline_params params = { a2Hinv.first, a2Hinv.second };
    F.function = &splineToIntegrate;
    F.params = &params;
    // input into integration function with relevant limits (1e-7, 1e-7 represent absolute and relative errors).
    gsl_integration_qags (&F, a_low, a_high, 0, 1e-4, 10000000, w, &result, &error);
    Result.first = result;
    Result.second = error;
    // free up workspace
    gsl_integration_workspace_free (w);
    return Result;
}

// Return a spline of chi to be used in Bessel functions
std::pair<gsl_spline *, gsl_interp_accel *> chiSpline(double* a, double* loga, double * a2Hinv, int a_length) {
    // make spline of a2Hinv (called a2H).
    std::pair<gsl_spline *, gsl_interp_accel*> aFunc;
    aFunc = spline1dSteffenData(a, a2Hinv, a_length);
    double * chi = new double[a_length];
    // calculate chi for multiple a lower limits, and use result to make spline for chi
    for (int i = 0; i < a_length; i++) {
        std::pair<double, double> result;
        result = integrateForChi(aFunc, loga[i], 0);
        chi[i] = result.first;
    }
    // free up a2Hinv spline (no longer needed).
    gsl_spline_free (aFunc.first);
    gsl_interp_accel_free (aFunc.second);
    std::pair<gsl_spline *, gsl_interp_accel *> Chi = spline1dSteffenData(loga, chi, a_length);
    delete[] chi;
    return Chi;
}

std::pair<double, double> renormalise(double z0, double a_low, double a_high) {
    double z_high = (1/a_low)-1;
    double z_low = (1/a_high)-1;
    // allocate max number of steps (100,000)
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000000);
    // initialise answer variables
    double result, error;
    std::pair<double, double> Result;
    // initialise function
    gsl_function F;
    struct normalisation_params params = { z0 };
    F.function = &normalisationToIntegrate;
    F.params = &params;
    // input into integration function with relevant limits (0, 1e-4 represent absolute and relative errors).
    gsl_integration_qags (&F, z_low, z_high, 0, 1e-4, 10000000, w, &result, &error);
    Result.first = result;
    Result.second = error;
    // free up workspace
    gsl_integration_workspace_free (w);
    return Result;
}

double numDenIntegration(double zp, double z0, double cb, double zb, double sigmab, double co, double zo, 
                        double sigmao, double fout, double z_low, double z_high, double* redshift, int a_length) {
    double* integrand = new double[a_length];
    for(int i = 0; i<a_length; i++) {
        double z = redshift[i];
        double exponent = -pow(z/z0, 1.5);
        double nz = (z*z/(z0*z0))*exp(exponent);
        double probfin = (1-fout)/(sigmab*(1+z)*pow(2.*PI, 0.5));
        double probfinexp = exp(-0.5*(pow(((z-(cb*zp)-zb)/(sigmab*(1+z))),2.)));
        double probfout = (fout)/(sigmao*(1+z)*pow(2.*PI, 0.5));
        double probfoutexp = exp(-0.5*(pow(((z-(co*zp)-zo)/(sigmao*(1+z))),2.)));
        integrand[i] = nz*((probfout*probfoutexp)+(probfin*probfinexp));
    }
    std::pair<gsl_spline *, gsl_interp_accel *> numDenInt = spline1dSteffenData(redshift, integrand, a_length);
    delete[] integrand;
    double numDenResult = gsl_spline_eval_integ(numDenInt.first, z_low, z_high, numDenInt.second);
    gsl_spline_free (numDenInt.first);
    gsl_interp_accel_free (numDenInt.second);
    return numDenResult;
}

double numDenIntegrationOriginal(double zz, double z0, double cb, double zb, double sigmab, double co, double zo, 
                        double sigmao, double fout, double z_low, double z_high, double* redshift, int a_length) {
    double* integrand = new double[a_length];
    for(int i = 0; i<a_length; i++) {
        double zp = redshift[i];
        double exponent = -pow(zz/z0, 1.5);
        double nz = (zz*zz/(z0*z0))*exp(exponent);
        double probfin = (1-fout)/(sigmab*(1+zz)*pow(2.*PI, 0.5));
        double probfinexp = exp(-0.5*(pow(((zz-(cb*zp)-zb)/(sigmab*(1+zz))),2.)));
        double probfout = (fout)/(sigmao*(1+zz)*pow(2.*PI, 0.5));
        double probfoutexp = exp(-0.5*(pow(((zz-(co*zp)-zo)/(sigmao*(1+zz))),2.)));
        integrand[i] = nz*((probfout*probfoutexp)+(probfin*probfinexp));
    }
    std::pair<gsl_spline *, gsl_interp_accel *> numDenInt = spline1dSteffenData(redshift, integrand, a_length);
    delete[] integrand;
    double numDenResult = gsl_spline_eval_integ(numDenInt.first, z_low, z_high, numDenInt.second);
    gsl_spline_free (numDenInt.first);
    gsl_interp_accel_free (numDenInt.second);
    return numDenResult;
}



////////////////////
// EXACT SOLUTION //
////////////////////


double integrateODE(double k, std::pair<gsl_spline *, gsl_interp_accel *> Bessel, 
                        std::pair<gsl_spline *, gsl_interp_accel *> function, double chiLow, double chiHigh, 
                        int a_length) {

    double step = (chiHigh-chiLow)/1000;
    if(k*step>0.1) {
        step = 0.01/k;
    }
    return solveSplineBesselODE(chiLow, chiHigh, step, function, Bessel, k).second;
}

void calcKernelIntegrals(parameterFile paramFile, double**& angularPowerKernels, double*** densityKernels, int* kernels,
                    double* k, double* kChi, double* chi, double chiLow, double chiHigh, double* lensingSources, 
                    double** galaxyBins, int k_length, int a_length, int lmax, int l, int BesselArray_length) {

    if(l<=lmax) {
        std::pair<gsl_spline *, gsl_interp_accel *> Bessel = BesselSpline(kChi, l, BesselArray_length);
        bool Integrate = true;
        double* kernelMax = new double[2+(2*paramFile.numGal)+(2*paramFile.numLens)]();
        int* kernelIntegrate = new int[2+(2*paramFile.numGal)+(2*paramFile.numLens)]();
        double* kernelLim = new double[2+(2*paramFile.numGal)+(2*paramFile.numLens)]();
        double integralCutoff = 1000000000.;

        // Initialise the results array
        angularPowerKernels = new double*[2+(2*paramFile.numGal)+(2*paramFile.numLens)]();
        int angKerIt = 0;
        for(int j=0; j<6; j++) {
            if(kernels[j]>0) {  // if calculating kernel header (T, CMBL, g, gL, L, IA)
                if(j==2 || j==3) {
                    for(int g = 0; g<paramFile.numGal; g++) {
                        angularPowerKernels[angKerIt] = new double[k_length]();
                        kernelMax[angKerIt] = 0;
                        kernelIntegrate[angKerIt] = 1;
                        angKerIt++;
                    }
                }
                else if(j==4 || j==5) {
                    for(int g = 0; g<paramFile.numLens; g++) {
                        angularPowerKernels[angKerIt] = new double[k_length]();
                        kernelMax[angKerIt] = 0;
                        kernelIntegrate[angKerIt] = 1;
                        angKerIt++;
                    }
                }
                else {
                    angularPowerKernels[angKerIt] = new double[k_length]();
                    kernelMax[angKerIt] = 0;
                    kernelIntegrate[angKerIt] = 1;
                    angKerIt++;
                }
            }
            else {  // if not calculating kernel header
                if(j==2 || j==3) {
                    for(int g = 0; g<paramFile.numGal; g++) {
                        angularPowerKernels[angKerIt] = new double[k_length]();
                        kernelMax[angKerIt] = 0;
                        kernelIntegrate[angKerIt] = 0;
                        angKerIt++;
                    }
                }
                else if(j==4 || j==5) {
                    for(int g = 0; g<paramFile.numLens; g++) {
                        angularPowerKernels[angKerIt] = new double[k_length]();
                        kernelMax[angKerIt] = 0;
                        kernelIntegrate[angKerIt] = 0;
                        angKerIt++;
                    }
                }
                else {
                    angularPowerKernels[angKerIt] = new double[k_length]();
                    kernelMax[angKerIt] = 0;
                    kernelIntegrate[angKerIt] = 0;
                    angKerIt++;
                }
            }
        }
        // Perform the integral
        double lowerBound, upperBound;
        for (int i = 0; i<k_length; i++) {
            double upperLim;
            if(chiHigh*k[i]>CUTOFF) {upperLim=CUTOFF/k[i];}
            else {upperLim=chiHigh;}
            if(Integrate) {  // If we are still integrating
                angKerIt = 0;
                for(int j = 0; j<6; j++) {  // iterating over kernel headers (T, CMBL, g, gL, L, IA)
                    if(kernels[j]>0) {  // if kernel is calculated at least once
                        if(j==2 || j==3) {  // if kernel header is "g" or "gL"
                            for(int g = 0; g<paramFile.numGal; g++) {
                                std::pair<gsl_spline *, gsl_interp_accel *> kernelFunc = spline1dSteffenData(chi, densityKernels[angKerIt][i], a_length);
                                lowerBound = std::max(galaxyBins[g][0], chiLow); upperBound = std::min(galaxyBins[g][1], upperLim);
                                if(upperBound<lowerBound) { 
                                    angularPowerKernels[angKerIt][i] = 0;
                                }
                                else {
                                    angularPowerKernels[angKerIt][i] = integrateODE(k[i], Bessel, kernelFunc, lowerBound, upperBound, a_length);
                                }
                                gsl_spline_free (kernelFunc.first);
                                gsl_interp_accel_free (kernelFunc.second);
                                double absValue = angularPowerKernels[angKerIt][i];
                                if(absValue<0) {absValue=-1.*absValue;}
                                if(paramFile.calcSpeedUp) {  // if calculation is being sped-up
                                    if(absValue>kernelMax[angKerIt]) {
                                        kernelMax[angKerIt] = absValue;  // record new highest kernel value
                                    }
                                    else if(absValue<kernelMax[angKerIt]/integralCutoff) {
                                        kernelIntegrate[angKerIt] = 0;  // Smallest value reached, set Integrate to F
                                    }
                                }
                                angKerIt++;
                            }
                        }
                        else if(j==4 || j==5) {  // if kernel header is "L" or "IA"
                            for(int g = 0; g<paramFile.numLens; g++) {
                                std::pair<gsl_spline *, gsl_interp_accel *> kernelFunc = spline1dSteffenData(chi, densityKernels[angKerIt][i], a_length);
                                lowerBound = chiLow; upperBound = std::min(lensingSources[g], upperLim);
                                angularPowerKernels[angKerIt][i] = integrateODE(k[i], Bessel, kernelFunc, lowerBound, upperBound, a_length);
                                gsl_spline_free (kernelFunc.first);
                                gsl_interp_accel_free (kernelFunc.second);
                                double absValue = angularPowerKernels[angKerIt][i];
                                if(absValue<0) {absValue=-1.*absValue;}
                                if(paramFile.calcSpeedUp) {  // if calculation is being sped-up
                                    if(absValue>kernelMax[angKerIt]) {
                                        kernelMax[angKerIt] = absValue;  // record new highest kernel value
                                    }
                                    else if(absValue<kernelMax[angKerIt]/integralCutoff) {
                                        kernelIntegrate[angKerIt] = 0;  // Smallest value reached, set Integrate to F
                                    }
                                }
                                angKerIt++;
                            }
                        }
                        else {  // if kernel header is "T" or "CMBL"
                            std::pair<gsl_spline *, gsl_interp_accel *> kernelFunc = spline1dSteffenData(chi, densityKernels[angKerIt][i], a_length);
                            lowerBound = chiLow; upperBound = upperLim;
                            angularPowerKernels[angKerIt][i] = integrateODE(k[i], Bessel, kernelFunc, lowerBound, upperBound, a_length);
                            gsl_spline_free (kernelFunc.first);
                            gsl_interp_accel_free (kernelFunc.second);
                            double absValue = angularPowerKernels[angKerIt][i];
                            if(absValue<0) {absValue=-1.*absValue;}
                            if(paramFile.calcSpeedUp) {  // if calculation is being sped-up
                                if(absValue>kernelMax[angKerIt]) {
                                    kernelMax[angKerIt] = absValue;
                                }
                                else if(absValue<kernelMax[angKerIt]/integralCutoff) {
                                    kernelIntegrate[angKerIt] = 0;
                                }
                            }
                            angKerIt++;
                        }
                    }
                    else {  // if kernel is never calculated
                        if(j==2 || j==3) {  // if kernel header is "g" or "gL"
                            for(int g = 0; g<paramFile.numGal; g++) {
                                angularPowerKernels[angKerIt][i] = 0;
                                angKerIt++;
                            }
                        }
                        else if(j==4 || j==5) {  // if kernel header is "L" or "IA"
                            for(int g = 0; g<paramFile.numLens; g++) {
                                angularPowerKernels[angKerIt][i] = 0;
                                angKerIt++;
                            }
                        }
                        else {  // if kernel header is "T" or "CMBL"
                            angularPowerKernels[angKerIt][i] = 0;
                            angKerIt++;
                        }
                    }

                }
            }
            int sum = 0;
            for(int sumIt=0; sumIt<2+(2*paramFile.numGal)+(2*paramFile.numLens); sumIt++) {  // sum kernel Integrate T/F
                sum=sum+kernelIntegrate[sumIt];
            }
            if(sum==0) {  // if all kernel Integrate entries are F
                Integrate = false;  // stop integrating
                for(int j = 0; j<2+(2*paramFile.numGal)+(2*paramFile.numLens); j++) {
                    double absValue = angularPowerKernels[j][i];
                    if(absValue<0) {absValue = -1.*absValue;}
                    kernelLim[j] = absValue;  // record smallest values for each kernel
                }
            }
        }
        for (int i = 0; i<k_length; i++) {
            for(int j = 0; j<2+(2*paramFile.numGal)+(2*paramFile.numLens); j++) {
                double absValue = angularPowerKernels[j][i];
                if(absValue<0) {absValue = -1.*absValue;}
                if(absValue>kernelLim[j]) {
                    angularPowerKernels[j][i]=angularPowerKernels[j][i];
                }
                else {
                    angularPowerKernels[j][i] = 0;
                }
            }
        }
        gsl_spline_free (Bessel.first);
        gsl_interp_accel_free (Bessel.second);
        delete[] kernelMax;
        delete[] kernelIntegrate;
        delete[] kernelLim;
    }

}

void angularPowerSpectrumCurvedSky(double** angularPowerSpectra, double*** densityKernels, parameterFile paramFile,
                                    double* P_primordial, int* kernels, int* kernelInstruments, int** obsIndices,
                                    double* k, double* kChi, double* chi, double chiLow, double chiHigh,
                                    double* lensingSources, double** galaxyBins, int k_length, int a_length, 
                                    int l_max, int multiplier, int offset, int BesselArray_length) {

    for(int i = 0; i < ((1.*l_max-offset+1.)/multiplier); i++) {
        if(chiHigh<chi[1]) { // Check upper limit is bigger than lower limit
            chiHigh=chi[1];
        }

        // Kernel integrals
        double** angularPowerKernels;
        calcKernelIntegrals(paramFile, angularPowerKernels, densityKernels, kernels, k, kChi, chi, chiLow, chiHigh, 
                            lensingSources, galaxyBins, k_length, a_length, l_max, i*multiplier+offset, BesselArray_length);
        int totalObs = 0;
        for(int obsType=0; obsType<paramFile.numObs; obsType++) {
            if(paramFile.calcObs[obsType]) {
                // Galaxy auto-spectra
                if(kernelInstruments[obsType]==2) {
                    for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                        for(int g2=0; g2<paramFile.numGal; g2++) {  // column galaxies
                            if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                                if(g1==g2) {
                                    double* obs = new double[k_length];
                                    for(int j = 0; j<k_length; j++) {
                                        obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]+g1][j]*angularPowerKernels[obsIndices[obsType][1]+g1][j]*P_primordial[j]/(pow(2.*PI, 3.));
                                    }
                                    std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                                    double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                                    angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                                    totalObs++;
                                    delete[] obs;
                                    gsl_spline_free (obsFunction.first);
                                    gsl_interp_accel_free (obsFunction.second);
                                }
                            }
                            else {  // if calculate full covar
                                double* obs = new double[k_length];
                                for(int j = 0; j<k_length; j++) {
                                    obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]+g1][j]*angularPowerKernels[obsIndices[obsType][1]+g2][j]*P_primordial[j]/(pow(2.*PI, 3.));
                                }
                                std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                                double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                                angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                                totalObs++;
                                delete[] obs;
                                gsl_spline_free (obsFunction.first);
                                gsl_interp_accel_free (obsFunction.second);

                            }
                        }
                    }
                }
                // Lensing auto-spectra
                else if(kernelInstruments[obsType]==6) {
                    for(int g1=0; g1<paramFile.numLens; g1++) {  // row galaxies
                        for(int g2=0; g2<paramFile.numLens; g2++) {  // column galaxies
                            if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                                if(g1==g2) {
                                    double* obs = new double[k_length];
                                    for(int j = 0; j<k_length; j++) {
                                        obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]+g1][j]*angularPowerKernels[obsIndices[obsType][1]+g1][j]*P_primordial[j]/(pow(2.*PI, 3.));
                                    }
                                    std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                                    double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                                    angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                                    totalObs++;
                                    delete[] obs;
                                    gsl_spline_free (obsFunction.first);
                                    gsl_interp_accel_free (obsFunction.second);
                                }
                            }
                            else {  // if calculate full covar
                                double* obs = new double[k_length];
                                for(int j = 0; j<k_length; j++) {
                                    obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]+g1][j]*angularPowerKernels[obsIndices[obsType][1]+g2][j]*P_primordial[j]/(pow(2.*PI, 3.));
                                }
                                std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                                double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                                angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                                totalObs++;
                                delete[] obs;
                                gsl_spline_free (obsFunction.first);
                                gsl_interp_accel_free (obsFunction.second);

                            }
                        }
                    }
                }
                // Lensing-Galaxy cross-spectra
                else if(kernelInstruments[obsType]==4) {
                    for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                        for(int g2=0; g2<paramFile.numLens; g2++) {  // column lenses
                            double* obs = new double[k_length];
                            for(int j = 0; j<k_length; j++) {
                                obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]+g1][j]*angularPowerKernels[obsIndices[obsType][1]+g2][j]*P_primordial[j]/(pow(2.*PI, 3.));
                            }
                            std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                            double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                            angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                            totalObs++;
                            delete[] obs;
                            gsl_spline_free (obsFunction.first);
                            gsl_interp_accel_free (obsFunction.second);
                        }
                    }
                }
                // Galaxy cross-spectra (no lensing)
                else if(kernelInstruments[obsType]==1) {
                    for(int g=0; g<paramFile.numGal; g++) {
                        double* obs = new double[k_length];
                        for(int j = 0; j<k_length; j++) {
                            obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]][j]*angularPowerKernels[obsIndices[obsType][1]+g][j]*P_primordial[j]/(pow(2.*PI, 3.));
                        }
                        std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                        double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                        angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                        totalObs++;
                        delete[] obs;
                        gsl_spline_free (obsFunction.first);
                        gsl_interp_accel_free (obsFunction.second);
                    }
                }
                // Lensing cross-spectra (no galaxy)
                else if(kernelInstruments[obsType]==3) {
                    for(int g=0; g<paramFile.numLens; g++) {
                        double* obs = new double[k_length];
                        for(int j = 0; j<k_length; j++) {
                            obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]][j]*angularPowerKernels[obsIndices[obsType][1]+g][j]*P_primordial[j]/(pow(2.*PI, 3.));
                        }
                        std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                        double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                        angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                        totalObs++;
                        delete[] obs;
                        gsl_spline_free (obsFunction.first);
                        gsl_interp_accel_free (obsFunction.second);
                    }
                }
                // Non-galaxy, non-lensing spectra
                else {
                    double* obs = new double[k_length];
                    for(int j = 0; j<k_length; j++) {
                        obs[j] = k[j]*k[j]*k[j]*2.*PI*PI*angularPowerKernels[obsIndices[obsType][0]][j]*angularPowerKernels[obsIndices[obsType][1]][j]*P_primordial[j]/(pow(2.*PI, 3.));
                    }
                    std::pair<gsl_spline *, gsl_interp_accel *> obsFunction = spline1dSteffenData(k, obs, k_length);
                    double obsResult = gsl_spline_eval_integ(obsFunction.first, k[0], 0.2, obsFunction.second);
                    angularPowerSpectra[totalObs][multiplier*i+offset] = obsResult;
                    totalObs++;
                    delete[] obs;
                    gsl_spline_free (obsFunction.first);
                    gsl_interp_accel_free (obsFunction.second);
                }
            }
        }
        for(int j = 0; j < 2+(2*paramFile.numGal)+(2*paramFile.numLens); j++) {
            delete[] angularPowerKernels[j];
        }
        delete[] angularPowerKernels;
        std::cout<<"ell = "<<multiplier*i+offset<<" computed."<<std::endl;
    }
}




//////////////////////////
// LIMBER APPROXIMATION //
//////////////////////////


void matterPowerSpectrum(double** powerSpectrum, double* P_primordial, double* k, 
                                int a_length, int k_length) {
    for (int i = 0; i < a_length; i++) {
        for (int j = 0; j < k_length; j++) {
            powerSpectrum[j][i] = 2. * PI * PI * k[j] * P_primordial[j];
        }
    }
    return;
}

void ISWPowerSpectrumCorr(double** ISWCorrFunc, double* fOfK, double* HubbleRate, double** densityISW,
                            int a_length, int k_length) {
    for (int i = 0; i < a_length; i++) {
        for (int j = 0; j < k_length; j++) {
            ISWCorrFunc[j][i] = 8.*PI*fOfK[j]*HubbleRate[i]*densityISW[j][i]/(pow(SPEED_OF_LIGHT, 3.));
        }
    }
    return;
}

void galPowerSpectrumCorr(double** galCorrFunc, double* bias, double* selectionFunc, double** densityGal,
                            int a_length, int k_length) {
    for (int i = 0; i < a_length; i++) {
        for (int j = 0; j < k_length; j++) {
            galCorrFunc[j][i] = 4.*PI*bias[i]*selectionFunc[i]*densityGal[j][i];
        }
    }
    return;
}

void CMBLensePowerSpectrumCorr(double** CMBLenseCorrFunc, double* lenseKernel, double** densityLense, int a_length, 
                            int k_length) {
    double* invertedLenseKernel = invertArray(lenseKernel, a_length);
    for (int i = 0; i < a_length; i++) {
        for (int j = 0; j < k_length; j++) {
            CMBLenseCorrFunc[j][i] = 8.*PI*invertedLenseKernel[i]*densityLense[j][i];
        }
    }
    delete[] invertedLenseKernel;
    return;
}

void LensePowerSpectrumCorr(double** LenseCorrFunc, double* lenseKernel, double** density, double* fOfK, double* k,
                            int a_length, int k_length) {
    double* invertedLenseKernel = invertArray(lenseKernel, a_length);
    for (int i = 0; i < a_length; i++) {
        for (int j = 0; j < k_length; j++) {
            LenseCorrFunc[j][i] = 4.*PI*fOfK[j]*k[j]*k[j]*invertedLenseKernel[i]*density[j][i];
        }
    }
    delete[] invertedLenseKernel;
    return;
}

void GalLensePowerSpectrumCorr(double** GalLenseCorrFunc, double* galNumDen, double* bias, double** densityGal,
                                int a_length, int k_length) {
    double* invertedGalNumDen = invertArray(galNumDen, a_length);
    for (int i = 0; i < a_length; i++) {
        for (int j = 0; j < k_length; j++) {
            GalLenseCorrFunc[j][i] = 4.*PI*invertedGalNumDen[i]*bias[i]*densityGal[j][i];
        }
    }
    delete[] invertedGalNumDen;
    return;
}

void IntAlignPowerSpectrumCorr(double** IntAlignCorrFunc, double* IAKernel, double* numDen, double** density, 
                                int a_length, int k_length) {
    double* invIAKernel = invertArray(IAKernel, a_length);
    double* invNumDen = invertArray(numDen, a_length);
    for (int i = 0; i < a_length; i++) {
        for (int j = 0; j < k_length; j++) {
            IntAlignCorrFunc[j][i] = 4.*PI*invIAKernel[i]*invNumDen[i]*density[j][i];
        }
    }
    delete[] invIAKernel;
    delete[] invNumDen;
    return;
}

std::pair<gsl_spline *, gsl_interp_accel *>* powerSpectrum(double** matterPowerSpectrum, double** corr1, 
                                                        double** corr2, double* k, int a_length, int k_length) {

    std::pair<gsl_spline *, gsl_interp_accel *>* powerSpectrumArray = new std::pair<gsl_spline *, gsl_interp_accel *>[a_length];
    for (int i = 0; i < a_length; i++) {
        double* powerSpectrum = new double[k_length];
        for (int j = 0; j < k_length; j++) {
            powerSpectrum[j] = corr1[j][i]*corr2[j][i]*matterPowerSpectrum[j][i]/(pow(2.*PI, 3.));
        }
        powerSpectrumArray[i] = spline1dSteffenData(k, powerSpectrum, k_length);
        delete[] powerSpectrum;
    }
    return powerSpectrumArray;
}


// Get spline of Limber approximation function
std::pair<gsl_spline *, gsl_interp_accel *> limberFunction(std::pair<gsl_spline *, gsl_interp_accel *>* powerSpectrum, 
                                            std::pair<gsl_spline *, gsl_interp_accel *> Chi, double* HubbleRate,
                                            double* a, double* loga, int a_length, double kmax, int l, int lmax) {

    if(l<=lmax) {
        double * LimberFunc = new double[a_length-1];
        double C = PI/2.;
        //for (int i = 0; i<a_length; i++) {
        for (int i = 0; i<a_length-1; i++) {
            double chi = gsl_spline_eval(Chi.first, loga[i], Chi.second);
            double x = (l+0.5)/chi;
            if(x<kmax) {  // For high ell calculations. else case isn't calculated, this is precautionary.
                double Power = gsl_spline_eval(powerSpectrum[i].first, x, powerSpectrum[i].second);
                LimberFunc[i] = C * Power * SPEED_OF_LIGHT / (a[i] * a[i] * chi * chi * HubbleRate[i]);
            }
            else {
                LimberFunc[i] = 0;
            }
        }
        std::pair<gsl_spline *, gsl_interp_accel *> limber;
        limber = spline1dSteffenData(a, LimberFunc, a_length-1);
        delete[] LimberFunc;
        return limber;
    }
    else {
        std::pair<gsl_spline *, gsl_interp_accel *> voidSpline;
        voidSpline.first = nullptr;
        voidSpline.second = nullptr;
        return voidSpline;
    }
}


// function to be used in multithreading for limber solution
void threadFunctionLimber(double* C_ell, std::pair<gsl_spline *, gsl_interp_accel *>* Power, 
                    double* scale, double* loga, double* HubbleRate, int a_length, 
                    std::pair<gsl_spline *, gsl_interp_accel *> chi, int l_lim, int l_start, int l_max, int multiplier, 
                    int offset, int k_length, double aLow, double aHigh, double kmax) {
    
    for(int i = 0; i < ((1.*l_max-l_start-offset+1.)/multiplier); i++) {
        if(aHigh>loga[a_length-2]) {
            aHigh=loga[a_length-2];
        }
        if(((multiplier*i)+offset+l_start)<l_lim) {
            std::pair<gsl_spline *, gsl_interp_accel *> LimberFunction;
            LimberFunction = limberFunction(Power, chi, HubbleRate, scale, loga, a_length, kmax, multiplier*i+offset+l_start, l_max);
            double result = gsl_spline_eval_integ(LimberFunction.first, exp(aLow), exp(aHigh), LimberFunction.second);
            C_ell[multiplier*i+offset + l_start] = result;
            gsl_spline_free (LimberFunction.first);
            gsl_interp_accel_free (LimberFunction.second);
        }
        else {
            C_ell[multiplier*i+offset + l_start] = 0;
        }
    }
}



// Check time dependent part of integrateAjl() for D(a)
std::pair<double, double> BesselLimitChecker(double* a, int l, double k,  double aLow,
                                                std::pair<gsl_spline *, gsl_interp_accel *> chi,
                                                std::pair<gsl_spline *, gsl_interp_accel *> growth) {
    
    // similar integration process as chi
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000000);
    double result, error;
    std::pair<double, double> Result;
    gsl_function F;
    struct A_Bessel_params params = { chi.first, chi.second, growth.first, growth.second, l, k};
    F.function = &BesselToIntegrate;
    F.params = &params;

    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    //switch off default error handler, store old error handler in old_handler
    double relerr=1e-3;   //initial error tolerance (relative error)
    int status=1;
    while(status) {
        status=gsl_integration_qags(&F, aLow, 1, 0., relerr, 10000000, w, &result, &error);
        relerr *= 5.;
        if(status) std::cout<<"Increased tolerance="<<relerr<<std::endl;
        if(relerr >= 1) {
            break;
        }
    }
    relerr = 1e-3;
    while(status) {
        status=gsl_integration_qags(&F, aLow, 1, 0., relerr, 10000000, w, &result, &error);
        relerr *= 1./5.;
        if(status) std::cout<<"Decreased tolerance="<<relerr<<std::endl;
        if(relerr <= 1e-7) {
            gsl_integration_qags(&F, aLow, 1, 0., relerr, 10000000, w, &result, &error);
            break;
        }
    }
    gsl_set_error_handler(old_handler); //reset error handler

    Result.first = result;
    Result.second = error;
    gsl_integration_workspace_free (w);
    return Result;
}

double timeFunction(std::pair<gsl_spline *, gsl_interp_accel *> chi, std::pair<gsl_spline *, gsl_interp_accel *>  A, 
                    double k, int l, double a) {

    double Chi = gsl_spline_eval(chi.first, a, chi.second);
    double Bessel = gsl_sf_bessel_jl(l, Chi*k);
    double I = Bessel*gsl_spline_eval(A.first, a, A.second);
    return I;
}

double** timeIntegralFunction(std::pair<gsl_spline *, gsl_interp_accel *> chi, double* A, double k, double* a,
                                int a_length, int l) {
    
    std::pair<gsl_spline *, gsl_interp_accel *> splineFunction;
    splineFunction = spline1dSteffenData(a, A, a_length);
    double** timeFunc = new double*[2];
    timeFunc[0] = new double[1000];
    timeFunc[1] = new double[1000];
    for(int i = 0; i < 1000; i++) {
        timeFunc[0][i] = 0.05 + (i/(19980./19.));
        timeFunc[1][i] = timeFunction(chi, splineFunction, k, l, timeFunc[0][i]);
    }
    gsl_spline_free (splineFunction.first);
    gsl_interp_accel_free (splineFunction.second);
    return timeFunc;
}



/////////////////////////////
// BESSEL FUNCTION SPLINER //
/////////////////////////////

// Retrieve length of kChi array
int kChiLength() {
    double x = 0;
    double kChiMax = CUTOFF;
    int i = 0;
    while(x<kChiMax) {
        i++;
        x=x+0.001;
    }
    i++;
    return i;
}

// Makes chi*k array for bessel func
double* kChiArray(int kChi_length) {
    double* kChi = new double[kChi_length];
    for(int i = 0; i < kChi_length; i++) {
        kChi[i] = i/1000.;
    }
    return kChi;
}

// Makes bessel array for given ell, using kChi array
double* BesselArray(double* kChi, int ell, int BesselArray_length) {
    double* Bessel = new double[BesselArray_length];
    for(int i = 0; i < BesselArray_length; i++) {
        Bessel[i] = gsl_sf_bessel_jl(ell, kChi[i]);
    }
    return Bessel;
}

// checks if bessel file exists for given ell
bool BesselFileChecker(int ell) {
    // Finds Bessel file for ell value. If true, BesselSpline opens file, else generates Bessel file.
    std::string fileName; fileName = FILEPATH + "/BesselFunctions/Bessel_"+std::to_string(ell)+".dat";
    std::ifstream f(fileName.c_str());
    return f.good();
}

// makes bessel array file
void BesselFileMaker(double* Bessel, int ell, int BesselArray_length) {
    std::string fileName; fileName = FILEPATH + "/BesselFunctions/Bessel_"+std::to_string(ell)+".dat";
    std::ofstream file(fileName);
    for (int i = 0; i < BesselArray_length; i++) {
        file << Bessel[i] << std::endl;
    }
    file.close();
}

// produces spline of bessel function, either using file or making from scratch if file does not exist
std::pair<gsl_spline *, gsl_interp_accel *> BesselSpline(double* kChi, int ell, int BesselArray_length) {
    std::pair<gsl_spline *, gsl_interp_accel *> BesselSpline;
    if(BesselFileChecker(ell)) {
        // Read in
        double * Bessel = new double[BesselArray_length];
        std::string fileName; fileName = FILEPATH + "/BesselFunctions/Bessel_"+std::to_string(ell)+".dat";
        std::ifstream file(fileName);
        int i = 0;
        while(!file.eof()) {
            std::string line;
            file >> line;
            if(!line.empty()) {
                double dat = std::stod(line);
                Bessel[i] = dat;
            }
            if(i>BesselArray_length) {
                break;
            }
            i++;
        }
        file.close();
        if(i<BesselArray_length) {
            delete[] Bessel;
            std::cout<<"Bessel file not long enough. Remaking for more values."<<std::endl;
            // Generate Bessel array
            double* Bessel = BesselArray(kChi, ell, BesselArray_length);
            // Make file
            BesselFileMaker(Bessel, ell, BesselArray_length);
        }
        // Spline the data
        BesselSpline = spline1dSteffenData(kChi, Bessel, BesselArray_length);
        delete[] Bessel;
    }
    else {
        // Generate Bessel array
        double* Bessel = BesselArray(kChi, ell, BesselArray_length);
        // Make file
        BesselFileMaker(Bessel, ell, BesselArray_length);
        // Spline the data
        BesselSpline = spline1dSteffenData(kChi, Bessel, BesselArray_length);
        delete[] Bessel;
    }
    return BesselSpline;
}

void BesselSplineThreadFunction(std::pair<gsl_spline *, gsl_interp_accel *>* BesselArray, double* kChi,
                                int BesselArray_length, int l_max, int multiplier, int offset) {

    for(int i = 0; i < ((1.*l_max-offset+1.)/multiplier); i++) {
        std::cout<<multiplier*i+offset<<std::endl;
        BesselArray[multiplier*i+offset] = BesselSpline(kChi, multiplier*i+offset, BesselArray_length);
    }
}


#endif // !ISWINTEGRATOR_HPP
