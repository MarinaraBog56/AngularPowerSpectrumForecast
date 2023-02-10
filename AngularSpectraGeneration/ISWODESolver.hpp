#pragma once

#ifndef ISWODESOLVER_HPP
#define ISWODESOLVER_HPP


// hpp file for integrating ISW Power spectrum for C_ells
// Elizabeth Brown

#include "stdafx.h"


struct spline_params {gsl_spline * spline_ptr; gsl_interp_accel * accel_ptr;};
struct splineBessel_params {gsl_spline * spline_ptr; gsl_interp_accel * spline_accel_ptr; 
                        gsl_spline * bessel_ptr; gsl_interp_accel * bessel_accel_ptr; double k;};

std::pair<double,double> RK4(double init, double fin, double step, double q1Init, double q2Init,
                                std::pair<gsl_spline *, gsl_interp_accel *> function, 
                                std::pair<gsl_spline *, gsl_interp_accel *> Bessel, double k);
void RK4Alg(double x, double q1, double q2, double step, double* q1Arr, double* q2Arr, int index,
            std::pair<gsl_spline *, gsl_interp_accel *> function, std::pair<gsl_spline *, gsl_interp_accel *> Bessel, 
            double k);
double Afunction(double q1, double q2, double x, std::pair<gsl_spline *, gsl_interp_accel *> function, 
                    std::pair<gsl_spline *, gsl_interp_accel *> Bessel, double k);
double Bfunction(double q1, double q2, double x);


// fi(t, y0, y1, ..., yn) = dyi/dt
int splineFunc(double t, const double y[], double f[], void *params) {
    struct spline_params * p = (struct spline_params *)params;
    gsl_spline * spline = (p->spline_ptr);
    gsl_interp_accel * accelerator = (p->accel_ptr);
    //double mu = *(double *)params;
    f[0] = gsl_spline_eval(spline, t, accelerator);
    return GSL_SUCCESS;
}

// J=(df0/dy0, df0/dy1, ...)
//   (df1/dy0, df1/dy1, ...)
//   (  ...      ...    ...)
// Plus dfi/dt vector
int splineJac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    //(void)(t); /* avoid unused parameter warning */
    struct spline_params * p = (struct spline_params *)params;
    gsl_spline * spline = (p->spline_ptr);
    gsl_interp_accel * accelerator = (p->accel_ptr);
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 1, 1);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    dfdt[0] = gsl_spline_eval_deriv(spline, t, accelerator);
    return GSL_SUCCESS;
}

std::pair<double, double> solveSplineODE(double t0, double tf, double y0, 
                                        std::pair<gsl_spline *, gsl_interp_accel *> function) {
    struct spline_params params = { function.first, function.second };
    //                dyi/dx func, Jacobian, dimensions, parameters
    gsl_odeiv2_system sys = {splineFunc, splineJac, 1, &params};
    //                                                 System, Stepper type, init step val, abserr, relerr
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 0.0, 1e-6);
    // initial conditions
    double y[1] = { y0 };
    int status = gsl_odeiv2_driver_apply (d, &t0, tf, y);

    if (status != GSL_SUCCESS) {
        printf ("error, return value=%d\n", status);
    }
    gsl_odeiv2_driver_free (d);
    std::pair<double, double> result;
    result.first=t0;
    result.second=y[0];
    return result;
}


// fi(t, y0, y1, ..., yn) = dyi/dt
int splineBesselFunc(double t, const double y[], double f[], void *params) {
    struct splineBessel_params * p = (struct splineBessel_params *)params;
    gsl_spline * spline = (p->spline_ptr);
    gsl_interp_accel * splineAccelerator = (p->spline_accel_ptr);
    gsl_spline * bessel = (p->bessel_ptr);
    gsl_interp_accel * besselAccelerator = (p->bessel_accel_ptr);
    double k = (p->k);
    //double mu = *(double *)params;
    double bess;
    if(k*t<CUTOFF) {
        bess = gsl_spline_eval (bessel, k*t, besselAccelerator);
    }
    else {bess=0;}
    f[0] = gsl_spline_eval(spline, t, splineAccelerator)*bess;
    return GSL_SUCCESS;
}

// J=(df0/dy0, df0/dy1, ...)
//   (df1/dy0, df1/dy1, ...)
//   (  ...      ...    ...)
// Plus dfi/dt vector
int splineBesselJac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    //(void)(t); /* avoid unused parameter warning */
    struct splineBessel_params * p = (struct splineBessel_params *)params;
    gsl_spline * spline = (p->spline_ptr);
    gsl_interp_accel * splineAccelerator = (p->spline_accel_ptr);
    gsl_spline * bessel = (p->bessel_ptr);
    gsl_interp_accel * besselAccelerator = (p->bessel_accel_ptr);
    double k = (p->k);
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 1, 1);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    double bess, bessDeriv;
    if(k*t<CUTOFF) {
        bess = gsl_spline_eval (bessel, k*t, besselAccelerator);
        bessDeriv = gsl_spline_eval_deriv(bessel, k*t, besselAccelerator);
    }
    else {bess=0; bessDeriv=0;}
    double splineDeriv = bess*gsl_spline_eval_deriv(spline, t, splineAccelerator);
    double besselDeriv = k*bessDeriv*gsl_spline_eval(spline, t, splineAccelerator);
    dfdt[0] = splineDeriv+besselDeriv;
    return GSL_SUCCESS;
}

std::pair<double, double> solveSplineBesselODE(double t0, double tf, double step, std::pair<gsl_spline *, gsl_interp_accel *> function, 
                                                std::pair<gsl_spline *, gsl_interp_accel *> Bessel, double k) {
    struct splineBessel_params params = { function.first, function.second , Bessel.first, Bessel.second, k};
    //                dyi/dx func, Jacobian, dimensions, parameters
    gsl_odeiv2_system sys = {splineBesselFunc, splineBesselJac, 1, &params};
    //                                                 System, Stepper type, init step val, abserr, relerr
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, step, 0.0, 1e-6);
    // initial conditions
    double y0 = gsl_spline_eval(function.first, t0, function.second)*gsl_spline_eval(Bessel.first, k*t0, Bessel.second);
    double y[1] = { y0 };
    int status = gsl_odeiv2_driver_apply (d, &t0, tf, y);

    if (status != GSL_SUCCESS) {
        printf ("error, return value=%d\n", status);
    }
    gsl_odeiv2_driver_free (d);
    std::pair<double, double> result;
    result.first=t0;
    result.second=y[0];
    return result;
}



// Runge-Kutta 4th order solver. Not used above.
std::pair<double,double> RK4(double init, double fin, double step, double q1Init, double q2Init, 
                                std::pair<gsl_spline *, gsl_interp_accel *> function, 
                                std::pair<gsl_spline *, gsl_interp_accel *> Bessel, double k) {
    int length = (int)round(1 + ((fin-init)/step));
    if(length<=0){
        std::cout<<"Step size not small enough. Returning 0."<<std::endl;
        std::pair<double,double> nullResult;
        nullResult.first=0; nullResult.second=0;
        return nullResult;
    }
	double* q1Arr = new double[length];
    double* q2Arr = new double[length];
	double* xArr = new double[length];
	double x = init;
	double q1 = q1Init;
    double q2 = q2Init;
    int i = 0;
    q1Arr[i] = q1;
    q2Arr[i] = q2;
    xArr[i] = init;
        i++;
    while(i<length) {
        RK4Alg(x, q1Arr[i-1], q2Arr[i-1], step, q1Arr, q2Arr, i, function, Bessel, k);
        x = x+step;
        xArr[i] = x;
        i++;
    }
    std::pair<double,double> result;
    result.first=q1Arr[length-1]; result.second=q2Arr[length-1];
    delete[] q1Arr; delete[] q2Arr; delete[] xArr;
    return result;
}

void RK4Alg(double x, double q1, double q2, double step, double* q1Arr, double* q2Arr, int index,
            std::pair<gsl_spline *, gsl_interp_accel *> function, std::pair<gsl_spline *, gsl_interp_accel *> Bessel, 
            double k) {
	double q10 = q1;
	double q20 = q2;
	double x1 = x+(step/2.);
	double x2 = x+step;
	// k^(0)
	double K11 = Afunction(q10, q20, x, function, Bessel, k);
	double q11 = q10 + (K11 * step/2.);
	double K21 = Bfunction(q10, q20, x);
	double q21 = q20 + (K21 * step/2.);
	// k^(2)
	double K12 = Afunction(q11, q21, x1, function, Bessel, k);
	double q12 = q10 + (K12 * step/2.);
	double K22 = Bfunction(q11, q21, x1);
	double q22 = q20 + (K22 * step/2.);
	// k^(3)
	double K13 = Afunction(q12, q22, x1, function, Bessel, k);
	double q13 = q10 + (K13 * step);
	double K23 = Bfunction(q12, q22, x1);
	double q23 = q20 + (K23 * step);
	// k^(4)
	double K14 = Afunction(q13, q23, x2, function, Bessel, k);
	double q14 = q10 + (K14 * step);
	double K24 = Bfunction(q13, q23, x2);
	double q24 = q20 + (K24 * step);
	// Estimates
	q1Arr[index] = q1 + (K11 + 2*K12 + 2*K13 + K14)*(step/6.);
	q2Arr[index] = q2 + (K21 + 2*K22 + 2*K23 + K24)*(step/6.);
	return;
}

double Afunction(double q1, double q2, double x, std::pair<gsl_spline *, gsl_interp_accel *> function, 
                    std::pair<gsl_spline *, gsl_interp_accel *> Bessel, double k) {
    // sort parameters here
    //struct ode_params params = p;

    // calculate numerical value of function
    double func = gsl_spline_eval (function.first, x, function.second);
    double bess;
    if(k*x<CUTOFF) {
        bess = gsl_spline_eval (Bessel.first, k*x, Bessel.second);
    }
    else {bess=0;}
    return func*bess;
}

double Bfunction(double q1, double q2, double x) {
    return 1;
}

std::pair<gsl_spline *, gsl_interp_accel *> invertSpline(std::pair<gsl_spline *, gsl_interp_accel *> spline, 
                                                            double* x, int length) {
    double* y = new double[length];
    for (int i = 0; i < length; i++) {
        y[i] = gsl_spline_eval (spline.first, x[i], spline.second);
    }
    gsl_interp_accel * my_accel_ptr;
    gsl_spline * spline_ptr;
    my_accel_ptr = gsl_interp_accel_alloc();
    spline_ptr = gsl_spline_alloc(gsl_interp_steffen, length);
    gsl_spline_init (spline_ptr, y, x, length);
    std::pair<gsl_spline*, gsl_interp_accel*> invertedSpline;
    invertedSpline.first = spline_ptr;
    invertedSpline.second = my_accel_ptr;
    delete[] y;
    return invertedSpline;
}

double* retrieveSplineArray(std::pair<gsl_spline *, gsl_interp_accel *> spline, double* x, int length) {
    double* y = new double[length];
    for (int i = 0; i < length; i++) {
        y[i] = gsl_spline_eval(spline.first, x[i], spline.second);
    }
    return y;
}

double* invertArray(double* array, int length) {
    double* invertedArray = new double[length];
    for (int i = 0; i < length; i++) {
        invertedArray[length - 1 - i] = array[i];
    }
    return invertedArray;
}


#endif // !ISWODESOLVER_HPP