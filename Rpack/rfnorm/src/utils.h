#ifndef __RFNORM_
#include <math.h>
#define __RFNORM_
#define _PI_ 3.14159265358979323846
#define _PI_T2_ 6.2831853071795853438
#define _SQRT_PI_T2_ 2.5066282746310002416
#define _INVSQRT_PI_T2_ 0.39894228040143270286

inline int rcpp_indicator(double x){return 1 - (x<0);}
inline double sq(double x) {return x*x;}
double density_double(double x, double u, double s);

/**
 * @param x value
 * @param u mean
 * @param s sigma
 * @return density
 */
double density_double(double x, double u, double s){
    return(
        _INVSQRT_PI_T2_/sqrt(s)*(exp(-sq(x - u)/(2.*s))+exp(-sq(x + u)/(2. * s))*rcpp_indicator(x))
    );
}
#endif
