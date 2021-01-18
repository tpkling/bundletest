//
// Created by Thomas Kling on 8/30/18.
//

#include "KerrRay.h"
#include <cmath>

double KerrRay::T1func(double r, double theta, double m, double a){ // Confirmed to work and have appropriate value
    return 0.5 * rhosqf(r, a, theta);
}

double KerrRay::T2func(double r, double theta, double m, double a, int dir){ // Confirmed to work and have appropriate value
    return -pow(-1.0, dir) * a * sin(theta) * sin(theta);
}

double KerrRay::T3func(double r, double theta, double m, double a, int dir){ // Confirmed to work and have appropriate value
    return 1.0*pow(-1.0, dir);
}

double KerrRay::T4func(double r, double theta, double m, double a){ // Confirmed to work and have appropriate value
    return 0.5 * sin(theta) * sin(theta) * ((r * r + a * a) * (r * r + a * a) -
                                            deltaf(r,m,a) * a * a * sin(theta) * sin(theta)) / rhosqf(r, a, theta);
}

double KerrRay::T5func(double r, double theta, double m, double a){ // Confirmed to work and have appropriate value
    return -2.0 * a * m * r * sin(theta) * sin(theta) / rhosqf(r, a, theta);
}

double KerrRay::T6func(double r, double theta, double m, double a){ // Confirmed to work and have appropriate value
    return -0.5 + m * r / rhosqf(r, a, theta);
}

// Rho Squared function
double KerrRay::rhosqf(double r, double a, double theta) {
    return (r * r + a * a * cos(theta) * cos(theta));
}

double KerrRay::deltaf(double r, double m, double a) {
    return (r * r - 2.0 * m * r + a * a);
}