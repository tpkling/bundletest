//
// Created by Thomas Kling on 8/30/18.
//

#include "KerrRay.h"
#include <cmath>

double KerrRay::T1r(double r, double theta, double m, double a){  // Confirmed to work
    return r; // cleaned up by tpk on July 22
}

double KerrRay::T1theta(double r, double theta, double m, double a){  // Confirmed to work
    return -a * a * sin(theta) * cos(theta); // cleaned up by tpk on July 22
}

double KerrRay::T2r(double r, double theta, double m, double a, int dir){  // Confirmed to work
    return 0; // cleaned up by tpk on July 22;
}

double KerrRay::T2theta(double r, double theta, double m, double a, int dir){  // Confirmed to work
    return -2.0 * pow(-1.0, dir) * a * sin(theta) * cos(theta); // cleaned up by tpk on July 22
}

double KerrRay::T3r(double r, double theta, double m, double a, int dir){  // Confirmed to work
    return 0; // cleaned up by tpk on July 22
}

double KerrRay::T3theta(double r, double theta, double m, double a, int dir){  // Confirmed to work
    return 0; // cleaned up by tpk on July 22
}

double KerrRay::T4r(double r, double theta, double m, double a){  // Confirmed to work
    double rhosq = rhosqf(r, a, theta);

    return ( (cos(theta) * cos(theta) * (2.0 * a * a * r * (a * a + r * r) + a * a * a * a * (m - r) * sin(theta) * sin(theta)))
             + (r * (-a * a * a * a + r * r * r * r + a * a * (a * a - m * r) * sin(theta) * sin(theta))) ) *
           sin(theta) * sin(theta) / rhosq / rhosq; // cleaned up by tpk on July 22 - confirmed in mathematica file.
}

double KerrRay::T4theta(double r, double theta, double m, double a){  // Confirmed to work
    double rhosq = rhosqf(r, a, theta);
    double delta = deltaf(r, m, a);

    return (-a * a * delta * cos(theta) * sin(theta) * sin(theta) * sin(theta) / rhosq) +
           (cos(theta) * sin(theta) * ((a * a + r * r) * (a * a + r * r) - a * a * delta * sin(theta) * sin(theta)) / rhosq) +
           (a * a * cos(theta) * sin(theta) * sin(theta) * sin(theta) * ((a * a + r * r) * (a * a + r * r) - a * a * delta * sin(theta) * sin(theta)) / rhosq / rhosq);
    // cleaned up by tpk on July 22 -- confirmed in mathematica file
}

double KerrRay::T5r(double r, double theta, double m, double a){  // Confirmed to work
    double rhosq = rhosqf(r, a, theta);

    return (4.0 * a * m * r * r * sin(theta) * sin(theta) / rhosq / rhosq) +
           (-2.0 * a * m * sin(theta) * sin(theta) / rhosq);
}

double KerrRay::T5theta(double r, double theta, double m, double a){  // Confirmed to work
    double rhosq = rhosqf(r, a, theta);

    return (-4.0 * a * m * r * cos(theta) * sin(theta) / rhosq) +
           (-4.0 * a * a * a * m * r * cos(theta) * sin(theta) * sin(theta) * sin(theta) / rhosq / rhosq);
}

double KerrRay::T6r(double r, double theta, double m, double a){  // Confirmed to work
    double rhosq = rhosqf(r, a, theta);

    //return 4.0 * m * r * r / rhosq / rhosq - 2 * m / rhosq;
    return m / rhosq - 2.0 * m * r * r / (rhosq * rhosq);
}

double KerrRay::T6theta(double r, double theta, double m, double a){  // Confirmed to work
    double rhosq = rhosqf(r, a, theta);

    return 2.0 * a * a * m * r * cos(theta) * sin(theta) / (rhosq * rhosq);
}
