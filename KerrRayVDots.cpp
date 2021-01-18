//
// Created by Thomas Kling on 8/30/18.
//

#include <cmath>
#include "KerrRay.h"

double KerrRay::Qu(Coords coords, double m, double a){
    double T5rval, T5thval, T6rval, T6thval;

    T5rval  = T5r(coords.r_, coords.theta_, m, a);
    T5thval = T5theta(coords.r_, coords.theta_, m, a);
    T6rval  = T6r(coords.r_, coords.theta_, m, a);
    T6thval = T6theta(coords.r_, coords.theta_, m, a);

    return (-T5rval * coords.vr_ * coords.vphi_ - T5thval * coords.vtheta_ * coords.vphi_ -
            2.0 * T6rval * coords.vr_ * coords.vu_ - 2.0 * T6thval * coords.vtheta_ * coords.vu_);
}

double KerrRay::Qphi(Coords coords, double m, double a, int dir){
    double T4rval, T4thval, T5rval, T5thval, T2thval;

    T4rval  = T4r(coords.r_, coords.theta_, m, a);
    T4thval = T4theta(coords.r_, coords.theta_, m, a);
    T5rval  = T5r(coords.r_, coords.theta_, m, a);
    T5thval = T5theta(coords.r_, coords.theta_, m, a);
    T2thval = T2theta(coords.r_, coords.theta_, m, a, dir);

    return (-2.0 * T4rval * coords.vr_ * coords.vphi_ - 2.0 * T4thval * coords.vtheta_ * coords.vphi_ - T5rval * coords.vr_ * coords.vu_ -
            T5thval * coords.vtheta_ * coords.vu_ - T2thval * coords.vr_ * coords.vtheta_);
}

double KerrRay::Qr(Coords coords, double m, double a, int dir){
    double T1rval, T4rval, T5rval, T6rval, T2thval;

    T1rval  = T1r(coords.r_, coords.theta_, m, a);
    T4rval  = T4r(coords.r_, coords.theta_, m, a);
    T5rval  = T5r(coords.r_, coords.theta_, m, a);
    T6rval  = T6r(coords.r_, coords.theta_, m, a);
    T2thval = T2theta(coords.r_, coords.theta_, m, a, dir);

    return (-T2thval * coords.vtheta_ * coords.vphi_ + T1rval * coords.vtheta_ * coords.vtheta_ + T4rval * coords.vphi_ * coords.vphi_ +
            T5rval * coords.vphi_ * coords.vu_ + T6rval * coords.vu_ * coords.vu_);
}

double KerrRay::R1(Coords coords, double m, double a, int dir){
    double T3val, T2val, T5val, T6val, T4val;

    T2val  = T2func(coords.r_, coords.theta_, m, a, dir);
    T3val  = T3func(coords.r_, coords.theta_, m, a, dir);
    T4val  = T4func(coords.r_, coords.theta_, m, a);
    T5val  = T5func(coords.r_, coords.theta_, m, a);
    T6val  = T6func(coords.r_, coords.theta_, m, a);

    return 2.0 * (T2val * T2val * T6val/ T3val / T3val  + T4val - T2val * T5val / T3val);
}

double KerrRay::R2(Coords coords, double m, double a, int dir){
    double T2val, T3val, T5val, T6val, Qphival, Qrval, Quval;

    T2val  = T2func(coords.r_, coords.theta_, m, a, dir);
    T3val  = T3func(coords.r_, coords.theta_, m, a, dir);
    T5val  = T5func(coords.r_, coords.theta_, m, a);
    T6val  = T6func(coords.r_, coords.theta_, m, a);
    Qphival = Qphi(coords, m, a, dir);
    Quval   = Qu(coords, m, a);
    Qrval   = Qr(coords, m, a, dir);

    return Qphival - T2val / T3val * Quval + 2.0 * T2val * T6val /T3val / T3val * Qrval - T5val /T3val * Qrval;
}

double KerrRay::VPhiDot(Coords coords, double m, double a, int dir){
    double R1val, R2val;

    R1val = R1(coords, m, a, dir);
    R2val = R2(coords, m, a, dir);

    return R2val / R1val;
}

double KerrRay::VUDot(Coords coords, double m, double a, int dir){
    double T2val, T3val, Qrval, VPhiDotVal;

    T2val  = T2func(coords.r_, coords.theta_, m, a, dir);
    T3val = T3func(coords.r_, coords.theta_, m, a, dir);
    Qrval   = Qr(coords, m, a, dir);
    VPhiDotVal = VPhiDot(coords, m, a, dir);

    return (Qrval - T2val * VPhiDotVal)/T3val;
}

double KerrRay::VRDot(Coords coords, double m, double a, int dir){
    double T3val, T5val, T6val, Quval, VUDotVal, VPhiDotVal;

    T3val = T3func(coords.r_, coords.theta_, m, a, dir);
    T5val  = T5func(coords.r_, coords.theta_, m, a);
    T6val  = T6func(coords.r_, coords.theta_, m, a);
    Quval  = Qu(coords, m, a);
    VUDotVal = VUDot(coords, m, a, dir);
    VPhiDotVal = VPhiDot(coords, m, a, dir);

    return (Quval - T5val * VPhiDotVal - 2 * T6val * VUDotVal)/T3val;
}

double KerrRay::VThetaDot(Coords coords, double m, double a, int dir){
    double T1val, T1rval, T1thval, T2thval, T4thval, T5thval, T6thval;

    T1val    = T1func(coords.r_, coords.theta_, m, a);
    T1rval   = T1r(coords.r_, coords.theta_, m, a);
    T1thval  = T1theta(coords.r_, coords.theta_, m, a);
    T2thval  = T2theta(coords.r_, coords.theta_, m, a, dir);
    T4thval  = T4theta(coords.r_, coords.theta_, m, a);
    T5thval  = T5theta(coords.r_, coords.theta_, m, a);
    T6thval  = T6theta(coords.r_, coords.theta_, m, a);

    return (-T1thval * coords.vtheta_ * coords.vtheta_ + T2thval * coords.vr_ * coords.vphi_ + T4thval * coords.vphi_ * coords.vphi_ +
            T5thval * coords.vphi_ * coords.vu_ + T6thval * coords.vu_ * coords.vu_ - 2 * T1rval * coords.vr_ * coords.vtheta_) / (2 * T1val);
}

double KerrRay::VPhiTrue(Coords coords, double m, double a, int dir) {
    return coords.vphi_ - pow(-1.0 , dir) * (a / deltaf(coords.r_, m, a)) * coords.vr_ ;
}

double KerrRay::lengthDot(Coords coords, double m, double a, int dir) {
    double vrTerm = rhosqf(coords.r_, a, coords.theta_) / deltaf(coords.r_, m, a);
    double vthetaTerm = rhosqf(coords.r_, a, coords.theta_);

    double vphiTerm = (coords.r_ * coords.r_ + a * a) * sin(coords.theta_) * sin(coords.theta_) +
                      (2.0 * m * coords.r_ * a * a * pow(sin(coords.theta_), 4.0))/rhosqf(coords.r_, a, coords.theta_) +
                      (4.0 * m * m * coords.r_ * coords.r_ * a * a * pow(sin(coords.theta_), 4.0)) /
                      (rhosqf(coords.r_, a, coords.theta_) * (rhosqf(coords.r_, a, coords.theta_) - 2.0 * m * coords.r_));

    double vphiTrue = VPhiTrue(coords, m, a, dir);

    return sqrt(vrTerm * coords.vr_ * coords.vr_ + vthetaTerm * coords.vtheta_ * coords.vtheta_ + vphiTerm * vphiTrue * vphiTrue);
}

double KerrRay::tDot(Coords coords, double m, double a, int dir){

    return coords.vu_ - pow(-1.0, dir) * (coords.r_*coords.r_ + a*a)/(coords.r_ * coords.r_ - 2.0*m*coords.r_ + a*a) * coords.vr_;

} // returns tdot from Hawking and Ellis