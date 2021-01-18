//
// Created by Thomas Kling on 8/30/18.
//

#include "InitialValues.h"
#include "KerrRay.h"
#include <math.h>


void InitialValues::computeValues(Coords coords) {
    t1_ = KerrRay::T1func(coords.r_, coords.theta_, m_, a_);
    t2_ = KerrRay::T2func(coords.r_, coords.theta_, m_, a_, direction_);
    t3_ = KerrRay::T3func(coords.r_, coords.theta_, m_, a_, direction_);
    t4_ = KerrRay::T4func(coords.r_, coords.theta_, m_, a_);
    t5_ = KerrRay::T5func(coords.r_, coords.theta_, m_, a_);
    t6_ = KerrRay::T6func(coords.r_, coords.theta_, m_, a_);

    tDot_ = 1.0*pow(-1.0, direction_); // changed for directon = 1 = past, 0 = future

    delta_ = coords.r_ * coords.r_ - 2 * m_ * coords.r_ + a_ * a_;
    rhosq_ = coords.r_ * coords.r_ + a_ * a_ * cos(coords.theta_)*  cos(coords.theta_);

    termD_ = (coords.r_ * coords.r_ + a_ * a_) / delta_;

    termA_ = (t3_ * termD_ + t6_ * termD_ * termD_);
    termB_ = (t2_ * coords.vphi_ + t3_ * tDot_ + t5_ * termD_ * coords.vphi_ + 2 * t6_ * termD_ * tDot_);
    termC_ = (t1_ * coords.vtheta_ * coords.vtheta_ + t4_ * coords.vphi_ * coords.vphi_ +
              t5_ * coords.vphi_ * tDot_ + t6_ * tDot_ * tDot_);
    termE_ = (t2_ * t2_ + t5_ * t5_ * termD_ * termD_ + 2 * t2_ * t5_ * termD_) / (4 * termA_) - t4_;
    termF_ = (2 * t2_ * t3_ * tDot_ + 4 * t2_ * t6_ * termD_ * tDot_ +  2 * t3_ * t5_ * termD_ * tDot_ +
              4 * t5_ * t6_ * termD_ * termD_ * tDot_) / (4 * termA_) - t5_ * tDot_;
    termG_ = (t3_ * t3_ * tDot_ * tDot_ + 4 * t6_ * t6_ * termD_ * termD_ * tDot_ * tDot_ +
              4 * t3_ * t6_ * termD_ * tDot_ * tDot_) / (4 * termA_) - t6_ * tDot_ * tDot_;

    g00_ = - 1.0 + 2.0 * m_ * coords.r_ / rhosq_;
    g03_ = - 2.0 * m_ * a_ * coords.r_ * sin(coords.theta_) * sin(coords.theta_) / rhosq_;
    g33_ = (coords.r_ * coords.r_ + a_*a_ ) * sin(coords.theta_) * sin(coords.theta_)
           + 2.0 * m_ * a_ * a_ * coords.r_ * pow(sin(coords.theta_), 4.0) / rhosq_;

}