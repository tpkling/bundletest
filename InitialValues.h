//
// Created by Thomas Kling on 8/30/18.
//

#ifndef KERRSINGLE_INITIALVALUES_H
#define KERRSINGLE_INITIALVALUES_H


#include "Coords.h"

class InitialValues {
public:
    double termA_, termB_, termC_, termD_, termE_, termF_, termG_;
    double g00_, g33_, g03_;
    double t1_, t2_, t3_, t4_, t5_, t6_;
    double a_, m_, delta_, tDot_, rhosq_, direction_;
    void computeValues(Coords coords);
    InitialValues(double m, double a, int direction) : m_(m), a_(a), direction_(direction) {};
};



#endif //KERRSINGLE_INITIALVALUES_H
