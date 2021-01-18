//
// Created by Thomas Kling on 8/30/18.
// Modified for bundle testing in Sept / Oct 2020
// Have cut out many of the functions from the original KerrRay Class.
//

#ifndef KERRSINGLE_KERRRAY_H
#define KERRSINGLE_KERRRAY_H

#include <fstream>
#include "Coords.h"
#include "InitialValues.h"

class KerrRay {

public:
    // KerrRay
    KerrRay(double vphi, double vtheta, double yfin, double zfin); // constructor when considering a group of rays in seeking minima
    KerrRay(double vphi, double vtheta); // constructor for just initial directions
    
    ~KerrRay();

    Coords coords_;
    double yfin_, zfin_;
    // Mandatory Element functions
    void step(int direction);

    static void computeNewCoordinates(Coords initCoords, Coords &newCoords, double err[], double h,
            double m, double a, int direction);

    double getError();

    // KerrRayVDots
    static double Qu(Coords coords, double m, double a);
    static double Qphi(Coords coords, double m, double a, int direction);
    static double Qr(Coords coords, double m, double a, int direction);
    static double R1(Coords coords, double m, double a, int direction);
    static double R2(Coords coords, double m, double a, int direction);
    static double VPhiDot(Coords coords, double m, double a, int direction);
    static double VUDot(Coords coords, double m, double a, int direction);
    static double VRDot(Coords coords, double m, double a, int direction);
    static double VThetaDot(Coords coords, double m, double a, int direction);
    static double VPhiTrue(Coords coords, double m, double a, int direction);
    static double lengthDot(Coords coords, double m, double a, int direction);
    static double tDot(Coords coords, double m, double a, int direction);
    // KerrRayTFuncs
    static double T1func(double r, double theta, double m, double a);
    static double T2func(double r, double theta, double m, double a, int dir);
    static double T3func(double r, double theta, double m, double a, int dir);
    static double T4func(double r, double theta, double m, double a);
    static double T5func(double r, double theta, double m, double a);
    static double T6func(double r, double theta, double m, double a);
    static double rhosqf(double r, double a, double theta);
    static double deltaf(double r, double m, double a);
    // KerrRayTDerivs
    static double T1r(double r, double theta, double m, double a);
    static double T1theta(double r, double theta, double m, double a);
    static double T2r(double r, double theta, double m, double a, int dir);
    static double T2theta(double r, double theta, double m, double a, int dir);
    static double T3r(double r, double theta, double m, double a, int dir);
    static double T3theta(double r, double theta, double m, double a, int dir);
    static double T4r(double r, double theta, double m, double a);
    static double T4theta(double r, double theta, double m, double a);
    static double T5r(double r, double theta, double m, double a);
    static double T5theta(double r, double theta, double m, double a);
    static double T6r(double r, double theta, double m, double a);
    static double T6theta(double r, double theta, double m, double a);

    double vThetaLimit(Coords coords, InitialValues iv);
    double posVPhiLimit(Coords coords, InitialValues iv);
    double negVPhiLimit(Coords coords, InitialValues iv);
    double initVR(Coords coords, InitialValues iv);
    
    double dfunc();

//    void runray(int direction, double timeStop);
    void runtox(int direction, double xstop);
    void runScreen(int direction);
//    int runset(int direction, int completed, double timeStop);

private:

    double rplus_, a_, m_, h_, totalError_;
    int goodStepCount_, badStepCount_;
    bool distOrTime_, phiPlusOrPhiTrue_;
   

};


#endif //KERRSINGLE_KERRRAY_H
