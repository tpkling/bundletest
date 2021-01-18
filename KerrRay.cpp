//
// Created by Thomas Kling on 8/30/18.
// Modified for testing bundle only in September / October 2020
//

#include "KerrRay.h"
#include "InitialValues.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <math.h>

using namespace std;

// Begin constants for cash-karp method of RKF
// a2 thru a6 unused??

const double a2  = 1.0/5.0;
const double a3  = 3.0/10.0;
const double a4  = 3.0/5.0;
const double a5  = 1.0;
const double a6  = 7.0/8.0;

const double b21 = 1.0/5.0;
const double b31 = 3.0/40.0;
const double b41 = 3.0/10.0;
const double b51 = -11.0/54.0;
const double b61 = 1631.0/55296.0;

const double b32 = 9.0/40.0;
const double b42 = -9.0/10.0;
const double b52 = 5.0/2.0;
const double b62 = 175.0/512.0;

const double b43 = 6.0/5.0;
const double b53 = -70.0/27.0;
const double b63 = 575.0/13824.0;

const double b54 = 35.0/27.0;
const double b64 = 44275.0/110592.0;

const double b65 = 253.0/4096.0;

const double c1  = 37.0/378.0;
const double c2  = 0.0;
const double c3  = 250.0/621.0;
const double c4  = 125.0/594.0;
const double c5  = 0.0;
const double c6  = 512.0/1771.0;

const double c1s = 2825.0/27648.0;
const double c2s = 0.0;
const double c3s = 18575.0/48384.0;
const double c4s = 13525.0/55296.0;
const double c5s = 277.0/14336.0;
const double c6s = 1.0/4.0;
// End constants for cash-karp method of RKF

const double EPS = 1.0e-10; // values were about 1e-3
const double SAFETY = 0.9;

const double PI = acos(-1);

KerrRay::KerrRay(double vphi_in, double vtheta_in,  double yfin, double zfin) {

    // version of the consrtuctor where vphi and vtheta values are central to the bundle, and we want to modify them somehow
    // constructor being used for series of rays where we are trying to find a ray that comes close to yfin, zfin
    
    a_ = 0.9;
    m_ = 1.0;
    h_ = 0.0010;//starting h value
    yfin_ = yfin;
    zfin_ = zfin; // passed in, desired final location
    int direction = +1;
    totalError_ = 0.0;
    goodStepCount_ = 0;
    badStepCount_ = 0;

    Coords initCoords;
    InitialValues initialValues(m_, a_, direction);

    rplus_ = m_ + sqrt(m_ * m_ - a_ * a_);  // we think this locates the event horizon

    initCoords.alpha_ = 0; // don't think we need this for this constructor
    initCoords.r_ = 20.0;
    initCoords.theta_ = PI/2.0;
    initCoords.phi_ = PI;
    initCoords.u_ = 10.0;
    initCoords.lambda_ = 0.0;
    initCoords.length_ = 0.0; // for spatial distance
    initCoords.t_      = 0.0; // for time
    initCoords.phitrue_ = initCoords.phi_;
    double initDelta = initCoords.r_ * initCoords.r_ - 2.0 * m_ * initCoords.r_ + a_ * a_;
    initCoords.vt_ = 1.0*pow(-1.0, direction); // direction = +1 is past, 0 is future

    // these next two lines are the ones we want to modify for the bundle integer.
//    double a = abs(0.1*vphi_in); // added absolute values
//    double b = abs(0.1*vtheta_in);
    
    // set a floor
//    if(a < 1e-6){
//        a = b;
//    }
//    if(b<1e-6){
//        b=a;
//    }
    initCoords.vphitrue_ = vphi_in;
    initCoords.vtheta_ = vtheta_in;
    
/*    double t = M_PI/8.0 + bundle_int*2.0*M_PI/8.0;
    initCoords.vphitrue_ =  vphi_in  + a*cos(t); // vphi_in + <something>
    initCoords.vtheta_ = vtheta_in  + b*sin(t); 
    
    if(bundle_int == 8){
        initCoords.vphitrue_ = vphi_in;
        initCoords.vtheta_ = vtheta_in;
    }
*/   
    initCoords.vr_ = initVR(initCoords, initialValues);

    initCoords.vu_= initCoords.vt_ + pow(-1.0, direction) * initCoords.vr_ * (initCoords.r_ * initCoords.r_ + a_ * a_) / initDelta;  // vu = vt + (r^2+a^2)/Delta vr
    initCoords.vphi_ = initCoords.vphitrue_ + pow(-1.0, direction) *  a_*initCoords.vr_/initDelta;

    coords_ = initCoords;
   
}

KerrRay::KerrRay(double vphi_in, double vtheta_in) {

    // version where we just set vphi and vtheta
    
    a_ = 0.9;
    m_ = 1.0;
    h_ = 0.0010;//starting h value
    int direction = +1;
    totalError_ = 0.0;
    goodStepCount_ = 0;
    badStepCount_ = 0;

    Coords initCoords;
    InitialValues initialValues(m_, a_, direction);

    rplus_ = m_ + sqrt(m_ * m_ - a_ * a_);  // we think this locates the event horizon

    initCoords.alpha_ = 0; // don't think we need this for this constructor
    initCoords.r_ = 20.0;
    initCoords.theta_ = PI/2.0;
    initCoords.phi_ = PI;
    initCoords.u_ = 10.0;
    initCoords.lambda_ = 0.0;
    initCoords.length_ = 0.0; // for spatial distance
    initCoords.t_      = 0.0; // for time
    initCoords.phitrue_ = initCoords.phi_;
    double initDelta = initCoords.r_ * initCoords.r_ - 2.0 * m_ * initCoords.r_ + a_ * a_;
    initCoords.vt_ = 1.0*pow(-1.0, direction); // direction = +1 is past, 0 is future
    initCoords.vphitrue_ = vphi_in;
    initCoords.vtheta_ = vtheta_in; 
    initCoords.vr_ = initVR(initCoords, initialValues);
    initCoords.vu_= initCoords.vt_ + pow(-1.0, direction) * initCoords.vr_ * (initCoords.r_ * initCoords.r_ + a_ * a_) / initDelta;  // vu = vt + (r^2+a^2)/Delta vr
    initCoords.vphi_ = initCoords.vphitrue_ + pow(-1.0, direction) *  a_*initCoords.vr_/initDelta;

    coords_ = initCoords;
}

KerrRay::~KerrRay() {
    //file_.close();
}

void KerrRay::computeNewCoordinates(Coords initCoords, Coords &newCoords, double err[], double h, double m, double a, int direction) {
    // Need to add phitrue calculated throughout as well
    double r_0, theta_0, phi_0, u_0, vr_0, vtheta_0, vphi_0, vu_0, length_0, t_0, phitrue_0;
    double k1r, k1theta, k1phi, k1u, k1vr, k1vtheta, k1vphi, k1vu, k1length, k1t, k1phitrue;
    double k2r, k2theta, k2phi, k2u, k2vr, k2vtheta, k2vphi, k2vu, k2length, k2t, k2phitrue;
    double k3r, k3theta, k3phi, k3u, k3vr, k3vtheta, k3vphi, k3vu, k3length, k3t, k3phitrue;
    double k4r, k4theta, k4phi, k4u, k4vr, k4vtheta, k4vphi, k4vu, k4length, k4t, k4phitrue;
    double k5r, k5theta, k5phi, k5u, k5vr, k5vtheta, k5vphi, k5vu, k5length, k5t, k5phitrue;
    double k6r, k6theta, k6phi, k6u, k6vr, k6vtheta, k6vphi, k6vu, k6length, k6t, k6phitrue;

    Coords tempCoords;

    r_0		  = initCoords.r_;
    theta_0   = initCoords.theta_;
    phi_0     = initCoords.phi_;
    u_0		  = initCoords.u_;
    vr_0	  = initCoords.vr_;
    vtheta_0  = initCoords.vtheta_;
    vphi_0    = initCoords.vphi_;
    vu_0	  = initCoords.vu_;
    length_0  = initCoords.length_;
    t_0       = initCoords.t_;
    phitrue_0 = initCoords.phitrue_;


    k1r       = h * vr_0;
    k1theta   = h * vtheta_0;
    k1phi     = h * vphi_0;
    k1u       = h * vu_0;
    k1vr      = h * VRDot(initCoords, m, a, direction);
    k1vtheta  = h * VThetaDot(initCoords, m, a, direction);
    k1vphi    = h * VPhiDot(initCoords, m, a, direction);
    k1vu      = h * VUDot(initCoords, m, a, direction);
    k1length  = h * lengthDot(initCoords, m, a, direction);
    k1t       = h * tDot(initCoords, m, a, direction);
    k1phitrue = h * VPhiTrue(initCoords, m, a, direction);

    tempCoords.r_ 		= 	r_0       + b21 * k1r;
    tempCoords.theta_ 	= 	theta_0   + b21 * k1theta;
    tempCoords.phi_ 	= 	phi_0     + b21 * k1phi;
    tempCoords.u_ 		= 	u_0       + b21 * k1u;
    tempCoords.vr_ 		= 	vr_0      + b21 * k1vr;
    tempCoords.vtheta_ 	= 	vtheta_0  + b21 * k1vtheta;
    tempCoords.vphi_ 	= 	vphi_0    + b21 * k1vphi;
    tempCoords.vu_ 		= 	vu_0      + b21 * k1vu;
    tempCoords.length_  =   length_0  + b21 * k1length;
    tempCoords.t_       =   t_0       + b21 * k1t;
    tempCoords.phitrue_ =   phitrue_0 + b21 * k1phitrue;

    k2r       = h * tempCoords.vr_;
    k2theta   = h * tempCoords.vtheta_;
    k2phi     = h * tempCoords.vphi_;
    k2u       = h * tempCoords.vu_;
    k2vr      = h * VRDot(tempCoords, m, a, direction);
    k2vtheta  = h * VThetaDot(tempCoords, m, a, direction);
    k2vphi    = h * VPhiDot(tempCoords, m, a, direction);
    k2vu      = h * VUDot(tempCoords, m, a, direction);
    k2length  = h * lengthDot(tempCoords, m, a, direction);
    k2t       = h * tDot(tempCoords, m, a, direction);
    k2phitrue = h * VPhiTrue(tempCoords, m, a, direction);

    tempCoords.r_ 		= r_0      	+ b31 * k1r      	+ b32 * k2r;
    tempCoords.theta_ 	= theta_0  	+ b31 * k1theta  	+ b32 * k2theta;
    tempCoords.phi_ 	= phi_0		+ b31 * k1phi	   	+ b32 * k2phi;
    tempCoords.u_ 		= u_0		+ b31 * k1u	   		+ b32 * k2u;
    tempCoords.vr_ 		= vr_0		+ b31 * k1vr     	+ b32 * k2vr;
    tempCoords.vtheta_ 	= vtheta_0 	+ b31 * k1vtheta 	+ b32 * k2vtheta;
    tempCoords.vphi_ 	= vphi_0	+ b31 * k1vphi   	+ b32 * k2vphi;
    tempCoords.vu_ 		= vu_0		+ b31 * k1vu	   	+ b32 * k2vu;
    tempCoords.length_  = length_0  + b31 * k1length    + b32 * k2length;
    tempCoords.t_       = t_0       + b31 * k1t         + b32 * k2t;
    tempCoords.phitrue_ = phitrue_0 + b31 * k1phitrue   + b32 * k2phitrue;

    k3r		  = h * tempCoords.vr_;
    k3theta   = h * tempCoords.vtheta_;
    k3phi	  = h * tempCoords.vphi_;
    k3u		  = h * tempCoords.vu_;
    k3vr	  = h * VRDot(tempCoords, m, a, direction);
    k3vtheta  = h * VThetaDot(tempCoords, m, a, direction);
    k3vphi	  = h * VPhiDot(tempCoords, m, a, direction);
    k3vu	  = h * VUDot(tempCoords, m, a, direction);
    k3length  = h * lengthDot(tempCoords, m, a, direction);
    k3t       = h * tDot(tempCoords, m, a, direction);
    k3phitrue = h * VPhiTrue(tempCoords, m, a, direction);

    tempCoords.r_       = r_0		+ b41 * k1r	   	  + b42 * k2r	    + b43 * k3r;
    tempCoords.theta_   = theta_0   + b41 * k1theta   + b42 * k2theta   + b43 * k3theta;
    tempCoords.phi_     = phi_0	    + b41 * k1phi	  + b42 * k2phi     + b43 * k3phi;
    tempCoords.u_       = u_0		+ b41 * k1u	   	  + b42 * k2u	    + b43 * k3u;
    tempCoords.vr_      = vr_0		+ b41 * k1vr	  + b42 * k2vr	    + b43 * k3vr;
    tempCoords.vtheta_  = vtheta_0  + b41 * k1vtheta  + b42 * k2vtheta  + b43 * k3vtheta;
    tempCoords.vphi_    = vphi_0	+ b41 * k1vphi    + b42 * k2vphi	+ b43 * k3vphi;
    tempCoords.vu_      = vu_0		+ b41 * k1vu      + b42 * k2vu	    + b43 * k3vu;
    tempCoords.length_  = length_0  + b41 * k1length  + b42 * k2length  + b43 * k3length;
    tempCoords.t_       = t_0       + b41 * k1t       + b42 * k2t       + b43 * k3t;
    tempCoords.phitrue_ = phitrue_0 + b41 * k1phitrue + b42 * k2phitrue + b43 * k3phitrue;

    k4r		  = h * tempCoords.vr_;
    k4theta   = h * tempCoords.vtheta_;
    k4phi	  = h * tempCoords.vphi_;
    k4u		  = h * tempCoords.vu_;
    k4vr	  = h * VRDot(tempCoords, m, a, direction);
    k4vtheta  = h * VThetaDot(tempCoords, m, a, direction);
    k4vphi	  = h * VPhiDot(tempCoords, m, a, direction);
    k4vu	  = h * VUDot(tempCoords, m, a, direction);
    k4length  = h * lengthDot(tempCoords, m, a, direction);
    k4t       = h * tDot(tempCoords, m, a, direction);
    k4phitrue = h * VPhiTrue(tempCoords, m, a, direction);

    tempCoords.r_       = r_0		+ b51 * k1r	      + b52 * k2r	    + b53 * k3r	      + b54 * k4r;
    tempCoords.theta_   = theta_0   + b51 * k1theta   + b52 * k2theta   + b53 * k3theta   + b54 * k4theta;
    tempCoords.phi_     = phi_0	    + b51 * k1phi     + b52 * k2phi	    + b53 * k3phi	  + b54 * k4phi;
    tempCoords.u_       = u_0		+ b51 * k1u	      + b52 * k2u	    + b53 * k3u	      + b54 * k4u;
    tempCoords.vr_      = vr_0		+ b51 * k1vr      + b52 * k2vr	    + b53 * k3vr	  + b54 * k4vr;
    tempCoords.vtheta_  = vtheta_0  + b51 * k1vtheta  + b52 * k2vtheta  + b53 * k3vtheta  + b54 * k4vtheta;
    tempCoords.vphi_    = vphi_0    + b51 * k1vphi    + b52 * k2vphi    + b53 * k3vphi    + b54 * k4vphi;
    tempCoords.vu_      = vu_0		+ b51 * k1vu      + b52 * k2vu      + b53 * k3vu	  + b54 * k4vu;
    tempCoords.length_  = length_0  + b51 * k1length  + b52 * k2length  + b53 * k3length  + b54 * k4length;
    tempCoords.t_       = t_0       + b51 * k1t       + b52 * k2t       + b53 * k3t       + b54 * k4t;
    tempCoords.phitrue_ = phitrue_0 + b51 * k1phitrue + b52 * k2phitrue + b53 * k3phitrue + b54 * k4phitrue;

    // Everything above this comment has been put into new metric structure from 2014

    k5r      = h * tempCoords.vr_;
    k5theta  = h * tempCoords.vtheta_;
    k5phi    = h * tempCoords.vphi_;
    k5u      = h * tempCoords.vu_;
    k5vr     = h * VRDot(tempCoords, m, a, direction);
    k5vtheta = h * VThetaDot(tempCoords, m, a, direction);
    k5vphi   = h * VPhiDot(tempCoords, m, a, direction);
    k5vu     = h * VUDot(tempCoords, m, a, direction);
    k5length = h * lengthDot(tempCoords, m, a, direction);
    k5t      = h * tDot(tempCoords, m, a, direction);
    k5phitrue = h * VPhiTrue(tempCoords, m, a, direction);

    tempCoords.r_       = r_0		+ b61 * k1r	      + b62 * k2r	    + b63 * k3r	      + b64 * k4r	    + b65 * k5r;
    tempCoords.theta_   = theta_0   + b61 * k1theta   + b62 * k2theta   + b63 * k3theta   + b64 * k4theta   + b65 * k5theta;
    tempCoords.phi_     = phi_0	    + b61 * k1phi     + b62 * k2phi	    + b63 * k3phi	  + b64 * k4phi     + b65 * k5phi;
    tempCoords.u_       = u_0	    + b61 * k1u	      + b62 * k2u	    + b63 * k3u	      + b64 * k4u	    + b65 * k5u;
    tempCoords.vr_      = vr_0		+ b61 * k1vr	  + b62 * k2vr	    + b63 * k3vr	  + b64 * k4vr      + b65 * k5vr;
    tempCoords.vtheta_  = vtheta_0  + b61 * k1vtheta  + b62 * k2vtheta  + b63 * k3vtheta  + b64 * k4vtheta  + b65 * k5vtheta;
    tempCoords.vphi_    = vphi_0    + b61 * k1vphi    + b62 * k2vphi    + b63 * k3vphi    + b64 * k4vphi    + b65 * k5vphi;
    tempCoords.vu_      = vu_0		+ b61 * k1vu	  + b62 * k2vu	    + b63 * k3vu	  + b64 * k4vu		+ b65 * k5vu;
    tempCoords.length_  = length_0  + b61 * k1length  + b62 * k2length  + b63 * k3length  + b64 * k4length  + b65 * k5length;
    tempCoords.t_       = t_0       + b61 * k1t       + b62 * k2t       + b63 * k4t       + b64 * k4t       + b65 * k5t;
    tempCoords.phitrue_ = phitrue_0 + b61 * k1phitrue + b62 * k2phitrue + b63 * k3phitrue + b64 * k4phitrue + b65 * k5phitrue;

    k6r      = h * tempCoords.vr_;
    k6theta  = h * tempCoords.vtheta_;
    k6phi    = h * tempCoords.vphi_;
    k6u      = h * tempCoords.vu_;
    k6vr     = h * VRDot(tempCoords, m, a, direction);
    k6vtheta = h * VThetaDot(tempCoords, m, a, direction);
    k6vphi   = h * VPhiDot(tempCoords, m, a, direction);
    k6vu     = h * VUDot(tempCoords, m, a, direction);
    k6length = h * lengthDot(tempCoords, m, a, direction);
    k6t      = h * tDot(tempCoords, m, a, direction);
    k6phitrue = h * VPhiTrue(tempCoords, m, a, direction);

    // r, theta, phi, u, vr, vtheta, vphi, vu

    newCoords.r_       = r_0 	   + c1 * k1r 	    + c2 * k2r 	     + c3 * k3r 	  + c4 * k4r 	   + c5 * k5r 	    + c6 * k6r;
    newCoords.theta_   = theta_0   + c1 * k1theta   + c2 * k2theta   + c3 * k3theta   + c4 * k4theta   + c5 * k5theta   + c6 * k6theta;
    newCoords.phi_     = phi_0 	   + c1 * k1phi 	+ c2 * k2phi 	 + c3 * k3phi 	  + c4 * k4phi 	   + c5 * k5phi 	+ c6 * k6phi;
    newCoords.u_       = u_0 	   + c1 * k1u 	    + c2 * k2u 	     + c3 * k3u 	  + c4 * k4u 	   + c5 * k5u 	    + c6 * k6u;
    newCoords.vr_      = vr_0 	   + c1 * k1vr 	    + c2 * k2vr 	 + c3 * k3vr 	  + c4 * k4vr 	   + c5 * k5vr 	    + c6 * k6vr;
    newCoords.vtheta_  = vtheta_0  + c1 * k1vtheta  + c2 * k2vtheta  + c3 * k3vtheta  + c4 * k4vtheta  + c5 * k5vtheta  + c6 * k6vtheta;
    newCoords.vphi_    = vphi_0    + c1 * k1vphi 	+ c2 * k2vphi 	 + c3 * k3vphi 	  + c4 * k4vphi    + c5 * k5vphi 	+ c6 * k6vphi;
    newCoords.vu_      = vu_0 	   + c1 * k1vu 	    + c2 * k2vu 	 + c3 * k3vu 	  + c4 * k4vu 	   + c5 * k5vu 	    + c6 * k6vu;
    newCoords.length_  = length_0  + c1 * k1length  + c2 * k2length  + c3 * k3length  + c4 * k4length  + c5 * k5length  + c6 * k6length;
    newCoords.t_       = t_0       + c1 * k1t       + c2 * k2t       + c3 * k3t       + c4 * k4t       + c5 * k5t       + c6 * k6t;
    newCoords.phitrue_ = phitrue_0 + c1 * k1phitrue + c2 * k2phitrue + c3 * k3phitrue + c4 * k4phitrue + c5 * k5phitrue + c6 * k6phitrue;
    newCoords.lambda_  = initCoords.lambda_  + h;

    err[0] = (c1 - c1s) * k1r 	 	+ (c2-c2s) * k2r      + (c3-c3s) * k3r      + (c4-c4s) * k4r      + (c5-c5s) * k5r      + (c6-c6s) * k6r;
    err[1] = (c1 - c1s) * k1theta 	+ (c2-c2s) * k2theta  + (c3-c3s) * k3theta  + (c4-c4s) * k4theta  + (c5-c5s) * k5theta  + (c6-c6s) * k6theta;
    err[2] = (c1 - c1s) * k1phi 	+ (c2-c2s) * k2phi    + (c3-c3s) * k3phi    + (c4-c4s) * k4phi    + (c5-c5s) * k5phi    + (c6-c6s) * k6phi;
    err[3] = (c1 - c1s) * k1u 		+ (c2-c2s) * k2u      + (c3-c3s) * k3u      + (c4-c4s) * k4u      + (c5-c5s) * k5u      + (c6-c6s) * k6u;
    err[4] = (c1 - c1s) * k1vr 	    + (c2-c2s) * k2vr     + (c3-c3s) * k3vr     + (c4-c4s) * k4vr     + (c5-c5s) * k5vr     + (c6-c6s) * k6vr;
    err[5] = (c1 - c1s) * k1vtheta 	+ (c2-c2s) * k2vtheta + (c3-c3s) * k3vtheta + (c4-c4s) * k4vtheta + (c5-c5s) * k5vtheta + (c6-c6s) * k6vtheta;
    err[6] = (c1 - c1s) * k1vphi 	+ (c2-c2s) * k2vphi	  + (c3-c3s) * k3vphi	+ (c4-c4s) * k4vphi	  + (c5-c5s) * k5vphi	+ (c6-c6s) * k6vphi;
    err[7] = (c1 - c1s) * k1vu 		+ (c2-c2s) * k2vu	  + (c3-c3s) * k3vu		+ (c4-c4s) * k4vu	  + (c5-c5s) * k5vu		+ (c6-c6s) * k6vu;
}

void KerrRay::step(int direction) {
    Coords newCoords;
    newCoords.alpha_ = coords_.alpha_;
    newCoords.vphiinit_ = coords_.vphiinit_;
    double err[8];

    KerrRay::computeNewCoordinates(coords_, newCoords, err, h_, m_, a_, direction);

    double errmax = 0.0;

    for (int j = 0; j<8; j++){ // this block picks out the biggest error value (in absolute value)
        if (fabs(err[j])>errmax){
            errmax = fabs(err[j]);
        }
    }

    double errRatio = EPS/errmax;

    if(errRatio > 1) {
        // GOOD STEP
        goodStepCount_++;

        // MAKE SURE TO COPY OVER ALL VALUES THAT ARE IMPORTANT

        coords_ = newCoords;
        totalError_ = totalError_ +  errmax;

        double h1 = SAFETY * h_ * pow(errRatio, 0.25); // from old NRC book
        double h2 = h_ * 1.1; // decreased from 10*h on Sept 9 to limit step size growth (Further on Sept 16th to debug)
        h_ = (h1 > h2) ? h2 : h1;
        if(coords_.r_ < 1.3*rplus_){
            h_ = 0.0005; // putting in an artificial floor to have slow steps at small radii
        }

    } else {
        //BAD STEP
        badStepCount_++;

        h_ = SAFETY * h_ * pow(errRatio, 0.2); // .2 from old NRC book

        if (fabs(h_) <= 5.0e-6) {
            h_ = 5.0e-6;
        }
    }
}

double KerrRay::getError() {
    return totalError_;
}

void KerrRay::runScreen(int direction) {
    // this is a version in bundle testing to just run a ray straight (no Kerr bending) so that we can look at 
    //     positions of a bundle of rays in the vtheta / vphi spread
    // since we want to run the ray out along a straight line (initially) we set m = a = 0 here
    
    a_ = 0.0;
    m_ = 0.0;
    double vphi_init = coords_.vphitrue_;
    double vtheta_init = coords_.vtheta_;
    
    // let's just run the steps 5 times
    for(int i = 0; i< 5; i++){
        step(direction);
    }
    double ystep = coords_.r_ * sin(coords_.theta_) * sin(coords_.phitrue_);
    double zstep = coords_.r_ * cos(coords_.theta_);
    double d = dfunc();
    ofstream fout;
    fout.open("output.csv", ios::app);
    fout<< "vphi = , " << vphi_init << " vtheta = , " << vtheta_init<< " y and z = ," << 1000*ystep << "," << 1000*zstep << "," << d << endl;
    fout.close();

}

double KerrRay::dfunc(){
    double ystep = coords_.r_ * sin(coords_.theta_) * sin(coords_.phitrue_);
    double zstep = coords_.r_ * cos(coords_.theta_);
    return pow(ystep - yfin_,2) + pow(zstep-zfin_,2);
}

void KerrRay::runtox(int direction, double xsource) {
    double xstep;

    while(goodStepCount_ < 50000){

        step(direction);
        xstep = coords_.r_ * sin(coords_.theta_) * cos(coords_.phitrue_);
 //       if(goodStepCount_ % 50 == 0)
 //           cout << coords_.t_ <<"  "<< coords_.r_ <<"  "<< xstep << endl;
        if (coords_.r_ < rplus_) {
            cout << "fell into black hole with t = " << coords_.t_ << endl;
            cout << "bad steps = " << badStepCount_ << " good steps = " << goodStepCount_ << endl;
            break;
        }
        else if (xstep > xsource) {
//            cout << "reached final x plane at t = " << coords_.t_ << " and x = " << xstep << endl;
//            double ystep = coords_.r_ * sin(coords_.theta_) * sin(coords_.phitrue_);
//            double zstep = coords_.r_ * cos(coords_.theta_);
//           cout << "at y = " << ystep <<" and z = " << zstep << endl;
            break;
        }
        else if (pow(-1, direction)*coords_.t_ > abs(250.0)) {
            cout << "reached time limit" << endl << coords_.t_ << " " << xstep << endl;
            break;
        }
    }

//    cout << "good steps = " << goodStepCount_ <<" and bad steps = " <<badStepCount_ << endl;
//    cout << "final error " << totalError_ << endl;

}


double KerrRay::vThetaLimit(Coords coords, InitialValues iv){
    iv.computeValues(coords);
    return sqrt((-iv.g00_ * pow(coords.vt_ * sin(coords.alpha_), 2.0)
                 -2.0 * iv.g03_ * coords.vphitrue_*coords.vt_ * pow(sin(coords.alpha_), 2.0)
                 - iv.g33_*coords.vphitrue_*coords.vphitrue_)/iv.rhosq_);
}

double KerrRay::posVPhiLimit(Coords coords, InitialValues iv){
    iv.computeValues(coords);

    return - coords.vt_/iv.g33_ *(iv.g03_ * sin(coords.alpha_) * sin(coords.alpha_)
                                  + sin(coords.alpha_) * sqrt(pow(iv.g03_*sin(coords.alpha_), 2.0) - iv.g00_*iv.g33_));
}

double KerrRay::negVPhiLimit(Coords coords, InitialValues iv){
    iv.computeValues(coords);
    return - coords.vt_/iv.g33_ *(iv.g03_ * sin(coords.alpha_) * sin(coords.alpha_)
                                  - sin(coords.alpha_) * sqrt(pow(iv.g03_*sin(coords.alpha_), 2.0) - iv.g00_*iv.g33_));
}

double KerrRay::initVR(Coords coords, InitialValues iv){
    iv.computeValues(coords);
    return -1.0 * sqrt((-1.0*iv.g00_*coords.vt_*coords.vt_ - 2.0*iv.g03_*coords.vt_*coords.vphitrue_
                        - iv.rhosq_*coords.vtheta_*coords.vtheta_ - iv.g33_*coords.vphitrue_*coords.vphitrue_)*iv.delta_/iv.rhosq_);
}


/* code not used below here */

/* 

KerrRay::KerrRay(double alpha, int direction) {
    a_ = 0.9;
    m_ = 1.0;
    h_ = 0.0010;//starting h value
    totalError_ = 0.0;
    goodStepCount_ = 0;
    badStepCount_ = 0;

    Coords initCoords;
    double vThetaLimitVal, vPhiPlusVal, vPhiMinusVal, initDelta;
    InitialValues initialValues(m_, a_, direction);

    rplus_ = m_ + sqrt(m_ * m_ - a_ * a_);  // we think this locates the event horizon

    initCoords.alpha_ = alpha;
    initCoords.r_ = 20.0;
    initCoords.theta_ = PI/2.0;
    initCoords.phi_ = PI;
    initCoords.u_ = 10.0;
    initCoords.lambda_ = 0.0;
    initCoords.length_ = 0.0; // for spatial distance
    initCoords.t_      = 0.0; // for time
    initCoords.phitrue_ = initCoords.phi_;
    initDelta = initCoords.r_ * initCoords.r_ - 2.0 * m_ * initCoords.r_ + a_ * a_;

    initCoords.vt_ = 1.0*pow(-1.0, direction); // direction = +1 is past, 0 is future

    vPhiPlusVal = negVPhiLimit(initCoords, initialValues);
    vPhiMinusVal = posVPhiLimit(initCoords, initialValues); // turns out that this one is less than the other one

    initCoords.vphitrue_ = vPhiPlusVal;
//    initCoords.vphitrue_ = vPhiMinusVal;

    initCoords.vtheta_ = 0.0 * vThetaLimitVal; // working in the plane.
    initCoords.vr_ = initVR(initCoords, initialValues);

    // these two lines added August  31, 2018, need to check directions
    initCoords.vu_= initCoords.vt_ + pow(-1.0, direction) * initCoords.vr_ * (initCoords.r_ * initCoords.r_ + a_ * a_) / initDelta;  // vu = vt + (r^2+a^2)/Delta vr
    initCoords.vphi_ = initCoords.vphitrue_ + pow(-1.0, direction) *  a_*initCoords.vr_/initDelta;

    coords_ = initCoords;

//    cout << coords_.vr_ <<" "<< coords_.vu_ << "  " << coords_.vphi_ << "  " << coords_.vphitrue_ << endl;

}

int KerrRay::runset(int direction, int completed, double timeStop) {

    ofstream fout;
    fout.open("phiplus.dat", ofstream::app);
    while(goodStepCount_ < 20000){

        step(direction);
//        if(goodStepCount_ % 5 == 0){
//          cout << coords_.t_ <<"  "<< coords_.r_ <<"  "<< coords_.phitrue_ << endl;
//        }
//        fout << coords_.t_ <<" "<< coords_.r_ * cos(coords_.phitrue_) <<" " << coords_.r_*sin(coords_.phitrue_)<<endl;

        if (coords_.r_ < rplus_) {
//            cout << "fell into black hole with t = " << coords_.t_ << endl;
//            cout << "bad steps = " << badStepCount_ << " good steps = " << goodStepCount_ << endl;
            break;
        }
        else if (coords_.r_ > 30) {
            //          cout << "escaped black hole with t = " << coords_.t_ << " and phi = " << coords_.phitrue_ << endl;
//            fout << coords_.alpha_ <<" " << coords_.phitrue_ << endl;
            completed++;
            cout  << completed << "  "<<  coords_.alpha_ <<"  " << coords_.phitrue_ << endl;
            break;
        }
        else if (pow(-1, direction)*coords_.t_ > abs(timeStop)) {
            cout << "reached time limit" << endl << coords_.t_ << " " << coords_.r_ << " " << coords_.phitrue_ << endl;
            break;
        }
    }

    cout << "good steps = " << goodStepCount_ <<" and bad steps = " <<badStepCount_ << endl;
    cout << "final error " << totalError_ << endl;
    fout.close();
    return completed;

}

void KerrRay::runray(int direction, double timeStop) {

ofstream fout;
fout.open("test8.csv", ofstream::app);
while(goodStepCount_ < 50000){

step(direction);
        if(goodStepCount_ % 5 == 0){
          cout << coords_.t_ <<"  "<< coords_.r_ <<"  "<< coords_.phitrue_ << endl;
        }
        fout << coords_.t_ <<" "<< coords_.r_ * cos(coords_.phitrue_) <<" " << coords_.r_*sin(coords_.phitrue_)<<endl;

    if (coords_.r_ < rplus_) {
        cout << "fell into black hole with t = " << coords_.t_ << endl;
        cout << "bad steps = " << badStepCount_ << " good steps = " << goodStepCount_ << endl;
        break;
    }
    else if (coords_.r_ > 30) {
          cout << "escaped black hole with t = " << coords_.t_ << " and phi = " << coords_.phitrue_ << endl;
//        fout << coords_.alpha_ <<" " << coords_.phitrue_ << endl;
        break;
    }
    else if (pow(-1, direction)*coords_.t_ > abs(timeStop)) {
        cout << "reached time limit" << endl << coords_.t_ << " " << coords_.r_ << " " << coords_.phitrue_ << endl;
        break;
    }
}

    cout << "good steps = " << goodStepCount_ <<" and bad steps = " <<badStepCount_ << endl;
    cout << "final error " << totalError_ << endl;
    fout.close();

}


*/


