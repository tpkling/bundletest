//
// Created by Thomas Kling on 8/30/18.
// Purpose is to examine individual rays / make plots outside walnut
// Working to test bundle ideas
//

#include <iostream>
#include <cmath>
#include "Coords.h"
#include "KerrRay.h"
#include <math.h>
#include <algorithm>

using namespace std;

int main()
{
    cout << "Hello World" << endl;
    double v = 0.015;
    int direction = +1; // for past directed rays
    double a, b, t, vphi_c, vtheta_c, xt,yt,zt;
    int N = 50; // number of times to run the refining loop
    double mins[N]; // array to store the smallest distance on each step
    double testB[13];
    double vphiB[13];
    double vthetaB[13];
    double xfin = 20.0;
    bool centerVal;
    int repeatC = 0;
    
//    KerrRay ray(v,v); // desired ray
//    ray.runtox(direction, xfin);
//    double y = ray.coords_.r_ * sin(ray.coords_.theta_) * sin(ray.coords_.phitrue_);
//    double z = ray.coords_.r_ * cos(ray.coords_.theta_);
//    double x = ray.coords_.r_ * sin(ray.coords_.theta_) * cos(ray.coords_.phitrue_);
//    cout<< x <<"  " << y << "  " << z << endl;
  
    double y = -2;
    double z = -4;
    double dummy = 1000.0;
    double fbig    = 0.2;
    double fsmall  = 0.0005;
    double f = 0.1;
    vphi_c = 0.01362;
    vtheta_c = -0.00574;
    a = abs(f*vphi_c); // added absolute values
    b = abs(f*vtheta_c);

    
    for(int j = 0; j<N; j++){
        for(int i = 0; i<13; i++){
            if(i<12){
                t = M_PI/12.0 + 2.0*i*M_PI/12.0;
                vphiB[i]   = vphi_c   + a*cos(t);
                vthetaB[i] = vtheta_c + b*sin(t); 
            }
            else{
                vphiB[i]   = vphi_c;
                vthetaB[i] = vtheta_c;
            }
            KerrRay bundleRay(vphiB[i], vthetaB[i], y, z);
            bundleRay.runtox(direction, xfin);
            yt = bundleRay.coords_.r_ * sin(bundleRay.coords_.theta_) * sin(bundleRay.coords_.phitrue_);
            zt = bundleRay.coords_.r_ * cos(bundleRay.coords_.theta_);
            xt = bundleRay.coords_.r_ * sin(bundleRay.coords_.theta_) * cos(bundleRay.coords_.phitrue_);
            testB[i] = bundleRay.dfunc();
//            cout << i << "  " <<  vphiB[i] <<"  " << vthetaB[i] <<"  " << testB[i] << "  " << xt << "  " << yt << "  " << zt << endl;
        } // closes loop on i
        dummy = 1000.0;
        for(int k=0; k<13; k++){
            if(testB[k] < dummy){
                dummy = testB[k];
                vphi_c   = vphiB[k];
                vtheta_c = vthetaB[k];
            }
        }

        if(testB[12] == dummy) {
            centerVal = true;
            f = 0.3*f;
            repeatC++;
        }
        else f = 2.5*f;
//        if(dummy > 1) f = fbig;


//        else f= 0.95*f;
//        if(f<fsmall) f = fsmall;
        mins[j] = dummy;
 /*       else{
            centerVal = false;
            repeatC = 0;
        }
        mins[j] = dummy;
        if(centerVal) { // make the a and b smaller if best value was central value
            if (repeatC > 4) {
                f =  f1;
                repeatC = 0;
            }
            else f = 0.5 * f;
        }
        else f = f1/double(j+1);
*/        a = abs(f*vphi_c); // make a new a and b
        b = abs(f*vtheta_c);
        if(a < 1e-6) a=b;
        if(b < 1e-6) b=a;
        
        cout << j << "  " << vphi_c <<"  " << vtheta_c << "  " << mins[j] <<"  " <<  a << "  " << b << "  " << f << endl ;
    } // closes loop on j

//    double val = *min_element(testB, testB+9);
//    cout << "minimum distance is " << val << " and we got that in the loop too " << dummy <<"  "<< bundle_min << endl;  // this works pretty good, but we really want the index of this value.


 /*   for(int j = 0; j<8; j++){
        double theta = j*2*M_PI/8.0; 
        double vtheta = v*sin(theta);
        double vphi = v*cos(theta);
        //cout<< theta << "  " << vphi << "  " << vtheta << endl;
        for(int i = 0; i<8; i++){
            KerrRay ray1(vphi, vtheta, i);
            ray1.runScreen(direction);
        }
    }
    */
    
    //will want to make a for loop to run over bundle_int below
   // KerrRay ray1(vphi, vtheta, bundle_int); // creates a kerr ray called ray1
   // ray1.runScreen(direction);

    return 0;
}
