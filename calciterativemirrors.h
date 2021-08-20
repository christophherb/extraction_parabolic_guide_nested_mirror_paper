#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double * get_r_at_z0(int number, double z_0, double r_0, double LStart, double LEnd, double lStart, double lEnd) {
    /*
    number: how many mirrors,
    z_0: z-position of returned points on mirror
    r_0: r-distance of returned points on mirror
    LStart: Position of the first focal point
    LEnd: Postion of the second focal point
    lStart: Beginning of the mirror
    lEnd: End of the mirror
    */
    int n = number;
    double *r_z0s = malloc(n*sizeof(double_t)); /* n is an array of 10 integers */
	r_z0s[0] = r_0;
    //helper variables as in conic_finite_mirror.h and explained in swissneutronics_überlegungen
    double k1;
    double k2;
    double k3;
    double c;
    double u;
    double a;
    double r_lEnd;
    double r_lStart;
    //initial mirror is calculated from the initial point z0, r0
    c = (LEnd - LStart)/2;
    u = (z_0 + c - LEnd);
    a = sqrt((u*u+c*c+r_0*r_0+sqrt(pow(u*u+c*c+r_0*r_0, 2)-4*c*c*u*u))/2);
    k3 = c*c/(a*a)-1;
    k2 = 2*k3*(c-LEnd);
    k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a;
    r_lEnd = sqrt(k1+ k2*lEnd + k3*lEnd*lEnd);
    r_lStart = r_lEnd*(lStart-LStart)/(lEnd-LStart);
    r_z0s[0] = r_0;
	//next mirror will be calculated with the point on the surface being lStart, r_lStart
	for( int k = 1; k < number;++k){
        c = (LEnd - LStart)/2;
        u = (lStart + c - LEnd);
        a = sqrt((u*u+c*c+r_lStart*r_lStart+sqrt(pow(u*u+c*c+r_lStart*r_lStart, 2)-4*c*c*u*u))/2);
        k3 = c*c/(a*a)-1;
        k2 = 2*k3*(c-LEnd);
        k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a;
        r_lEnd = sqrt(k1+ k2*lEnd + k3*lEnd*lEnd);
        r_lStart = r_lEnd*(lStart-LStart)/(lEnd-LStart);
		if (r_z0s[k-1] > 0)
        {
		r_z0s[k] = sqrt(k1+ k2*z_0 + k3*z_0*z_0);
        }
        else
        {
        r_z0s[k] = -sqrt(k1+ k2*z_0 + k3*z_0*z_0);
        }
	};
   return r_z0s;
}



