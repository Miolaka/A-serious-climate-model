#include <stdio.h>
#include <math.h>



//double sum(double x, double y){}
double B_plank(double h_planck, double c_light, double lambda, double k_boltzman, double T){
    double B_plank = (2.0*h_planck*c_light*c_light)/pow(lambda,5)/(exp(h_planck*c_light/lambda/k_boltzman/T)-1);
    return B_plank;
}


//temperature to pot. temp
double T2theta (double T, double p_ground, double p_middle){
    double theta = T*pow(p_ground/p_middle,2.0/7.0);
    return theta;
}

//pot.temp to temp, 1000hPa=100000Pa
double theta2T (double theta, double p_ground, double p_middle){
    double T = theta*pow(p_ground/p_middle,-1*(2.0/7.0));
    return T;
}

//void(){}
    //print is a void

//stability check
int stability(double theta_i, double theta_j){
    if(theta_i>theta_j){
        return 1;
    }
    else{
        return 0;
    }
}

//radiation, E_em, epsilon, sigma_b, alpha
//need E_up, E_dw
//float E_up[];
//float E_dw[];
//float E_delta;


//E_up[nlayer]=sigma_b*pow(T[nlayer-1],4);
//E_up[nlayer-1]=E_up[nlayer]*(1-alpha)+epsilon*sigma_b*T[nlayer-2];
//E_up[i]=E_up[i+1]*(1-alpha)+epsilon*sigma_b*T[?]

//E_dw[0]=0; //emission of space
//E_delta = E_dw[i]+E_up[i+1] - E_dw[i+1] - E_up[i];



