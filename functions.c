#include <stdio.h>
#include <math.h>


//float sum(float x, float y){}

//temperature to pot. temp
float T2theta (float T, float p_ground, float p_middle){
    float theta = T*pow(p_ground/p_middle,2.0/7.0);
    return theta;
}

//pot.temp to temp, 1000hPa=100000Pa
float theta2T (float theta, float p_ground, float p_middle){
    float T = theta*pow(p_ground/p_middle,-1*(2.0/7.0));
    return T;
}

//void(){}
    //print is a void

//stability check
int stability(float theta_i, float theta_j){
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



