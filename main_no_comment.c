#include <stdio.h>
#include <math.h>
#include "header.h"

int main(int argc, char **argv)
{
    int nlayer=100;  
    double p_ground = 100000.0;  
    int nlevel = nlayer+1;  
    double p[nlevel];  
    double theta[nlayer];  
    double T[nlayer];  
    double adiabatic_constant = 2.0/7.0;  
    int unstable = 1;  
    double p_delta = p_ground/nlayer;  
    double p_middle;
    int ntime = 20000;
    double sigma_b = 5.67*pow(10.0,-8.0);  
    double E_up[nlevel];
    double E_down[nlevel];
    double delta_E[nlevel];
    double delta_time = 60.0*60.0*1.0*1.0;  
    double delta_Temp[nlayer];
    double C_p = 1004.0;  
    double g = 9.81;  
    double L_up[nlevel];
    double L_down[nlevel]; 
    int imu_max = 10;  
    double delta_mu=1.0/imu_max;
    double absorb;
    double emi;
    double pi = 3.14159265;
    double B[nlayer];
    double lambda;
    double h_planck = 6.6260628E-34;
    double c_light = 299792458;
    double k_boltzman = 1.3806503E-23;
    double tau_co2 = 1.0;
    double tau_h2o = 1.0;
    double tau = tau_co2+tau_h2o;
    double lambda_start = 1.0*pow(10.0,-6.0); 
    double lambda_end = 28.0*pow(10.0,-6.0);
    int gridnumber = 20;
    double delta_lambda = (lambda_end-lambda_start)/gridnumber;

    printf("test B is %f \n\n", B_plank(h_planck, c_light, 0.0001, k_boltzman, 300.0));
    for(int i=0; i<nlevel; i++){
        p[i]=p_ground-(nlevel-i)*p_delta;  
         
    }
    printf("--------------------------------\nInitiating temperature \n");
    for(int i=0; i<nlayer; i++){ 
        T[i]= 180 + (i)*5.0;  
        theta[i]= T2theta(T[i],p_ground,(p[i]+p[i+1])/2);
        printf("%d, theta=%.2f, T=%.2f, p_middle=%.2f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
    }
    printf("done \n--------------------------------\n");
    for(int it=0; it<ntime; it++){
        printf("Iteration: %d, \n--------------------------------\n recalculating Theta temperature \n", it);
        for (int i=0; i<nlayer; i++){  
            theta[i] = T2theta(T[i],p_ground,(p[i]+p[i+1])/2.0);
        }             
        for(int i=0; i<nlayer; i++){                        
        }
        printf("done \n--------------------------------\n");     
        unstable = 1;
        while(unstable == 1){
            unstable = 0;  
            for(int i=1; i<nlayer; i++){  
                if(theta[i-1]<theta[i]){   
                    double temp_cache = theta[i-1];
                    theta[i-1] = theta[i];
                    theta[i] = temp_cache; 
                    temp_cache = 0.0; 
                    unstable = 1;
                }
            }
        }
        for (int i=0; i<nlayer; i++){                
            T[i] = theta2T(theta[i],p_ground, (p[i]+p[i+1])/2.0);

        }    
            printf("\n");
        for (double lambda = lambda_start; lambda<lambda_end; lambda = lambda+delta_lambda){
            if(lambda>0&& lambda< 8.0 *pow(10.0,-6.0)){
                tau_co2 = 1.0;
                tau_h2o = 1.0;
                tau = tau_co2+tau_h2o;
            }
            else if(lambda > 8.0*pow(10,-6) && lambda<12*pow(10.0,-6.0)){
                tau_co2 = 0.01;
                tau_h2o = 0.01;
                tau = tau_co2+tau_h2o;
            }
            else if(lambda>12.0 *pow(10.0,-6.0)){
                tau_co2 = 1.0;
                tau_h2o = 1.0;
                tau = tau_co2+tau_h2o;
            }
            else{
                tau_co2 = 1.0;
                tau_h2o = 1.0;
                tau = tau_co2+tau_h2o;
            }
            double L_B[nlayer];
            for(int i=0; i<nlayer; i++){
                B[i] = B_plank(h_planck, c_light, lambda, k_boltzman, T[i]);
                L_B[i] = B[i]*delta_lambda;
            }
            for(int k=0; k<nlevel; k++){
                E_up[k] = 0.0;
                E_down[k] = 0.0;
            }
            for(int imu=0; imu<imu_max; imu++){
                double mu= (1.0+imu)/imu_max - delta_mu/2;  
                absorb=1.0-exp(-tau/mu);  
                emi=absorb;
                L_up[nlevel-1]=L_B[nlayer-1];  
                for(int i=nlevel-2; i>-1; i--){
                    L_up[i]=L_up[i+1]*(1-absorb)+(1/pi)*L_B[i]*emi;
                    E_up[i]=E_up[i]+2*pi*L_up[i]*mu*delta_mu;
                }
                L_down[0]=0;
                for(int i=1; i<nlevel; i++){
                    L_down[i]=L_down[i-1]*(1-absorb)+(1.0/pi)*L_B[i-1]*emi;
                    E_down[i]=E_down[i]+2.0 *pi*L_down[i]*mu*delta_mu;
                }
            }
        }
        for(int i=0; i<nlayer; i++){
            delta_E[i] = E_down[i]+E_up[i+1]-E_up[i]-E_down[i+1];
        }  
        for(int i=0; i<nlayer; i++){
            delta_Temp[i]=delta_time*delta_E[i]*g/((p[i]-p[i+1])*C_p);
            T[i]=T[i]+delta_Temp[i];
        }
        printf("--------------------------------\n Profile: \n");
        for(int i=0; i<nlayer; i++){ 
                printf("%d, theta=%.2f, T=%.2f, p_middle=%.2f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
            }
            printf("--------------------------------\n \n");
        if(theta[nlayer-1]<0){
            printf("%f", theta[nlayer-1]);
            break;
        }
    }
return 0; 
}
/*
edit note: 17.Okt.2023
edit note: 30.Okt.2023
edit date: 26.Nov.2023s

last time: CO2 and gas
this time: real radiative transfer
last week Wedsday CO2 

comment paul: Woher kommt aktuell die Energie in deinem Modell?
Mir scheint, es gibt keinen Input von der Sonne. Daher wird es immer kälter.
Am besten, du gibst in der Energiebilanz der untersten Schicht noch einen Term E0/4(1-A) dazu.
E0/4(1-A) 
*/
