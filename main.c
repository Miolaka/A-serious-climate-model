#include <stdio.h>
#include <math.h>
#include "header.h"
//#include "functions.c" //maybe not required

int main(int argc, char **argv)
{
    int nlayer=100; //10 elemente 0-9, p[9]<p[nlayer]
    double p_ground = 100000.0; //1000hPa at ground
    int nlevel = nlayer+1; //level is 1 more than layers
    double p[nlevel]; //pressure
    double theta[nlayer]; //potential temperature
    double T[nlayer]; //temperature
    double adiabatic_constant = 2.0/7.0; //R_a/C_p
    int unstable = 1; // stability condition 1 or 0
    double p_delta = p_ground/nlayer; //100hPa per level
    double p_middle;
    int ntime = 200;
    //double heat_temp = 1;
    //legacy:
    //double epsilon = 0.33; //thickness???
    //double alpha = epsilon;
    double sigma_b = 5.67*pow(10.0,-8.0); //5.67*10^-8
    double E_up[nlevel];
    double E_down[nlevel];
    double delta_E[nlevel];
    //int delta_time = 60*60*1; //12 stunden
    double delta_time = 60.0*60.0*1.0*1.0; //1 stunden
    double delta_Temp[nlayer];
    double C_p = 1004.0; // specific constant
    double g = 9.81; //gravitational constant
    double L_up[nlevel];
    double L_down[nlevel]; 
    int imu_max = 10; //split in 10 steps
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
    //double tau2d[lambda][p]; //not working, non-int
                

    //printf("class for convetion\n");
    //printf("nlayer=%d, nlevel=%d\n", nlayer, nlevel);

    //1. Initial variables
    //1.1 Create pressure gradient, loop time = nlevel
    printf("test B is %f \n\n", B_plank(h_planck, c_light, 0.0001, k_boltzman, 300.0));

    for(int i=0; i<nlayer; i++){
        p[i]=p_ground-(nlayer-i)*p_delta; //1000hPa - n*100hPa
        //printf("%dpressure=%f\n",i,p[i]);
    }
    //1.2 Create temperature gradient, loop time = nlayer
    printf("--------------------------------\nInitiating temperature \n");
    for(int i=0; i<nlayer; i++){ 
        T[i]= 180 + (i)*5.0; //-5K each layer,to ground temp 293K
        //printf("T=%f, i=%d \n",T[i],i);
        theta[i]= T2theta(T[i],p_ground,(p[i]+p[i+1])/2);
        //theta[i]=T[i]*pow((p_ground/p_middle),adiabatic_constant); //gamma?
        

        printf("%d, theta=%.2f, T=%.2f, p_middle=%.2f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
    }
    printf("done \n--------------------------------\n");

    

    //2. Time loop
    for(int it=0; it<ntime; it++){
        //2.1 calculate Theta
        printf("Iteration: %d, \n--------------------------------\n recalculating Theta temperature \n", it);
        for (int i=0; i<nlayer; i++){ //runs nlevel times
            theta[i] = T2theta(T[i],p_ground,(p[i]+p[i+1])/2.0);
            //printf("hi2\n");
            //printf("i is %d\n", i);
        }
            //printf("i is", i);
        for(int i=0; i<nlayer; i++){ 
            //!!!
            //printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
        }
        printf("done \n--------------------------------\n");
    
        //printf("2.1 done \n");

        //2.2 Convection //unstable if 0
        unstable = 1;
        while(unstable == 1){
            unstable = 0; //stable assumption
            for(int i=1; i<nlayer; i++){ //till which i?
                //print info
                //for(int i=0; i<nlayer+1; i++){ 
                //    printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
                //}
                //printf("\n");
                if(theta[i-1]<theta[i]){ //looking find a problem
                    //exchange
                    double temp_cache = theta[i-1];
                    theta[i-1] = theta[i];
                    theta[i] = temp_cache; 
                    temp_cache = 0.0; 
                    //print info
                    //for(int i=0; i<nlayer+1; i++){ 
                    //    printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
                    //}/webmail/webmail/
                    //printf("\n");
                    //if(i=nlayer){
                    //    i=1;
                    //}
                    unstable = 1;
                }
            }
        }
        //printf("2.2 done \n");
        //2.3 new Theta to temperature
        for (int i=0; i<nlayer; i++){                
            T[i] = theta2T(theta[i],p_ground, (p[i]+p[i+1])/2.0);
            //print info
            
            //printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
        }    
            printf("\n");
                    
        //printf("2.3 done \n");
        //2.4 Change phase ####### (Radiation)
        //radiation
        //2.4.1 Calculate E (naive cloud radiation)
        
        //for(int i=0; i<nlevel-1; i++){
        //    E_up[i]=E_up[i+1]*(1+alpha)+epsilon*sigma_b*pow(T[i],4);
        //    E_down[i+1]=epsilon*sigma_b*pow(T[i+1],4)+(1-alpha)*E_down[i];
        //}
  
        //2.4.2 Calculate better E (better cloud radiation)
        //Radiative Transfer
        //for wavelength-> for directions-> for ilayer
        //lambda from 1 micro to 26 micro
        //add radiative transfer window
        for (double lambda = lambda_start; lambda<lambda_end; lambda = lambda+delta_lambda){
            // printf("lambda is %f \n", lambda);
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
            //printf("2.4 lambda done \n");

            
            //zz=zz+ B[i]*delta_lambda;



            //layer loop, calculate Radiance B, irradiance E_up_B
            double L_B[nlayer];
            for(int i=0; i<nlayer; i++){
                //calculate plack law -> B(T) -> L units
                //a function to integrate B
                //B[i]=(2*h_planck*c_light*c_light)/(((pow(lambda,5)))*((exp(h_planck*c_light/(lambda*k_boltzman*T[i])))-1));
                
                B[i] = B_plank(h_planck, c_light, lambda, k_boltzman, T[i]);
                //printf("B_i %d = %.9f\n", i, B[i]);
                L_B[i] = B[i]*delta_lambda;
                //printf("L_B %d = %.9f\n", i, L_B[i]);
            }

            //loop for angle mu
            //set initial array element value = 0
            for(int k=0; k<nlevel; k++){
                E_up[k] = 0.0;
                E_down[k] = 0.0;
            }
            for(int imu=0; imu<imu_max; imu++){
                double mu= (1.0+imu)/imu_max - delta_mu/2; //verschiebung integral punkte
                absorb=1.0-exp(-tau/mu); //works
                //printf("mu=%f   ", mu);
                //printf("mu/2=%f  ", delta_mu/2);
                //printf("absorb=%f \n", absorb);
                emi=absorb;

                //upward L_up
                //(p,T,tau -> E_up, E_dw)
                //L_up[nlevel]=(1/pi)*sigma_b*pow(T[nlevel-1],4);
                L_up[nlevel-1]=L_B[nlayer-1]; //the actually highest defined elements

                for(int i=nlevel-2; i>-1; i--){
                    //calculate L
                    L_up[i]=L_up[i+1]*(1-absorb)+(1/pi)*L_B[i]*emi;
                    //printf("L_up %d = %f\n", i, L_up[i]);
                    //L_up[i]=L_up[i+1]*(1-absorb)+(1/pi)*sigma_b*pow(T[i],4)*emi;
                    E_up[i]=E_up[i]+2*pi*L_up[i]*mu*delta_mu;
                    //E_up[i]=E_up[i]+E_up_B;

                }
                //printf("E_up[0] = %f\n", E_up[0]);

                //downward L_down
                //
                L_down[0]=0;
                //ERROR! : L_B[i] defined with nlayer not nlevel! L_B[nlevel] overshooting!
                for(int i=1; i<nlevel; i++){
                    L_down[i]=L_down[i-1]*(1-absorb)+(1.0/pi)*L_B[i-1]*emi;
                    //printf("L_B %d = %.9f \n", i, L_B[i]);
                    //printf("L_down %d = %.9f \n", i, L_down[i]);
                    //L_down[i]=L_down[i-1]*(1-absorb)+(1/pi)*sigma_b*pow(T[i-1],4)*emi;
                    E_down[i]=E_down[i]+2.0 *pi*L_down[i]*mu*delta_mu;
                    //printf("E_dw %d = %.9f\n", i, E_down[i]);
                    
                }

                //Radiance to Irradiance

            }
        }
        //printf("2.4 done \n");
        //???
        
        //2.5 Calculate delta_E (at each level) (Radiance)
        //Problem Paul: what about i=nlayer-1?
        for(int i=0; i<nlayer; i++){
            delta_E[i] = E_down[i]+E_up[i+1]-E_up[i]-E_down[i+1];
            //there's a problem at ground layer.
                    //Solar radiation
                    if (i==nlayer-1){
                    //E=2*pi*(L(mu)*mu[i]*delta_mu;)
                        double E_0 = 1361.0; //solar irradiance E_0 [W/m2]
                        printf("absorb = %f\n", absorb);
                        printf("E_0 = %f\n", E_0/4.0/(1-absorb));
                        E_down[i] = E_down[i] + E_0/4.0/(1-absorb); //yet A is now calculated with window
                    }
            
            //printf("E_down[i]=%f    E_up[i+1]%f E_up[i]%f   E_down[i+1]%f   \n  ",E_down[i],E_up[i+1],E_up[i],E_down[i+1]);
            //printf("delta_E by i= %d = %.9f\n", i, delta_E[i]);
        } //delta_E[nlayer]

        //printf("2.5 done \n");
        //for (int i=0; i<nlayer; i++){
        //    double delta[nlayer-1]
        //    delta[i]=
        //}

        //heat with 1K
        //T[nlevel] = T[nlevel]+heat_temp; 
        //stability = 1; //make unstable

        //2.6 From delta_E to temperature
        for(int i=0; i<nlayer; i++){
            delta_Temp[i]=delta_time*delta_E[i]*g/((p[i]-p[i+1])*C_p);
            T[i]=T[i]+delta_Temp[i];
            //printf("delta_T by i=%d is %.9f \n", i, delta_Temp[i]);
        }

        //2.7 Tempeature to Theta (loop 2.1)
        //7. loop
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
//edit note: 17.Okt.2023
//edit note: 30.Okt.2023
//edit date: 26.Nov.2023s

//last time: CO2 and gas
//this time: real radiative transfer
//last week Wedsday CO2 
//comment paul: Woher kommt aktuell die Energie in deinem Modell?
//Mir scheint, es gibt keinen Input von der Sonne. Daher wird es immer kälter.
//Am besten, du gibst in der Energiebilanz der untersten Schicht noch einen Term E0/4(1-A) dazu.
//E0/4(1-A)