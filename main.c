#include <stdio.h>
#include <math.h>
#include "header.h"
//#include "functions.c" //maybe not required

int main(int argc, char **argv)
{
    int nlayer=10; //p[10] = T[9]
    float p_ground = 100000.0; //1000hPa at ground
    int nlevel = nlayer+1; //level is 1 more than layers
    float p[nlevel]; //pressure
    float theta[nlayer]; //potential temperature
    float T[nlayer]; //temperature
    float rho[nlayer]; //density
    float altitude[nlevel]; //Altitude
    float adiabatic_constant = 2.0/7.0; //R_a/C_p
    int unstable=1; // stability condition 1 or 0
    float p_delta = p_ground/nlevel; //100hPa per level
    float p_middle;
    int heat_time = 15;
    float heat_temp = 1;
    float epsilon=0.33; //thickness???
    float alpha = epsilon;
    float sigma_b=5.67*pow(10,-8); //5.67*10^-8
    float E_up[nlevel];
    float E_down[nlevel];
    float delta_E[nlevel];
    int delta_time = 60*60*12; //12 stunden
    float delta_Temp[nlevel];
    float C_p = 1004; // specific constant
    float g= 9.81; //gravitational constant
    float L_up[nlayer];
    float L_down[nlayer];
    float delta_mu;
    float mu_max=180;
    float absorb;
    float emi;
    float tau = 0.1;
    float pi = 3.14159265;
    float B[nlayer];
    float lambda;
    float h_planck;
    float c_light;
    float k_boltzman;
    float delta_lambda;
    
                

    //printf("class for convetion\n");
    //printf("nlayer=%d, nlevel=%d\n", nlayer, nlevel);

    //1. Initial variables
    //1.1 Create pressure gradient, loop time = nlevel
    for(int i=0; i<nlevel+1; i++){
        p[i]=p_ground-(nlevel-i)*p_delta; //1000hPa - n*100hPa
        //printf("%dpressure=%f\n",i,p[i]);
    }
    //1.2 Create temperature gradient, loop time = nlayer
    for(int i=0; i<nlayer+1; i++){ 
        T[i]= 293.0-(i)*5.0; //-5K each layer,to ground temp 293K
        //printf("T=%f, i=%d \n",T[i],i);
        theta[i]= T2theta(T[i],p_ground,(p[i]+p[i+1])/2);
        //theta[i]=T[i]*pow((p_ground/p_middle),adiabatic_constant); //gamma?
        

        printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
    }
        //test
        T[2]=T[3]-100;
        T[5]=T[6]-100.0;
        printf("\n");
    for(int i=0; i<nlayer+1; i++){ 
        printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
        printf("--------------------------------\n");
    
    }

    //2. Time loop
    for(int it=0; it<heat_time+1; it++){
        
        //printf("timeloop,it=%d \n",it);
        //2.1 calculate Theta
        for (int i=0; i<nlevel; i++){ //runs nlevel times
            theta[i] = T2theta(T[i],p_ground,(p[i]+p[i+1])/2);
            //printf("hi2\n");
            //printf("i is %d\n", i);
        }
            //printf("i is", i);
        //2.2 Convection
        //unstable if 0

        while(unstable == 1){
            unstable = 0; //stable assumption
            for(int i=1; i<nlevel; i++){ //till which i?
                //print info
                //for(int i=0; i<nlayer+1; i++){ 
                //    printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
                //}
                //printf("\n");
                if(theta[i-1]<theta[i]){ //looking find a problem
                    //exchange
                    float temp_cache = 0.0; 
                    temp_cache = theta[i-1];
                    theta[i-1] = theta[i];
                    theta[i] = temp_cache; 
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
        //2.3 new Theta to temperature
        for (int i=0; i<nlevel+1; i++){                
            T[i] = theta2T(theta[i],p_ground, (p[i]+p[i+1])/2);
            //print info
            //for(int i=0; i<nlayer+1; i++){ 
            //    printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
            //}
            //printf("\n");
                    
        }
        
        //2.4 Change phase ####### (Radiation)
        //radiation
        //2.4.1 Calculate E (naive cloud radiation)
        for(int i=0; i<nlevel+1; i++){
            E_up[i]=E_up[i+1]*(1+alpha)+epsilon*sigma_b*pow(T[i],4);
            E_down[i+1]=epsilon*sigma_b*pow(T[i+1],4)+(1-alpha)*E_down[i];
        }

        //2.4.2 Calculate better E (better cloud radiation)
        //for wavelength-> for directions-> for ilayer
        for(int ilambda=0; ilambda<3; ilambda++){
                delta_lambda = 2*pow(10,-6);
            if(ilambda>0&& ilambda<2){
                tau = 0;
                lambda = 10*pow(10,-6); //800nm in meter
            }
            else if(ilambda<1){
                tau = 0.1;
                lambda = 8*pow(10,-6);
            }
            else{
                tau = 0.1;
                lambda = 12*pow(10,-6);
            }
            for(int i=0; i<nlayer+1; i++){
                //calculate plack law -> B(T) -> L units
                B[i]=(2*h_planck*c_light*c_light)/(pow(lambda,5)*((exp(h_planck*c_light/(lambda*k_boltzman*T[i])))-1));
                //a function to integrate B
                for(int i=0; i<N; i++){
                    float zz; //geschickter?
                    zz=B[i]*delta_lambda;
                }

            }

            //loop for angle mu
            int imu_max = 10; //split in 10 steps
            for(int imu=0; imu<imu_max; imu++){
                float mu=(1+imu)/mu_max;
                absorb=1.0-exp(-tau/mu);
                emi=absorb;

                //upward L_up
                //(p,T,tau -> E_up, E_dw)
                L_up[nlayer]=(1/pi)*sigma_b*pow(T[nlevel-1],4);
                for(int i; i<nlayer+1; i++){
                    //calculate L
                    
                    //L_up[i]=L_up[i+1]*(1-absorb)+(1/pi)*sigma_b*pow(T[i],4)*emi;
                    L_up[i]=L_up[i+1]*(1-absorb)+B[i];  //
                    E_up[i]=E_up[i]+2*pi*L_up[i]*mu*delta_mu;
                }

                //downward L_down
                //
                L_down[0]=0;
                for(int i=1; i<nlayer+1; i++){
                    L_down[i]=L_down[i-1]*(1-absorb)+(1/pi)*sigma_b*pow(T[i-1],4);
                    E_down[i]=E_down[i]+2*pi*L_down[i]*mu*delta_mu;
                }

                //Radiance to Irradiance
                //E=2*pi*(L(mu)*mu[i]*delta_mu;)
            }
        }

        
        //???
        
        //2.5 Calculate delta_E (at each level) (Radiance)
        for(int i=0; i<nlayer+1; i++){
            delta_E[i] = E_down[i]+E_up[i+1]-E_up[i]-E_down[i+1];
        }
        
        //for (int i=0; i<nlayer; i++){
        //    float delta[nlayer-1]
        //    delta[i]=
        //}

        //heat with 1K
        //T[nlevel] = T[nlevel]+heat_temp; 
        //stability = 1; //make unstable

        //2.6 From delta_E to temperature
        for(int i=0; i<nlayer+1; i++){
            delta_Temp[i]=delta_time*delta_E[i]*g/((p[i]-p[i+1])*C_p);
            T[i]=T[i]+delta_Temp[i];
        }

        //2.7 Tempeature to Theta (loop 2.1)
        //7. loop
        for(int i=0; i<nlayer+1; i++){ 
                printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
            }
            printf("\n");
    }
        

return 0;
    
}
//edit note: 17.Okt.2023
//edit note: 30.Okt.2023
