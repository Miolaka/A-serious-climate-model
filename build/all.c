#include <stdio.h>
#include <math.h>

//#include "functions.c" //maybe not required
//header
double T2theta (double T, double p_ground, double p_middle);
double theta2T (double theta, double p_ground, double p_middle);
int stability(double theta_i, double theta_j);
double B_plank(double h_planck, double c_light, double lambda, double k_boltzman, double T);


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

//#include "functions.c" //maybe not required

int main(int argc, char **argv)
{   
    double E_0 = 1361.0; //solar irradiance E_0 [W/m2]
    double albedo = 0.30;   //albedo
    int nlayer=30; //10 elemente 0-9, p[9]<p[nlayer]
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
    double sigma_b = 5.67*pow(10.0,-8.0); //5.67*10^-8
    double E_up[nlevel];
    double E_down[nlevel];
    double delta_E[nlevel];
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

    //1. Initial variables, create pressure gradient, loop time = nlevel
    printf("test B is %f \n\n", B_plank(h_planck, c_light, 0.0001, k_boltzman, 300.0));

    for(int i=0; i<nlevel; i++){
        p[i]=p_ground-(nlayer-i)*p_delta;   //1000hPa - n*100hPa
    }

    //1.2 Create temperature gradient, loop time = nlayer
    printf("--------------------------------\nInitiating temperature \n");
    for(int i=0; i<nlayer; i++){ 
        T[i]= 180 + (i)*5.0;    //-5K each layer,to ground temp 293K
        theta[i]= T2theta(T[i],p_ground,(p[i]+p[i+1])/2);
        printf("%d, theta=%.2f, T=%.2f, p_middle=%.2f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
    }
    printf("done \n--------------------------------\n");

    //2. Time loop
    for(int it=0; it<ntime; it++){
        //2.1 calculate Theta
        printf("Iteration: %d, \n--------------------------------\n recalculating Theta temperature \n", it);
        for (int i=0; i<nlayer; i++){   //runs nlevel times
            theta[i] = T2theta(T[i],p_ground,(p[i]+p[i+1])/2.0);
        }
        printf("done \n--------------------------------\n");
    
        //2.2 Convection //unstable if 0
        unstable = 1;
        while(unstable == 1){
            unstable = 0;   //stable assumption
            for(int i=1; i<nlayer; i++){
                if(theta[i-1]<theta[i]){
                    //layer exchange
                    double temp_cache = theta[i-1];
                    theta[i-1] = theta[i];
                    theta[i] = temp_cache; 
                    temp_cache = 0.0; 
                    unstable = 1;
                }
            }
        }

        //2.3 new Theta to temperature
        for (int i=0; i<nlayer; i++){                
            T[i] = theta2T(theta[i],p_ground, (p[i]+p[i+1])/2.0);
        }    
        printf("\n");

        //2.4 Radiative Transfer
        //2.4.1 Calculate E (naive cloud radiation)  
        //2.4.2 Calculate better E (better cloud radiation)

        //for wavelength-> for directions-> for ilayer
        //lambda from 1 micro to 28 micro
        //add radiative transfer window

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
            //layer loop, calculate Radiance B, irradiance E_up_B
            double L_B[nlayer];
            for(int i=0; i<nlayer; i++){
                //calculate plack law -> B(T) -> L units
                //a function to integrate B
                B[i] = B_plank(h_planck, c_light, lambda, k_boltzman, T[i]);
                L_B[i] = B[i]*delta_lambda;
            }
            //loop for angle mu
            //set initial array element value = 0
            for(int k=0; k<nlevel; k++){
                E_up[k] = 0.0;
                E_down[k] = 0.0;
            }
            for(int imu=0; imu<imu_max; imu++){
                double mu= (1.0+imu)/imu_max - delta_mu/2;  //verschiebung integral punkte
                absorb=1.0-exp(-tau/mu); //works
                emi=absorb;

                //upward L_up
                L_up[nlevel-1]=L_B[nlayer-1]; //the actually highest defined elements
                for(int i=nlevel-2; i>-1; i--){
                    L_up[i]=L_up[i+1]*(1-absorb)+(1/pi)*L_B[i]*emi;
                    E_up[i]=E_up[i]+2*pi*L_up[i]*mu*delta_mu;
                }

                //downward L_down
                L_down[0]=0;
                for(int i=1; i<nlevel; i++){
                    L_down[i]=L_down[i-1]*(1-absorb)+(1.0/pi)*L_B[i-1]*emi;
                    E_down[i]=E_down[i]+2.0 *pi*L_down[i]*mu*delta_mu;
                }
                //Radiance to Irradiance
            }
        }
        
        //2.5 Calculate delta_E (at each level) (Radiance)
        for(int i=0; i<nlayer; i++){
            delta_E[i] = E_down[i]+E_up[i+1]-E_up[i]-E_down[i+1];
            //there's a problem at ground layer?
                    //Solar radiation
                    if (i==nlayer-1){
                        printf("E_0 = %f\n", E_0/4.0/(1-albedo));
                        delta_E[i] = delta_E[i] + E_0/4.0/(1-absorb); //yet A is now calculated with window
                    }
        }

        //2.6 From delta_E to temperature
        for(int i=0; i<nlayer; i++){
            delta_Temp[i]=delta_time*delta_E[i]*g/((p[i+1]-p[i])*C_p);
            T[i]=T[i]+delta_Temp[i];
            //printf("-----------------------\n delta T[0] = %f\n", delta_Temp[0]);
            //printf("delta E[0] = %f\n", delta_E[0]);
            printf("p[i+1] = %f\n", p[i+1]);
            printf("p[i]-p[i+1]) = %f\n-----------------------\n", p[i+1]-p[i]);
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