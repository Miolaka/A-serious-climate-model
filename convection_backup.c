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
    float epsilon[nlayer]; //thickness???
    float rho[nlayer]; //density
    float altitude[nlevel]; //Altitude
    float adiabatic_constant = 2.0/7.0; //R_a/C_p
    int unstable=1; // stability condition 1 or 0
    float p_delta = p_ground/nlevel; //100hPa per level
    float p_middle;
    int heat_time = 5;
    float heat_temp = 1;

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
                    for(int i=0; i<nlayer+1; i++){ 
                        printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
                    }
                    printf("\n");
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
            for(int i=0; i<nlayer+1; i++){ 
                printf("%d, theta=%f, T=%f, p_middle=%f\n", i, theta[i], T[i], (p[i]+p[i+1])/2);
            }
            printf("\n");
                    
        }
        
        //2.4 Change phase ####### (Radiation)
        //radiation
        //2.4.1 Calculate E (naive cloud radiation)
        //2.4.2 Calculate better E (better cloud radiation)

        
        //2.5 Calculate delta_E (at each level) (Radiance)
        
        //for (int i=0; i<nlayer; i++){
        //    float delta[nlayer-1]
        //    delta[i]=
        //}

        //heat with 1K
        //T[nlevel] = T[nlevel]+heat_temp; 
        //stability = 1; //make unstable

        //2.6 From delta_E to temperature

        //2.7 Tempeature to Theta (loop 2.1)
        //7. loop
    }
        

return 0;
    
}
//edit note: 17.Okt.2023
//edit note: 30.Okt.2023
