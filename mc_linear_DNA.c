#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main() {
    int MC = 10000000;  // no of monte carlo moves per simulation
    int sim = 10000;  // No of independent simulations
    int u = 1;
    double k_on = 0.01;
    double k_off = 1;
    double s1 = u + u + k_off;
    double s2 = k_on;
    double Time;
    int NB = 15;
    int TS = 8;

    srand(time(NULL));  // Initialize random seed

    //for (int k=2;k<=1002;k+=50) {
    //NB=k;TS=NB/2;Time=0;
    
    for (int j = 0; j < sim; j++) {
        int pos = 0;  // position of protein, pos = 0 indicates protein is in the solution and not on the DNA
        double time = 0;

        for (int i = 0; i < MC; i++) {
            double r1 = (double)rand() / RAND_MAX;
            double S = r1 * s1;


            if (pos == TS) {  // target site is at 30th position
	      break;
            }
	    
            if (pos == 0) {
	      pos = rand() % NB + 1;
	      double r2 = (double)rand() / RAND_MAX;
	      time += (1 / s2)*log(1/r2);
	    }
	    
	    else if (pos>0) {
	      
	      if (S <= u && pos < NB) {
                pos += 1;
	      } else if (S <= u && pos == NB) {
                pos = NB-1;
	      }
	      
	      else if (S > u && S <= (u+u) && pos > 1) {
                pos -= 1;
	      } else if (S > u && S <= (u+u) && pos == 1) {
                pos += 1;
	      }
	      
	      else if (S > (u+u) && S <= (u + u + k_off)) {
                pos = 0;
	      }

	      time += (1 / s1)*log(1/r1);
	    }
	    
	}
	    Time += time;
	
    }
    
    double av_Time = Time / sim;
    printf("%d  %lf\n", NB,av_Time);
    //}
    return 0;
}

