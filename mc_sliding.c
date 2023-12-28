#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main() {
    long int MC = 10000000;  // no of monte carlo moves per simulation
    int sim = 10000;  // No of independent simulations
    int u = 1;
    double u_1 = 10;
    double k_s = 0.0;
    double k_s_1 = 0.0;   
    double k_on = 0.01;
    double k_on_1 = 0.0001;
    double k_off = 1;
    double k_off_1 = 0.01;
    double s_on = k_on + k_on_1;
    double s_l = u + u + k_off + k_s + k_s_1;
    double s_n = u_1 + u_1 + k_off_1 + k_s + k_s_1;
    int NB1 , NB2 , NB3 , NB4;
    int l = 5;
    NB1 =1 ;
    NB2 = 5;
    NB3 = 20;
    NB4 = 25;
  //  double s1 = u + u + k_off;
 // double s2 = k_on;
    double Time;
    //double Time = 0;
    int NB = NB4;
    int TS;
    //int TS = 10;

    srand(time(NULL));  // Initialize random seed

    for (int k=1; k<=NB; k++) {
      TS=10; Time=0;
      
      for (int j = 0; j < sim; j++) {
        int pos = 0;  // position of protein, pos = 0 indicates protein is in the solution and not on the DNA
        double time = 0;

        for (int i = 0; i < MC; i++) {
            double r1 = (double)rand() / RAND_MAX;
            double S_on = r1 * s_on;
            double S_n = r1 * s_n;
            double S_l = r1 * s_l;


            if (pos == TS) {  // target site is at 30th position
	      //printf("reached");
	      break;
            }
	    
            if (pos == 0) {
	      if (S_on <= k_on ) {
		double r2 = rand() % (NB2+(NB4-NB3)) + 1;
		if (r2<=NB2) {pos = r2;}
		else if (r2>NB2) {pos = r2+(NB3-NB2);}
	      }
	      else if (S_on > k_on) {
		double r2 = rand() % (NB3-NB2) + 1;
		pos = r2+NB2; }
	      time += (1 / s_on)*log(1/r1);
	    }
	    
	    else if (pos>0) {
	      if (pos <= NB2 || pos > NB3) {
		
		if (S_l <= u && pos < NB4) {pos += 1;}
		else if (S_l <= u && pos == NB4) {pos = pos-1;}
		else if (S_l > u && S_l <= (u+u) && pos > 1) {pos-=1;}
		else if (S_l > u && S_l <= (u+u) && pos == 1) {pos+=1;}
		else if (S_l > (u+u) && S_l <= (u + u + k_off)) {pos=0;}

		else if (S_l > (u + u + k_off) && S_l <= ( u + u + k_off + k_s )) {
		  if ((NB4-NB3) >= l) {NB2+=l; NB3+=l;}
		  else if ((NB4-NB3) < l) {NB2-=l; NB3-=l;}}

		else if (S_l > ( u + u + k_off + k_s ) && S_l <= s_l) {
		  if (NB2 >= l) {NB2-=l; NB3-=l;}
		  else if (NB2 < l) {NB2+=l; NB3+=l;}}
		
		
		time += (1 / s_l)*log(1/r1);
	      }
	      
	      
	      else if (pos > NB2 && pos <= NB3) {
		if (S_n <= u_1) {pos+=1;}
		else if (S_n > u_1 && S_n <= (u_1+u_1)){pos-=1;}
		else if (S_n > (u_1+u_1) && S_n <= (u_1 + u_1 + k_off_1)) {pos = 0;}
		
		else if (S_n > (u_1 + u_1 + k_off_1) && S_n <= (u_1 + u_1 + k_off_1 + k_s)) {
		  if ((NB4-NB3) >= l) {NB2+=l; NB3+=l;}
		  else if ((NB4-NB3) < l) {NB2-=l; NB3-=l;}}
		
		else if (S_n > (u_1 + u_1 + k_off_1 + k_s) && S_n <= s_n) {
		  if (NB2 >= l) {NB2-=l; NB3-=l;}
		  else if (NB2 < l) {NB2+=l; NB3+=l;}}
		
		time += (1 / s_n)*log(1/r1);
	      }
	    }
	    
	}
	
	Time += time;
      }
      
      double av_Time = Time / sim;
      printf("%d %lf\n",TS, av_Time);
    }
    return 0;
}

