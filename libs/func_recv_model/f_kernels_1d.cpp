// -*- coding: utf-8 -*-

#include <stdio.h>
#include <math.h>

//#include "def_const.h"
#include "def_types.h"

#include "useful_tools.h"

#include "f_kernels_1d.h"


float psi_2_1d(float t){

  float psi = 0.;
  float abs_t = fabs(t);

  if(abs_t<=1){
    psi = 1.-abs_t;
  }

  return psi;
}

float psi_4_1d(float t){

  float psi = 0.;
  float abs_t = fabs(t);

  if (abs_t<=1){
    psi = (3*abs_t*abs_t*abs_t - 6*abs_t*abs_t+4)/6.;
  }
  if(abs_t>1 && abs_t<=2){
    psi = ((-1)*abs_t*abs_t*abs_t+6*abs_t*abs_t-12*abs_t+8)/6.;
  }

  return psi;
}

float psi_1d(int r, float t){

  float psi = 0.;
  float abs_t = fabs(t);

  float r_2 = r/2.;

  int k;

  float arg;
  int sign = 1;

  if(abs_t <= r_2){
    for(k=0; k<(abs_t+r_2); k++){
      arg = abs_t+r_2-k;
      psi += sign*binomial_coeff(r,k)*power(arg, r-1);
      sign*=-1; 
    }
    psi*=(1./factorial(r-1));
  }

  return psi;
}
