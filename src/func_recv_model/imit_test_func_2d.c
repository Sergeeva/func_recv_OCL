// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: imit_test_func_2d.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                  
//******************************************************************************

#include <stdio.h>
#include <math.h>

//#include "def_const.h"
#include "def_types.h"

#include "useful_tools.h"
#include "param_imit_2d.h"

#include "imit_test_func_2d.h"

void const_2d(float* node_values, int n_1, int n_2, t_vector_2d_f h, float coeff){

  int k1, k2;

  for(k1 = 0; k1 < n_1; k1++) 
	  for(k2 = 0; k2 < n_2; k2++){
		  *node_values = coeff;
		  node_values++;
	  }

  return;

}

void line_func_2d(float* node_values, int n_1, int n_2, t_vector_2d_f h, t_vector_2d_f coeff[2]){

  int k1, k2;

  for(k1=0; k1<n_1; k1++) 
	  for(k2 = 0; k2<n_2; k2++){
		  *node_values = coeff[1]._1*k1*h._1 + coeff[1]._2*k2*h._2 + coeff[0]._1 + coeff[0]._2;
		  node_values++;
	  }
  return;
}

// a_1*x_1^2 + a_2*x_2^2 + b
void polinomial_func_2_2d_v1(float* node_values, int n_1, int n_2, t_vector_2d_f h, t_vector_2d_f coeff[2]){

  int k1, k2;

  float x1, x2;

  for(k1 = 0; k1 < n_1; k1++)
    for(k2 = 0; k2<n_2; k2++){
      x1 = k1*h._1;
      x2 = k2*h._2;
      *node_values = coeff[1]._1*x1*x1 + coeff[1]._2*x2*x2 + coeff[0]._1 + coeff[0]._2;
	  node_values++;
  }

  return;
}

void reverse_polinomial_func_2_2d_v1(float* node_values, int n_1, int n_2, t_vector_2d_f h, t_vector_2d_f coeff[2]){

  int k1, k2;

  float x1, x2, tmp_val;

  for(k1 = 0; k1 < n_1; k1++)
    for(k2 = 0; k2 < n_2; k2++){
      x1 = k1*h._1;
      x2 = k2*h._2;
      tmp_val = coeff[1]._1*x1*x1 + coeff[1]._2*x2*x2 + coeff[0]._1 + coeff[0]._2;;
      *node_values = 1./tmp_val;
	  node_values++;
  }

  return;
}
