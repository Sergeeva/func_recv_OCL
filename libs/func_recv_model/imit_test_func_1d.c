// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: imit_test_func_1d.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                   
//******************************************************************************

#include <stdio.h>
#include <math.h>

//#include "def_const.h"
#include "def_types.h"

#include "useful_tools.h"
#include "param_imit_1d.h"

#include "imit_test_func_1d.h"

void const_1d(float* node_values, int n_nodes, float h, float coeff){

  int k;

  for(k=0; k<n_nodes; k++) node_values[k] = coeff;

  return;

}

void line_func_1d(float* node_values, int n_nodes, float h, float coeff[2]){

  int k;

  for(k=0; k<n_nodes; k++) node_values[k] = coeff[1]*k*h + coeff[0];

  return;
}

void polinomial_func_2_1d(float* node_values, int n_nodes, float h, float coeff[3]){

  int k;

  float x;

  for(k = 0; k < n_nodes; k++){
    x = (k - n_nodes/4)*h;
    node_values[k] = coeff[2]*x*x + coeff[1]*x + coeff[0];
  }

  return;
}

void polinomial_func_3_1d(float* node_values, int n_nodes, float h, float coeff[4]){

  int k;

  float x;

  for(k = 0; k < n_nodes; k++){
    x = (k)*h;
    node_values[k] = coeff[3]*x*x*x + coeff[2]*x*x + coeff[1]*x + coeff[0];
  }

  return;
}

void polinomial_func_1d(float* node_values, int n_nodes, float h, float coeff[P+1], int polinom_size){

  int k,p;

  float x;

  for(k = 0; k < n_nodes; k++){
    x = (k)*h;
    node_values[k] = coeff[0];
    for(p = polinom_size; p>0; p--)
      node_values[k] += coeff[p]*power(x,p);
  }

  return;
}

void reverse_polinomial_func_2_1d(float* node_values, int n_nodes, float h, float coeff[3]){

  int k;

  float x, tmp_val;

  for(k = 0; k < n_nodes; k++){
    x = (k-n_nodes/4)*h;
    tmp_val = coeff[2]*x*x + coeff[1]*x + coeff[0];
    node_values[k] = 1./tmp_val;
  }

  return;
}

void reverse_polinomial_func_3_1d(float* node_values, int n_nodes, float h, float coeff[4]){

  int k;

  float x, tmp_val;

  for(k = 0; k < n_nodes; k++){
    x = (k-n_nodes/2)*h;
    tmp_val = coeff[3]*x*x*x + coeff[2]*x*x + coeff[1]*x + coeff[0];
    node_values[k]  = 1./tmp_val;
  }

  return;
}

void reverse_polinomial_func_1d(float* node_values, int n_nodes, float h, float coeff[P+1], int polinom_size){

  int k,p;

  float x, tmp_val;

  for(k = 0; k < n_nodes; k++){
    x = (k-n_nodes/2)*h;
    node_values[k] = coeff[0];
    for(p = polinom_size; p>0; p--)
      tmp_val += coeff[p]*power(x,p);
    node_values[k] = 1./tmp_val;
  }

  return;
}

