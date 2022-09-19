// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_approximation_2d.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                 
//******************************************************************************

#include <stdio.h>
#include <math.h>

//#include "def_const.h"
#include "def_types.h"

#include "useful_tools.h"

#include "f_kernels_2d.h"
#include "f_approximation_2d.h"

#include "param_imit_2d.h"

float approximation_2d(t_vector_2d_f x, float* node_values, t_vector_2d_f l, t_vector_2d_f h){

  float g = 0.;

  int n_1 = l._1/h._1;
  int n_2 = l._2/h._2;
 
  int k1,k2;

  t_vector_2d_f t;

  for(k1 = 0; k1 < n_1; k1++){
	  for(k2 = 0; k2 < n_2; k2++){
		  t._1 = x._1/h._1 - k1;
		  t._2 = x._2/h._2 - k2;
		  g+=(*node_values)*psi_2_2d(t);
		  node_values++;
	  }
  }

  return g;
}

float approximation_2d_S_smooth(t_vector_2d_f x, float* node_values, t_vector_2d_f l, t_vector_2d_f h, t_vector_2d_int r){

  float g = 0.;

  int n_1 = l._1/h._1;
  int n_2 = l._2/h._2;
 
  int k1,k2;

  t_vector_2d_f t;

  for(k1 = 0; k1 < n_1; k1++){
	  for(k2 = 0; k2 < n_2; k2++){
		  t._1 = x._1/h._1 - k1;
		  t._2 = x._2/h._2 - k2;
		  g+=(*node_values)*psi_2d(r,t);
		  node_values++;
	  }
  }


  return g;
}
