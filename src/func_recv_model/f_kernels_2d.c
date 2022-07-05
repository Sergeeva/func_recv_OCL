// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_kernels_2d.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                  
//******************************************************************************

#include <stdio.h>
#include <math.h>

//#include "def_const.h"
#include "def_types.h"

#include "useful_tools.h"

#include "f_kernels_1d.h"
#include "f_kernels_2d.h"


float psi_2_2d(t_vector_2d_f t){
  float psi = psi_2_1d(t._1)*psi_2_1d(t._2);
  return psi;
}


float psi_4_2d(t_vector_2d_f t){
  float psi = psi_4_1d(t._1)*psi_4_1d(t._2);
  return psi;
}


float psi_2d(t_vector_2d_int r, t_vector_2d_f t){
  float psi = psi_1d(r._1,t._1)*psi_1d(r._2,t._2);
  return psi;
}

