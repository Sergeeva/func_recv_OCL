// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: imit_test_func_2d.h 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                  
//******************************************************************************

#ifndef IMIT_TEST_FUNC_2D_H
#define IMIT_TEST_FUNC_2D_H

//#include "def_const.h"
#include "def_types.h"

#include "param_imit_2d.h"

void const_2d(float* node_values, int n_1, int n_2, t_vector_2d_f h, float coeff);
void line_func_2d(float* node_values, int n_1, int n_2, t_vector_2d_f h, t_vector_2d_f coeff[2]);

void polinomial_func_2_2d_v1(float* node_values, int n_1, int n_2, t_vector_2d_f h, t_vector_2d_f coeff[3]);
void reverse_polinomial_func_2_2d_v1(float* node_values, int n_1, int n_2, t_vector_2d_f h, t_vector_2d_f coeff[3]);

#endif 
