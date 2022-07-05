// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: imit_test_func_1d.h 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                  
//******************************************************************************

#ifndef IMIT_TEST_FUNC_1D_H
#define IMIT_TEST_FUNC_1D_H

//#include "def_const.h"
#include "def_types.h"

#include "param_imit_1d.h"

void const_1d(float* node_values, int n_nodes, float h, float coeff);
void line_func_1d(float* node_values, int n_nodes, float h, float coeff[2]);

void polinomial_func_2_1d(float* node_values, int n_nodes, float h, float coeff[3]);
void polinomial_func_3_1d(float* node_values, int n_nodes, float h, float coeff[4]);
void polinomial_func_1d(float* node_values, int n_nodes, float h, float coeff[P+1], int polinom_size);

void reverse_polinomial_func_2_1d(float* node_values, int n_nodes, float h, float coeff[3]);
void reverse_polinomial_func_3_1d(float* node_values, int n_nodes, float h, float coeff[4]);
void reverse_polinomial_func_1d(float* node_values, int n_nodes, float h, float coeff[P+1], int polinom_size);

#endif 
