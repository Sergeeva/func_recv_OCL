// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_approximation_1d.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.
//******************************************************************************

#include <stdio.h>
#include <math.h>

//#include "def_const.h"
#include "def_types.h"

#include "param_imit_1d.h"

#include "useful_tools.h"

#include "f_kernels_1d.h"
#include "f_approximation_1d_etal.h"

float approximation_1d(float x, float* node_values, int l, float h){

  float g = 0.;

//  int n_nodes = (l/h) + 1;
  int n_nodes = 3*l;

  int k;
  for(k=0; k<n_nodes; k++){
    g+=node_values[k]*psi_2_1d(x/h - k);
  }

  return g;
}

float approximation_1d_S_smooth(float x, float* node_values, int l, float h, int r){

  float g = 0.;

//  int n_nodes = (l/h) + 1;
  int n_nodes = 3*l;

  int k;
  for(k=0; k<n_nodes; k++){
    g+=node_values[k]*psi_1d(r+2, x/h - k);
  }

  return g;
}

//повышение частоты дискретизации входного сигнала в целое число раз
void interpolate_signal_1d(float* in_sig, //отсчеты входного сигнала (предыдущая+текущая+ новая реализации)
	                       float* out_sig, // отсчеты выходного сигнала (текущая реализация)
                           int n_in_nodes, //кол-во отсчетов входного сигнала (в одной входной реализации)
						   int n_out_nodes,  // кол-во отсчетов (в одной выходной реализации)
                           int interpolate_coeff //коэф-т повышения fd в целое число раз
                          ){

	int i;
	float step = 1./interpolate_coeff;

	for(i=0; i<n_out_nodes; i++){
		out_sig[i] = approximation_1d((n_in_nodes+i)*step, in_sig, n_in_nodes, 1.);
	}

	return;
}

//повышение частоты дискретизации входного сигнала в целое число раз
//со сглаживанием
void interpolate_signal_1d_smooth(float* in_sig, //отсчеты входного сигнала (предыдущая+текущая+ новая реализации)
	                              float* out_sig, // отсчеты выходного сигнала (текущая реализация)
                                  int n_in_nodes, //кол-во отсчетов входного сигнала (в одной входной реализации)
						          int n_out_nodes,  // кол-во отсчетов (в одной выходной реализации)
                                  int interpolate_coeff, //коэф-т повышения fd в целое число раз
                                  float smooth // парамерт сглаживания
                          ){
	int i;
	float step = 1./interpolate_coeff;

	for(i=0; i<n_out_nodes; i++){
		out_sig[i] = approximation_1d_S_smooth((n_in_nodes+i)*step, in_sig, n_in_nodes, 1., smooth);
	}


	return;
}
