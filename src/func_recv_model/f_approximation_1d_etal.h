// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_approximation_1d.h 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                 
//******************************************************************************

#ifndef F_APPROXIMATION_1D_H
#define F_APPROXIMATION_1D_H

#include "def_types.h"

float approximation_1d(float x, float* node_values, int l, float h);

float approximation_1d_S_smooth(float x, float* node_values, int l, float h, int r);

void interpolate_signal_1d(float* in_sig, //отсчеты входного сигнала (предыдущая+текущая+ новая реализации)
	                       float* out_sig, // отсчеты выходного сигнала (текущая реализация)
                           int n_in_nodes, //кол-во отсчетов входного сигнала (в одной входной реализации)
						   int n_out_nodes,  // кол-во отсчетов (в одной выходной реализации)
                           int interpolate_coeff //коэф-т повышения fd в целое число раз
	);

//повышение частоты дискретизации входного сигнала в целое число раз
//со сглаживанием
void interpolate_signal_1d_smooth(float* in_sig, //отсчеты входного сигнала (предыдущая+текущая+ новая реализации)
	                              float* out_sig, // отсчеты выходного сигнала (текущая реализация)
                                  int n_in_nodes, //кол-во отсчетов входного сигнала (в одной входной реализации)
						          int n_out_nodesy,  // кол-во отсчетов (в одной выходной реализации)
                                  int interpolate_coeff, //коэф-т повышения fd в целое число раз
                                  float smooth // парамерт сглаживания
	);

#endif 
