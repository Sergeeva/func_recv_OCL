// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_approximation_2d.h 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                  
//******************************************************************************

#ifndef F_APPROXIMATION_2D_H
#define F_APPROXIMATION_2D_H

#include "def_types.h"

float approximation_2d(t_vector_2d_f x, float* node_values, t_vector_2d_f l, t_vector_2d_f h);

float approximation_2d_S_smooth(t_vector_2d_f x, float* node_values, t_vector_2d_f l, t_vector_2d_f h, t_vector_2d_int r);

#endif 
