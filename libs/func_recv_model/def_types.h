// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_def_types.h 				       
// Назначение:			
// Разработчик: Сергеева Е.И.
//******************************************************************************

#ifndef _DEF_TYPES_H_
#define _DEF_TYPES_H_



typedef struct {
  float _1;
  float _2;
} t_vector_2d_f;

typedef struct {
  int _1;
  int _2;
} t_vector_2d_int;



typedef struct {
	float x1;
	float x2;
} t_point;

typedef struct {
	t_point x;
	float val;
} t_node;

typedef struct {

} t_lattice;


typedef struct {

} t_rectangle;

//------------------------------------
// Комплексное число
//------------------------------------
typedef struct 
{
	float re;
	float im;
}
  complex_float;


//------------------------------------
// Комплексное число
//------------------------------------
typedef struct 
{
	float r;
	float i;
}
  cfloat;


#endif  /* _DEF_TYPES_H_ */
