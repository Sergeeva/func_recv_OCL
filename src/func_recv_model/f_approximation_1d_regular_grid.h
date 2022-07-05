// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_approximation_1d.h 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                 
//******************************************************************************

#ifndef F_APPROXIMATION_1D_REGULAR_GRID_H
#define F_APPROXIMATION_1D_REGULAR_GRID_H

#include "def_types.h"

/*
    Инициализация 
    N значений ядер ПОРЯДКА 2 
    в узлах мелкой сетки 
    (на единичном интервале крупной)
    
    Функция возвращает:
    массив kernVal размером N
 */
void kernVals_2_1d(float* kernVals, int N);

/*
    Инициализация 
    N значений ядер ПОРДКА 4 
    в узлах мелкой сетки 
    (на единичном интервале крупной)
    
    Функция возвращает:
    массив kernVal размером N+1
 */
void kernVals_4_1d(float* kernVals, int N);

/*
    Вычисление произведений 
    данного значения функции 
    в узле крупной сетки на 
    значения ядер  
    в узлах мелкой сетки 
    (на единичном интервале крупной)
    
    Функция возвращает:
    массив произведений kernValsxNode 
    размером 2N 
    (значения)
 */

void kernVals_Node_mul(float* kernVals, float Node, float* kernValsxNode, int Nn);

/*
    А1: восстановление данных без сглаживания
    Вычислительная процедура 1 
 */
void F1_Mul_2_1d(float* Nodes, //исходные значения в узлах крупной сетки
		 float* kernVals2, //значения ядер Стеклова порядка 2
		 float* prodVals_buffer2,
		 unsigned short flag_extra_node,
		 float* prodVals_extra,
		 unsigned int part_K,
		 unsigned int N);


/*
    А1: восстановление данных без сглаживания
    Вычислительная процедура 2
 */
void F2_MulSum_2_1d(float* Nodes, //исходные значения в узлах крупной сетки
		    float* kernVals2, //значения ядер Стеклова порядка 2
		    float* prodVals_next,
		    unsigned short flag_last,   
		    float* prodVals_buffer2, //иниц. значениями prodVals из F1
		    float* gridVals, //результат - значения в узлах мелкой сетки
		    unsigned int part_K,
		    unsigned int N);


/*
    А2: восстановление данных co сглаживанием
    Вычислительная процедура 1 
 */
void F1_Mul_4_1d(float* Nodes, //исходные значения в узлах крупной сетки
		 float* kernVals4, //значения ядер Стеклова порядка 4
		 float* prodVals_buffer4, //произведения с первыми 2 узлами
		 float* prodVals_last, //произведения с последним узлом
		 unsigned short flag_extra_node,
		 float* prodVals_extra_left,
		 float* prodVals_extra_right,
		 unsigned int part_K,
		 unsigned int N);


/*
    А2: восстановление данных co сглаживанием
    Вычислительная процедура 2 
 */
void F2_MulSum_4_1d(float* Nodes, //исходные значения в узлах крупной сетки
		    float* kernVals4, //значения ядер Стеклова порядка 4
		    float* prodVals_prev,
		    float* prodVals_next,
		    unsigned short flag_last,  
		    float* prodVals_last,
		    float* prodVals_buffer4, //иниц. значениями prodVals_first из F1
		    float* gridVals, //результат - значения в узлах мелкой сетки
		    unsigned int part_K,
		    unsigned int N);




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// UNUSED!!

/*
    Супершаг 2: 

    Вычисление приближённых значения  
    в узлах мелкой сетки со сглаживанием 
   
!!
ТРЕБУЕТ ХРАНЕНИЯ 
ДВУХ МЕЛКИХ СЕТОК
!!

Используются ядра порядка 4

    Входные данные: 
    kernValsxNode 
    kernValsxNode_prev --- размером N
    kernValsxNode_next --- размером N
    
    Функция возвращает:
    двумерный массив приближённых
    gridVals  
    размером partK на N

    где partK --- количество обрабатывемых узлов крупной сетки, 
        partK = K --- обрабатывается вся исходная сетка;

        N --- количество точек мелкой сетки в одном интервале крупной
 */

void SS2_sumResult_4_1d(float* kernValsxNode,
			float* kernValsxNode_prev, float* kernValsxNode_next,
			float* gridVals, int part_K, int N);



#endif 
