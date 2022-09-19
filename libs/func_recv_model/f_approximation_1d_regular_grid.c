// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: f_approximation_1d_regular_grid.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.
//******************************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//#include "def_const.h"
#include "def_types.h"

#include "param_imit_1d.h"

#include "useful_tools.h"

#include "f_kernels_1d.h"
#include "f_approximation_1d_regular_grid.h"

//****************************************************************************
// Обозначения:
// N --- кол-во узлов мелкой сетки в единичном интервале крупной
// K --- кол-во узлов крупной сетки (в которых известны значения ф-ии)
// шаг крупной сетки принят единичный (h=1)
//****************************************************************************

/*
    Инициализация 
    N значений ядер ПОРЯДКА 2 
    в узлах мелкой сетки 
    (на единичном интервале крупной)
    
    Функция возвращает:
    массив kernVal размером N
 */
void kernVals_2_1d(float* kernVals, int N){

  char* fname = "kernVals_2_1d";
  
  for(int i=0; i<(N+1); i++){
    kernVals[i] = psi_2_1d(i*(1./N));
  }
  
  printf("%s end \n", fname);
  
  return;
}
 

/*
    Инициализация 
    N значений ядер ПОРДКА 4 
    в узлах мелкой сетки 
    (на 2х интервалах крупной)
    
    Функция возвращает:
    массив kernVal размером N+1
 */
void kernVals_4_1d(float* kernVals, int N){

  char* fname = "kernVals_4_1d";
  
  for(int i=0; i<2*(N+1); i++){
    kernVals[i] = psi_4_1d(i*(1./N));
  }

  printf("%s end \n", fname);

  return;
}



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

void kernVals_Node_mul(float* kernVals, float Node, float* kernValsxNode, int Nn){

  char* fname = "kernVals_Node_mul";
  
  for(int i=0; i<Nn; i++){
    kernValsxNode[i] = kernVals[i]*Node;
  }

  //  printf("%s end \n", fname);
  
  return;
}


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
		 unsigned int N){

  char* fname = "F1_Mul_2_1d";
  
  kernVals_Node_mul (&kernVals2[0], Nodes[0], &prodVals_buffer2[0], (N+1));

  printf("%s end \n", fname);
  
  return;
}

void swap(int* a, int* b){
  int tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

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
		    unsigned int N){
  char* fname = "F2_MulSum_2_1d";
  
  int ind_k0 = 0;
  int ind_k1 = 1;
  int N1 = N+1;
  int tmp;
  
  kernVals_Node_mul (&kernVals2[0], Nodes[1], &prodVals_buffer2[ind_k1*N1], N1);

  for(int n = 0; n<N1; n++){
    gridVals[n] = prodVals_buffer2[n] +  prodVals_buffer2[N1 + (N-n)];
    //    printf("[%i] %f, %f \n", n, prodVals_buffer2[ind_k0*N1 + n], prodVals_buffer2[ind_k1*N1 + (N-n)]);
  }
  
  for (int k=1; k<(part_K-1); k++){
    swap(&ind_k0, &ind_k1);
    
    //    printf("%i, %i \n", ind_k0, ind_k1);
    
    kernVals_Node_mul (&kernVals2[0], Nodes[k+1], &prodVals_buffer2[ind_k1*N1], N1);
    for(int n = 0; n<N1; n++){
      //      printf("[%i] %f, %f \n", n, prodVals_buffer2[ind_k0*N1 + n], prodVals_buffer2[ind_k1*N1 + (N-n)]);
      gridVals[k*N1 + n] = prodVals_buffer2[ind_k0*N1 + n] + prodVals_buffer2[ind_k1*N1 + (N-n)];
    }
  }

  //интервал ПОСЛЕ последней точки --- если есть данные next
  if(flag_last!=0){
    swap(&ind_k0, &ind_k1);

    for(int n = 0; n<=N; n++)
      gridVals[part_K*N1 + n] = prodVals_buffer2[ind_k0*N1 + n] + prodVals_next[(N-n)];
  }
  
  printf("%s end \n", fname);
  
  return;
}


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
		 unsigned int N){

  char* fname = "F1_Mul_4_1d";

  int N1 = N+1;
  
  for (int k = 0; k<2; k++){
    kernVals_Node_mul (&kernVals4[0], Nodes[k], &prodVals_buffer4[(k+1)*2*N1], 2*N1);
  }
  
  kernVals_Node_mul (&kernVals4[0], Nodes[part_K-1], &prodVals_last[0], 2*N1);
  
  
  return;
}


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
		    unsigned int N){
  
  char* fname = "F2_MulSum_4_1d";

  printf("%s start\n", fname);

  int ind_k0, ind_k1, ind_k2, ind_k3, ind_head = 3;

  int N1 = N+1;

  //инициализируем старт кольцевого буфера
  memcpy(&prodVals_buffer4[N1], &prodVals_prev[0], N1*sizeof(float));
  kernVals_Node_mul(&kernVals4[0], Nodes[2], &prodVals_buffer4[3*2*N1], 2*N1);
  
  for (int k=0; k<(part_K-1); k++){
    
    // пересчитываем индексы в кольцевом буфере
    ind_k0 = (ind_head + 4 - 3) % 4;
    ind_k1 = (ind_head + 4 - 2) % 4;
    ind_k2 = (ind_head + 4 - 1) % 4;
    ind_k3 = ind_head;

    //    printf("ind: %i, %i, %i, %i \n", ind_k0, ind_k1, ind_k2, ind_k3);
    //    printf("[%i] head %i:  (%i, %i, %i, %i)\n", k, ind_head, ind_k0, ind_k1, ind_k2, ind_k3);
    
    gridVals[k*N1 + 0] = prodVals_buffer4[ind_k0*2*N1 + N1] +
                         prodVals_buffer4[ind_k1*2*N1] +
                         prodVals_buffer4[ind_k2*2*N1 + N1];
    /*
    printf("[%i] %f, %f, %f \n", 0, prodVals_buffer4[ind_k0*2*N1 + N1],
	                            prodVals_buffer4[ind_k1*2*N1],
	                            prodVals_buffer4[ind_k2*2*N1 + N1]);
    */
    for(int n = 1; n<N1; n++){
      gridVals[k*N1 + n] = prodVals_buffer4[ind_k0*2*N1 + N1 + n] +
	                   prodVals_buffer4[ind_k1*2*N1 + n] +
	                   prodVals_buffer4[ind_k2*2*N1 + (N1 - n)] +
	                   prodVals_buffer4[ind_k3*2*N1 + (2*N1 - n)];
      /*
      printf("[%i] %f, %f, %f, %f \n", n, prodVals_buffer4[ind_k0*2*N1 + N1 + n],
	                            prodVals_buffer4[ind_k1*2*N1 + n],
	                            prodVals_buffer4[ind_k2*2*N1 + (N1 - n)],
	                            prodVals_buffer4[ind_k3*2*N1 + (2*N1 - n)]);
      */
    }
    
    //считаем массив произведений для следующей точки
    ind_head = (ind_head+1) % 4;
    if (k == (part_K - 2)){
      memcpy(&prodVals_buffer4[ind_head*2*N1], &prodVals_next[0], 2*N1*sizeof(float));
    }else if (k == (part_K - 3)){
      memcpy(&prodVals_buffer4[ind_head*2*N1], &prodVals_last[0], 2*N1*sizeof(float));
    }else{
      kernVals_Node_mul(&kernVals4[0], Nodes[k+3], &prodVals_buffer4[ind_head*2*N1], 2*N1);
    }    
  }

  //интервал ПОСЛЕ последней точки --- если есть данные next
  if(flag_last!=0){
    // пересчитываем индексы в кольцевом буфере
    ind_k0 = (ind_head + 4 - 3) % 4;
    ind_k1 = (ind_head + 4 - 2) % 4;
    ind_k2 = (ind_head + 4 - 1) % 4;
    ind_k3 = ind_head;
    
    gridVals[part_K*N1 + 0] = prodVals_buffer4[ind_k0*2*N1 + N1] +
                              prodVals_buffer4[ind_k1*2*N1] +
                              prodVals_buffer4[ind_k2*2*N1 + N1];
    
    for(int n = 1; n<N1; n++)
      gridVals[part_K*N1 + n] = prodVals_buffer4[ind_k0*2*N1 + N1 + n] +
	                   prodVals_buffer4[ind_k1*2*N1 + n] +
	                   prodVals_buffer4[ind_k2*2*N1 + (N1 - n)] +
	                   prodVals_buffer4[ind_k3*2*N1 + (2*N1 - n)];
  }
      
  printf("%s end \n", fname);
  
  return;
}



 
