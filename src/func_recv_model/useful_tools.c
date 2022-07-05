// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: useful_tools.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                  
//******************************************************************************

#include <stdio.h>
#include <math.h>

#include "def_types.h"

#include "useful_tools.h"

// вычисление факториала по определению
int factorial(int x){

  int i;
  int f = 1;
 
  for(i=1; i<=x; i++) f*=i;

  return f;
 
 }

//вычисление биномиального коэф-та по определению
int binomial_coeff(int n, int k){

  int c;

  c = factorial(n)/(factorial(k)*factorial(n-k));

  return c;
}

//возведение в степень
float power(float x, int n){
  float p = 1.;
  int i;

  if(n>=0)
    for (i=1; i<=n; i++) p *= x;
    
  return p;
}
