// -*- coding: utf-8 -*-

//******************************************************************************
// Обозначение: print_file_utils.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.                  
//******************************************************************************

#include <stdio.h>
#include <math.h>

#include "def_types.h"

#include "print_file_utils.h"

void print_file(FILE* f, float* data, int n){

  int i;

  for(i=0; i<n; i++){
	  fprintf(f, "%f\n", (*data));
	  data++;
  }

  return;
}


