//******************************************************************************
// Обозначение: main.c 				       
// Назначение:			
// Разработчик: Сергеева Е.И.  
//******************************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "def_types.h"
//#include "def_const.h"

#include "param_imit_1d.h"
#include "param_imit_2d.h"

#include "control.h"

#include "useful_tools.h"
#include "print_file_utils.h"

#include "f_kernels_1d.h"
#include "f_kernels_2d.h"

#include "f_approximation_1d_etal.h"
#include "f_approximation_1d_regular_grid.h"

#include "f_approximation_2d.h"

#include "imit_test_func_1d.h"
#include "imit_test_func_2d.h"

/*
#include "param_imit.h"
#include "param_scrk.h"

#include "imit_ton_signal.h"

#include "param_imit.h"
#include "def_types_imit.h"
#include "f_init_imit_data.h"
#include "f_print_imit_data.h"
#include "f_imit_main.h"
*/

/*
Имитируемый сигнал



float sig_in_chd3_cur [SIZE_L][SIZE_FFT/2];
float sig_in_chd3_prev[SIZE_L][SIZE_FFT/2];
 */

/*
Структуры для имитатора

t_str_common_imit_data imit_data_chd3;
 */

/*
float node_values_2d_const          [n_nodes_1][n_nodes_2];
float node_values_2d_line           [n_nodes_1][n_nodes_2];
float node_values_2d_polinom        [n_nodes_1][n_nodes_2];
float node_values_2d_reverse_polinom[n_nodes_1][n_nodes_2];

float etalon_2d_const          [n_test_nodes_1][n_test_nodes_2];
float etalon_2d_line           [n_test_nodes_1][n_test_nodes_2];
float etalon_2d_polinom        [n_test_nodes_1][n_test_nodes_2];
float etalon_2d_reverse_polinom[n_test_nodes_1][n_test_nodes_2];

float result1_2d        [n_test_nodes_1][n_test_nodes_2];
float result1_2d_smooth1[n_test_nodes_1][n_test_nodes_2];
float result1_2d_smooth2[n_test_nodes_1][n_test_nodes_2];
float result1_2d_smooth3[n_test_nodes_1][n_test_nodes_2];

float result2_2d        [n_test_nodes_1][n_test_nodes_2];
float result2_2d_smooth1[n_test_nodes_1][n_test_nodes_2];
float result2_2d_smooth2[n_test_nodes_1][n_test_nodes_2];
float result2_2d_smooth3[n_test_nodes_1][n_test_nodes_2];

float result3_2d        [n_test_nodes_1][n_test_nodes_2];
float result3_2d_smooth1[n_test_nodes_1][n_test_nodes_2];
float result3_2d_smooth2[n_test_nodes_1][n_test_nodes_2];
float result3_2d_smooth3[n_test_nodes_1][n_test_nodes_2];

float result4_2d        [n_test_nodes_1][n_test_nodes_2];
float result4_2d_smooth1[n_test_nodes_1][n_test_nodes_2];
float result4_2d_smooth2[n_test_nodes_1][n_test_nodes_2];
float result4_2d_smooth3[n_test_nodes_1][n_test_nodes_2];

float diff1_2d[n_test_nodes_1][n_test_nodes_2];
float diff2_2d[n_test_nodes_1][n_test_nodes_2];
float diff3_2d[n_test_nodes_1][n_test_nodes_2];
float diff4_2d[n_test_nodes_1][n_test_nodes_2];
*/

int main() {

int i;
int step;
//int  smooth_factor_2, smooth_factor_4, smooth_factor_8;

  //testing kernels_1d and kernels_2d  

#if KERNELS_TEST
  float t = 0.5;
  int r = 8;

  printf("t = %f\t r = %i\n", t, r);

  printf("psi_2_1d = %f\n", psi_2_1d(t));
  printf("psi_4_1d = %f\n", psi_4_1d(t));
  printf("psi_1d = %f\n", psi_1d(r,t));


  t_vector_2d_f t_2;
  t_vector_2d_int r_2;

  t_2._1 = 0.5;
  t_2._2 = 0.5;

  r_2._1 = 8;
  r_2._2 = 8;

  printf("t_2 = (%f,%f)\t r_2 = (%i,%i)\n", t_2._1, t_2._2, r_2._1, r_2._2);

  printf("psi_2_2d = %f\n", psi_2_2d(t_2));
  printf("psi_4_2d = %f\n", psi_4_2d(t_2));
  printf("psi_2d = %f\n", psi_2d(r_2,t_2));

#endif

  //testing approximation_1d

#if APPROXIMATION_1D_TEST

  printf("test approximation 1D\n");

  ///////////////////////////////////////////////////////////
  int K = N_NODES;
  int N = 7; // кол-во интервалов мелкой сетки в крупной
  float test_step = H/N;
  int n_test_nodes = (K-1)*(N+1); //(L/test_step) - L/H;
  ///////////////////////////////////////////////////////////

  float kernVals_2[N+1];
  float kernVals_4[2*(N+1)];
  memset(kernVals_2, 0, (N+1) * sizeof(float));
  memset(kernVals_4, 0, 2*(N+1) * sizeof(float));
  
  float prodVals_next[N+1];
  memset(prodVals_next, 0, (N+1) * sizeof(float));
  
  float prodVals_last[2*(N+1)]; 
  memset(prodVals_last, 0, 2*(N+1) * sizeof(float));
  
  float kernValsxNode_prev[N+1];
  float kernValsxNode_next[2*(N+1)];

  memset(kernValsxNode_prev, 0, (N+1) * sizeof(float));
  memset(kernValsxNode_next, 0, 2*(N+1) * sizeof(float));
  
  //кольцевой буффер
  float buff2[2*(N+1)]; //для хранения
  float buff4[4*2*(N+1)]; //для хранения
  memset(buff2, 0, 2*(N+1) * sizeof(float));
  memset(buff4, 0, 4*2*(N+1) * sizeof(float));
  
  //тестовая функция 1
  float node_values_1d_1[N_NODES]; // известные значения в узлах
  float result_1d_1[n_test_nodes]; //[n_test_nodes];
  float result_1d_1_smooth[n_test_nodes]; //[n_test_nodes];

  memset(node_values_1d_1, 0, N_NODES * sizeof(float));
  memset(result_1d_1, 0, n_test_nodes * sizeof(float));
  memset(result_1d_1_smooth, 0, n_test_nodes * sizeof(float));
  
  //для сравнения
  //  float etalon_1d_1[n_test_nodes];
  //  float approx_etalon_1d_1[n_test_nodes];
  //  float approx_etalon_1d_1_smooth[n_test_nodes];
  
  //тестовая функция 2
  float node_values_1d_2[N_NODES];
  float result_1d_2[n_test_nodes]; //[n_test_nodes];
  float result_1d_2_smooth[n_test_nodes]; //[n_test_nodes];

  memset(node_values_1d_2, 0, N_NODES * sizeof(float));
  memset(result_1d_2, 0, n_test_nodes * sizeof(float));
  memset(result_1d_2_smooth, 0, n_test_nodes * sizeof(float));

  
  //для сравнения 
  //  float etalon_1d_2[n_test_nodes];
  //  float approx_etalon_1d_2[n_test_nodes];
  //  float approx_etalon_1d_2_smooth[n_test_nodes];
  
  //тестовая функция 3
  float node_values_1d_3[N_NODES];
  float result_1d_3[n_test_nodes]; //[n_test_nodes];
  float result_1d_3_smooth[n_test_nodes]; //[n_test_nodes];

  memset(node_values_1d_3, 0, N_NODES * sizeof(float));
  memset(result_1d_3, 0, n_test_nodes * sizeof(float));
  memset(result_1d_3_smooth, 0, n_test_nodes * sizeof(float));

  
  //для сравнения
  //  float etalon_1d_3[n_test_nodes];
  //  float approx_etalon_1d_3[n_test_nodes];
  //  float approx_etalon_1d_3_smooth[n_test_nodes];

  
  //тестовая функция 4
  float node_values_1d_4[N_NODES];
  float result_1d_4[n_test_nodes]; //[n_test_nodes];
  float result_1d_4_smooth[n_test_nodes]; //[n_test_nodes];

  memset(node_values_1d_4, 0, N_NODES * sizeof(float));
  memset(result_1d_4, 0, n_test_nodes * sizeof(float));
  memset(result_1d_4_smooth, 0, n_test_nodes * sizeof(float));

  
  //для сравнения
  //  float etalon_1d_4[n_test_nodes];
  //  float approx_etalon_1d_4[n_test_nodes];
  //  float approx_etalon_1d_4_smooth[n_test_nodes];


  ///////////////////////////////////////////////////
  //  float node_values_1d_ton_sig[N_NODES];
  //  float result_1d_ton_sig[n_test_nodes];
  //  float result_1d_ton_sig_smooth_2[n_test_nodes];

  float coeff1 = 10.;
  float coeff2[2] = {2.,3.};
  float coeff3[3] = {-2., 5., 12.};

  float coeff4[P+1];
  
  for(i=0; i<P+1; i++){
    coeff4[i] =  i*10.+(i+P);
  }
  
  const_1d(node_values_1d_1, N_NODES, H, coeff1);
  //  const_1d(etalon_1d_1, n_test_nodes, test_step, coeff1);

  line_func_1d(node_values_1d_2, N_NODES, H, coeff2);
  //  line_func_1d(etalon_1d_2, n_test_nodes, test_step, coeff2);

  polinomial_func_2_1d(node_values_1d_3, N_NODES, H, coeff3);
  //  polinomial_func_2_1d(etalon_1d_3, n_test_nodes, test_step, coeff3);
  //polinomial_func_3_1d(node_values_1d_3, N_NODES, H, coeff3);
  //polinomial_func_1d(node_values_1d_3, N_NODES, H, coeff3, P);

  reverse_polinomial_func_2_1d(node_values_1d_4, N_NODES, H, coeff3);
  //  reverse_polinomial_func_2_1d(etalon_1d_4, n_test_nodes, test_step, coeff3);
  //reverse_polinomial_func_3_1d(node_values_1d_3, N_NODES, H, coeff3);
  //reverse_polinomial_func_1d(node_values_1d_3, N_NODES, H, coeff3, P);

 
  //////////////////////////////////////////////////////////////////////////////
  //Проверка функций для параллельных алгоритмов

  //Инициализация ядер
  kernVals_2_1d(&kernVals_2[0], N);
  kernVals_4_1d(&kernVals_4[0], N);

  for(i = 0; i<N+1; i++) printf("[%i] %f \n", i, kernVals_2[i]);
  printf("------\n");
  
  for(i = 0; i<2*(N+1); i++) printf("[%i] %f \n", i, kernVals_4[i]);
  printf("------\n");
  
  // Тестовая ф-ия 1 : const
  //
  
  ///////////////////////////////////////  
  //1. Без сглаживания
  ///////////////////////////////////////

  //Обработка по известным значениям в узлах

  // 1. Умножение 
  F1_Mul_2_1d(&node_values_1d_1[0], &kernVals_2[0],
	      &buff2[0],
	      0, NULL, 
	      K, N);
  
  // 2. Суммирование
  F2_MulSum_2_1d(&node_values_1d_1[0], &kernVals_2[0],
		 &prodVals_next[0],
		 0,
		 &buff2[0],
		 &result_1d_1[0], K, N);

  ///////////////////////////////////////
  //2. Со сглаживанием
  ///////////////////////////////////////

  //Обработка по известным значениям в узлах
    
  // 1. Умножение 
  F1_Mul_4_1d(&node_values_1d_1[0], &kernVals_4[0],
		 &buff4[0], &prodVals_last[0],
		 0, NULL, NULL, 
		 K, N);
  
  // 2. Суммирование
  F2_MulSum_4_1d(&node_values_1d_1[0], &kernVals_4[0],
		 &kernValsxNode_prev[0], &kernValsxNode_next[0],
		 0,
		 &prodVals_last[0],
		 &buff4[0],
		 &result_1d_1_smooth[0], K, N);

  
  // Тестовая ф-ия 2 : line
  //

  ///////////////////////////////////////  
  //1. Без сглаживания
  ///////////////////////////////////////
  
  //Обработка по известным значениям в узлах

  
  // 1. Умножение 
  F1_Mul_2_1d(&node_values_1d_2[0], &kernVals_2[0],
	      &buff2[0],
	      0, NULL, 
	      K, N);
  
  // 2. Суммирование
  F2_MulSum_2_1d(&node_values_1d_2[0], &kernVals_2[0],
		 &prodVals_next[0],
		 0,
		 &buff2[0],
		 &result_1d_2[0], K, N);

  ///////////////////////////////////////
  //2. Со сглаживанием
  ///////////////////////////////////////

  //Обработка по известным значениям в узлах
  
  // 1. Умножение 
  F1_Mul_4_1d(&node_values_1d_2[0], &kernVals_4[0],
		 &buff4[0], &prodVals_last[0],
		 0, NULL, NULL, 
		 K, N);
  
  // 2. Суммирование
  F2_MulSum_4_1d(&node_values_1d_2[0], &kernVals_4[0],
		 &kernValsxNode_prev[0], &kernValsxNode_next[0],
		 0,
		 &prodVals_last[0],
		 &buff4[0],
		 &result_1d_2_smooth[0], K, N);
  
  
  
  // Тестовая ф-ия 3 : polinom
  //
  
  ///////////////////////////////////////  
  //1. Без сглаживания
  ///////////////////////////////////////

  //Обработка по известным значениям в узлах

    // 1. Умножение 
  F1_Mul_2_1d(&node_values_1d_3[0], &kernVals_2[0],
	      &buff2[0],
	      0, NULL, 
	      K, N);
  
  // 2. Суммирование
  F2_MulSum_2_1d(&node_values_1d_3[0], &kernVals_2[0],
		 &prodVals_next[0],
		 0,
		 &buff2[0],
		 &result_1d_3[0], K, N);

  ///////////////////////////////////////
  //2. Со сглаживанием
  ///////////////////////////////////////

  //Обработка по известным значениям в узлах
  
  // 1. Умножение 
  F1_Mul_4_1d(&node_values_1d_3[0], &kernVals_4[0],
		 &buff4[0], &prodVals_last[0],
		 0, NULL, NULL, 
		 K, N);
  
  // 2. Суммирование
  F2_MulSum_4_1d(&node_values_1d_3[0], &kernVals_4[0],
		 &kernValsxNode_prev[0], &kernValsxNode_next[0],
		 0,
		 &prodVals_last[0],
		 &buff4[0],
		 &result_1d_3_smooth[0], K, N);



  // Тестовая ф-ия 4 : reverse polinom
  //
  
  ///////////////////////////////////////  
  //1. Без сглаживания
  ///////////////////////////////////////
  
  //Обработка по известным значениям в узлах

      // 1. Умножение 
  F1_Mul_2_1d(&node_values_1d_4[0], &kernVals_2[0],
	      &buff2[0],
	      0, NULL, 
	      K, N);
  
  // 2. Суммирование
  F2_MulSum_2_1d(&node_values_1d_4[0], &kernVals_2[0],
		 &prodVals_next[0],
		 0,
		 &buff2[0],
		 &result_1d_4[0], K, N);

  ///////////////////////////////////////
  //2. Со сглаживанием
  ///////////////////////////////////////

  //Обработка по известным значениям в узлах
  
  // 1. Умножение 
  F1_Mul_4_1d(&node_values_1d_4[0], &kernVals_4[0],
		 &buff4[0], &prodVals_last[0],
		 0, NULL, NULL, 
		 K, N);
  
  // 2. Суммирование
  F2_MulSum_4_1d(&node_values_1d_4[0], &kernVals_4[0],
		 &kernValsxNode_prev[0], &kernValsxNode_next[0],
		 0,
		 &prodVals_last[0],
		 &buff4[0],
		 &result_1d_4_smooth[0], K, N);

  
  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#if PRINT_FILE_1D
  
  FILE *f1_in_data, *f2_in_data, *f3_in_data, *f4_in_data;
  FILE *f1_out_data, *f1_out_data_sm, *f1_out_data_sm_v2;
  FILE *f2_out_data, *f2_out_data_sm, *f2_out_data_sm_v2;
  FILE *f3_out_data, *f3_out_data_sm, *f3_out_data_sm_v2;
  FILE *f4_out_data, *f4_out_data_sm, *f4_out_data_sm_v2;
  FILE *f_etalon_3, *f_etalon_4;

  /*
  FILE *f1_approx_etalon_sm, *f2_approx_etalon_sm, *f3_approx_etalon_sm, *f4_approx_etalon_sm;

  f1_approx_etalon_sm = fopen("approx_etalon_const_1d.dat", "w");
  f2_approx_etalon_sm = fopen("approx_etalon_line_1d.dat", "w");
  f3_approx_etalon_sm = fopen("approx_etalon_polinom_1d.dat", "w");
  f4_approx_etalon_sm = fopen("approx_etalon_reverse_polinom_1d.dat", "w");
  
  print_file(f1_approx_etalon_sm, approx_etalon_1d_1_smooth, n_test_nodes);
  print_file(f2_approx_etalon_sm, approx_etalon_1d_2_smooth, n_test_nodes);
  print_file(f3_approx_etalon_sm, approx_etalon_1d_3_smooth, n_test_nodes);
  print_file(f4_approx_etalon_sm, approx_etalon_1d_4_smooth, n_test_nodes);

  fclose(f1_approx_etalon_sm);
  fclose(f2_approx_etalon_sm);
  fclose(f3_approx_etalon_sm);
  fclose(f4_approx_etalon_sm);
  */
  
  f1_in_data = fopen("node_values_const_1d.dat", "w");
  f2_in_data = fopen("node_values_line_1d.dat", "w");
  f3_in_data = fopen("node_values_polinom_1d.dat", "w");
  f4_in_data = fopen("node_values_reverse_polinom_1d.dat", "w");
  
  print_file(f1_in_data, node_values_1d_1, N_NODES);
  print_file(f2_in_data, node_values_1d_2, N_NODES);
  print_file(f3_in_data, node_values_1d_3, N_NODES);
  print_file(f4_in_data, node_values_1d_4, N_NODES);

  fclose(f1_in_data);
  fclose(f2_in_data);
  fclose(f3_in_data);
  fclose(f4_in_data);

  f1_out_data = fopen("result_const_1d.dat", "w");
  f1_out_data_sm = fopen("result_const_1d_smooth.dat", "w");

  f2_out_data = fopen("result_line_1d.dat", "w");
  f2_out_data_sm = fopen("result_line_1d_smooth.dat", "w");

  f3_out_data = fopen("result_polinom_1d.dat", "w");
  f3_out_data_sm = fopen("result_polinom_1d_smooth.dat", "w");

  f4_out_data = fopen("result_reverse_polinom_1d.dat", "w");
  f4_out_data_sm = fopen("result_reverse_polinom_1d_smooth.dat", "w");

  print_file(f1_out_data, result_1d_1, n_test_nodes);
  print_file(f1_out_data_sm, result_1d_1_smooth, n_test_nodes);

  
  print_file(f2_out_data, result_1d_2, n_test_nodes);
  print_file(f2_out_data_sm, result_1d_2_smooth, n_test_nodes);

  
  print_file(f3_out_data, result_1d_3, n_test_nodes);
  print_file(f3_out_data_sm, result_1d_3_smooth, n_test_nodes);

  
  print_file(f4_out_data, result_1d_4, n_test_nodes);
  print_file(f4_out_data_sm, result_1d_4_smooth, n_test_nodes);


  fclose(f1_out_data);
  fclose(f1_out_data_sm);

  fclose(f2_out_data); 
  fclose(f2_out_data_sm);

  fclose(f3_out_data);
  fclose(f3_out_data_sm);

  fclose(f4_out_data);
  fclose(f4_out_data_sm);
  
#endif
  
  /*
  float x = 1./sqrt(2);
  float test_1 = approximation_1d(x, node_values_1d_4, N_NODES, H);
  float test_2 = approximation_1d_S_smooth(x, node_values_1d_4, N_NODES, H, 2);

   float tmp_val_etalon = coeff3[2]*x*x + coeff3[1]*x + coeff3[0];
   float etalon_val = 1./tmp_val_etalon;

   printf("%f\t%f\t%f\n", test_1, test_2, etalon_val);
  */
  
#endif

   //testing approximation 1d for imit ton and harmonic signals

#if APPROXIMATION_1D_TON_SIGNALS

   printf("test approximation 1D for ton signal\n");

   float node_values_1d_ton_signal[3*N_NODES];

   int interpolate_coeff_1 = 64;


//   int n_test_nodes = (Length/test_step);
   n_test_nodes = interpolate_coeff_1*N_NODES;

   float result_1d_ton_signal[n_test_nodes];
   float result_1d_ton_signal_smooth2[n_test_nodes];
   float result_1d_ton_signal_smooth4[n_test_nodes];
   float result_1d_ton_signal_smooth8[n_test_nodes];


   float etalon_1d_ton_signal[n_test_nodes];
  
  smooth_factor_2 = 2; // параметр сглаживания
  smooth_factor_4 = 4; // параметр сглаживания
  smooth_factor_8 = 8; // параметр сглаживания


  FILE *f1_in_data;
  FILE *f1_out_data, *f1_out_data_sm2, *f1_out_data_sm4, *f1_out_data_sm8;

  FILE *f_etalon;

  f1_in_data = fopen("node_values_ton_signal.dat", "w");

  f1_out_data = fopen("result_ton_signal_1d.dat", "w");
  f1_out_data_sm2 = fopen("result_ton_signal_1d_smooth2.dat", "w");
  f1_out_data_sm4 = fopen("result_ton_signal_1d_smooth4.dat", "w");
  f1_out_data_sm8 = fopen("result_ton_signal_1d_smooth8.dat", "w");

  f_etalon = fopen("etalon_ton_signal_1d.dat", "w");

  ton_signal(&node_values_1d_ton_signal[0],         N_NODES, N_NODES, A1, f1, 6000.0, 0, 0.0, 0);
  ton_signal(&node_values_1d_ton_signal[N_NODES],   N_NODES, N_NODES, A1, f1, 6000.0, 0, 0.0, 1);
  ton_signal(&node_values_1d_ton_signal[2*N_NODES], N_NODES, N_NODES, A1, f1, 6000.0, 0, 0.0, 2);

	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 4*16*6000., 0, 0.0, 0);
//	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 16*6000.,  0,  0.0, 0);
//	   print_file(f_etalon, etalon_1d_ton_signal, n_test_nodes);

	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 4*16*6000., 0, 0.0, 1);
//	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 16*6000.,  0,  0.0, 1);
//       print_file(f_etalon, etalon_1d_ton_signal, n_test_nodes);

	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 4*16*6000., 0, 0.0, 2);
//	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 16*6000.,  0,  0.0, 2);
//       print_file(f_etalon, etalon_1d_ton_signal, n_test_nodes);

   for(step = 3; step < 10; step++){

	   printf("---- %i ----\n", step);

	   memcpy(&node_values_1d_ton_signal[0], &node_values_1d_ton_signal[N_NODES], N_NODES*sizeof(float));
	   memcpy(&node_values_1d_ton_signal[N_NODES], &node_values_1d_ton_signal[2*N_NODES], N_NODES*sizeof(float));
	   ton_signal(&node_values_1d_ton_signal[2*N_NODES], N_NODES, N_NODES, A1, f1, 6000.0, 0, 0.0, step);

	   interpolate_signal_1d(node_values_1d_ton_signal, result_1d_ton_signal, N_NODES, n_test_nodes, interpolate_coeff_1);
	   interpolate_signal_1d_smooth(node_values_1d_ton_signal, result_1d_ton_signal_smooth2, 
                                    N_NODES, n_test_nodes, interpolate_coeff_1, smooth_factor_2);
	   interpolate_signal_1d_smooth(node_values_1d_ton_signal, result_1d_ton_signal_smooth4, 
                                    N_NODES, n_test_nodes, interpolate_coeff_1, smooth_factor_4);
	   interpolate_signal_1d_smooth(node_values_1d_ton_signal, result_1d_ton_signal_smooth8, 
                                    N_NODES, n_test_nodes, interpolate_coeff_1, smooth_factor_8);



	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 4*16*6000., 0, 0.0, step);
//	   ton_signal(etalon_1d_ton_signal, n_test_nodes, n_test_nodes, A1, f1, 16*6000.,  0,  0.0, step);

  print_file(f1_in_data, &node_values_1d_ton_signal[2*N_NODES], N_NODES);

  print_file(f1_out_data, result_1d_ton_signal, n_test_nodes);
  print_file(f1_out_data_sm2, result_1d_ton_signal_smooth2, n_test_nodes);
  print_file(f1_out_data_sm4, result_1d_ton_signal_smooth4, n_test_nodes);
  print_file(f1_out_data_sm8, result_1d_ton_signal_smooth8, n_test_nodes);

  print_file(f_etalon, etalon_1d_ton_signal, n_test_nodes);
  }

  fclose(f1_in_data);

  fclose(f1_out_data);
  fclose(f1_out_data_sm2);
  fclose(f1_out_data_sm4);
  fclose(f1_out_data_sm8);

  fclose(f_etalon);

#endif

#if APPROXIMATION_1D_HARMONIC_SIGNALS

   printf("test approximation 1D for harmonic signal\n");


  FILE *f2_in_data;
  FILE *f2_out_data, *f2_out_data_sm2, *f2_out_data_sm4, *f2_out_data_sm8;

  f2_in_data = fopen("node_values_harmonic_signal.dat", "w");

  f2_out_data = fopen("result_harmonic_signal_1d.dat", "w");
  f2_out_data_sm2 = fopen("result_harmonic_signal_1d_smooth2.dat", "w");
  f2_out_data_sm4 = fopen("result_harmonic_signal_1d_smooth4.dat", "w");
  f2_out_data_sm8 = fopen("result_harmonic_signal_1d_smooth8.dat", "w");



   float node_values_1d_harmonic_signal[3*N_NODES_1024];

   // int n_test_nodes = (Length/test_step);

   int interpolate_coeff_2 = 16;
   n_test_nodes = interpolate_coeff_2*N_NODES_1024;

   float result_1d_harmonic_signal[n_test_nodes];
   float result_1d_harmonic_signal_smooth2[n_test_nodes];
   float result_1d_harmonic_signal_smooth4[n_test_nodes];
   float result_1d_harmonic_signal_smooth8[n_test_nodes];

  smooth_factor_2 = 2; // параметр сглаживания
  smooth_factor_4 = 4; // параметр сглаживания
  smooth_factor_8 = 8; // параметр сглаживания

/*
Инициализация параметров для имитатора
 */
	f_init_imit_data(&imit_data_chd3, 6000.0, F_delta, SIZE_FFT/2, SIZE_L, c0_sound, N_TARGETS);
	f_print_common_imit_data(imit_data_chd3);

	f_imit_main(&sig_in_chd3_cur[0][0], imit_data_chd3, 0);
	memcpy(&node_values_1d_harmonic_signal[0], &sig_in_chd3_cur[0][0], N_NODES_1024*sizeof(float));

	f_imit_main(&sig_in_chd3_cur[0][0], imit_data_chd3, 1);
    memcpy(&node_values_1d_harmonic_signal[N_NODES_1024], &sig_in_chd3_cur[0][0], N_NODES_1024*sizeof(float));

	f_imit_main(&sig_in_chd3_cur[0][0], imit_data_chd3, 2);
    memcpy(&node_values_1d_harmonic_signal[2*N_NODES_1024], &sig_in_chd3_cur[0][0], N_NODES_1024*sizeof(float));


   for(step = 3; step<10; step++){

	   printf("---- %i ----\n", step);

//	   memcpy(sig_in_chd3_prev, sig_in_chd3_cur, SIZE_L*(SIZE_FFT/2)*sizeof(float));

	   f_imit_main(&sig_in_chd3_cur[0][0], imit_data_chd3, step);

	   memcpy(&node_values_1d_harmonic_signal[0], &node_values_1d_ton_signal[N_NODES_1024], N_NODES_1024*sizeof(float));
	   memcpy(&node_values_1d_ton_signal[N_NODES_1024], &node_values_1d_ton_signal[2*N_NODES_1024], N_NODES_1024*sizeof(float));
	   memcpy(&node_values_1d_harmonic_signal[2*N_NODES_1024], &sig_in_chd3_cur[0][0], (SIZE_FFT/2)*sizeof(float));

	   interpolate_signal_1d(node_values_1d_harmonic_signal, result_1d_harmonic_signal, N_NODES_1024, n_test_nodes, interpolate_coeff_2);
	   interpolate_signal_1d_smooth(node_values_1d_harmonic_signal, result_1d_harmonic_signal_smooth2, 
                                    N_NODES_1024, n_test_nodes, interpolate_coeff_2, smooth_factor_2);
	   interpolate_signal_1d_smooth(node_values_1d_harmonic_signal, result_1d_harmonic_signal_smooth4, 
                                    N_NODES_1024, n_test_nodes, interpolate_coeff_2, smooth_factor_4);
	   interpolate_signal_1d_smooth(node_values_1d_harmonic_signal, result_1d_harmonic_signal_smooth8, 
                                    N_NODES_1024, n_test_nodes, interpolate_coeff_2, smooth_factor_8);


  print_file(f2_in_data, node_values_1d_harmonic_signal, N_NODES_1024);

  print_file(f2_out_data, result_1d_harmonic_signal, n_test_nodes);
  print_file(f2_out_data_sm2, result_1d_harmonic_signal_smooth2, n_test_nodes);
  print_file(f2_out_data_sm4, result_1d_harmonic_signal_smooth4, n_test_nodes);
  print_file(f2_out_data_sm8, result_1d_harmonic_signal_smooth8, n_test_nodes);

  }

  fclose(f2_in_data);

  fclose(f2_out_data);
  fclose(f2_out_data_sm2);
  fclose(f2_out_data_sm4);
  fclose(f2_out_data_sm8);

#endif



   //testing approximation_2d

#if APPROXIMATION_2D_TEST

   printf("test approximation 2D \n");

   //define variables

   // шаг исходной равномерной сетки
   t_vector_2d_f h_2d; 
   h_2d._1 = H_1;
   h_2d._2 = H_2;

   // размеры исходного прямоугольника
   t_vector_2d_f l_2d;
   l_2d._1 = L_1;
   l_2d._2 = L_2;

   //кол-во входных узлов
//   int n_nodes_1 = l_2d._1/h_2d._1;
//   int n_nodes_2 = l_2d._2/h_2d._2;

   // шаг сетки для тестирования
   t_vector_2d_f h_test_2d;
   h_test_2d._1 = 0.001;
   h_test_2d._2 = 0.001;

   // размеры выходных данных
//   int n_test_nodes_1 = l_2d._1/h_test_2d._1 - n_nodes_1;
//   int n_test_nodes_2 = l_2d._2/h_test_2d._2 - n_nodes_2;

   //константа
//   float coeff1_2d = 10.;

   //коэф-ты полинома
   t_vector_2d_f coeff2_2d[2];
   coeff2_2d[1]._1 = -10.;
   coeff2_2d[1]._2 = 20.;
   coeff2_2d[0]._1 = 3.;
   coeff2_2d[0]._2 = -7.;

   //параметры сглаживания
   t_vector_2d_int r1_2d;
   r1_2d._1 = 2;
   r1_2d._2 = 2;

   t_vector_2d_int r2_2d;
   r2_2d._1 = 4;
   r2_2d._2 = 4;

   t_vector_2d_int r3_2d;
   r3_2d._1 = 8;
   r3_2d._2 = 8;

   // имитация известных значений в узлах 

   printf("imit node values \n");

   // const_2d(&node_values_2d_const[0][0], n_nodes_1, n_nodes_2, h_2d, coeff1_2d);
   //line_func_2d(&node_values_2d_line[0][0], n_nodes_1, n_nodes_2, h_2d, coeff2_2d);
   polinomial_func_2_2d_v1(&node_values_2d_polinom[0][0], n_nodes_1, n_nodes_2, h_2d, coeff2_2d);
   //reverse_polinomial_func_2_2d_v1(&node_values_2d_reverse_polinom[0][0], n_nodes_1, n_nodes_2, h_2d, coeff2_2d);

   printf("write to file node values \n");

   FILE *f1_in_2d, *f2_in_2d, *f3_in_2d, *f4_in_2d;
   //f1_in_2d = fopen("node_values_const_2d.dat", "w");
   //f2_in_2d = fopen("node_values_line_2d.dat", "w");
   f3_in_2d = fopen("node_values_polinom_2d.dat", "w");
   //f4_in_2d = fopen("node_values_reverse_polinom_2d.dat", "w");

   //print_file(f1_in_2d, &node_values_2d_const[0][0],           n_nodes_1*n_nodes_2);
   //print_file(f2_in_2d, &node_values_2d_line[0][0],            n_nodes_1*n_nodes_2);
   print_file(f3_in_2d, &node_values_2d_polinom[0][0],         n_nodes_1*n_nodes_2);
   //print_file(f4_in_2d, &node_values_2d_reverse_polinom[0][0], n_nodes_1*n_nodes_2);

   //  fclose(f1_in_2d);
   //fclose(f2_in_2d);
   fclose(f3_in_2d);
   //fclose(f4_in_2d);

   //имитация эталонной функции для сравнения. 
   //Размерность совпадает с размерностью результата

   printf("imit etalon func \n");

//   const_2d                       (&etalon_2d_const[0][0],           n_test_nodes_1, n_test_nodes_2, h_test_2d, coeff1_2d);
//   line_func_2d                   (&etalon_2d_line[0][0],            n_test_nodes_1, n_test_nodes_2, h_test_2d, coeff2_2d);
   polinomial_func_2_2d_v1        (&etalon_2d_polinom[0][0],         n_test_nodes_1, n_test_nodes_2, h_test_2d, coeff2_2d);
//   reverse_polinomial_func_2_2d_v1(&etalon_2d_reverse_polinom[0][0], n_test_nodes_1, n_test_nodes_2, h_test_2d, coeff2_2d);

   printf("write to file etalon func\n");

   FILE *f1_et_2d, *f2_et_2d, *f3_et_2d, *f4_et_2d;

   // if ((f1_et_2d = fopen("etalon_const_2d.dat", "w"))          == NULL) printf("open file error\n");
   //if ((f2_et_2d = fopen("etalon_line_2d.dat", "w"))           == NULL) printf("open file error\n");
   if ((f3_et_2d = fopen("etalon_polinom_2d.dat", "w"))        == NULL) printf("open file error\n");
   //if ((f4_et_2d = fopen("etalon_reverse_polinom_2d.dat", "w"))== NULL) printf("open file error\n");

   //print_file(f1_et_2d, &etalon_2d_const[0][0],           n_test_nodes_1*n_test_nodes_2);
   //print_file(f2_et_2d, &etalon_2d_line[0][0],            n_test_nodes_1*n_test_nodes_2);
   print_file(f3_et_2d, &etalon_2d_polinom[0][0],         n_test_nodes_1*n_test_nodes_2);
   //print_file(f4_et_2d, &etalon_2d_reverse_polinom[0][0], n_test_nodes_1*n_test_nodes_2);

   //fclose(f1_et_2d);
   //fclose(f2_et_2d);
   fclose(f3_et_2d);
   //fclose(f4_et_2d);


   //тестирование 2d аппроксимации без сглаживания и со сглаживанием

   int k1, k2;

   t_vector_2d_f x_test_2d;

   printf("test approximation 2d \n");

   for(k1 = 0; k1 < n_test_nodes_1; k1++)
     for(k2 = 0; k2 < n_test_nodes_2; k2++){

       x_test_2d._1 = k1*h_test_2d._1;
       x_test_2d._2 = k2*h_test_2d._2;

	   // result1_2d[k1][k2]         = approximation_2d         (x_test_2d, &node_values_2d_const[0][0], l_2d, h_2d);
       //result1_2d_smooth1[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_const[0][0], l_2d, h_2d, r1_2d); 
       //result1_2d_smooth2[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_const[0][0], l_2d, h_2d, r2_2d);
       //result1_2d_smooth3[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_const[0][0], l_2d, h_2d, r3_2d);  

       //result2_2d[k1][k2]         = approximation_2d         (x_test_2d, &node_values_2d_line[0][0], l_2d, h_2d);
       //result2_2d_smooth1[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_line[0][0], l_2d, h_2d, r1_2d); 
       //result2_2d_smooth2[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_line[0][0], l_2d, h_2d, r2_2d);
       //result2_2d_smooth3[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_line[0][0], l_2d, h_2d, r3_2d);  

       result3_2d[k1][k2]         = approximation_2d         (x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d);
       result3_2d_smooth1[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d, r1_2d); 
       result3_2d_smooth2[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d, r2_2d);
       result3_2d_smooth3[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d, r3_2d);  

	   diff1_2d[k1][k2] = fabs(etalon_2d_polinom[k1][k2] - result3_2d[k1][k2]);
	   diff2_2d[k1][k2] = fabs(etalon_2d_polinom[k1][k2] - result3_2d_smooth1[k1][k2]);
	   diff3_2d[k1][k2] = fabs(etalon_2d_polinom[k1][k2] - result3_2d_smooth2[k1][k2]);
	   diff4_2d[k1][k2] = fabs(etalon_2d_polinom[k1][k2] - result3_2d_smooth3[k1][k2]);

	   // result4_2d[k1][k2]         = approximation_2d         (x_test_2d, &node_values_2d_reverse_polinom[0][0], l_2d, h_2d);
       //result4_2d_smooth1[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_reverse_polinom[0][0], l_2d, h_2d, r1_2d); 
       //result4_2d_smooth2[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_reverse_polinom[0][0], l_2d, h_2d, r2_2d);
       //result4_2d_smooth3[k1][k2] = approximation_2d_S_smooth(x_test_2d, &node_values_2d_reverse_polinom[0][0], l_2d, h_2d, r3_2d);  
     }

   printf("write to file results of 2D approximation \n");

   FILE *f_res1_2d, *f_res1_2d_smooth1, *f_res1_2d_smooth2, *f_res1_2d_smooth3;
   FILE *f_res2_2d, *f_res2_2d_smooth1, *f_res2_2d_smooth2, *f_res2_2d_smooth3;
   FILE *f_res3_2d, *f_res3_2d_smooth1, *f_res3_2d_smooth2, *f_res3_2d_smooth3;
   FILE *f_res4_2d, *f_res4_2d_smooth1, *f_res4_2d_smooth2, *f_res4_2d_smooth3;

   FILE *f_diff1_2d, *f_diff2_2d, *f_diff3_2d, *f_diff4_2d;

//   f_res1_2d = fopen("result_2d_const.dat", "w");
//   f_res1_2d_smooth1 = fopen("result_2d_const_smooth1.dat", "w");
//   f_res1_2d_smooth2 = fopen("result_2d_const_smooth2.dat", "w");
//   f_res1_2d_smooth3 = fopen("result_2d_const_smooth3.dat", "w");

//   f_res2_2d = fopen("result_2d_line.dat", "w");
//   f_res2_2d_smooth1 = fopen("result_2d_line_smooth1.dat", "w");
//   f_res2_2d_smooth2 = fopen("result_2d_line_smooth2.dat", "w");
//   f_res2_2d_smooth3 = fopen("result_2d_line_smooth3.dat", "w");

   f_res3_2d = fopen("result_2d_polinom.dat", "w");
   f_res3_2d_smooth1 = fopen("result_2d_polinom_smooth1.dat", "w");
   f_res3_2d_smooth2 = fopen("result_2d_polinom_smooth2.dat", "w");
   f_res3_2d_smooth3 = fopen("result_2d_polinom_smooth3.dat", "w");

   f_diff1_2d = fopen("diff1_2d.dat", "w");
   f_diff2_2d = fopen("diff2_2d.dat", "w");
   f_diff3_2d = fopen("diff3_2d.dat", "w");
   f_diff4_2d = fopen("diff4_2d.dat", "w");


//   f_res4_2d = fopen("result_2d_reverse_polinom.dat", "w");
//   f_res4_2d_smooth1 = fopen("result_2d_reverse_polinom_smooth1.dat", "w");
//   f_res4_2d_smooth2 = fopen("result_2d_reverse_polinom_smooth2.dat", "w");
//   f_res4_2d_smooth3 = fopen("result_2d_reverse_polinom_smooth3.dat", "w");
   
//   print_file(f_res1_2d,         &result1_2d[0][0],         n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res1_2d_smooth1, &result1_2d_smooth1[0][0], n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res1_2d_smooth2, &result1_2d_smooth2[0][0], n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res1_2d_smooth3, &result1_2d_smooth3[0][0], n_test_nodes_1*n_test_nodes_2);

//   print_file(f_res2_2d,         &result2_2d[0][0],         n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res2_2d_smooth1, &result2_2d_smooth1[0][0], n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res2_2d_smooth2, &result2_2d_smooth2[0][0], n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res2_2d_smooth3, &result2_2d_smooth3[0][0], n_test_nodes_1*n_test_nodes_2);

   print_file(f_res3_2d,         &result3_2d[0][0],         n_test_nodes_1*n_test_nodes_2);
   print_file(f_res3_2d_smooth1, &result3_2d_smooth1[0][0], n_test_nodes_1*n_test_nodes_2);
   print_file(f_res3_2d_smooth2, &result3_2d_smooth2[0][0], n_test_nodes_1*n_test_nodes_2);
   print_file(f_res3_2d_smooth3, &result3_2d_smooth3[0][0], n_test_nodes_1*n_test_nodes_2);

   print_file(f_diff1_2d, &diff1_2d[0][0], n_test_nodes_1*n_test_nodes_2);
   print_file(f_diff2_2d, &diff2_2d[0][0], n_test_nodes_1*n_test_nodes_2);
   print_file(f_diff3_2d, &diff3_2d[0][0], n_test_nodes_1*n_test_nodes_2);
   print_file(f_diff4_2d, &diff4_2d[0][0], n_test_nodes_1*n_test_nodes_2);

//   print_file(f_res4_2d,         &result4_2d[0][0],         n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res4_2d_smooth1, &result4_2d_smooth1[0][0], n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res4_2d_smooth2, &result4_2d_smooth2[0][0], n_test_nodes_1*n_test_nodes_2);
//   print_file(f_res4_2d_smooth3, &result4_2d_smooth3[0][0], n_test_nodes_1*n_test_nodes_2);

//   fclose(f_res1_2d);
//   fclose(f_res1_2d_smooth1);
//   fclose(f_res1_2d_smooth2);
//   fclose(f_res1_2d_smooth3);
//   fclose(f_res2_2d);
//   fclose(f_res2_2d_smooth1);
//   fclose(f_res2_2d_smooth2);
//   fclose(f_res2_2d_smooth3);
   fclose(f_res3_2d);
   fclose(f_res3_2d_smooth1);
   fclose(f_res3_2d_smooth2);
   fclose(f_res3_2d_smooth3);

   fclose(f_diff1_2d);
   fclose(f_diff2_2d);
   fclose(f_diff3_2d);
   fclose(f_diff4_2d);

//   fclose(f_res4_2d);
//   fclose(f_res4_2d_smooth1);
//   fclose(f_res4_2d_smooth2);
//   fclose(f_res4_2d_smooth3);

//

   x_test_2d._1 = 0.2;
   x_test_2d._2 = 0.2;

float x_etalon = coeff2_2d[1]._1*x_test_2d._1*x_test_2d._1 + coeff2_2d[1]._2*x_test_2d._1*x_test_2d._2 + coeff2_2d[0]._1 + coeff2_2d[0]._2;

printf("etalion value:\t%f\n", x_etalon);

float test1_2d = approximation_2d         (x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d);
float test2_2d = approximation_2d_S_smooth(x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d, r1_2d); 
float test3_2d = approximation_2d_S_smooth(x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d, r2_2d);
float test4_2d = approximation_2d_S_smooth(x_test_2d, &node_values_2d_polinom[0][0], l_2d, h_2d, r3_2d);  

printf("2d approximation:\t%f\n", test1_2d);
printf("2d approximation smooth 2:\t%f\n", test2_2d);
printf("2d approximation smooth 4:\t%f\n", test3_2d);
printf("2d approximation smooth 8:\t%f\n", test4_2d);

#endif

  return 0;
}
