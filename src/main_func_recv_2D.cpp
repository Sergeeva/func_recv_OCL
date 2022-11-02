#include <libutils/misc.h>
#include <libutils/timer.h>
#include <libutils/fast_random.h>

#include <libgpu/context.h>
#include <libgpu/shared_device_buffer.h>


#include <f_kernels_2d.h>
#include <def_types.h>
#include <imit_test_func_2d.h>

// Этот файл будет сгенерирован автоматически в момент сборки - см. convertIntoHeader в CMakeLists.txt:22
//#include "cl/sum_cl.h"
#include "cl/func_recv_2D_prod_cl.h"
#include "cl/func_recv_2D_sum_cl.h"


#include <vector>
#include <fstream>

//////////////////////////////////////////
//Последовательная Эмуляция выч. ядер
//////////////////////////////////////////

void prod_kern_emul(    const float* in_data, 
                        int K, // кол-во узлов крупной сетки
                        const float* kern_vals,
                        int Nr, 
                        float* prod_vals,
                        int globalID){

//std::cout<<"prod_kern_emul STARTED"<<std::endl;

    for (int n_y = 0; n_y<Nr; n_y++){
        for(int n_x = 0; n_x<Nr; n_x++){
            prod_vals[globalID*Nr*Nr + n_y*Nr + n_x]=
                kern_vals[n_y*Nr + n_x] * in_data[globalID];
        };
    };  

//std::cout<<"prod_kern_emul ENDED"<<std::endl;
}
        
void sum_kern_emul(const float* prod_vals, 
                    int K_x, int K_y, 
                    int N, int r, int Nr, 
                    float* result, 
                    int globalID_x, int globalID_y){


//std::cout<<"sum_kern_emul STARTED"<<std::endl;

    //индекс ближайшей точки крупной сетки слева-сверху
    //т.е. округление вниз!!
    int k_x = globalID_x / N;
    int k_y = globalID_y / N;

    //индекс в ячейке мелкой сетки 
    int m_x = globalID_x % N;
    int m_y = globalID_y % N;

    //индекс в мелкой сетке
    int x = globalID_x;
    int y = globalID_y;

    result[y*(N*K_x) + x] = 0.0;

//std::cout<<"sum_kern_emul [x dec]"<<std::endl;                             
    // "убывающая" половина по оси x    
    for (int i_x = 0; i_x < r/2; i_x++){
 //       std::cout<<"sum_kern_emul [y dec]"<<std::endl;  
        // "убывающая" половина по оси y
        for (int i_y = 0; i_y < r/2; i_y++){ 
            if((k_x - i_x)>=0 && ((k_y-i_y)>=0)){ //проверка выхода за границу сетки
                result[y*(N*K_x) + x] += prod_vals[((k_x-i_x)+(k_y-i_y)*K_x)*Nr*Nr + (i_y*N + m_y)*Nr + (i_x*N + m_x)];
            }
        }
//        std::cout<<"sum_kern_emul [y inc]"<<std::endl;  
        //"возрастающая" половина по оси y
        for (int i_y = 1; i_y <= r/2; i_y++){  
        if((k_x-i_x)>=0 && (k_y+i_y)<K_y) {//проверка выхода за границу сетки
            result[y*(N*K_x) + x] += prod_vals[((k_x-i_x)+(k_y+i_y)*K_x)*Nr*Nr + (i_y*N  - m_y)*Nr + (i_x*N + m_x)];
        }
        }
    }

 //   std::cout<<"sum_kern_emul [x inc]"<<std::endl;  
    // "возрастающая" половина по оси x 
    for (int i_x = 1; i_x <= r/2; i_x++){    
    // "убывающая" половина по оси y
 //       std::cout<<"sum_kern_emul [y dec]"<<std::endl;  
        for (int i_y = 0; i_y < r/2; i_y++){ 
        if((k_x + i_x)<K_x && (k_y-i_y)>=0){ //проверка выхода за границу сетки
            result[y*(N*K_x) + x] += prod_vals[((k_x+i_x)+(k_y-i_y)*K_x)*Nr*Nr + (i_y*N + m_y)*Nr + (i_x*N  - m_x)];
        }
        }

//        std::cout<<"sum_kern_emul [y inc]"<<std::endl;  
    //"возрастающая" половина по оси y
        for (int i_y = 1; i_y <= r/2; i_y++){  
        if((k_x + i_x)<K_x && (k_y+i_y)<K_y){ //проверка выхода за границу сетки
            result[y*(N*K_x) + x] += prod_vals[((k_x+i_x)+(k_y+i_y)*K_x)*Nr*Nr + (i_y*N  - m_y)*Nr + (i_x*N  - m_x)];
        }
        }

    }
//std::cout<<"sum_kern_emul ENDED"<<std::endl;
        
}

///////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
void raiseFail(const T &a, const T &b, std::string message, std::string filename, int line)
{
    if (a != b) {
        std::cerr << message << " But " << a << " != " << b << ", " << filename << ":" << line << std::endl;
        throw std::runtime_error(message);
    }
}

#define EXPECT_THE_SAME(a, b, message) raiseFail(a, b, message, __FILE__, __LINE__)


int main(int argc, char **argv)
{
    int benchmarkingIters = 10;

// Определяем размер задачи 

// Кол-во узлов крупной сетки
const int K_x = 16;
const int K_y = 16; // а если K_x!=K_y ??
const int K = K_x*K_y;

std::cout << "Coarse net size  " << K_x << "x" << K_y << std::endl;

//Количество узлов мелкой сетки в интервале крупной
//без учёта левой границы
//одинаковое по каждой координате
const int N = 3; //==N_x ==N_y ??

//Количество узлов мелкой сетки
const int M_x = (N+1)*K_x;
const int M_y = (N+1)*K_y;
const int M = M_x*M_y;

std::cout << "Fine net size  " << M_x << "x" << M_y << std::endl;

//Порядок ядра (чётное значение)
int r = 4; //4

// Массивы (матрицы) в DRAM
// (host memory)

std::vector<float> in_data(K); //NB: в общем случае значения мб не только float, но и вектора

//Инициализируем ТЕСТОВЫЕ 
//входные данные
// см. imit_test_func_2d.h
t_vector_2d_f h;
h._1 = 1.0f;
h._2 = 1.0f;

//t_vector_2d_f coeffs[2];
//coeffs[0]._1 = -2.0f;
//coeffs[0]._2 = 2.0f;

//coeffs[1]._1 = 5.0f;
//coeffs[1]._2 = -3.0f;


//const_2d(in_data.data(), K_y, K_x, h, 5.0);

//line_func_2d(in_data.data(), K_y, K_x, h, coeffs);


t_vector_2d_f coeffs3[3];
coeffs3[0]._1 = 10.0f;
coeffs3[0]._2 = 10.0f;

coeffs3[1]._1 = 5.0f;
coeffs3[1]._2 = -2.0f;

coeffs3[2]._1 = 5.0f;
coeffs3[2]._2 = -2.0f;

//reverse_polinomial_func_2_2d_v1(in_data.data(), K_y, K_x, h, coeffs3);
test_surface_2(in_data.data(), K_y, K_x, h, coeffs3);

//polinomial_func_2_2d_v1(in_data.data(), K_y, K_x, h, coeffs3);
//////////////////////////////////////////////////////////////////////

std::vector<float> out_data(M);
std::vector<float> out_data_etal(M);

//////////////////////////////////////////////////////////////////////

// Размер !четверти! ядра
int Nr = (N+1)*(r/2) + 1; 

// значения ядра 
//храним четверть (правую нижнюю) с учётом симметричности по обоим координатам
std::vector<float> kern_vals(Nr*Nr);
//Инициализируем значения ядер
std::cout<<"init kernels r="<<r<<std::endl;
const int Nr1 = (N+1)*(r/2);
if(r == 2){
    for (int y = 0, i = 0; y >= -Nr1; y--, i++){
        for(int x = 0, j = 0; x <= Nr1; x++, j++){
            t_vector_2d_f t;
            t._1 = x/(float)Nr1;
            t._2 = y/(float)Nr1;
            kern_vals[i*Nr+j] = psi_2_2d(t);
//            std::cout<<kern_vals[i*Nr+j]<<"  ";
        }
//        std::cout << std::endl;
    }
}else if (r == 4){
    for (int y = 0, i = 0; y >= -Nr1; y--, i++){
        for(int x = 0, j = 0; x <= Nr1; x++, j++){
            t_vector_2d_f t;
            t._1 = x/(float)(N+1);
            t._2 = y/(float)(N+1);
            kern_vals[i*Nr+j] = psi_4_2d(t);
//            std::cout<<kern_vals[i*Nr+j]<<"  ";
        }
//       std::cout << std::endl;
    }
}

/*
//for print and draw
//r = 4
std::ofstream kern4_file;
kern4_file.open("../data_save/kern4.dat");
 for (int y = Nr, i = 0; y >-Nr; y--, i++){
        for(int x = -Nr, j = 0; x < Nr; x++, j++){
            t_vector_2d_f t;
            t._1 = x/(float)(N+1);
            t._2 = y/(float)(N+1);
            kern4_file<<psi_4_2d(t)<<"  ";
        }
    kern4_file<<std::endl;
 }
*/


// значения произведений
// для каждого элемента крупной сетки хранится вычисляется матрица произведений
std::vector<float> prod_vals((Nr*Nr)*K); //(((N+1)*r/2)*((N+1)*r/2)*K); //результат работы kern_prod


{
// Считаем эталон 
//замеряем производительность

timer t;
for (int iter = 0; iter < benchmarkingIters; ++iter) {

for(int k = 0; k<K; k++){
    prod_kern_emul( in_data.data(), K, // кол-во узлов крупной сетки
                    kern_vals.data(), Nr, 
                    prod_vals.data(),
                    k); //итерация по globalID
}

for (int m_y = 0; m_y < M_y; m_y++) {
    for (int m_x = 0; m_x < M_x; m_x++) {

        sum_kern_emul(prod_vals.data(),
            K_x, K_y,
            (N + 1), r, Nr,
            out_data_etal.data(),
            m_x, m_y);
    }
}
///
t.nextLap();
}

std::cout << "CPU:     " << t.lapAvg() << "+-" << t.lapStd() << " s" << std::endl;
//std::cout << "CPU:     " << (n / 1000.0 / 1000.0) / t.lapAvg() << " millions/s" << std::endl;
}



//
#define OCL 0
#if OCL
{
    // chooseGPUDevice:
    // - Если не доступо ни одного устройства - кинет ошибку
    // - Если доступно ровно одно устройство - вернет это устройство
    // - Если доступно N>1 устройства:
    //   - Если аргументов запуска нет или переданное число не находится в диапазоне от 0 до N-1 - кинет ошибку
    //   - Если аргумент запуска есть и он от 0 до N-1 - вернет устройство под указанным номером
        gpu::Device device = gpu::chooseGPUDevice(argc, argv);

    // Этот контекст после активации будет прозрачно использоваться при всех вызовах в libgpu библиотеке
    // это достигается использованием thread-local переменных, т.е. на самом деле контекст будет активирован для текущего потока исполнения
        gpu::Context context;
        context.init(device.device_id_opencl);
        context.activate();

    // Создаем буфер в видеопамяти
        gpu::gpu_mem_32f in_data_gpu;
        in_data_gpu.resizeN(K);

        gpu::gpu_mem_32f out_data_gpu;
        out_data_gpu.resizeN(M);

        gpu::gpu_mem_32f kern_vals_gpu;
        kern_vals_gpu.resizeN(Nr*Nr);

        gpu::gpu_mem_32f prod_vals_gpu;
        prod_vals_gpu.resizeN((Nr*Nr)*K);
    
    // Исходники кернела написаны в src/cl/aplusb.cl
    // Но благодаря convertIntoHeader(src/cl/aplusb.cl src/cl/aplusb_cl.h aplusb_kernel) (см. CMakeLists.txt:18)
    // при компиляции автоматически появится файл src/cl/aplusb_cl.h с массивом aplusb_kernel состоящим из байт исходника
    // т.о. программе не будет нужно в runtime читать файл с диска, т.к. исходник кернелов теперь хранится в массиве данных основной программы
        //ocl::Kernel kern_sum(sum_kernel, sum_kernel_length, "sum");
        //kern_sum.compile(false);

    ocl::Kernel kern_prod(prod_kernel, prod_kernel_length, "prod");
    kern_prod.compile(true);

    ocl::Kernel kern_sum(sum_kernel, sum_kernel_length, "sum");
    kern_sum.compile(true);
    

// 2D размер задачи для kern_prod
        unsigned int workGroupSize_1 = 128;
        unsigned int global_work_size_1 = (K + workGroupSize_1 - 1) / workGroupSize_1 * workGroupSize_1; //for kern_prod
        
        std::cout<<"Global Work Size1  " << global_work_size_1 << std::endl;


// 2D размер задачи для kern_sum
        unsigned int workGroupSize_2x = 128;
        unsigned int workGroupSize_2y = 1;
        unsigned int global_work_size_2x = M_x; //(M_x + workGroupSize_2x - 1) / workGroupSize_2x * workGroupSize_2x;
        unsigned int global_work_size_2y = M_y;//(M_y + workGroupSize_2y - 1) / workGroupSize_2y * workGroupSize_2y;

         std::cout<<"Global Work Size2  " << global_work_size_2x << " x " << global_work_size_2y << std::endl;

        timer t;  
        for (int i = 0; i < benchmarkingIters; ++i) {

            // Прогружаем данные из векторов as
            // DRAM --> VRAM
            // (есть нетипизированный метод write для которого количество измеряется в байтах,
            // и типизированный writeN, для которого количество измеряется в количестве float-элементов, т.к. gpu::gpu_mem_32f - это shared_device_buffer_typed<float>)
            in_data_gpu.writeN(in_data.data(), K);
            out_data_gpu.writeN(out_data.data(), M);
            kern_vals_gpu.writeN(kern_vals.data(), (Nr * Nr));
            prod_vals_gpu.writeN(prod_vals.data(), (Nr * Nr) * K);

            kern_prod.exec(gpu::WorkSize(workGroupSize_1, global_work_size_1),
                in_data_gpu, K,
                kern_vals_gpu, Nr,
                prod_vals_gpu);
            //std::cout<<"[kern_prod] exec done"<<std::endl;

            kern_sum.exec(gpu::WorkSize(workGroupSize_2x, workGroupSize_2y,
                global_work_size_2x, global_work_size_2y),
                prod_vals_gpu,
                K_x, K_y,
                (N + 1), r, Nr,
                out_data_gpu);

            //std::cout<<"[kern_sum] exec done"<<std::endl;

            out_data_gpu.readN(out_data.data(), M); //VRAM -> DRAM

            t.nextLap();
        }

        std::cout << "GPU <OCL kern>: " << t.lapAvg() << "+-" << t.lapStd() << " s" << std::endl;
        //std::cout << "GPU <OCL kern>: " << (n / 1000.0 / 1000.0) / t.lapAvg() << " millions/s" << std::endl;
    }
#endif

// EXPECT_THE_SAME(reference_sum, sum, "CPU result should be consistent!");


/////////////////////////////////////////////////
//Запись в файл результата из RAM
/////////////////////////////////////////////////
std::ofstream out_file;
out_file.open("../data_save/func_recv_2d__result.dat");

if (out_file.is_open()){
    std::cout<<"WRITING file with result"<<std::endl;
        for(int i = 0; i < M_y; i++){
            for(int j = 0; j < M_x; j++){
                out_file<<out_data_etal[i*M_x + j]<<"  ";
            }
            out_file<<std::endl;
        }            
        out_file.close();
    }


/////////////////////////////////////////////////
//Запись в файл входа для проверки
/////////////////////////////////////////////////
std::ofstream in_file;
in_file.open("../data_save/func_recv_2d__input.dat");

if (in_file.is_open()){
    std::cout<<"WRITING file with input"<<std::endl;
        for(int i = 0; i < K_y; i++){
            for(int j = 0; j < K_x; j++){
                in_file<<in_data[i*K_x + j]<<"  ";
            }
            in_file<<std::endl;
        }            
        in_file.close();
    }


    return 0;
}
