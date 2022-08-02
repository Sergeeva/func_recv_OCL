#include <libutils/misc.h>
#include <libutils/timer.h>
#include <libutils/fast_random.h>

#include <libgpu/context.h>
#include <libgpu/shared_device_buffer.h>

// Этот файл будет сгенерирован автоматически в момент сборки - см. convertIntoHeader в CMakeLists.txt:22
//#include "cl/sum_cl.h"
#include "cl/func_recv_2D_prod_cl.h"
#include "cl/func_recv_2D_sum_cl.h"


#include <vector>


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
int K_x = 128;
int K_y = 128; // а если K_x!=K_y ??
int K = K_x*K_y;

//Количество узлов мелкой сетки в интервале крупной
//без учёта левой границы
//одинаковое по каждой координате
int N = 5; //==N_x ==N_y ??

//Количество узлов мелкой сетки
int M_x = (N+1)*K_x;
int M_y = (N+1)*K_y;
int M = M_x*M_y;

//Порядок ядра (чётное значение)
int r = 2; //4

// Размер !четверти! ядра
int Nr = (N+1)*r/4; 

// Массивы (матрицы) в DRAM
// (host memory)

std::vector<float> in_data(K); //NB: в общем случае значения мб не только float, но и вектора
std::vector<float> out_data(M);
std::vector<float> out_data_etal(M);

// значения ядра 
//храним четверть (правую нижнюю) с учётом симметричности по обоим координатам
std::vector<float> kern_vals(Nr*Nr); 

// значения произведений
// для каждого элемента крупной сетки хранится вычисляется матрица произведений
std::vector<float> prod_vals((Nr*Nr)*K); //(((N+1)*r/2)*((N+1)*r/2)*K); //результат работы kern_prod


// Имитируем входные данные :: крупную сетку 
// TODO::использовать imit_test_fun_2d
{

}

// Считаем эталон 
{

}

//
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


// 2D размер задачи для kern_sum
        unsigned int workGroupSize_2x = 32;
        unsigned int workGroupSize_2y = 32;
        unsigned int global_work_size_2x = (M_x + workGroupSize_2x - 1) / workGroupSize_2x * workGroupSize_2x;
        unsigned int global_work_size_2y = (M_y + workGroupSize_2y - 1) / workGroupSize_2y * workGroupSize_2y;
//

        //timer t;  
        //for (int i=0; i<benchmarkingIters; ++i){

        // Прогружаем данные из векторов as
        // DRAM --> VRAM
        // (есть нетипизированный метод write для которого количество измеряется в байтах,
        // и типизированный writeN, для которого количество измеряется в количестве float-элементов, т.к. gpu::gpu_mem_32f - это shared_device_buffer_typed<float>)
            in_data_gpu.writeN(in_data.data(), K);
            out_data_gpu.writeN(out_data.data(), M);
            kern_vals_gpu.writeN(kern_vals.data(), (Nr*Nr));
            prod_vals_gpu.writeN(prod_vals.data(), (Nr*Nr)*K);

            kern_prod.exec( gpu::WorkSize(workGroupSize_1, global_work_size_1),
                            in_data_gpu, K, 
                            kern_vals_gpu, Nr, 
                            prod_vals_gpu);
            std::cout<<"[kern_prod] exec done"<<std::endl;

            kern_sum.exec( gpu::WorkSize(   workGroupSize_2x, workGroupSize_2y, 
                                            global_work_size_2x, global_work_size_2y),
                            prod_vals_gpu, 
                            K_x, K_y, 
                            (N+1), r, 
                            out_data_gpu);
        
            std::cout<<"[kern_sum] exec done"<<std::endl;

            //res_gpu.readN(&sum, 1); //VRAM -> DRAM

            //EXPECT_THE_SAME(reference_sum, sum, "GPU <OCL kern> result should be consistent!");

            //std::cout<<"reference_sum = "<< reference_sum <<" :: sum = "<< sum <<std::endl;

            //t.nextLap();
        //}
        //std::cout << "GPU <OCL kern>: " << t.lapAvg() << "+-" << t.lapStd() << " s" << std::endl;
        //std::cout << "GPU <OCL kern>: " << (n/1000.0/1000.0) / t.lapAvg() << " millions/s" << std::endl;

    }
    return 0;
}
