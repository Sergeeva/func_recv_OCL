#define WORK_GROUP_SIZE 128 //кол-во выч. потоков (=элементов) в WG

//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable


__kernel void prod(    __global const float* in_data, 
                        int K, // кол-во узлов крупной сетки
                        __global const float* kern_vals,
                        int Nr, //размер 
                        volatile __global float* prod_vals){

    int localID = get_local_id(0);
    int globalID = get_global_id(0);

    //__local int in_local[WORK_GROUP_SIZE]; //буфер в локальной памяти workgroup (WG)
    //in_local[localID] = in_data[globalID]; //каждый выч. элемент WG берёт свой узел крупной сетки

    //TODO:: в локальной памяти WG должны также оказаться инициализированные ячейки значений ядра
    //мб скопировать "ведущим" выч. элементом WG

    //TODO:: выделить буфер в локальной памяти WG для результата
    //ячеек с произведениями в одной точке
    //Выгружается в глобальную VRAM по globalID  каждым выч. элементом WG

   // barrier(CLK_LOCAL_MEM_FENCE); //синхр по чтению из VRAM

    for (int n_y = 0; n_y<Nr; n_y++){
        for(int n_x = 0; n_x<Nr; n_x++){
            prod_vals[globalID*Nr*Nr + n_y*Nr + n_x]=
                kern_vals[n_y*Nr + n_x] * in_data[globalID];//in_local[localID]
        };
    };  
    
    //TODO::выгрузить результат для этой точки в VRAM по globalID

}
        