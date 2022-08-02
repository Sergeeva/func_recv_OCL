#define WORK_GROUP_SIZE 128

//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable


__kernel void prod(    __global const float* in_data, 
                        int K, // кол-во узлов крупной сетки
                        __global const float* kern_vals,
                        int Nr, //размер 
                        volatile __global float* prod_vals){

    int localID = get_local_id(0);
    int globalID = get_global_id(0);

    //__local int input_local[WORK_GROUP_SIZE];
    //input_local[localID] = input[globalID];

    //barrier(CLK_LOCAL_MEM_FENCE); //синхр по чтению из VRAM

    for (int n_y = 0; n_y<Nr; n_y++){
        for(int n_x = 0; n_x<Nr; n_x++){
            prod_vals[globalID*Nr*Nr + n_y*Nr + n_x]=
                kern_vals[n_y*Nr + n_x] * in_data[globalID];
        }
        
    }
    
    
}