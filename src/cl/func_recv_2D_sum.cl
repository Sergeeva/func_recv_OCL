#define WORK_GROUP_SIZE 128

//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable


__kernel void sum(__global const float* prod_vals, 
                  int K_x, int K_y, 
                  int N, int r,
                  volatile __global float* result){

    int localID_x = get_local_id(0);
    int globalID_x = get_global_id(0);

    int localID_y = get_local_id(1);
    int globalID_y = get_global_id(1);

    int Nr = N*r/4;

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

    //__local int input_local[WORK_GROUP_SIZE];
    //input_local[localID] = input[globalID];

    //barrier(CLK_LOCAL_MEM_FENCE); //барьер по чтению из VRAM

    //инициализация суммы в точке мелной сетки
    result[y*(N*K_x) + x] = prod_vals[(k_x+k_y)*Nr*Nr + m_y*Nr + m_x];

    //суммирование по оси х
    // "убывающая" половина, исключая начальную точку
    for(int i = 1; i < r/2; i++){
        if((k_x - i)>=0)
            result[y*(N*K_x) + x] += prod_vals[((k_x-i)+k_y*K_x)*Nr*Nr + m_y*Nr + (i*N + m_x)];
    }
    // "возрастающая" половина
    for(int i = 0; i < r/2; i++){
        if((k_x + i)<=K_x)
            result[y*(N*K_x) + x] += prod_vals[((k_x+i)+k_y*K_x)*Nr*Nr + m_y*Nr + (i*N - m_x)];
    }
    

    //суммирование по оси y  --- аналогично
    // "убывающая" половина, исключая начальную точку
    for(int i = 1; i < r/2; i++){
        if((k_y - i)>=0)
            result[y*(N*K_x) + x] += prod_vals[(k_x+(k_y-i)*K_x)*Nr*Nr + (i*N + m_y)*Nr + m_x];
    }
    // "возрастающая" половина
    for(int i = 0; i < r/2; i++){
        if((k_y + i)<=K_y)
            result[y*(N*K_x) + x] += prod_vals[(k_x+(k_y+i)*K_x)*Nr*Nr + (i*N - m_y)*Nr + m_x];
    }

    /*
    //atomics example
    if(localID==0){
        atom_add(res, input_local[localID]);
    }
    */
    
}