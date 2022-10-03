#define WORK_GROUP_SIZE 128

//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable


__kernel void sum(__global const float* prod_vals, 
                  int K_x, int K_y, 
                  int N, int r, int Nr, 
                  volatile __global float* result){

    int localID_x = get_local_id(0);
    int globalID_x = get_global_id(0);

    int localID_y = get_local_id(1);
    int globalID_y = get_global_id(1);

    //индекс ближайшей точки крупной сетки слева-сверху
    //т.е. округление вниз!!
    int k_x = globalID_x / N;
    int k_y = globalID_y / N;

    //индекс в ячейке мелкой сетки 
    int m_x = globalID_x % N;
    int m_y = globalID_y % N;

    //индекс в четвертинке ядра
    //int mm_x = globalID_x % Nr;
    //int mm_y = globalID_y % Nr;

    //индекс в мелкой сетке
    int x = globalID_x;
    int y = globalID_y;


    
    //__local int input_local[WORK_GROUP_SIZE];
    //input_local[localID] = input[globalID];

    //barrier(CLK_LOCAL_MEM_FENCE); //барьер по чтению из VRAM

    // Суммировать ВСЕ точки вокруг!
    
    
    result[y*(N*K_x) + x] = 0.0;
                            

            
    // "убывающая" половина по оси x    
    for (int i_x = 0; i_x < r/2; i_x++){
        // "убывающая" половина по оси y
        for (int i_y = 0; i_y < r/2; i_y++){ 
        if((k_x - i_x)>=0 && ((k_y-i_y)>=0)){ //проверка выхода за границу сетки
            result[y*(N*K_x) + x] += prod_vals[((k_x-i_x)+(k_y-i_y)*K_x)*Nr*Nr + (i_y*N + m_y)*Nr + (i_x*N + m_x)];
        }
        }
        //"возрастающая" половина по оси y
        for (int i_y = 1; i_y <= r/2; i_y++){  
        if((k_x-i_x)>=0 && (k_y+i_y)<K_y) {//проверка выхода за границу сетки
            result[y*(N*K_x) + x] += prod_vals[((k_x-i_x)+(k_y+i_y)*K_x)*Nr*Nr + (i_y*N  - m_y)*Nr + (i_x*N + m_x)];
        }
        }
    }

    // "возрастающая" половина по оси x 
    for (int i_x = 1; i_x <= r/2; i_x++){    
    // "убывающая" половина по оси y
        for (int i_y = 0; i_y < r/2; i_y++){ 
        if((k_x + i_x)<K_x && (k_y-i_y)>=0){ //проверка выхода за границу сетки
            result[y*(N*K_x) + x] += prod_vals[((k_x+i_x)+(k_y-i_y)*K_x)*Nr*Nr + (i_y*N + m_y)*Nr + (i_x*N  - m_x)];
        }
        }
    //"возрастающая" половина по оси y
        for (int i_y = 1; i_y <= r/2; i_y++){  
        if((k_x + i_x)<K_x && (k_y+i_y)<K_y){ //проверка выхода за границу сетки
            result[y*(N*K_x) + x] += prod_vals[((k_x+i_x)+(k_y+i_y)*K_x)*Nr*Nr + (i_y*N  - m_y)*Nr + (i_x*N  - m_x)];
        }
        }

    }
        
}