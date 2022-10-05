cmake_minimum_required(VERSION 3.1)

project(func_recv_model)

set(HEADERS
        def_types.h
        f_kernels_2d.h
        f_kernels_1d.h
	useful_tools.h
	imit_test_func_2d.h
	)

set(SOURCES
        f_kernels_2d.cpp
        f_kernels_1d.cpp
	useful_tools.cpp
	imit_test_func_2d.cpp
        )

set(CMAKE_CXX_STANDARD 11)
add_library(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries (${PROJECT_NAME} ${CMAKE_DL_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

