cmake_minimum_required(VERSION 3.10)

project(colm_share LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/home/kong/github/CUG-hydro/CoLM202X/cmake")

# 添加 include 目录
add_compile_options(-fallow-argument-mismatch)

# 搜索当前目录下所有的Fortran源文件
file(GLOB SRC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
    "../main/MOD_Vars_Global.F90"
)

# 使用NetCDF的链接配置 (如果需要NetCDF库)
find_package(NetCDF REQUIRED)

include_directories("/home/kong/github/CUG-hydro/CoLM202X/include" ${NetCDF_F90_INCLUDEDIR})

# 编译为静态链接库
add_library(${PROJECT_NAME} STATIC ${SRC_FILES})

# if(NetCDF_FOUND)
    message(STATUS "NetCDF include directories: ${NetCDF_F90_INCLUDEDIR}")
    message(STATUS "NetCDF libraries: ${NetCDF_F90_LIBRARIES}")
    # target_include_directories(${PROJECT_NAME} PRIVATE ${NetCDF_F90_INCLUDEDIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${NetCDF_F90_LIBRARIES})
# endif()


# # 可以根据需要添加其他库的配置，比如 MPI, LAPACK 等
# find_package(MPI REQUIRED)
# if(MPI_FOUND)
#     target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran)
# endif()
