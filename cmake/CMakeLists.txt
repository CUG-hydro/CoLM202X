cmake_minimum_required(VERSION 3.10)

project(colm_share LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/home/kong/github/CUG-hydro/CoLM202X/hello")

# 添加 include 目录
add_compile_options(-fallow-argument-mismatch)

# 搜索当前目录下所有的Fortran源文件
# file(GLOB SRC_FILES
#     "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
#     "../main/MOD_Vars_Global.F90"
# )

# 使用NetCDF的链接配置 (如果需要NetCDF库)
find_package(NetCDF REQUIRED)

# include_directories("/home/kong/github/CUG-hydro/CoLM202X/include" ${NetCDF_INCLUDE_DIR})

# # 编译为静态链接库
# add_library(${PROJECT_NAME} STATIC ${SRC_FILES})

# # if(NetCDF_FOUND)
#     message(STATUS "NetCDF include directories: ${NetCDF_INCLUDE_DIRS}")
#     message(STATUS "NetCDF libraries: ${NetCDF_LIBRARIES}")
#     target_include_directories(${PROJECT_NAME} PRIVATE ${NetCDF_INCLUDE_DIRS})
#     target_link_libraries(${PROJECT_NAME} PRIVATE ${NetCDF_LIBRARIES})
# # endif()
