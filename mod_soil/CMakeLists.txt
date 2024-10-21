cmake_minimum_required(VERSION 3.10)

project(soil LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

add_compile_options(-fPIC)

set(SRC_FILES 
  ../main/MOD_Runoff.F90
  ../main/HYDRO/MOD_Hydro_VIC.F90
  ../main/HYDRO/MOD_Hydro_VIC_Variables.F90
  ../main/HYDRO/MOD_Hydro_SoilWater.F90
  ../main/HYDRO/MOD_Hydro_SoilFunction.F90
  ../main/MOD_Const_Physical.F90
  ../main/MOD_Vars_TimeInvariants.F90
  ../main/MOD_Vars_1DFluxes.F90
  ../main/MOD_Vars_1DPFTFluxes.F90
  ../mksrfdata/MOD_LandPFT.F90
  MOD_SoilSnowHydrology.F90
)
# message(STATUS "SRC_FILES: ${SRC_FILES}")

find_package(NetCDF REQUIRED)
include_directories("../include" ${NetCDF_F90_INCLUDEDIR} 
  ${CMAKE_BINARY_DIR}/share
  ${CMAKE_BINARY_DIR}/mksrfdata
)

# add_subdirectory(../share ${CMAKE_BINARY_DIR}/share)
add_subdirectory(../mksrfdata ${CMAKE_BINARY_DIR}/mksrfdata)

# add_library(${PROJECT_NAME} STATIC ${SRC_FILES})
add_library(${PROJECT_NAME} SHARED ${SRC_FILES})

# add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE mksrfdata colm_share netcdff z)
