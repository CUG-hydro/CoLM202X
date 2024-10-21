cmake_minimum_required(VERSION 3.10)

project(mkinidata LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(SRC_FILES 
  CoLMINI.F90
  # ../main/MOD_Vars_TimeInvariants.F90
  # ../main/MOD_Vars_TimeVariables.F90
  ../main/MOD_Const_PFT.F90
  ../main/MOD_Const_Physical.F90
  ../main/MOD_LAIEmpirical.F90
  ../main/MOD_LAIReadin.F90
  ../main/MOD_OrbCoszen.F90
  ../main/HYDRO/MOD_Hydro_SoilFunction.F90
  ../main/MOD_Albedo.F90
  ../main/MOD_Aerosol.F90
  ../main/MOD_Vars_1DForcing.F90
  ../main/MOD_SnowSnicar.F90
  
  # MOD_DBedrockReadin.F90
  MOD_HtopReadin.F90
  MOD_Initialize.F90
  MOD_IniTimeVariable.F90
  MOD_LakeDepthReadin.F90
  MOD_PercentagesPFTReadin.F90
  MOD_SoilColorRefl.F90
  MOD_SoilParametersReadin.F90
  MOD_UrbanIniTimeVariable.F90
  MOD_UrbanReadin.F90
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

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE mksrfdata colm_share netcdff z)
