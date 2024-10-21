cmake_minimum_required(VERSION 3.10)

project(mksrfdata LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# 定义选项
option(LULC_IGBP_PFT "Enable LULC_IGBP_PFT" TRUE)
option(LULC_IGBP_PC "Enable LULC_IGBP_PC" OFF)

# 搜索当前目录下所有的Fortran源文件
file(GLOB SRC_FILES
    # ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
    ../main/MOD_Const_LC.F90
    MOD_AggregationRequestData.F90
    # MOD_ElmVector.F90
    # MOD_HRUVector.F90
    # MOD_LandPFT.F90
    MOD_LandCrop.F90
    MOD_LandElm.F90
    MOD_LandHRU.F90
    MOD_LandPatch.F90
    MOD_LandUrban.F90
    MOD_MeshFilter.F90
    MOD_PixelsetShared.F90
    MOD_RegionClip.F90
    MOD_SingleSrfdata.F90
    MOD_SrfdataDiag.F90
    MOD_SrfdataRestart.F90
)

if(LULC_IGBP_PFT OR LULC_IGBP_PC)
  list(APPEND SRC_FILES
    MOD_LandPFT.F90
  )
endif()

set(exe_FILES 
  MKSRFDATA.F90
  Aggregation_DBedrock.F90
  Aggregation_ForestHeight.F90
  Aggregation_LAI.F90
  Aggregation_LakeDepth.F90
  Aggregation_PercentagesPFT.F90
  Aggregation_SoilBrightness.F90
  Aggregation_SoilParameters.F90
  Aggregation_Topography.F90
  Aggregation_TopographyFactors.F90
  Aggregation_Urban.F90
)

set(exe mksrfdata.exe)

find_package(NetCDF REQUIRED)
include_directories("../include" ${NetCDF_F90_INCLUDEDIR} ${CMAKE_BINARY_DIR}/share)

add_subdirectory(../share ${CMAKE_BINARY_DIR}/share)

# 编译为静态链接库
add_library(${PROJECT_NAME} STATIC ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE colm_share ${NetCDF_F90_LIBRARIES})

add_executable(${exe} ${exe_FILES})
target_link_libraries(${exe} ${PROJECT_NAME} ${NetCDF_F90_LIBRARIES} z)
