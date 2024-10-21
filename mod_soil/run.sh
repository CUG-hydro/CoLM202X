export NETCDF_DIR="/opt/netcdf_v4.9.2_openmpi"
export CMAKE_MODULE_PATH="/home/kong/github/CUG-hydro/CoLM202X/cmake"

rm bld -rf
mkdir bld
cd bld
cmake -DCMAKE_MODULE_PATH="$CMAKE_MODULE_PATH" ..
# cmake ..
