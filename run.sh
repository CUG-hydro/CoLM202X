export NETCDF="/opt/netcdf_v4.9.2_openmpi"
export LD_LIBRARY_PATH="$NETCDF/lib:$NETCDF/hdf5/lib" # :$NETCDF/zlib/lib
export NETCDF_LIB=$NETCDF/lib
export NETCDF_INC=$NETCDF/include

# export LD_LIBRARY_PATH="$NETCDF/lib:$NETCDF/hdf5/lib:$NETCDF/zlib/lib"
# export NETCDF="/home/kong/github/CUG-hydro/CoLM202X/netcdf"
echo $LD_LIBRARY_PATH

export FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
           -ffpe-trap=invalid,zero,overflow -fbounds-check \
           -mcmodel=medium -fbacktrace -fdump-core -cpp \
           -ffree-line-length-0 -fallow-argument-mismatch 


# sudo ln -s /opt/netcdf_v4.9.2_openmpi/zlib/lib/libcurl.so.4.8.0 /usr/lib/x86_64-linux-gnu/libcurl.so.4 -f
# ln -s -f /opt/netcdf_v4.9.2_openmpi/zlib/lib/libcurl.so.4.8.0 /opt/netcdf_v4.9.2_openmpi/zlib/lib/libcurl.so.4


export NETCDF_DIR="/opt/netcdf_v4.9.2_openmpi"
export CMAKE_MODULE_PATH="/home/kong/github/CUG-hydro/CoLM202X/cmake"

rm bld -rf
mkdir bld
cd bld
cmake -DCMAKE_MODULE_PATH="$CMAKE_MODULE_PATH" ..
# cmake ..
