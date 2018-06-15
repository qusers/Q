#!/bin/sh

#
# Build GCC with support for offloading to NVIDIA GPUs.
# After installation add $install_dir/lib64 to path:
# export LD_LIBRARY_PATH=/usr/local/gcc-7.3/lib64/
# export LIBRARY_PATH=/usr/local/gcc-7.3/lib64/
#
# For the latest version of gcc, run
# svn ls svn://gcc.gnu.org/svn/gcc/tags/ | grep '^gcc[_-][0-9]'\
#     | sort | tail -n1
#

work_dir=$HOME/gcc-7_offload/wrk
install_dir=/usr/local/gcc-7.3

# Location of the installed CUDA toolkit
cuda=/usr/local/cuda

# Number of cores to use for building
cores=8


# Build assembler and linking tools
mkdir -p $work_dir
cd $work_dir
git clone https://github.com/MentorEmbedded/nvptx-tools
cd nvptx-tools
./configure \
    --with-cuda-driver-include=$cuda/include \
    --with-cuda-driver-lib=$cuda/lib64 \
    --prefix=$install_dir
make
make install
cd ..

# Set up the GCC source tree
git clone https://github.com/MentorEmbedded/nvptx-newlib
svn co svn://gcc.gnu.org/svn/gcc/tags/gcc_7_3_0_release gcc
cd gcc
contrib/download_prerequisites
ln -s ../nvptx-newlib/newlib newlib
cd ..
target=$(gcc/config.guess)

# Build nvptx GCC
mkdir build-nvptx-gcc
cd build-nvptx-gcc
../gcc/configure \
    --target=nvptx-none --with-build-time-tools=$install_dir/nvptx-none/bin \
    --enable-as-accelerator-for=$target \
    --disable-sjlj-exceptions \
    --enable-newlib-io-long-long \
    --enable-languages="c,c++,fortran,lto" \
    --prefix=$install_dir
make -j$cores
make install
cd ..

# Build host GCC
mkdir build-host-gcc
cd  build-host-gcc
../gcc/configure \
    --enable-offload-targets=nvptx-none \
    --with-cuda-driver-include=$cuda/include \
    --with-cuda-driver-lib=$cuda/lib64 \
    --disable-bootstrap \
    --disable-multilib \
    --enable-languages="c,c++,fortran,lto" \
    --prefix=$install_dir
make -j$cores
make install
cd ..
