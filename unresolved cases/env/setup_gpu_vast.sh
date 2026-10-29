#!/bin/bash
# UMID-1 GPU instance setup (vastai/pytorch image). Same recipe as runs/127_exact64_gpu/vast_setup.sh with the
# five fixes from its NOTES.md folded in, plus the `sigma` make target (the evGW W ladder needs sigma on the GPU).
# 1) NVIDIA HPC SDK  2) QE 7.4.1 CPU (gfortran/OpenMPI, patched pw2bgw)  3) HDF5 parallel+Fortran with nvhpc
# 4) ScaLAPACK with nvhpc  5) BerkeleyGW 4.0 GPU (cplx): epsilon sigma bse meanfield
set -uo pipefail
exec > >(tee -a /workspace/setup.log) 2>&1
export DEBIAN_FRONTEND=noninteractive
cd /workspace
echo "== $(date -u) apt"
apt-get update -qq
apt-get install -y -qq build-essential gfortran cmake pkg-config libopenblas-dev libfftw3-dev libfftw3-mpi-dev \
  libopenmpi-dev openmpi-bin git wget curl jq python3 zlib1g-dev >/dev/null
echo "== $(date -u) nvhpc"
curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' > /etc/apt/sources.list.d/nvhpc.list
apt-get update -qq
NVPKG=$(apt-cache search --names-only '^nvhpc-2[45]-[0-9]+$' | awk '{print $1}' | sort -V | tail -1)
echo "installing $NVPKG"; apt-get install -y -qq "$NVPKG" >/dev/null
NVROOT=$(ls -d /opt/nvidia/hpc_sdk/Linux_x86_64/*/ | sort -V | tail -1)
echo "NVROOT=$NVROOT" > /workspace/nv.env
echo "export PATH=$NVROOT/compilers/bin:$NVROOT/comm_libs/mpi/bin:\$PATH" >> /workspace/nv.env
echo "export LD_LIBRARY_PATH=$NVROOT/compilers/lib:$NVROOT/comm_libs/mpi/lib:$NVROOT/math_libs/lib64:$NVROOT/cuda/lib64:\${LD_LIBRARY_PATH:-}" >> /workspace/nv.env
echo "export NVROOT=$NVROOT" >> /workspace/nv.env
source /workspace/nv.env; nvfortran --version | head -2; which mpif90
NJ=$(nproc); [ "$NJ" -gt 48 ] && NJ=48
# ---- QE 7.4.1 (CPU, gfortran) in the background ----
( echo "== $(date -u) QE"; cd /workspace
  [ -d qe-7.4.1 ] || { wget -q https://gitlab.com/QEF/q-e/-/archive/qe-7.4.1/q-e-qe-7.4.1.tar.gz -O qe.tgz && tar xzf qe.tgz && mv q-e-qe-7.4.1 qe-7.4.1; }
  cp /workspace/inputs/pw2bgw.f90 qe-7.4.1/PP/src/pw2bgw.f90; cd qe-7.4.1 && ./configure MPIF90=/usr/bin/mpif90 F90=gfortran CC=gcc FFLAGS="-O3 -fallow-argument-mismatch" \
     BLAS_LIBS="-lopenblas" LAPACK_LIBS="" --with-scalapack=no --enable-openmp=no >/dev/null
  make -j$((NJ/2)) pw pp >/dev/null 2>&1 && ls -la bin/pw.x bin/pw2bgw.x && echo "QE OK" ) > /workspace/qe_build.log 2>&1 &
# ---- HDF5 with nvhpc ----
echo "== $(date -u) HDF5"; cd /workspace
wget -q https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_5/downloads/hdf5-1.14.5.tar.gz -O hdf5.tgz || wget -q https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz -O hdf5.tgz
tar xzf hdf5.tgz && cd hdf5-1.14.5
CC=mpicc FC=mpifort CFLAGS="-O2" FCFLAGS="-O2" ./configure --prefix=/workspace/hdf5 --enable-fortran --enable-parallel --disable-shared >/dev/null
make -j$((NJ/2)) >/dev/null 2>&1 && make install >/dev/null && echo "HDF5 OK"
# ---- ScaLAPACK with nvhpc ----
echo "== $(date -u) ScaLAPACK"; cd /workspace
wget -q https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v2.2.2.tar.gz -O scalapack.tgz && tar xzf scalapack.tgz
mkdir -p scalapack-2.2.2/build && cd scalapack-2.2.2/build
cmake .. -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=/workspace/scalapack \
  -DBLAS_LIBRARIES="$NVROOT/compilers/lib/libblas.so" -DLAPACK_LIBRARIES="$NVROOT/compilers/lib/liblapack.so" -DCMAKE_Fortran_FLAGS="-O2" >/dev/null
make -j$((NJ/2)) >/dev/null 2>&1 && make install >/dev/null && echo "ScaLAPACK OK"
# ---- BerkeleyGW GPU ----
echo "== $(date -u) BGW"; cd /workspace
[ -d BerkeleyGW-4.0 ] || tar xzf /workspace/inputs/BerkeleyGW-4.0-src.tgz
cd /workspace/BerkeleyGW-4.0
find . -name "*.o" -delete; find . -name "*.mod" -delete; find . -name "*.a" -delete; rm -rf bin; mkdir -p bin   # bin must exist before the per-target "ln -sf ../bin" links, else the first link becomes a FILE named bin
cp /workspace/inputs/arch_nvhpc_h100.mk arch.mk; cp flavor_cplx.mk flavor.mk
make -j$((NJ/2)) epsilon sigma bse meanfield > /workspace/bgw_build.log 2>&1 || { echo "BGW make rc=$?"; tail -40 /workspace/bgw_build.log; }; while pgrep -x make >/dev/null; do sleep 15; done
make install INSTDIR=/workspace/bgw >/dev/null 2>&1 || true
ls -la bin/ 2>/dev/null | grep -E "epsilon|sigma|kernel|absorption|wfn2hdf"
OK=1; for b in epsilon.cplx.x sigma.cplx.x kernel.cplx.x absorption.cplx.x wfn2hdf.x; do [ -x bin/$b ] || { echo "MISSING bin/$b"; OK=0; }; done
[ $OK = 1 ] && echo "BGW OK"
wait; echo "== $(date -u) SETUP DONE"; tail -3 /workspace/qe_build.log
