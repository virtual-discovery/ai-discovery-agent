# BerkeleyGW 4.0 arch.mk: NVIDIA HPC SDK (nvfortran) + its CUDA-aware OpenMPI, H100 (cc80), OpenACC GPU build.
COMPFLAG  = -DNVHPC -DNVHPC_API -DNVIDIA_GPU
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DOMP_TARGET -DOPENACC
NVCC    = nvcc
NVCCOPT = -O3 -use_fast_math
CUDALIB = -lcufft -lcublasLt -lcublas -lcudart -lcuda
FCPP    = /usr/bin/cpp -C -nostdinc
F90free = mpif90 -Mfree -acc -mp=multicore,gpu -gpu=cc80 -cudalib=cublas,cufft -traceback -Minfo=accel -gopt
LINK    = mpif90 -acc -mp=multicore,gpu -gpu=cc80 -cudalib=cublas,cufft
FOPTS   = -fast -Mfree -Mlarge_arrays
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I
C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpic++
C_COMP  = mpicc
C_LINK  = mpic++ -lstdc++
C_OPTS  = -fast -mp
C_DEBUGFLAG =
REMOVE  = /bin/rm -f
FFTWLIB      = -L/usr/lib/x86_64-linux-gnu -lfftw3_omp -lfftw3_threads -lfftw3 $(CUDALIB) -lstdc++
FFTWINCLUDE  = /usr/include
PERFORMANCE  =
SCALAPACKLIB = /workspace/scalapack/lib/libscalapack.a
LAPACKLIB    = -llapack -lblas
HDF5_LDIR    = /workspace/hdf5/lib
HDF5LIB      = $(HDF5_LDIR)/libhdf5hl_fortran.a $(HDF5_LDIR)/libhdf5_hl.a $(HDF5_LDIR)/libhdf5_fortran.a $(HDF5_LDIR)/libhdf5.a -lz -ldl
HDF5INCLUDE  = /workspace/hdf5/include
TESTSCRIPT = MPIRUN
