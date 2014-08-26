#-----------------------------------------------------------------
#
# You might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh 
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the 
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the 
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=GIZMO
# 
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype 
#     This file has priority over your shell variable.:
#
#    (1) Copy the file "Template-Makefile.systype"  to  "Makefile.systype"
#
#        cp Template-Makefile.systype Makefile.systype 
#
#    (2) Uncomment your system in  "Makefile.systype".
#
# If you add an ifeq for a new system below, also add that systype to
# Template-Makefile.systype
#
###########
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel (volker.springel@h-its.org). The code has been modified
#   slighty by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly 
#   dealing with new files and filename conventions)
#
#############

ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif


CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat GIZMO_config.h)


CC       = mpicc        # sets the C-compiler (default)
CXX       = mpiCC       # sets the C++-compiler (default)

FC 	 = mpif90

OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES =  -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES =  -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES =  -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif

# we only need fftw if PMGRID is turned on
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
# or if POWERSPEC_GRID is activated
ifeq (POWERSPEC_GRID, $(findstring POWERSPEC_GRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
  FFTW_LIB =
endif

endif




#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Stampede")
CC       =  mpicc
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -O3 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -g -Wall # compiler warnings
#OPTIMIZE += -parallel -openmp  # openmp (comment out this line if OPENMP not used)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp  # openmp required compiler flags
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(TACC_MKL_INC)
MKL_LIBS = -L$(TACC_MKL_LIB) -mkl=sequential
##MKL_LIBS = -L$(TACC_MKL_LIB) -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
#MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
## modules to load: 
## module load intel mvapich2 gsl hdf5 fftw2
##  -- performance is very similar with impi (intel-mpi) instead of mpavich2, 
##   if preferred use that with MPICHLIB line uncommented
## newest version of code needed for compatibility with calls in MPI-2 libraries
##
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Zwicky")
CC       =  mpicc
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -O3 -funroll-loops
OPTIMIZE += -g -Wall # compiler warnings
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKL_HOME)/include
MKL_LIBS = -L$(MKL_HOME)/lib/em64t -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += # -DUSE_MPI_IN_PLACE
## modules to load: 
## module load intel/2011.4.191 impi/4.0.2.003 gsl/1.15-gcc HDF5 
##  -- the machine is quite picky, impi seems to be the only working mpi option right now
##  --  currently fftw2 isnt pre-installed, built library in my directory, with config flags:
##       ./configure --prefix=/home/phopkins/fftw --enable-mpi --enable-type-prefix --enable-float --with-gcc
##      linked via the above FFTW2_HOME=/home/phopkins/fftw (where the libraries are installed)
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"SciNet")
CC       =  mpicc     # sets the C-compiler
OPTIMIZE =  -O1 -xHost -funroll-loops -no-prec-div -fast-transcendentals -fp-model fast=2 -ipo  # speed
#OPTIMIZE += -openmp -parallel -par-num-threads=4  # for openmp mode
OPTIMIZE += -g -debug parallel -Wall  # warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE +=-openmp -parallel  # openmp required compiler flags
endif
FC       =  $(CC)
GSL_INCL =  -I${SCINET_GSL_INC}
GSL_LIBS =  -L${SCINET_GSL_LIB} #-limf
FFTW_INCL=  -I${SCINET_FFTW_INC}
FFTW_LIBS=  -L${SCINET_FFTW_LIB}
MPICHLIB =
HDF5INCL =  -I${SCINET_HDF5_INC} -DH5_USE_16_API
HDF5LIB  =  -L${SCINET_HDF5_LIB} -lhdf5 -lz
MPICHLIB =
##
## Notes:
## 
### benchmarking suggests these optimizations, 256 cores with omp=4 or 2, no DOUBLE, multidomain=16 or 32, altogether gives best performance in
###   simple galaxy merger experiment (6x speedup over 16-core run with old-but-highly-optimized code).
##
## module load intel use.experimental openmpi/intel/1.6.0 gsl fftw/2.1.5-intel-openmpi hdf5/intel-openmpi/1.8.9
## NOTYPEPREFIX_FFTW should not be set on this machine
## 
## flags: 
## OPT      += -DNOCALLSOFSYSTEM -DMPICH_IGNORE_CXX_SEEK -DNO_ISEND_IRECV_IN_DOMAIN
##   -- these used to be recommended, with new compiler settings they don't seem necessary, but may help
## If memory problems crash the code, recommend small-scale chunking: MPISENDRECV_SIZELIMIT=10-100 (but this costs speed!)
##
## old options with intelmpi (not as good as openmpi):
##    module load intel intelmpi gsl fftw/2.1.5-intel-intelmpi4 use.experimental hdf5/intelmpi/1.8.9
##    OPTIMIZE =  -O2 -m64 -mt_mpi -openmp -xhost -g -debug parallel -mcmodel=medium -funroll-loops -Wall
##
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ranger_intel")
CC       =  mpicc
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -O3 -xO -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -parallel -openmp  # openmp
OPTIMIZE += -g -Wall -debug parallel # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB =      # must be empty if using openmpi
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
##
## Notes:
## 
## include the following in your .bashrc file (there is no default fftw2 module):
## module load intel/10.1 openmpi/1.2.4 gmp gsl hdf5 #now have to add fftw2 manually
## export TACC_FFTW2_INC=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/include
## export TACC_FFTW2_LIB=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/lib
## export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
##
## Options
## OPT += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
##   are not necessary, but may improve stability in some cases
##
endif


ifeq ($(SYSTYPE),"Ranger_pgi")
CC       =  mpicc 
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -tp barcelona-64 -fast -Mipa=fast,inline -Munroll -Mvect -O4
OPTIMIZE += -mp -Mconcur  # openmp
OPTIMIZE += -Wall  # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPT     += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
## 
## Notes:
##
## include the following in your .bashrc file:
##   module load pgi mvapich gmp gsl fftw2 hdf5
##   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
## 
endif


ifeq ($(SYSTYPE),"Octopus")
CC       =  mpicc
CXX      =  mpiCC
OPTIMIZE = -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS= -lfftw3
MPICHLIB = -lmpi
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"odyssey_comparison")
CC       =  mpicc
#gnu
OPT      +=  -DNOCALLSOFSYSTEM  -DMPICH_IGNORE_CXX_SEEK  -DNOTYPEPREFIX_FFTW
OPTIMIZE =   -O2 -g -Wall -m64
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -I/n/sw/hdf5-1.8.5_gnu/include
HDF5LIB  =  -L/n/sw/hdf5-1.8.5_gnu/lib -lhdf5
#intel
#OPT      +=  -DMPICH_IGNORE_CXX_SEEK  -DNOTYPEPREFIX_FFTW
#OPTIMIZE =   -parallel -O3 -ipo -funroll-loops -Wall
#GSL_INCL =
#GSL_LIBS =
#FFTW_INCL=
#FFTW_LIBS=
#MPICHLIB =
#HDF5INCL =
#HDF5LIB  =  -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc     # sets the C-compiler
FCC      =  $(CC)
OPT      +=  -DNOCALLSOFSYSTEM  -DMPICH_IGNORE_CXX_SEEK  -DNOTYPEPREFIX_FFTW
#OPTIMIZE =   -O0 -g -Wall -m64 -Wno-unused-but-set-variable -Wno-uninitialized
OPTIMIZE =  -O2 -mcmodel=medium -funroll-loops -Wall
GSL_INCL =
GSL_LIBS =
FFTW_INCL=  -I/n/home12/mvogelsberger/opt/include/
FFTW_LIBS=  -L/n/home12/mvogelsberger/opt/lib/
MPICHLIB =
HDF5INCL = -DH5_USE_16_API
HDF5LIB  =  -lhdf5
endif


ifeq ($(SYSTYPE),"odyssey_old")
CC       =  mpicc
CXX      =  mpiCC 
OPT      +=  -DMPICH_IGNORE_CXX_SEEK  
OPTIMIZE =   -O3 -g -Wall -m64
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -I/n/sw/hdf5-1.8.5_gnu/include
HDF5LIB  =  -L/n/sw/hdf5-1.8.5_gnu/lib -lhdf5
#intel
#OPT      +=  -DMPICH_IGNORE_CXX_SEEK  -DNOTYPEPREFIX_FFTW
#OPTIMIZE =   -parallel -O3 -ipo -funroll-loops -Wall
#GSL_INCL =
#GSL_LIBS =
#FFTW_INCL=
#FFTW_LIBS=
#MPICHLIB =
#HDF5INCL =
#HDF5LIB  =  -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"Darwin")
CC       =  mpicc -g -O2 -Wall
CXX      =  mpiCC -g -O2 -Wall
OPTIMIZE =  
GSL_INCL = -I/sw/include 
GSL_LIBS = -L/sw/lib 
FFTW_INCL= -I/sw/include 
FFTW_LIBS= -L/sw/lib 
MPICHLIB = -lmpi 
HDF5INCL =
HDF5LIB  =
endif


# modules for Genius
# module load hdf5/1.8.4p1

ifeq ($(SYSTYPE),"Genius")
CC       =  mpixlc_r
CXX      =  mpixlcxx_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450
GSL_INCL =  -I/u/vrs/BlueGene/Libs/include
GSL_LIBS =  -L/u/vrs/BlueGene/Libs/lib -static
FFTW_INCL=  -I/u/vrs/BlueGene/Libs/include
FFTW_LIBS=  -L/u/vrs/BlueGene/Libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"BGP_toy")
CC       =  mpixlc_r
CXX      =  mpixlcxx_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450
GSL_INCL =  -I/gpfs/scratch/userinternal/cin8263a/lib/gsl-1.9/include/
GSL_LIBS =  -L/gpfs/scratch/userinternal/cin8263a/lib/gsl-1.9/lib/ -static
FFTW_INCL=  -I/gpfs/scratch/userinternal/cin8263a/lib/fftw-2.1.5/include
FFTW_LIBS=  -L/gpfs/scratch/userinternal/cin8263a/lib/fftw-2.1.5/lib
MPICHLIB =
HDF5INCL =  -I/gpfs/scratch/userinternal/cin0753a/hdf5_2/hdf5/hdf5-1.8.5-patch1/hdf5/include -I/gpfs/scratch/userinternal/cin0753a/zlib/zlib-1.2.5/include
HDF5LIB  =  -L/gpfs/scratch/userinternal/cin0753a/hdf5_2/hdf5/hdf5-1.8.5-patch1/hdf5/lib -static -lhdf5 -L/gpfs/scratch/userinternal/cin0753a/zlib/zlib-1.2.5/lib -lz
endif

# modules for Jugene
# module load gsl/1.14
# module load fftw/v2.1.5
# module load hdf5/1.8.4_450

ifeq ($(SYSTYPE),"Jugene")
CC       =  mpixlc_r
CXX      =  mpixlcxx_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450 -qsmp=omp:noauto -qreport
GSL_INCL =  -I$(GSL_DIR)/include
GSL_LIBS =  -L$(GSL_DIR)/lib
FFTW_INCL=  -I$(FFTW_DIR)/include
FFTW_LIBS=  -L$(FFTW_DIR)/lib
MPICHLIB =
HDF5INCL =  $(HDF5_INCLUDE)
HDF5LIB  =  -L/bgsys/local/hdf5/v1.8.4/lib_450 -lhdf5 -lhdf5_hl -L/bgsys/local/zlib/lib/ -lz -L/bgsys/local/szip/lib/ -lsz # $(HDF5_LIB)
ZLIB     =  -L/bgsys/local/zlib/lib/ -lz
SZIP     =  -L/bgsys/local/szip/lib/ -lsz
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"Solaris_amd64")
CC       = /opt/SUNWhpc/bin/mpicc
CXX      = /opt/SUNWhpc/bin/mpiCC
OPTIMIZE = -xO4 -xchip=native -xarch=amd64a
#OPTIMIZE =   -g -xarch=amd64a
GSL_INCL = $(GSL64INCL)
GSL_LIBS = $(GSL64LIB)
FFTW_INCL= $(FFTW64INCL)
FFTW_LIBS= $(FFTW64LIB)
FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
R_PATH   = $(R64LIB)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL)  $(GSL_LIBS)
HDF5INCL =
HDF5LIB  =
LIBS     =   -lm  $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB) $(R64LIB)
endif


ifeq ($(SYSTYPE),"MBM") 
CC       =   mpicc
CXX      =   mpicxx
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/sw/tap/include
GSL_LIBS =  -L/sw/tap/lib  -Xlinker -R -Xlinker /sw/tap/lib
FFTW_INCL=  -I/sw/tap/include
FFTW_LIBS=  -L/sw/tap/lib  -Xlinker -R -Xlinker /sw/tap/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif

# modules for Magny
# module add mvapich2/gcc/64/1.4.1-qlc

ifeq ($(SYSTYPE),"Magny") 
CC       =   mpicc
CXX      =   mpicxx
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/hits/tap/sw/libs/include
GSL_LIBS =  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
FFTW_INCL=  -I/hits/tap/sw/libs/include
FFTW_LIBS=  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
GMP_INCL =  -I/hits/tap/sw/libs/include
GMP_LIBS =  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
MPICHLIB =
HDF5INCL =  -I/hits/tap/sw/libs/include -DH5_USE_16_API
HDF5LIB  =  -L/hits/tap/sw/libs/lib -lhdf5 
OPT      +=  -DNOCALLSOFSYSTEM
OPT      +=  -DIMPOSE_PINNING
endif

ifeq ($(SYSTYPE),"OpenSuse")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=
MPICHLIB = -L/usr/lib/mpi/gcc/openmpi/lib -Xlinker -R -Xlinker /usr/lib/mpi/gcc/openmpi/lib -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"OpenSuse64")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=  
MPICHLIB = -L/usr/lib64/mpi/gcc/openmpi/lib64 -Xlinker -R -Xlinker /usr/lib64/mpi/gcc/openmpi/lib64 -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"SUNNYVALE")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpicxx
OPTIMIZE =  -g -O3 #-Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib64
FFTW_INCL= -I${MCKENZIE_FFTW_INC_PATH}
FFTW_LIBS= -L${MCKENZIE_FFTW_LIB_PATH}
MPICHLIB =  -L/opt/lam-7.1.2-intel/lib -lmpi
HDF5INCL = -I${MCKENZIE_HDF_INC_PATH}
HDF5LIB  = -L${MCKENZIE_HDF_LIB_PATH} -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"McKenzie")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpicxx
OPTIMIZE =  -g -O3 #-Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib
FFTW_INCL= -I/opt/fftw/intel_8.1/2.1.5/include
FFTW_LIBS= -L/opt/fftw/intel_8.1/2.1.5/lib #-ldrfftw_mpi
MPICHLIB = -L/opt/lam-7.1.2b24-g77/lib -lmpi
HDF5INCL = -I/opt/hdf5-oscar-1.6.4/include
HDF5LIB  = -L/opt/hdf5-oscar-1.6.4/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"CITA")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/include/gsl
GSL_LIBS =  -L/usr/lib/libgsl
FFTW_INCL=  -I/opt/fftw-2.1.5/include
FFTW_LIBS=  -L/opt/fftw-2.1.5/lib
MPICHLIB =  -L/usr/lib/libmpi
HDF5INCL =  -I/usr/include
HDF5LIB  =  -L/usr/lib/libhdf5 -static -lhdf5 -lz
endif 


ifeq ($(SYSTYPE),"Stella")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/home/schaye/libs/include
GSL_LIBS =  -L/home/schaye/libs/lib -static
FFTW_INCL=  -I/home/schaye/libs/include
FFTW_LIBS=  -L/home/schaye/libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif


# modules for OPA  
# module load intel
# module load impi
# module load fftw/2.1.5
# module load gsl/1.14
# module load hdf5-serial

ifeq ($(SYSTYPE),"OPA-Cluster64-Intel")  
CC       =   mpiicc
CXX      =   mpicxx
OPTIMIZE =  -O3 -g -Wall -m64 -wd981 -wd2259 -wd1572 -ipo4 -funroll-loops -finline -finline-functions -finline-limit=200
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR) -Xlinker -R -Xlinker $(GSL_LIBDIR)
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib -Xlinker -R -Xlinker $(FFTW_HOME)/lib
MPICHLIB =  
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB  = -L$(HDF5_HOME)/lib -Xlinker -R -Xlinker $(HDF5_HOME)/lib -lhdf5 -lz
#OPT      +=  -DOLD_HDF5
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"opus")
CC       =   mpicc
CXX      =   mpicxx
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/snfs1/cgwork/dgad0043/include
GSL_LIBS =  -L/snfs1/cgwork/dgad0043/lib
FFTW_INCL=  -I/snfs1/cgwork/dgad0043/include
FFTW_LIBS=  -L/snfs1/cgwork/dgad0043/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"Sauron-gcc")
CC       =   mpicc.gcc   # sets the C-compiler
OPTIMIZE =   -O3 -funroll-loops -march=k8 -msse2 -static
GSL_INCL =   -I/usr/local/gsl.gcc/include
GSL_LIBS =   -L/usr/local/gsl.gcc/lib -static -lgsl -lgslcblas
FFTW_INCL=   -I/usr/local/fftw.gcc/include
FFTW_LIBS=   -L/usr/local/fftw.gcc/lib -static -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Sauron")
CC       =  mpicc  -m64 # sets the C-compiler
CXX      =  mpiCC  -m64
OPTIMIZE =   -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


ifeq ($(SYSTYPE),"OpteronMPA-Intel")
CC       =  /usr/local/appl/openmpi-intel/bin/mpicc  -m64 # sets the C-compiler
CXX      =  /usr/local/appl/openmpi-intel/bin/mpiCC  -m64
OPTIMIZE =   -O3 -Wall
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"OpteronMPA-Gnu")
CC       =  /usr/local/appl/openmpi-gnu/bin/mpicc  -m64 # sets the C-compiler
CXX      =  /usr/local/appl/openmpi-gnu/bin/mpiCC  -m64
OPTIMIZE =   -O3 -Wall
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpiCC
OPTIMIZE =   -g -Wall -fopenmp
# GSL_INCL =  -I/usr/common/pdsoft/include
# GSL_LIBS =  -L/usr/common/pdsoft/lib
GSL_INCL =  -I/afs/mpa/home/volker/Libs/include
GSL_LIBS =  -L/afs/mpa/home/volker/Libs/lib
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -lhdf5 -lz 
OPT     +=  -DOLD_HDF5
endif


ifeq ($(SYSTYPE),"Mako")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpiCC
OPTIMIZE =   -O3 -march=athlon-mp  -mfpmath=sse
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


# modules for VIP
# module load fftw/2.1.5
# module load gsl/1.14
# module load hdf5

ifeq ($(SYSTYPE),"VIP")
CC       =  mpcc_r -g -qcpluscmt # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
CCX      =  mpCC_r -g
FC       =  mpxlf_r
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(CONFIGVARS)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPT      += -DVIP
FOPT     += -WF,-DVIP
OPTIMIZE =  -q64 -O2 -qarch=auto -qtune=auto -qinline -qstrict -qipa -qhot -qipa=level=2 -qlanglvl=extc99 -qalias=ansi -qinline=200 # -qfloat=hdflt
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR)
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
OPT     +=   -DPOWER6 -DDO_NOT_BRACH_IF -DPAD_STRUCTURES
endif

# modules for SARA
# module load c/ibm/11.1
# module load fftw2/2.1.5
# module load gsl/1.11
# module load hdf5

ifeq ($(SYSTYPE),"SARA")
CC       =  mpcc -g -qcpluscmt # -qsmp -qreport # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
CXX      =  mpCC -g
FC       =  mpxlf_r
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(CONFIGVARS)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPT      += -DVIP
FOPT     += -WF,-DVIP
OPTIMIZE =  -q64 -O5 -qarch=auto -qtune=auto -qinline -qstrict -qhot -qipa=level=2 -qlanglvl=extc99 -qalias=ansi -qinline=200 # -qfloat=hsflt
GSL_INCL = -I$(SARA_GSL_INCLUDE)
GSL_LIBS = -L$(SARA_GSL_LIB)
FFTW_INCL= -I$(SARA_FFTW2_INCLUDE)
FFTW_LIBS= -L$(SARA_FFTW2_LIB)
MPICHLIB =
HDF5INCL = -I$(SARA_HDF5_INCLUDE)
HDF5LIB  = -L$(SARA_HDF5_LIB) -lhdf5 -lz -static
OPT     +=   -DPOWER6 -DDO_NOT_BRACH_IF -DPAD_STRUCTURES
endif


ifeq ($(SYSTYPE),"CINECA_SP6")
CC       =   mpcc # sets the gcc C-compiler
CXX      =   mpCC
#to be used for deugging: -g -qcheck=all -qdbxextra -qheapdebug -qkeepparm -qoptdebug
OPTIMIZE =  -O5 -q64  -qhot -qstrict -qarch=pwr6 -qtune=pwr6 -qcpluscmt -qlanglvl=extc99 -qalias=ansi -qinline=200 -qipa=level=2 # :inline=noauto:list

GSL_INCL = -I/cineca/prod/libraries/gsl/1.9/xl--10.1/include
GSL_LIBS = -L/cineca/prod/libraries/gsl/1.9/xl--10.1/lib
FFTW_INCL= -I/cineca/prod/libraries/fftw/2.1.5/xl--10.1/include
FFTW_LIBS= -L/cineca/prod/libraries/fftw/2.1.5/xl--10.1/lib
HDF5INCL = -I/sp6/userexternal/sborgani/code/HDF5/include
HDF5LIB  = -L/sp6/userexternal/sborgani/code/HDF5/lib -lhdf5
MPICHLIB =
OPT     +=   -DPOWER6 -DDO_NOT_BRACH_IF -DPAD_STRUCTURES
endif


ifeq ($(SYSTYPE),"Ingeld_LinuxCluster")
CC       =   mpicc  # sets the C-compiler
CXX      =   mpiCC
OPTIMIZE =   -O3 -Wall

GSL_INCL = -I/home/patricia/my-gsl/include
GSL_LIBS = -L/home/patricia/my-gsl/lib -static
FFTW_INCL= -I/home/patricia/my-fftw/include
FFTW_LIBS= -L/home/patricia/my-fftw/lib
endif


ifeq ($(SYSTYPE),"Warp")
CC       =   /opt/intel/impi/3.1/bin64/mpiicc
CXX      =   /opt/intel/impi/3.1/bin64/mpiiCC
OPTIMIZE =  -xW
GSL_INCL = -I/opt/gsl/impi/include/gsl
GSL_LIBS = -L/opt/gsl/impi/lib64 -Wl,-rpath,/opt/gsl/impi/lib64
FFTW_INCL= -I/opt/fftw/impi_intel-3.1/include
FFTW_LIBS= -L/opt/fftw/impi_intel-3.1/lib64 -Wl,-rpath,/opt/fftw/impi_intel-3.1/lib64
MPICHLIB = -L/opt/intel/impi/3.1/lib64 -Wl,-rpath,/opt/intel/impi/3.1/lib64
endif


ifeq ($(SYSTYPE),"Kraken")
OPT     +=   -DNOCALLSOFSYSTEM
CC       =   cc # sets the C-compiler
OPTIMIZE =   -fast
GSL_INCL =  -I/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/include
GSL_LIBS =   -L/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/lib
FFTW_INCL=  -I/opt/fftw/2.1.5/cnos/include
FFTW_LIBS=  -L/opt/fftw/2.1.5/cnos/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
endif

ifeq ($(SYSTYPE),"Kraken_FFTW3")
OPT     +=   -DNOCALLSOFSYSTEM
CC       =   cc # sets the C-compiler
OPTIMIZE =   -fast
GSL_INCL =  -I/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/include
GSL_LIBS =   -L/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/lib
FFTW_INCL=  -I/sw/xt/fftw/3.3_alpha/cnl2.2_pgi9.0.3/include
FFTW_LIBS=  -L/sw/xt/fftw/3.3_alpha/cnl2.2_pgi9.0.3/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))    # fftw installed with type prefix?
  FFTW_LIBNAMES = $(FFTW_LIBS) -lfftw3_mpi -lfftw3_threads -lfftw3
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))
  FFTW_LIBNAMES =  -lfftw3_mpi -lfftw3_threads -lfftw3
else
  FFTW_LIBNAMES =  -lfftw3f_mpi -lfftw3f_threads -lfftw3f
endif
endif

endif

ifeq ($(SYSTYPE),"INAF_IRA")
OPT     += -DNOCALLSOFSYSTEM
CC       = mpicc # sets the C-compiler
OPTIMIZE = -O3 -Wall
GSL_INCL = -I/homes/donnert/Libs/include
GSL_LIBS = -L/homes/donnert/Libs/lib
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB =
HDF5INCL =
HDF5LIB  = -lhdf5 
endif


ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5INCL =
HDF5LIB  =
endif

SYSTEM_OBJS =   system/system.o system/allocate.o system/mymalloc.o system/parallel_sort.o \
                system/peano.o system/parallel_sort_special.o system/mpi_util.o

GRAVITY_OBJS  = gravity/forcetree.o gravity/cosmology.o gravity/pm_periodic.o gravity/potential.o \
                gravity/gravtree.o gravity/forcetree_update.o gravity/pm_nonperiodic.o gravity/longrange.o \
                gravity/ags_hsml.o

HYDRO_OBJS = hydro/hydra_master.o hydro/density.o hydro/gradients.o

STRUCTURE_OBJS = structure/twopoint.o

L3_OBJS =


OPTIONS = $(OPTIMIZE) $(OPT) 

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = GIZMO

OBJS  =  main.o accel.o  timestep.o init.o restart.o io.o \
         predict.o global.o begrun.o run.o allvars.o read_ic.o \
         domain.o driftfac.o kicks.o ngb.o compile_time_info.o merge_split.o

OBJS	+= $(GRAVITY_OBJS) $(HYDRO_OBJS) $(SYSTEM_OBJS) $(STRUCTURE_OBJS)
OBJS	+= $(L3_OBJS)

INCL    += allvars.h proto.h gravity/forcetree.h domain.h system/myqsort.h kernel.h Makefile \


ifeq (GRAIN_FLUID,$(findstring GRAIN_FLUID,$(CONFIGVARS)))
OBJS    += solids/grain_physics.o
endif

ifeq (GALSF,$(findstring GALSF,$(CONFIGVARS)))
OBJS    += galaxy_sf/sfr_eff.o
endif

ifeq (GALSF_FB_HII_HEATING,$(findstring GALSF_FB_HII_HEATING,$(CONFIGVARS)))
OBJS    += galaxy_sf/hII_heating.o
endif

ifeq (GALSF_FB_SNE_HEATING,$(findstring GALSF_FB_SNE_HEATING,$(CONFIGVARS)))
OBJS    += galaxy_sf/mechanical_fb.o
endif

ifeq (GALSF_FB_RPWIND_FROMSTARS,$(findstring GALSF_FB_RPWIND_FROMSTARS,$(CONFIGVARS)))
OBJS    += galaxy_sf/rp_localwinds.o
endif

ifeq (BLACK_HOLES,$(findstring BLACK_HOLES,$(CONFIGVARS)))
OBJS	+= galaxy_sf/blackhole.o
endif

ifeq (SINKS,$(findstring SINKS,$(CONFIGVARS)))
OBJS    += galaxy_sf/sinks.o
endif

ifeq (SCFPOTENTIAL,$(findstring SCFPOTENTIAL,$(CONFIGVARS)))
OBJS    += modules/potentials/scf.o modules/potentials/scf_util.o
endif

ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS    += structure/fof.o
INCL	+= structure/fof.h
endif

ifeq (OUTPUTLINEOFSIGHT,$(findstring OUTPUTLINEOFSIGHT,$(CONFIGVARS)))
OBJS    += structure/lineofsight.o
endif

ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS    += cooling/cooling.o
INCL	+= cooling/cooling.h
endif

ifeq (BUBBLES,$(findstring BUBBLES,$(CONFIGVARS)))
OBJS    += modules/bubbles/bubbles.o
endif

ifeq (EOS_DEGENERATE,$(findstring EOS_DEGENERATE,$(CONFIGVARS)))
OBJS	+= nuclear/helm_eos.o 
INCL	+= nuclear/helm_eos.h 
endif

ifeq (IMPOSE_PINNING,$(findstring IMPOSE_PINNING,$(CONFIGVARS)))
OBJS	+= system/pinning.o
endif

ifeq (JD_DPP,$(findstring JD_DPP,$(CONFIGVARS)))
OBJS	+= modules/cosmic_rays/cr_electrons.o
INCL	+= modules/cosmic_rays/cr_electrons.h
endif

ifeq (DISTORTIONTENSORPS,$(findstring DISTORTIONTENSORPS,$(CONFIGVARS)))
OBJS	+= modules/phasespace/phasespace.o modules/phasespace/phasespace_math.o
endif

ifeq (MACHNUM,$(findstring MACHNUM,$(CONFIGVARS)))
OBJS	+= machfinder.o
INCL	+= machfinder.h
endif

ifeq (RAD_TRANSFER,$(findstring RAD_TRANSFER,$(CONFIGVARS)))
OBJS	+= modules/rt/rt_chem.o modules/rt/rt_bh_lum.o modules/rt/rt_sfr_lum.o modules/rt/rt_cooling.o \
    modules/rt/rt_eddington.o modules/rt/rt_n.o modules/rt/rt_CGmethod.o modules/rt/rt_stars_lum.o modules/rt/rt_gas_lum.o
endif

ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS	+= subfind/subfind.o subfind/subfind_vars.o subfind/subfind_collective.o subfind/subfind_serial.o subfind/subfind_so.o subfind/subfind_cont.o \
	subfind/subfind_distribute.o subfind/subfind_findlinkngb.o subfind/subfind_nearesttwo.o subfind/subfind_loctree.o subfind/subfind_alternative_collective.o subfind/subfind_reshuffle.o \
	subfind/subfind_potential.o subfind/subfind_density.o
INCL	+= subfind/subfind.h
endif

ifeq (COSMIC_RAYS,$(findstring COSMIC_RAYS,$(CONFIGVARS)))
OBJS	+= modules/cosmic_rays/cosmic_rays.o modules/cosmic_rays/cosmic_rays_diffusion.o modules/cosmic_rays/greenf_diffusion.o
INCL	+= modules/cosmic_rays/cosmic_rays.h
endif

ifeq (SIDM,$(findstring SIDM,$(CONFIGVARS)))
OBJS    +=  sidm/sidm_core.o sidm/sidm_allvars.o
INCL    +=
endif

ifeq (NUCLEAR_NETWORK,$(findstring NUCLEAR_NETWORK,$(CONFIGVARS)))
OBJS	+=  nuclear/utilities.o nuclear/integrate.o nuclear/network_solver.o nuclear/network.o nuclear/swap.o
INCL	+=  nuclear/utilities.h nuclear/integrate.h nuclear/network_solver.h nuclear/network.h nuclear/swap.h
endif

ifeq (TURB_DRIVING,$(findstring TURB_DRIVING,$(CONFIGVARS)))
OBJS	+= turb/turb_driving.o turb/turb_powerspectra.o
endif

ifeq (BP_REAL_CRs,$(findstring BP_REAL_CRs,$(CONFIGVARS))) # add bp cr part
OBJS += bp_cosmic_rays/bp_cosmic_rays.o
INCL += bp_cosmic_rays/bp_cosmic_rays.h
endif

ifeq (ADJ_BOX_POWERSPEC,$(findstring ADJ_BOX_POWERSPEC,$(CONFIGVARS)))
OBJS += modules/power_spec/adj_box_powerspec.o 
INCL += modules/power_spec/adj_box_powerspec_proto.h
endif

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(GMP_INCL)

ifeq (VIP,$(findstring VIP,$(CONFIGVARS)))
FFLAGS = $(FOPTIONS)
else
FFLAGS = $(OPTIONS)
endif


ifeq (ALTERNATIVE_PSORT,$(findstring ALTERNATIVE_PSORT,$(CONFIGVARS)))
OBJS  += fof_alt_psort.o modules/psort-1.0/error_handling.o
CXXFLAGS = $(CFLAGS)
FC    = $(CXX)
endif

FFTW = $(FFTW_LIBS)  $(FFTW_LIBNAMES) 


LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW)

ifeq (OMP_NUM_THREADS,$(findstring OMP_NUM_THREADS,$(CONFIGVARS))) 
LIBS   +=  -lpthread
endif

$(EXEC): $(OBJS) $(FOBJS)  
	$(FC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)  $(CONFIG)  compile_time_info.c


$(FOBJS): $(FINCL)

compile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) *.oo *.c~ compile_time_info.c GIZMO_config.h


