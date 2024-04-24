# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
source /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/init/csh
module purge 
module load ncarenv/1.2 intel/17.0.1 esmf_libs mkl esmf-7.0.0-defio-mpi-O mpt/2.16 netcdf-mpi/4.5.0 pnetcdf/1.9.0 ncarcompilers/0.4.1
setenv OMP_STACKSIZE 256M
setenv TMPDIR /glade/scratch/tykukla
setenv MPI_TYPE_DEPTH 16
setenv MPI_IB_CONGESTED 1
setenv MPI_USE_ARRAY None
setenv TMPDIR /glade/scratch/tykukla
setenv MPI_USE_ARRAY false