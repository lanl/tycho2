# Setup to get lmod modules
module purge
unset dracomodules
unset NoModules
unset _LMFILES_
unset MODULEPATH
unset LOADEDMODULES
unset MODULESHOME

export MODULE_HOME=/scratch/vendors/spack.20170502/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/lmod-7.4.8-oytncsoih2sa4jdogz2ojvwly6mwle4n
source $MODULE_HOME/lmod/lmod/init/bash || die "Can't find /mod/init/bash"
module use /scratch/vendors/spack.20170502/share/spack/lmod/linux-rhel7-x86_64/Core
module use --append /scratch/vendors/Modules.lmod



# Load modules
module load gcc/6.3.0
module load openmpi/1.10.5
module load metis/5.1.0
module load moab/5.0.0
module load petsc/3.7.5-netlib

# List modules
module list
