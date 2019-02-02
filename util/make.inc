# Settings to change for local environment
# Copy this file to make.inc and fill in your local environment settings


# Compile command
# Note: the c++11 standard is required so this must be in the compile
#       command for some compilers.

### compiler
CPP = mpicxx -Wall -std=c++11


# Library directories
MOAB_DIR   = 
CGNS_DIR   = 
HDF5_DIR   = 

#LIB_METIS = -lmetis