# Makefile settings for CCS servers


# Assert level: 0 = no asserts, 1 = all asserts ################################
ASSERT_ON = 1


# PETSC library ################################################################
USE_PETSC = 1   # 0 = no, 1 = yes
PETSC_INC = 
PETSC_LIB = -lpetsc



# Extra libraries ##############################################################
INC = 
LIBS = 


# Compile command ##############################################################
# Needs (1) MPI C++ compiler wrapper
#       (2) openmp enabled
#       (3) c++11 standard enabled  (-std=c++11 for gcc, intel)
MPICXX = mpicxx -Wall -Wextra -g -std=c++11 -O3 -fopenmp
MPICC = mpicc -Wall -Wextra -g -O3 -fopenmp
