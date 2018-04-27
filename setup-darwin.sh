module load openmpi/2.1.2-gcc_5.4.0 cuda/9.0
echo "Make sure to change kokkos/bin/nvcc_wrapper line: host_compiler=mpicxx"
echo "salloc -N 1 --constraint=\"gpu1_model:Tesla_K40m\""
