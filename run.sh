#/bin/bash -x 

export OMP_NUM_THREADS=1

# export OMP_PROC_BIND=spread
# export OMP_PLACE=threads

# Run one MPI rank

./util/PartitionMetis.x $1 ~/tycho2/util/cube-208.smesh temp.pmesh
mpirun -n $1 ./src/sweep.x temp.pmesh input.kokkos
