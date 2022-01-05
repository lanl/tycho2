#/bin/bash -x 

export OMP_NUM_THREADS=1

# export OMP_PROC_BIND=spread
# export OMP_PLACE=threads

# Run $2 MPI ranks

#./util/PartitionMetis.x $2 $1/util/cube-208.smesh temp.pmesh
#jsrun --np $2 -k ty.err -o ty.out ./src/tycho2.x temp.pmesh $1/input.kokkos

./util/PartitionMetis.x $2 $1/util/cube-12.smesh temp.pmesh
jsrun --np $2 -k ty.err -o ty.out ./src/tycho2.x temp.pmesh $1/input.sn2-cells12
