export OMP_NUM_THREADS=1
# export OMP_PROC_BIND=spread
# export OMP_PLACE=threads

# Run one MPI rank
mpirun -n 1 ./sweep.x test/cube-67249-1.pmesh test/input-big.deck --kokkos-device=1
mv test.psi test-1.psi

# Print space
echo "  "
echo "  "

