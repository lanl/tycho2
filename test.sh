# Run one MPI rank
mpirun -n 1 ./sweep.x test/cube-208-1.pmesh test/input-kokkos.deck
mpirun -n 1 ./sweep.x test/cube-208-1.pmesh test/input-no-kokkos.deck
mv kokkos.psi kokkos-1.psi
mv no-kokkos.psi no-kokkos-1.psi

# Run two MPI ranks
mpirun -n 2 ./sweep.x test/cube-208-2.pmesh test/input-kokkos.deck
mpirun -n 2 ./sweep.x test/cube-208-2.pmesh test/input-no-kokkos.deck
mv kokkos.psi kokkos-2.psi
mv no-kokkos.psi no-kokkos-2.psi

# Print space
echo "  "
echo "  "

# Diff runs from gold
python test/diff.py test/gold.psi kokkos-1.psi
python test/diff.py test/gold.psi no-kokkos-1.psi
python test/diff.py test/gold.psi kokkos-2.psi
python test/diff.py test/gold.psi no-kokkos-2.psi

# Diff runs against each other
python test/diff.py kokkos-1.psi no-kokkos-1.psi
python test/diff.py kokkos-2.psi no-kokkos-2.psi

# Remove output files
rm kokkos-1.psi
rm no-kokkos-1.psi
rm kokkos-2.psi
rm no-kokkos-2.psi

