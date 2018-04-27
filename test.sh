export OMP_NUM_THREADS=1
# export OMP_PROC_BIND=spread
# export OMP_PLACE=threads

# Run one MPI rank
mpirun -n 1 ./sweep.x test/cube-208-1.pmesh test/input.deck
mv test.psi test-1.psi

# Run two MPI ranks
# mpirun -n 2 ./sweep.x test/cube-208-2.pmesh test/input.deck
# mv test.psi test-2.psi

# Print space
echo "  "
echo "  "

# Diff runs from gold
python test/diff.py test/gold.psi test-1.psi
# python test/diff.py test/gold.psi test-2.psi

# Remove output files
# rm test-1.psi
# rm test-2.psi

