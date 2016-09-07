

NX=2
NY=2
NUM_PARTS=$((NX*NY))
IN_FILE="cube-208.smesh"
OUT_FILE="temp.pmesh"
INPUT_DECK="regression/input-intra4.deck"
export OMP_NUM_THREADS=1

./PartitionColumns.x $NX $NY $IN_FILE $OUT_FILE
mpirun -n $NUM_PARTS ./sweep.x $OUT_FILE $INPUT_DECK
rm $OUT_FILE
