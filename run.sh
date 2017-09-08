# This is an example run script to automate creating a parallel mesh and
# running the sweeper on it.


NUM_PARTS=24
IN_FILE="./util/cube-4128.smesh"
OUT_FILE="temp.pmesh"
INPUT_DECK="input.deck"
export OMP_NUM_THREADS=1

./util/PartitionMetis.x $NUM_PARTS $IN_FILE $OUT_FILE
mpirun -n $NUM_PARTS ./sweep.x $OUT_FILE $INPUT_DECK
