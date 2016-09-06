# Perform regression tests on all configurations

# Move binaries to this folder
cp ../sweep.x ./
cp ../util/PartitionColumns.x ./
cp ../util/cube-208.smesh ./






### TEST 1 ###
sh regression/run1.sh > out1.regression.txt
#diff -q regression/gold.psi out.psi
python diff.py regression/gold.psi out.psi 1e-10
if [ $? -eq 0 ]
    then echo "Test1: pass"
    else echo "Test1: fail"
fi






# Remove binaries from this folder
rm sweep.x
rm PartitionColumns.x
rm cube-208.smesh
