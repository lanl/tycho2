# Perform regression tests on all configurations

# Move binaries to this folder
cp ../sweep.x ./
cp ../util/PartitionColumns.x ./
cp ../util/cube-208.smesh ./






### TESTS ###
echo " "
echo " "
echo "Running Regression Tests"

# gold
sh regression/run-gold.sh > out-gold.regression.txt
python diff.py regression/gold.psi out.psi 1e-10
if [ $? -eq 0 ]
    then echo "   Test gold: pass"
    else echo "   Test gold: fail"
fi


# gold-mpi2
sh regression/run-gold-mpi2.sh > out-gold-mpi2.regression.txt
python diff.py regression/gold.psi out.psi 1e-10
if [ $? -eq 0 ]
    then echo "   Test gold-mpi2: pass"
    else echo "   Test gold-mpi2: fail"
fi


# gold-mpi3
sh regression/run-gold-mpi3.sh > out-gold-mpi3.regression.txt
python diff.py regression/gold.psi out.psi 1e-10
if [ $? -eq 0 ]
    then echo "   Test gold-mpi3: pass"
    else echo "   Test gold-mpi3: fail"
fi


# gold-mpi4
sh regression/run-gold-mpi4.sh > out-gold-mpi4.regression.txt
python diff.py regression/gold.psi out.psi 1e-10
if [ $? -eq 0 ]
    then echo "   Test gold-mpi4: pass"
    else echo "   Test gold-mpi4: fail"
fi







echo " "
echo " "


# Remove binaries from this folder
rm sweep.x
rm PartitionColumns.x
rm cube-208.smesh
