import subprocess


tolerance = "1e-10"


# Print what we're doing
print " "
print "--- Running Regression Tests ---"


# Move necessary files to this folder
subprocess.call(["cp", "../sweep.x", "./"])
subprocess.call(["cp", "../util/PartitionColumns.x", "./"])
subprocess.call(["cp", "../util/cube-208.smesh", "./"])


# Get all run*.sh files
output = subprocess.check_output(["ls", "regression"])
output = output.split()
output1 = []
for s in output:
    if s.endswith(".sh") > 0:
        output1.append(s)


# Do regressions
numFail = 0
numPass = 0
for s in output1:
    name = "regression/" + s
    name2 = s + ".regression.txt"
    
    print "Test", s
    subprocess.check_output(["sh", name, ">", name2])
    status = subprocess.call(["python", "diff.py", "regression/gold.psi", "out.psi", tolerance])
    if status == 0:
        print "                                                       Pass"
        print " "
        numPass = numPass + 1
    else:
        print "                                                       Fail"
        print " "
        numFail = numFail + 1
    
    subprocess.call(["rm", "out.psi"])



# Print overall stats
print "Num Pass/Fail", numPass, numFail
print " "


# Remove un-necessary files
subprocess.call(["rm", "sweep.x"])
subprocess.call(["rm", "PartitionColumns.x"])
subprocess.call(["rm", "cube-208.smesh"])
