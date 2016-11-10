# Calculates the relative error between two .psi files

import numpy
import sys


# Print usage info if not right number of args
if len(sys.argv) < 3:
    print "Calculates the L1 and Linf relative error between two .psi files."
    print "Usage 1: python diff.py <file1> <file2>"
    print "   Prints differences."
    print "Usage 2: python diff.py <file1> <file2> <tolerance>"
    print "   Prints differences."
    print "   Returns 1 to OS if errors not within tolerance."
    sys.exit(1)


# Get the vectors from file
v1 = numpy.fromfile(sys.argv[1])
v2 = numpy.fromfile(sys.argv[2])


# Get rid of header data
v1 = v1[8:]
v2 = v2[8:]


# Check vector lengths
if len(v1) != len(v2):
    print "Vector lengths not the same."
    sys.exit(1)


# Get the tolerance
tolerance = -1.0
if len(sys.argv) == 4:
    tolerance = float(sys.argv[3])


# Find the errors
v1 = v1 - v2
errL1   = numpy.linalg.norm(v1, 1) / numpy.linalg.norm(v2, 1)
errLinf = max(abs(v1)) / max(abs(v2))


# Print the errors if either:
# tolerance not specified or errors don't meet tolerance
print "   L1   relative error:", "%.2e" % errL1
print "   Linf relative error:", "%.2e" % errLinf


# Check tolerance
if errL1 < tolerance and errLinf < tolerance and tolerance > 0.0:
    sys.exit(0)

sys.exit(1)

