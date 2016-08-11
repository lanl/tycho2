#This compares two files line by line to determine if they are within a given tolerance

tol = float(0.00001)


#prompts for filenames
file1 = raw_input("Enter file 1: ")
file2 = raw_input("Enter file 2: ")
#out = raw_input("Enter file to write results to: ")

#open files
f1=open(file1, 'r')
f2=open(file2, 'r')
#outfile = open(out, 'w')

#read first line of files
f1line = f1.readline()
f2line = f2.readline()

#initialize counters
counter = 1
numzeros = 0
numones = 0

#check end of file
while f1line != '' or f2line != '':
	#remove whitespace
	f1line=f1line.rstrip()
	f2line=f2line.rstrip()
	#compare	
	if abs(float(f1line)-float(f2line))<=tol:
	#	outfile.write('0\n')
		numzeros += 1
	elif abs(float(f1line)-float(f2line))>tol:
	#	outfile.write('1\n')
		numones += 1
	else:
	#	outfile.write('ERROR!')
		pass	 

	#read the next line
	f1line=f1.readline()
	f2line=f2.readline()

	#increment
	counter += 1

#outfile.close()
f1.close()
f2.close()

percent = numzeros/(numzeros+numones)

print("Number of Correct Results: " + str(numzeros) + '\n')
print("Number of Incorrect Results: " + str(numones) + '\n')
	
		

