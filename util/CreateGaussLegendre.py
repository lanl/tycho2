import numpy as np

MAX_DEGREE = 30
file = open('GaussLegendre.dat', 'w')

for deg in range(2,MAX_DEGREE+1,2):
    x, w = np.polynomial.legendre.leggauss(deg)
    
    # Write nodes
    file.write('static double s_n' + str(deg) + '[] = { ' + format(x[0], '.16f'))
    for k in range(1,deg):
        file.write(', ' + format(x[k], '.16f'))
    file.write(' };\n')
    
    # Write weights
    file.write('static double s_w' + str(deg) + '[] = { ' + format(w[0], '.16f'))
    for k in range(1,deg):
        file.write(', ' + format(w[k], '.16f'))
    file.write(' };\n')
    
