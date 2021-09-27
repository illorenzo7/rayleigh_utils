import numpy as np
import sys, os
sys.path.append(os.environ['raco'])

f = open('polytrope_10000.txt', 'r')

lines = f.readlines()
xi = []
theta = []
dtheta = []

for i in range(len(lines)):
    line = lines[i].strip().split()
    xi.append(float(line[0]))
    theta.append(float(line[1]))
    dtheta.append(float(line[2]))

f.close()

xi = np.array(xi)
theta = np.array(theta)
dtheta = np.array(dtheta)
