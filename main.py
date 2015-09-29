# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 14:58:58 2015
main.py: Monte-carlo simulation of bacteria in fixed gradient of L-Aspertate.
Attempting to replicate results from "Microfluidic study of the chemotactic
response of Escherichia coli to amino acids, signaling molecules and
secondary metabolites." Fig 4a.
@author: jkk3
"""

import numpy as np
#import matplotlib.pyplot as plt

# Define parameters of simulation
x1 = -600.0 # Left boundary (source) in µm
x2 = 600.0 # Right boundary (sink) in µm
c1 = 0.1 # Concentration of L-aspartate at source in mM
c2 = 0.0 # Concentration of L-adpartate at sink in mM
a = -600.0 # Left boundary of initial position sampling in µm
b = 600.0 # Right boundary of initial position sampling in µm
wb = 2 # Wall behavior flag. 0) None 1) Bounce 2) Stick
params = (x1, x2, c1, c2)

D = 590 # Diffusion coefficient of E-coli in µm^2/s
dt = 1 # seconds
t_final = 901 # seconds
n_trials = 1000 # number of simulations
sigma = np.sqrt(2*D)

def calc_v(pos, params):
    x1, x2, c1, c2 = params
    v = 22. # µm/s
    chi = 50000. # µm^2/s
    k = 0.125 # mM
    if pos < x1: # in source
        dc = 0.0
        c = 0.1
    elif pos > x2: # in sink
        dc = 0.0
        c = 0.0
    else:
        dc = (c2 - c1)/(x2 - x1)
        c = dc*(pos - x1) + c1
    return ((8*v)/(3*np.pi))*np.tanh(((chi*np.pi)/(8*v))*(k/np.power((k+c),2))*dc)

data = np.zeros((t_final,n_trials))
for i in range(n_trials):
    # Set start position of bacteria and reset start time
    data[0,i] = np.random.uniform(a,b)
    t = 1
    while(t < t_final):
        data[t,i] = data[t-1,i] + np.random.normal(0.0, sigma, 1) +\
        calc_v(data[t-1,i], params)*dt
        if wb == 0:
            t += 1
        elif wb == 1:
            if data[t,i] < x1:
                data[t,i] = x1 + abs(data[t-1,i] - data[t,i])
            else:
                t += 1
        elif wb == 2:
            if data[t,i] < x1:
                data[t:,i] = x1
                break
            else:
                t += 1
        else:
            t += 1
        
#a = 500
#fig = plt.figure("Histogram", (9,6))
#ax = fig.add_subplot(111)
#n = ax.hist(data[a,:], bins = 120, range = (x1,x2))[0]
#ax.axis([x1,x2,0,1.5*n.max()])
#
#plt.show()

#plt.figure("Kymograph", (8,8))
graph = []
for k in np.arange(0,900,10):
    graph.append(np.histogram(data[k,:], bins = 120, range = (x1,x2))[0])
graph = np.array(graph)
#plt.imshow(graph, cmap = 'brg', interpolation = 'none')
#plt.show()

while(True):
    ans = raw_input("Save dataset to file? (Y/N): ")
    if ans.upper() == 'Y':
        fname = raw_input("Enter filename: ")
        np.savetxt(fname, graph, delimiter = ',')
        break
    elif ans.upper() == 'N':
        break
    else:
        print "Enter Y or N."