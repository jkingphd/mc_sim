# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 14:58:58 2015
main.py: Monte-carlo simulation of bacteria in a gradient of L-Aspertate.
Import appropriate lookup tables calculate c(x,t) and dc(x,t).
@author: jkk3
"""

import numpy as np
import matplotlib.pyplot as plt

# Define parameters of simulation
L = 1000. # Length of mixing region in µm
a = 7010. # Left boundary of initial position sampling in µm
b = a + L # Right boundary of initial position sampling in µm
D = 590 # Diffusion coefficient of E-coli in µm^2/s
dt = 1 # seconds
t_final = 3600 # seconds
n_trials = 1000 # number of simulations
sigma = np.sqrt(2*D)

# Import lookup tables
scale = 0.1
x_table = np.load("tables\L1_x.npy")
c_table = np.load("tables\L1_c.npy")*scale
dc_table = np.load("tables\L1_dc.npy")*1E-6*scale # Fix this in import tables later.

def int_tables(pos, t, x_table, c_table, dc_table):
    '''Interpolate values for c(pos,t) and dc(pos,t)'''
    idx = (np.abs(x_table - pos)).argmin()
    if pos == x_table[0]:
        c = c_table[0,t]
        dc = dc_table[0,t]
    elif pos == x_table[-1]:
        c = c_table[-1,t]
        dc = dc_table[-1,t]
    else:
        if x_table[idx] > pos:
            x1 = x_table[idx-1]; x2 = x_table[idx];
            c1 = c_table[idx-1,t]; c2 = c_table[idx,t];
            dc1 = dc_table[idx-1,t]; dc2 = dc_table[idx,t];
        else:
            x1 = x_table[idx]; x2 = x_table[idx+1];
            c1 = c_table[idx,t]; c2 = c_table[idx+1,t];
            dc1 = dc_table[idx-1,t]; dc2 = dc_table[idx,t];
        c = ((c2-c1)/(x2-x1))*(pos-x1) + c1
        dc = ((dc2-dc1)/(x2-x1))*(pos-x1) + dc1
    return c, dc
    
def calc_v(pos, t, x_table, c_table, dc_table):
    '''Calculate drift velocity as a function of c(pos,t) and dc(pos,t)'''
    v = 22. # µm/s
    chi = 50000. # µm^2/s
    k = 0.125 # mM
    if pos < x_table[0]:
        pos = x_table[0]
        c, dc = int_tables(pos, t, x_table, c_table, dc_table)
    elif pos > x_table[-1]:
        pos = x_table[-1]
        c, dc = int_tables(pos, t, x_table, c_table, dc_table)
    else:
        c, dc = int_tables(pos, t, x_table, c_table, dc_table)
    return ((8*v)/(3*np.pi))*np.tanh(((chi*np.pi)/(8*v))*(k/np.power((k+c),2))*dc)

data = np.zeros((t_final,n_trials))
for i in range(n_trials):
    # Set start position of bacteria and reset start time
    data[0,i] = np.random.uniform(a,b)
    t = 1
    while(t < t_final):
        data[t,i] = data[t-1,i] + np.random.normal(0.0, sigma, 1) +\
        calc_v(data[t-1,i], t-1, x_table, c_table, dc_table)*dt
#        if wb == 0:
#            t += 1
#        elif wb == 1:
#            if data[t,i] < x1:
#                data[t,i] = x1 + abs(data[t-1,i] - data[t,i])
#            else:
#                t += 1
#        elif wb == 2:
#            if data[t,i] < x1:
#                data[t:,i] = x1
#                break
#            else:
#                t += 1
#        else:
        t += 1
        
#t = 500
#fig = plt.figure("Histogram", (9,6))
#ax = fig.add_subplot(111)
#n = ax.hist(data[t,:], bins = 100, range = (x_table[0],x_table[-1]))[0]
#ax.axis([x_table[0],x_table[-1],0,1.5*n.max()])
#ax.grid(True)
#ax.axvspan(a, b, color = 'gray', alpha = 0.3, lw = 0)
#plt.show()

times = np.array([1,2,5,10,20,100,200,500,1000,2000])

for t in times:
    fig = plt.figure('Test', (9,6))
    ax1 = fig.add_subplot(111)
    ax1.axvspan(a, b, color = 'black', alpha = 0.1, lw = 0)
    ax1.hist(data[t,:], bins = 100, range = (x_table[0], x_table[-1]), lw = 0, color = 'blue')
    ax1.axis([x_table[0],x_table[-1],0,n_trials/2])
    ax1.set_ylabel('Population', color = 'blue', size = 16)
    ax1.set_title('t = %05d s' % t, size = 16)
    for tl in ax1.get_yticklabels():
        tl.set_color('blue')
    
    ax2 = ax1.twinx()
    ax2.plot(x_table, c_table[:,t], color = 'red', lw = 1)
    ax2.scatter(x_table, c_table[:,t], color = 'red', s = 10, lw = 0)
    ax2.axis([x_table[0],x_table[-1],0,c_table.max()])
    ax2.grid(True)
    ax2.set_ylabel('Concentration (mM)', color = 'red', size = 16)
    for t2 in ax2.get_yticklabels():
        t2.set_color('red')
    
    plt.savefig('Test_%d' % t)
    plt.close()

#plt.savefig('%05d.png' % t)
#plt.close()

##plt.figure("Kymograph", (8,8))
#graph = []
#for k in np.arange(0,900,10):
#    graph.append(np.histogram(data[k,:], bins = 100, range = (x_table[0],x_table[-1]))[0])
#graph = np.array(graph)
#np.savetxt('test', graph, delimiter = ',')
##plt.imshow(graph, cmap = 'brg', interpolation = 'none')
##plt.show()
#
#while(True):
#    ans = raw_input("Save dataset to file? (Y/N): ")
#    if ans.upper() == 'Y':
#        fname = raw_input("Enter filename: ")
#        np.savetxt(fname, graph, delimiter = ',')
#        break
#    elif ans.upper() == 'N':
#        break
#    else:
#        print "Enter Y or N."