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
L = 6000. # Length of mixing region in µm
a = 7010. # Left boundary of initial position sampling in µm
b = a + L # Right boundary of initial position sampling in µm
D = 590 # Diffusion coefficient of E-coli in µm^2/s
dt = 1 # seconds
t_final = 3600 # seconds
n_trials = 1000 # number of simulations
sigma = np.sqrt(2*D)

# Import lookup tables
s_attract = 0.1
s_repel = 0.25
x_table = np.load("tables\L6_x.npy")
c_attract = np.load("tables\L6_c.npy")*s_attract
dc_attract = np.load("tables\L6_dc.npy")*s_attract # Fix this in import tables later.
c_repel = np.flipud(np.load("tables\L6_c.npy")*s_repel)
dc_repel = np.flipud(np.load("tables\L6_dc.npy")*s_repel)

def int_tables(pos, t, x_table, c_attract, dc_attract):
    '''Interpolate values for c(pos,t) and dc(pos,t)'''
    idx = (np.abs(x_table - pos)).argmin()
    if pos == x_table[0]:
        c = c_attract[0,t]
        dc = dc_attract[0,t]
    elif pos == x_table[-1]:
        c = c_attract[-1,t]
        dc = dc_attract[-1,t]
    else:
        if x_table[idx] > pos:
            x1 = x_table[idx-1]; x2 = x_table[idx];
            c1 = c_attract[idx-1,t]; c2 = c_attract[idx,t];
            dc1 = dc_attract[idx-1,t]; dc2 = dc_attract[idx,t];
        else:
            x1 = x_table[idx]; x2 = x_table[idx+1];
            c1 = c_attract[idx,t]; c2 = c_attract[idx+1,t];
            dc1 = dc_attract[idx-1,t]; dc2 = dc_attract[idx,t];
        c = ((c2-c1)/(x2-x1))*(pos-x1) + c1
        dc = ((dc2-dc1)/(x2-x1))*(pos-x1) + dc1
    return c, dc
    
def calc_v(pos, t, x_table, c_attract, dc_attract, c_repel, dc_repel):
    '''Calculate drift velocity as a function of c(pos,t) and dc(pos,t)'''
    v = 22. # µm/s
    chi = 50000. # µm^2/s
    chi2 = 50000. # µm^2/s
    k = 0.125 # mM
    k2 = 0.125 # mM
    if pos < x_table[0]:
        pos = x_table[0]
        c, dc = int_tables(pos, t, x_table, c_attract, dc_attract)
        c2, dc2 = int_tables(pos, t, x_table, c_repel, dc_repel)
    elif pos > x_table[-1]:
        pos = x_table[-1]
        c, dc = int_tables(pos, t, x_table, c_attract, dc_attract)
        c2, dc2 = int_tables(pos, t, x_table, c_repel, dc_repel)
    else:
        c, dc = int_tables(pos, t, x_table, c_attract, dc_attract)
        c2, dc2 = int_tables(pos, t, x_table, c_repel, dc_repel)
    v1 = ((8*v)/(3*np.pi))*np.tanh(((chi*np.pi)/(8*v))*(k/np.power((k+c),2))*dc)
    v2 = ((8*v)/(3*np.pi))*np.tanh(((chi2*np.pi)/(8*v))*(k2/np.power((k+c2),2))*dc2)
    #print v1, v2
    if np.abs(v1) > np.abs(v2):
        return v1
    else:
        return v2

def pop_plot(data, thresh = a + L + 10., save = False, fname = '0.png'):
    '''Plot percentage of population in target zone.'''
    t = np.arange(0,np.shape(data)[0])
    frac = 100*np.sum(data >= thresh, axis = 1)/np.float(np.shape(data)[1])
    fig = plt.figure("Population Plot", (9,6))
    ax = fig.add_subplot(111)
    #ax.scatter(t, frac, c = 'black', lw = 0, s = 5)
    ax.plot(t, frac, c = 'red', lw = 1)
    idx = (np.abs(frac - 75.0)).argmin()
    ax.hlines(75.0, t[0], t[-1], color = 'black', linestyle = 'dashed')
    ax.text(idx, 70.0, '$t_{75} = %d$' % idx, fontsize = 16)
    idx = (np.abs(frac - 50.0)).argmin()
    ax.hlines(50.0, t[0], t[-1], color = 'black', linestyle = 'dashed')
    ax.text(idx, 45.0, '$t_{50} = %d$' % idx, fontsize = 16)
    ax.grid(True)
    ax.axis([t[0],t[-1],0.0,100.0])
    ax.set_ylabel("Population (%)", size = 16)
    ax.set_xlabel("Time (s)", size = 16)
    ax.set_title("Population vs. Time")
    # Choose whether or not to save the file.
    if save == False:
        plt.show()
    else:
        plt.savefig(fname)
        plt.close()
    return None
    
def k_plot(data, save = False, fname = '0.png'):
    '''Generate Kymograph of population in time.'''
    fig = plt.figure("Kymograph", (8,8))
    ax = fig.add_subplot(111)
    graph = []
    for k in np.arange(0,3600,10):
        graph.append(np.histogram(data[k,:], bins = np.floor(x_table[-1]/100), range = (x_table[0],x_table[-1]))[0])
    graph = np.array(graph)
    ax.imshow(graph, cmap = 'hot', interpolation = 'none', aspect = 'auto', extent = (x_table[0],x_table[-1],3600,0))
    ax.set_xlabel('Position ($\mu m$)', color = 'black', size = 16)
    ax.set_ylabel('Time (s)', color = 'black', size = 16)
    # Choose whether or not to save the file.
    if save == False:
        plt.show()
    else:
        plt.savefig(fname)
        plt.close()
    return None
    
def hist_plot(data, t, save = False, fname = '0.png'):
    '''Plot population distribution and concentraction at time t.'''
    fig = plt.figure('Histogram', (9,6))
    ax1 = fig.add_subplot(111)
    ax1.axvspan(a, b, color = 'black', alpha = 0.1, lw = 0)
    ax1.hist(data[t,:], bins = 100, range = (x_table[0], x_table[-1]), lw = 0, color = 'blue')
    ax1.axis([x_table[0],x_table[-1],0,n_trials/2])
    ax1.set_xlabel('Position ($\mu m$)', color = 'black', size = 16)
    ax1.set_ylabel('Population', color = 'blue', size = 16)
    ax1.set_title('t = %05d s' % t, size = 16)
    for tl in ax1.get_yticklabels():
        tl.set_color('blue')
        
    ax2 = ax1.twinx()
    ax2.plot(x_table, c_attract[:,t], color = 'green', lw = 1)
    ax2.scatter(x_table, c_attract[:,t], color = 'green', s = 10, lw = 0)
    ax2.plot(x_table, c_repel[:,t], color = 'red', lw = 1)
    ax2.scatter(x_table, c_repel[:,t], color = 'red', s = 10, lw = 0)
    c_max = np.max([c_attract.max(), c_repel.max()])
    ax2.axis([x_table[0],x_table[-1],0,c_max])
    ax2.grid(True)
    ax2.set_ylabel('Concentration (mM)', color = 'black', size = 16)
    for t2 in ax2.get_yticklabels():
        t2.set_color('red')
    # Choose whether or not to save the file.
    if save == False:
        plt.show()
    else:
        plt.savefig(fname)
        plt.close()
    return None
    
def gen_frames(data, t_o, t_f, dt = 1):
    for t in np.arange(t_o, t_f, dt):
        hist_plot(data, t, save = True, fname = 'frames\%05d.png' % t)
    return None
    
data = np.zeros((t_final,n_trials))
for i in range(n_trials):
    # Set start position of bacteria and reset start time
    data[0,i] = np.random.uniform(a,b)
    t = 1
    while(t < t_final):
        data[t,i] = data[t-1,i] + np.random.normal(0.0, sigma, 1) +\
        calc_v(data[t-1,i], t-1, x_table, c_attract, dc_attract, c_repel, dc_repel)*dt
        t += 1
        
k_plot(data, True, 'Kymograph.png')
pop_plot(data, L + a + 10., True, 'pop_plot.png')
gen_frames(data, 0, 3600, 10)