#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 22:01:12 2022

@author: konstantinos
"""
import numpy as np
import pickle as pk
import matplotlib.pyplot as plt
phi_real = -30 * np.pi/180
width_real = 4.5 * np.pi/180
phi_range = np.linspace(phi_real-width_real , phi_real+width_real, num = 6)
fig, ax = plt.subplots(1, figsize=(6,6), dpi=300)
for phi in phi_range:
    name = str(np.round(phi,3))
    plotdata = pk.load(open(name+'_run_plotdata.pk', 'rb'))
    beet_z  = np.zeros(len(plotdata))
    sun_z  = np.zeros(len(plotdata))
    beet_x  = np.zeros(len(plotdata))
    sun_x  = np.zeros(len(plotdata))
    beet_y  = np.zeros(len(plotdata))
    sun_y  = np.zeros(len(plotdata))
    dist = np.zeros(len(plotdata))
    
    for i, plot_data in enumerate(plotdata):
        tit, beetpos, cloudpos, sunpos = plot_data
        beet_z[i] = beetpos[2]
        sun_z[i] = sunpos[2]
        beet_x[i] = beetpos[0]
        sun_x[i] = sunpos[0]
        beet_y[i] = beetpos[1]
        sun_y[i] = sunpos[1]
        dist[i] = np.linalg.norm( [sunpos[0]-beetpos[0],
                                   sunpos[1] - beetpos[1],
                                   sunpos[2] - beetpos[2]])
        
        
    
    print(min(dist))
    ax.plot(dist, label=name)
plt.legend()
#%%
fig, ax = plt.subplots(1, figsize=(6,6), dpi=300)
phi = -20 * np.pi/180
name = str(np.round(phi,3))
plotdata = pk.load(open(name+'_run_plotdata.pk', 'rb'))
beet_z  = np.zeros(len(plotdata))
sun_z  = np.zeros(len(plotdata))
beet_x  = np.zeros(len(plotdata))
sun_x  = np.zeros(len(plotdata))
beet_y  = np.zeros(len(plotdata))
sun_y  = np.zeros(len(plotdata))
dist = np.zeros(len(plotdata))

for i, plot_data in enumerate(plotdata):
    tit, beetpos, cloudpos, sunpos = plot_data
    beet_z[i] = beetpos[2]
    sun_z[i] = sunpos[2]
    beet_x[i] = beetpos[0]
    sun_x[i] = sunpos[0]
    beet_y[i] = beetpos[1]
    sun_y[i] = sunpos[1]
    dist[i] = np.linalg.norm( [sunpos[0]-beetpos[0],
                               sunpos[1] - beetpos[1],
                               sunpos[2] - beetpos[2]])
    
    

print(min(dist))
ax.plot(dist, label=name)
plt.legend()