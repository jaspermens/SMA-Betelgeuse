#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 15:52:27 2022

@author: konstantinos
"""
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["figure.figsize"] = (6,3)
plt.rcParams['figure.dpi'] = 300

timesteps = [1, 1.5,0.75,0.5,0.1, 2.5, 0.05]
detections_1 = [26,28,25,25,29, 33, 30]
plt.scatter(timesteps, detections_1,  color='purple', s=100,
            label='Detection: 0.1 Myr')
detections_05 = [27,28,25,25,29, 33, 30]

plt.scatter(timesteps, detections_05, s=50, marker='h', color='goldenrod', 
            label='Detection: 0.05 Myr')

plt.grid()
plt.legend()
plt.xlabel('Timestep Size [Myr]')
plt.ylabel('No. of Detections')
plt.title('Convergence Check for MWG timestep, n = 10_000')
#%%

timesteps = [0.1, 0.2, 0.5, 0.75, 1,1.2, 1.5, 2, 2.5]
detections_1 = [278,279, 286, 255 ,279, 288, 255, 250, 232]
plt.scatter(timesteps, detections_1,  color='purple', s=100,
            label='Detection: 0.1 Myr')

plt.grid()
#plt.legend()
plt.xlabel('Timestep Size [Myr]')
plt.ylabel('No. of Detections')
plt.title('Convergence Check for MWG timestep, n = 100_000')

