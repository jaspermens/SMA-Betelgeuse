#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:26:52 2022

@author: konstantinos

postproccssing
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from amuse.lab import units
import random

#%% Create dummy data
# Dataframe [key, time, x, y, z, vx, vy, vz]
#           [ str, rest are amuse obj]
keys = []
letter_bank = ['a','n','t','h','e','e']
times = []
xs = []
ys = []
zs = []
vxs = []
vys = []
vzs = []
# Init while loop
dummy_len = 10_000
i = 0
while (i<dummy_len):
    # string
    # I could enforce these being uniques but I am not gonna
    keys.append(''.join(random.choices(letter_bank, k=6)))
    # 80-120 [Myr]
    times.append( np.random.randint(80, 120))
    # 0-1 [pc]
    xs.append(np.random.randint(2,5))
    ys.append(np.random.rand())
    zs.append(np.random.rand())
    # 12-17 [km/s]
    vxs.append(np.random.randint(12,17))
    vys.append(np.random.randint(12,17))
    vzs.append(np.random.randint(12,17))
    i += 1
# Make the DF
names = ['key', 'time', 'x', 'y', 'z', 'vx' , 'vy',
         'vz']
data = [ keys, times, xs, ys, zs, vxs, vys, vzs]
df = pd.DataFrame(data, index=names).T
#%% Stats
# kde which works
xy = np.vstack( [df['x'].to_numpy(dtype=np.float32),
                 df['y'].to_numpy(dtype=np.float32)])
# instantiate and fit the KDE model
kde = stats.gaussian_kde(xy)
xrange = np.linspace(0,10, num=10_000)
yrange = np.linspace(0,1,num=10_000)
ranges = np.vstack([xrange, yrange])
results = np.reshape(kde(ranges) , (100,100))
plt.imshow(results)




