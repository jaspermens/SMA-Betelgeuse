#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 14:48:44 2022

@author: konstantinos

Test run, just building a test oort cloud around a fucker beet
"""
# Python Modules, gotta have 'em
import numpy as np
import matplotlib.pyplot as plt
# General AMUSE Stuff
from amuse.lab import units
from amuse.ic import isotropic_cloud as ic # Has generate oort cloud
from amuse.lab import Particles
from amuse.units import nbody_system # Nbody units
# Community Modules
from amuse.community.ph4.interface import ph4
#%% Make beet and his cloud

beet = Particles(1)
beet.mass = 20 | units.MSun

# Make Oort Cloud
# Linear extrapolation from the mass of beet gets us 
# a_min = 60_000
# a_max = 200_000
# NOTE: For some reason its faster if α is bigger.

beet_cloud = ic.new_isotropic_cloud(number_of_particles=10_000,
                       m_star=beet.mass,
                       a_min = 60_000 | units.AU,
                       a_max = 200_000 | units.AU)

# Add the particles to beet
beet.add_particles(beet_cloud)

# Renaming
beet_system = beet

#%% Initialize gravity solver


# Caclulate max distance for converter
# NOTE: Definetly can be coded neater through AMUSE functions
# This is okay-ish quick. The distances_squared method runs loong.

max_dist = 0 | units.AU
for obj in beet_system:
    dist = np.sqrt(obj.x**2 + obj.y**2 + obj.z**2)
    if dist.value_in(units.AU) > max_dist.value_in(units.AU):
        print('Maximum Distance Updated')
        max_dist = dist
    
# Init with beet's mass and the distance of the furthest Oort Obj.
converter = nbody_system.nbody_to_si(beet_system[0].mass, max_dist)
gravity = ph4(converter)
gravity.particles.add_particles(beet_system)

# Gets information FROM the solver
channel_out = gravity.particles.new_channel_to(beet_system)

# Puts information TO the solver | Needed cause we need to reduce beet's mass
channel_in = beet_system.new_channel_to(gravity.particles)

#%% Evolve a bit
# Just to see what happens

end_time = 1_000 | units.yr
model_time = 0 | units.yr
j = 0 # counter for figure names

while(model_time < end_time):
    # Arbitrary time step
    dt = 1 | units.yr
    model_time += dt
    
    # Do the thing.
    gravity.evolve_model(model_time)
    channel_out.copy()
    
    # Check what happens, every 10 years
    # NOTE: Make this into trails because this looks like shit currently.
    if model_time.value_in(units.yr) % 50 == 0:
        time_percentage = model_time/end_time 
        print(str( int(100*time_percentage)) + ' percent done')
        plt.figure()
        plt.scatter(beet_system.x.value_in(units.AU), 
                    beet_system.y.value_in(units.AU),
                    s=1)
        plt.xlim(-max_dist.value_in(units.AU),max_dist.value_in(units.AU))
        plt.ylim(-max_dist.value_in(units.AU),max_dist.value_in(units.AU))
        plt.savefig('beet_test'+str(j)+'.png')
        j=j+1

#%% Make into a .gif
from PIL import Image

i=0
image_paths=[None]*j
im=[]
for i in range(len(image_paths)):
    image_paths[i]='beet_test'+str(i)+'.png' #load up the image names
    new_im= Image.open(image_paths[i]) #open one
    im.append(new_im)

im[0].save('beet_test.gif', format='GIF',
            append_images=im[1:],
            save_all=True,
            duration=500, #in ms
            loop=1) #for zero repeats forever

