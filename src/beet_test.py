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
from amuse.community.hermite.interface import Hermite
from galactic_potential import MilkyWay_galaxy
from amuse.couple import bridge
from datetime import datetime# for timing

#%% Make beet and his cloud

beet_system = Particles(2)
beet = beet_system[0]
beet.mass = 20 | units.MSun
sun = beet_system[1]
sun.mass = 1 | units.MSun

print(beet_system[0].mass)
print(beet.mass)

# Make Oort Cloud
# Linear extrapolation from the mass of beet gets us 
# a_min = 60_000
# a_max = 200_000
# NOTE: For some reason its faster if α is bigger.

beet_cloud = ic.new_isotropic_cloud(number_of_particles=100,
                       m_star=beet.mass,
                       a_min = 60_000 | units.AU,
                       a_max = 500_000 | units.AU)

beet.position = (1.0, 0, 0) * ( 26660 | units.lightyear)
beet.velocity = (0.0, 1.0, 0) * (250 | units.kms)

# distance between sun and beet is approx 2% of sun orbit radius
# small angle approximation or some shit: 
sun.position = (1.02, 0, 0) * (26660 | units.lightyear)
sun.velocity = (0, 0.98, 0) * (250 | units.kms)

beet_cloud.position += beet.position
beet_cloud.velocity += beet.velocity

# Add the particles to beet
beet_system.add_particles(beet_cloud)
# Everybody knows beet is red
beet_system[0].color = 'maroon'
# Everybody also knows sun is yellow
beet_system[1].color = 'yellow'
# and oort objects are blue
beet_system[2:].color = 'steelblue'



#%% Initialize gravity solver

# Init with beet's mass and the distance of the furthest Oort Obj.
converter = nbody_system.nbody_to_si(beet_system[0].mass, 0.01 |units.pc)
gravity_code = Hermite(converter)
gravity_code.particles.add_particles(beet_system)

# Gets information FROM the solver
channel_out = gravity_code.particles.new_channel_to(beet_system)

# Puts information TO the solver | Needed cause we need to reduce beet's mass
channel_in = beet_system.new_channel_to(gravity_code.particles)

#instance of the Milky Way galaxy
MWG = MilkyWay_galaxy()

# bridging the potential and gravity code
bridge_timestep = 1 | units.Myr
gravity = bridge.Bridge(use_threading=True)
gravity.add_system(gravity_code, (MWG,))
gravity.timestep = bridge_timestep 

#%% Evolve a bit
# Just to see what happens
start_time = datetime.now()
end_time = 1_000 | units.Myr
model_time = 0 | units.yr
j = 0 # counter for figure names
explode_time = 20 | units.Myr
while(model_time < end_time):
    # Arbitrary time step
    dt = .1 | units.Myr
    model_time += dt
    
    # Do the thing.
    gravity.evolve_model(model_time)
    
    # Implementing the mass loss near the end 
    if model_time == explode_time:
        beet_system[0].mass = 0 | units.MSun
        channel_in.copy()
        
    channel_out.copy()
    
    # Check what happens, every 10 years
    # NOTE: Make this into trails because this looks like shit currently.
    if model_time.value_in(units.Myr) % 100 == 0:
        time_percentage = model_time/end_time 
        print(str( int(100*time_percentage)) + ' percent done')
        plt.figure()
        plt.title(str(model_time.value_in(units.Myr)) + ' Myr')
        plt.scatter(0,0,c='k')
        plt.scatter(beet_system.x.value_in(units.AU), 
                    beet_system.y.value_in(units.AU),
                    c = beet_system.color,
                    s=1)
        plt.scatter(beet_system[1].x.value_in(units.AU), 
                    beet_system[1].y.value_in(units.AU),
                    c = 'xkcd:piss yellow',
                    s=20,
                    zorder = 10)
        if model_time == explode_time:
            plt.text(1e9,0,s='SUPERNOVA IS NOW')
        #plt.xlim((beet_system[0].position.x - (10 | units.pc)).value_in(units.AU),(beet_system[0].position.x + (10 |units.pc)).value_in(units.AU))
        #plt.ylim((beet_system[0].position.y - (10 | units.pc)).value_in(units.AU),(beet_system[0].position.y + (10 |units.pc)).value_in(units.AU))
        plt.xlim(-3e9, 3e9)
        plt.ylim(-3e9, 3e9)
        plt.savefig('beet_test'+str(j)+'.png')
        j=j+1
        plt.close()
        
end_time = datetime.now()

print('Duration: {}'.format(end_time - start_time))
#%% Make into a .gif
# from PIL import Image

# i=0
# image_paths=[None]*j
# im=[]
# for i in range(len(image_paths)):
#     image_paths[i]='beet_test'+str(i)+'.png' #load up the image names
#     new_im= Image.open(image_paths[i]) #open one
#     im.append(new_im)

# im[0].save('beet_test.gif', format='GIF',
#             append_images=im[1:],
#             save_all=True,
#             duration=50, #in ms
#             loop=1) #for zero repeats forever

