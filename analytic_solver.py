#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:32:25 2022

@author: konstantinos
"""
# Python Modules, gotta have 'em
import numpy as np
import matplotlib.pyplot as plt
# General AMUSE Stuff
from amuse.units import units, constants
from amuse.lab import Particles
from amuse.couple import bridge
from amuse.ic import isotropic_cloud as ic # Has generate oort cloud
from datetime import datetime # for timing
from galactic_potential import MilkyWay_galaxy
#%%
class star:
    def __init__(self, timestep, particles):
        self.particles = particles
        self.model_time = 0 | units.yr
        self.timestep = timestep
        
    def get_gravity_at_point(self, eps, x, y, z): 
        x -= self.particles.x
        y -= self.particles.y
        z -= self.particles.z
        
        r = (x**2+y**2+z**2).sqrt()
        a = - constants.G * self.particles.mass/r**2
          
        ax = a*x/r
        ay = a*y/r
        az = a*z/r
        
        return ax, ay, az
    
    def evolve_model(self, time):
        # self.M = self.M - 0.5 | units.MSun
        if self.particles.mass > 0.5 | units.MSun:
            self.particles.mass -= 0.5|units.MSun
        self.model_time += self.timestep
        
class test_particles:
    def __init__(self,
                 particles,
                 timestep):
        self.particles = particles
        self.timestep = timestep
        self.model_time = 0 | units.yr
        
    def update_pos(self):
        
        self.particles.x += self.particles.vx * self.timestep
        self.particles.y += self.particles.vy * self.timestep
        self.particles.z += self.particles.vz * self.timestep
        self.model_time += self.timestep
        
    def get_gravity_at_point(self, eps, x, y, z):
        return (0,0,0 )| (units.m * units.s**(-2))
        
    def evolve_model(self, time):
        self.update_pos()
        
#%% Testing sim

beet_cloud = ic.new_isotropic_cloud(number_of_particles=100,
                        m_star = 20 | units.MSun,
                        a_min = 60_000 | units.AU,
                        a_max = 200_000 | units.AU)

beet = Particles(1)
beet.mass = 20 | units.MSun
beet.position = (26660, 0, 0) | units.lightyear
beet.velocity = (0, 250, 0) | units.kms

beet_cloud.position += beet.position
beet_cloud.velocity += beet.velocity

colors = ['purple' , 'goldenrod', 'forestgreen', 'steelblue', 'teal']
colors = colors * (len(beet_cloud) // len(colors))
# EXPLICIT TIMESTEP
timestep = 1_000 | units.yr

# Make test particles and sun
TP = test_particles(beet_cloud, timestep)
BEET = star(timestep, beet)

# Bridge 'em
# Bridge Oort and Beet
gravity = bridge.Bridge(use_threading=True)
gravity.add_system(TP, (BEET,))
gravity.add_system(BEET, (TP,))
# Bridge Milky Way and Beet System
MWG = MilkyWay_galaxy()
gravity2 = bridge.Bridge(use_threading=True)
gravity2.add_system(gravity, (MWG,) )
gravity2.timestep = 0.01 | units.Myr

channel_out = gravity2.particles.new_channel_to(beet_cloud)

# Do the thing
end_time = 100 | units.Myr
model_time = 0 | units.yr
plot_times = np.linspace(0,1,num=101) * end_time
start_time = datetime.now()
while(model_time < end_time):
    # Arbitrary time step
    model_time += timestep
    # Do the thing.
    gravity2.evolve_model(model_time)
    channel_out.copy()

    #  Plotting
    if int(model_time.value_in(units.yr)) in plot_times.value_in(units.yr):
        # Progress check
        time_percentage = model_time/end_time
        # print(str( int(100*time_percentage)) + ' percent done')
        # print(model_time)
        plt.figure()
        plt.title(str(model_time.value_in(units.Myr)) + ' Myr')
        plt.scatter(0,0,c='maroon', zorder=2)
        plt.scatter(BEET.particles.x.value_in(units.kpc),
                    BEET.particles.y.value_in(units.kpc),
                    c = 'red', zorder=2)
        plt.scatter(beet_cloud.x.value_in(units.kpc), 
                    beet_cloud.y.value_in(units.kpc),
                    c = colors,
                    s = 7
                    )
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        plt.grid()
        plt.show()
        
end_time = datetime.now() - start_time
print('Duration: {}'.format(end_time))