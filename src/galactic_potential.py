#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 11:52:21 2022

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt 
from amuse.units import units, constants
from amuse.lab import Particles
from amuse.units import nbody_system # Nbody units
from amuse.couple import bridge
# Community Modules
from amuse.community.hermite.interface import Hermite
from amuse.community.ph4.interface import ph4

#%% Milky Way galaxy

class MilkyWay_galaxy(object):
    def __init__(self, 
                 Mb = 1.40592e10 | units.MSun, # Bulge Mass
                 Md = 8.5608e10 | units.MSun, # Disk Mass
                 Mh = 1.07068e11 | units.MSun): # Halo Mass
        self.Mb = Mb
        self.Md = Md
        self.Mh = Mh

    def get_potential_at_point(self, eps, x, y, z):
        # Allen, Santillan '91 
        r = (x**2+y**2+z**2)**0.5 # radial coord. of point
        R = (x**2+y**2)**0.5 # 2-D radial coord. in disk
        
        # Bulge Eq. (1)
        b1 = 0.3873 | units.kpc # Scale length bulge
        pot_bulge = -constants.G * self.Mb / (r**2+b1**2)**0.5 
        
        # Disk Eq. (3)
        a2 = 5.31 | units.kpc # Scale length disk 1
        b2 = 0.25 | units.kpc # Scale length disk 2
        ari = -constants.G * self.Md # Numerator 
        par = (R**2 + (a2+ (z**2+ b2**2)**0.5 )**2 )**0.5 # Denominator
        pot_disk = ari / par
        
        # Halo Eq. (5)
        a3 = 12.0 | units.kpc # Scale length halo
        cut_off = 100 | units.kpc
        d1 = r/a3  
        c = 1 + (cut_off/a3)**1.02
        pot_halo = -constants.G*(self.Mh/a3)*d1**1.02/(1+ d1**1.02) \
                  - (constants.G*self.Mh/(1.02*a3))\
                      * (-1.02/c +np.log(c) + 1.02/(1+d1**1.02) \
                           - np.log(1.0 +d1**1.02) )
                          
        # Multiply by 2 because it is a rigid potential
        potential = 2*(pot_bulge + pot_disk + pot_halo) 
        return potential 
                
    
    def get_gravity_at_point(self, eps, x, y, z): 
        r =  (x**2+y**2+z**2)**0.5
        R =  (x**2+y**2)**0.5
        #print('AAA MASTER I HAVE BEN SUMMONED')
        #bulge
        b1 = 0.3873 |units.kpc
        force_bulge = -constants.G*self.Mb/(r**2+b1**2)**1.5 
        
        #disk
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        d = a2 + (z**2 + b2**2)**0.5
        force_disk = -constants.G*self.Md/(R**2 + d**2)**1.5
        
        #halo
        a3 = 12.0 |units.kpc
        d1 = r/a3
        force_halo = -constants.G*self.Mh*d1**0.02/(a3**2*(1+d1**1.02))
       
        ax = force_bulge*x + force_disk*x  + force_halo*x/r
        ay = force_bulge*y + force_disk*y  + force_halo*y/r
        az = force_bulge*z + force_disk*d*z/(z**2 + b2**2)**0.5 + force_halo*z/r 

        return ax, ay, az
    
    def vel_circ(self, x , y ):
        z = 0 | units.kpc 
        r = (x**2+y**2+z**2)**0.5 # radial coord. of point
        b1 = 0.3873 |units.kpc
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        a3 = 12.0 |units.kpc
    
        rdphi_b = constants.G*self.Mb*r**2/(r**2+b1**2)**1.5
        rdphi_d = constants.G*self.Md*r**2/(r**2+(a2+(z**2+b2**2)**0.5)**2 )**1.5
        rdphi_h = constants.G*self.Mh*(r/a3)**0.02*r/(a3**2*(1+(r/a3)**1.02))
    
        vel_circb = rdphi_b
        vel_circd = rdphi_d
        vel_circh = rdphi_h

        return (vel_circb + vel_circd + vel_circh)**0.5 
    
#%% small simulation with beet and the sun

#%%

def test_run(end_time, step, bridge_timestep = 1 | units.Myr):
    # Builds the system
    beet_and_sun = Particles(2)
    beet = beet_and_sun[0]
    beet.mass = 20 | units.MSun
    beet.position = (1.0, 0, 0) * ( 8 | units.kpc)
    beet.velocity = (1.0, 0, 0) * (250 | units.kms)
    sun = beet_and_sun[1]
    sun.mass = 1 | units.MSun
    sun.position = (1.0,0,0) * (7 | units.kpc)
    sun.velocity = (0,1.0,0) * (200 | units.kms)
    
    # Converter
    converter = nbody_system.nbody_to_si(21 | units.MSun, 
                                         1 | units.kpc)
    
    # Initialising the gravity solver
    gravity_code = Hermite(converter)
    gravity_code.particles.add_particles(beet_and_sun)
    channel_out = gravity_code.particles.new_channel_to(beet_and_sun)
    
    #instance of the Milky Way galaxy
    MWG = MilkyWay_galaxy()
    
    # bridging the potential and gravity code
    gravity = bridge.Bridge(use_threading= False)
    gravity.add_system(gravity_code, (MWG,))
    gravity.timestep = bridge_timestep 
    
    model_time = 0 | units.yr
    
    while(model_time < end_time):
        # Arbitrary time step
        model_time += step
        
        # Do the thing.
        gravity.evolve_model(model_time)
        channel_out.copy()
        
        if model_time.value_in(units.yr) % 1000 == 0:
            time_percentage = model_time/end_time 
            print(str( int(100*time_percentage)) + ' percent done')
            plt.figure()
            # Galactic centre as an anchor
            plt.scatter(0,0, c='k')
            # Sun is yellow
            plt.scatter(beet_and_sun[1].x.value_in(units.kpc), 
                        beet_and_sun[1].y.value_in(units.kpc), 
                        c = 'goldenrod')
            # Beet is red
            plt.scatter(beet_and_sun[0].x.value_in(units.kpc), 
                        beet_and_sun[0].y.value_in(units.kpc), 
                        c = 'red')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)

if __name__ in '__main__':
    test_run(end_time = 10 | units.Myr, step = 1 | units.Myr)



