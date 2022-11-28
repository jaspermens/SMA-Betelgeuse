
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:32:25 2022
@author: Konstantinos, Rahul, Jasper
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
import pickle as pk

#%%
class star:
    def __init__(self, 
                particles, 
                timestep, 
                mdot=(0.1|units.MSun / units.yr)
                ):
        self.particles = particles
        self.model_time = 0 | units.yr
        self.timestep = timestep
        self.mdot = mdot
        
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

    def update_mass(self):
        # self.M = self.M - 0.5 | units.MSun
        
        if self.particles.mass > (0 | units.MSun):
            self.particles.mass -= self.mdot * self.timestep

    def update_pos(self):
        self.particles.x += self.particles.vx * self.timestep
        self.particles.y += self.particles.vy * self.timestep
        self.particles.z += self.particles.vz * self.timestep

    def evolve_model(self, time):
        while self.model_time < time:
            self.update_pos()
            self.update_mass()
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
        
    def get_gravity_at_point(self, eps, x, y, z):
        return (0, 0, 0) | (units.m * units.s**(-2))
        
    def evolve_model(self, time):
        while self.model_time < time:
            self.update_pos()
            self.model_time += self.timestep

#%% progress bar function
def bar_x_out_of_y(x, y, text: str) -> None:
    maxbars = 20
    nbars = int(x / y * maxbars)
    print('\rProcessing: | ' + '█' * nbars + '-' * (maxbars - nbars) + ' |', end=' ')
    print(f'~ {x}/{y} {text}', end='')
    if x == y:
        print('')

#%% Testing sim

def detect_encounters(cloud, sun):
    detections = []
    for particle in cloud:
        # AMUSE bullshit dictates that cloud.position[0] - sun.position[0] is wrong
        relx = particle.position.x - sun.position.x
        rely = particle.position.y - sun.position.y 
        relz = particle.position.z - sun.position.z
        # relx, rely, relz = sun.position[0] - particle.position[0] 
        if (relx*relx + rely*rely + relz*relz).sqrt() < (20 | units.pc):
            detections.append(particle)
            
    return detections
            
    
def run_simulation(end_time=(100|units.Myr), 
                    timestep_oort=(1e-3|units.Myr), 
                    timestep_mw=(1e-2|units.Myr),
                    n_oort_objects=100,
                    do_plot: bool=True):

    """ Performs a simple test run"""
    
    # beet
    beet = Particles(1)
    beet.mass = 20 | units.MSun
    beet.position = (26660, 0, 0) | units.lightyear
    beet.velocity = (0, 250, 0) | units.kms
    
    # sun
    sun = Particles(1)
    sun.mass = 1 | units.MSun
    sun.position = (26660*1.02, 0, 0) | units.lightyear
    sun.velocity = (0, 250*.98, 0) | units.kms

    # Oort cloud
    beet_cloud = ic.new_isotropic_cloud(number_of_particles=n_oort_objects,
                            m_star = 18 | units.MSun,
                            a_min = 60_000 | units.AU,
                            a_max = 200_000 | units.AU,
                            seed=42)

    beet_cloud.position += beet.position
    beet_cloud.velocity += beet.velocity

    # Make test particles and beet and sun
    OORT = test_particles(beet_cloud, timestep_oort)
    BEET = star(beet, timestep_oort)
    SUN = star(sun, timestep_mw)    

    # Bridge 'em
    # Bridge Oort and Beet
    oort_cloud = bridge.Bridge(use_threading=True)
    oort_cloud.add_system(OORT, (BEET,))
    oort_cloud.add_system(BEET, (OORT,))
    oort_cloud.timestep = timestep_oort

    # Bridge Milky Way and Beet System
    MWG = MilkyWay_galaxy()
    milky_way = bridge.Bridge(use_threading=True)
    milky_way.add_system(oort_cloud, (MWG,) )
    milky_way.add_system(SUN, (MWG,))
    milky_way.timestep = timestep_mw

    # Make channel
    channel_out = milky_way.particles.new_channel_to(beet_cloud)
    channel_out_sun = milky_way.particles.new_channel_to(sun)
    
    # Do the thing
    model_time = 0 | units.yr
    plot_times = np.linspace(0,1,num=101) * end_time
    plot_interval = 1 | units.Myr
    plotdata = []

    start_time = datetime.now()
    while(model_time < end_time):
        # Evolve in milky way timesteps because those should be the largest
        channel_out.copy()
        model_time += timestep_mw
        
        # Do the thing.
        milky_way.evolve_model(model_time)

        # collision detection
        detections = detect_encounters(beet_cloud, sun)
        if len(detections) > 0:
            print(len(detections), 'detection(s)!')
        # Progress check
        bar_x_out_of_y(model_time, end_time, '')

        #  Saving data for future plotting
        if do_plot:
            if (model_time % plot_interval).value_in(units.Myr) == 0.0:

                # print(str( int(100*time_percentage)) + ' percent done')
                # print(model_time)
                plotdata.append((str(model_time.value_in(units.Myr)) + ' Myr', 
                                (BEET.particles.x.value_in(units.kpc), 
                                    BEET.particles.y.value_in(units.kpc)), 
                                (beet_cloud.x.value_in(units.kpc), 
                                    beet_cloud.y.value_in(units.kpc)),
                                (sun.x.value_in(units.kpc),
                                     sun.y.value_in(units.kpc))
                ))
            # if int(model_time.value_in(units.yr)) in plot_times.value_in(units.yr):
            #     channel_out.copy()
            #     # print(str( int(100*time_percentage)) + ' percent done')
            #     # print(model_time)
            #     plotdata.append((str(model_time.value_in(units.Myr)) + ' Myr', 
            #                     (BEET.particles.x.value_in(units.kpc), 
            #                         BEET.particles.y.value_in(units.kpc)), 
            #                     (beet_cloud.x.value_in(units.kpc), 
            #                         beet_cloud.y.value_in(units.kpc))
            #     ))
    
    end_time = datetime.now() - start_time
    print(f'Duration: {end_time}')
            
    if do_plot:
        pk.dump(plotdata, open('last_run_plotdata.pk', 'wb'))
        make_plots(plotdata)      


def make_plots(plotdata=None):
    import matplotlib.patches as patches
    if not plotdata:
        plotdata = pk.load(open('last_run_plotdata.pk', 'rb'))
    print(len(plotdata))
    fignum = 0
    for plot_data in plotdata:
        tit, beetpos, cloudpos, sunpos = plot_data
        
        colors = ['purple' , 'goldenrod', 'forestgreen', 'steelblue', 'teal']
        colors = colors * (len(cloudpos[0]) // len(colors))

        fig, ax = plt.subplots(1, figsize=(4,4), dpi=200)
        ax.set_title(tit)
        ax.scatter(0, 0, c='maroon', zorder=2)
        ax.scatter(beetpos[0],
                    beetpos[1],
                    c = 'red', zorder=2)
        ax.scatter(sunpos[0],
                    sunpos[1],
                    c = 'gold', zorder=2)
        circ = patches.Circle((sunpos[0], sunpos[1]), radius=0.02, transform=ax.transData, fill=False, color='purple', linestyle='--')
        ax.add_patch(circ)
        ax.scatter(cloudpos[0], 
                    cloudpos[1],
                    c = colors,
                    s = 1
                    )
        ax.set_xlim(sunpos[0] - 0.1, sunpos[0] + 0.1)
        ax.set_ylim(sunpos[1] - 0.1, sunpos[1] + 0.1)
        # plt.grid()
        plt.savefig(f'../figures/fig_{fignum:03d}.png')
        plt.close()
        fignum +=1


# %%
if __name__ in '__main__':
    run_simulation(end_time=(100|units.Myr), timestep_oort=(0.01|units.Myr), timestep_mw=(0.01|units.Myr), n_oort_objects=100)
    # make_plots()