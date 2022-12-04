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
import os
import pandas as pd
from amuse.community.seba.interface import SeBa

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

        self.stellar = SeBa()
        self.stellar.particles.add_particles(self.particles)
        self.channel_mass = self.stellar.particles.new_channel_to(self.particles)
        
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

    def update_mass_linear(self):
        # self.M = self.M - 0.5 | units.MSun
        
        if self.particles.mass > (0 | units.MSun):
            self.particles.mass -= self.mdot * self.timestep

    def update_mass(self):
        self.stellar.evolve_model(self.model_time + self.timestep)
        self.channel_mass.copy()

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
def bar_x_out_of_y(x, y, text: str='') -> None:
    maxbars = 20
    nbars = int(x / y * maxbars)
    print('\rProcessing: | ' + '█' * nbars + '-' * (maxbars - nbars) + ' |', end=' ')
    print(f'~ {x}/{y} {text}', end='')
    if x == y:
        print('')

#%% Testing sim

def detect_encounters_slow(cloud, sun, model_time, detections, detection_keys, detection_radius):
    for particle in cloud:
        if particle.key in detection_keys:
            continue

        relx, rely, relz = particle.position - sun[0].position

        if (relx*relx + rely*rely + relz*relz).sqrt() < detection_radius:            
            relvx, relvy, relvz = particle.velocity - sun[0].velocity
            print('Detection!')            
            detections.append([particle.key, model_time, relx, rely, relz, relvx, relvy, relvz])
            detection_keys.append(particle.key)
            
    return detections


def detect_pseudo_encounters_slow(cloud, sun, model_time, detections, detection_keys, detection_radius):
    for particle in cloud:
        if particle.key in detection_keys:
            continue

        relx, rely, relz = particle.position - sun[0].position
        if (relx*relx + rely*rely).sqrt() < detection_radius:
            relvx, relvy, relvz = particle.velocity - sun[0].velocity
            print('Projected detection!')            
            detections.append([particle.key, model_time, relx, rely, relz, relvx, relvy, relvz])
            detection_keys.append(particle.key)
            
    return detections


def detect_encounters(cloud, sun, model_time, detections, detection_keys, detection_radius):
    relpos = cloud.position - sun.position
    distances = relpos.lengths()

    if distances.min() < detection_radius:
        print('detection!!')
        cpi = distances.argmin()  # closest particle index
        closest = cloud[cpi]
        relvx, relvy, relvz = closest.velocity
        relx, rely, relz = closest.position
        detections.append([closest.key, model_time, relx, rely, relz, relvx, relvy, relvz])
        detection_keys.append(closest.key)
        cloud.remove_particle(cloud[cpi])
        detect_encounters(cloud, sun, model_time, detections, detection_keys, detection_radius)


def run_simulation(end_time=(100|units.Myr), 
                    timestep=(1e-3|units.Myr),
                    n_oort_objects=100,
                    detection_radius=(20|units.pc)):

    """ Performs a simple test run"""
    
    run_params = {'t_end': end_time, 
                    'timestep': timestep,
                    'n_objects': n_oort_objects,
                    'det_rad': detection_radius}

    pk.dump(run_params, open('run_params.pk', 'wb'))

    # beet
    beet = Particles(1)
    beet.mass = 18 | units.MSun
    beet.position = (26660, 0, 0) | units.lightyear
    beet.velocity = (0, 250, 0) | units.kms
    
    # sun
    sun = Particles(1)
    sun.mass = 1 | units.MSun
    sun.position = (26660*1.02, 0, 0) | units.lightyear
    sun.velocity = (0, 250*.98, 0) | units.kms

    # Oort cloud
    beet_cloud = ic.new_isotropic_cloud(number_of_particles=n_oort_objects,
                            m_star = beet.mass,
                            a_min = 60_000 | units.AU,
                            a_max = 300_000 | units.AU,
                            q_min = 20_000 | units.AU,
                            seed=42069)

    beet_cloud.position += beet.position
    beet_cloud.velocity += beet.velocity

    # Make test particles and beet and sun
    OORT = test_particles(beet_cloud, timestep)
    BEET = star(beet, timestep)
    SUN = test_particles(sun, timestep)    

    # Bridge 'em
    # Bridge Oort and Beet
    MWG = MilkyWay_galaxy()
    milky_way = bridge.Bridge(use_threading=True)
    milky_way.add_system(OORT, (BEET, MWG))
    milky_way.add_system(BEET, (MWG,))
    milky_way.add_system(SUN, (MWG,))
    milky_way.timestep = timestep

    # Bridge Milky Way and Beet System
    # MWG = MilkyWay_galaxy()``
    # milky_way = bridge.Bridge(use_threading=True)
    # milky_way.add_system(oort_cloud, (MWG,) )
    # milky_way.add_system(SUN, (MWG,))
    # milky_way.timestep = timestep_mw

    # Make channel
    channel_out = milky_way.particles.new_channel_to(beet_cloud)
    channel_in = beet_cloud.new_channel_to(milky_way.particles)
    channel_out_sun = milky_way.particles.new_channel_to(sun)
    channel_out_beet = milky_way.particles.new_channel_to(beet)
    
    # Do the thing
    model_time = 0 | units.yr
    plot_interval = 1 | units.Myr
    plotdata = []
    detections = []
    detection_keys = []

    start_time = datetime.now()
    while(model_time < end_time):
        # Evolve in milky way timesteps because those should be the largest
        model_time += timestep
        
        # Do the thing.
        milky_way.evolve_model(model_time)

        # collision detection
        channel_out.copy()
        channel_out_sun.copy()          # I don't know if these two do anything
        channel_out_beet.copy()         # It runs just as well without them two channels
        detect_encounters(beet_cloud, sun, model_time, detections, detection_keys, detection_radius)
        channel_in.copy()

        #  Saving data for future plotting
        if (model_time % plot_interval).value_in(units.Myr) == 0.0:
            # Progress check
            bar_x_out_of_y(model_time.in_(end_time.unit), end_time, '')

            # print(str( int(100*time_percentage)) + ' percent done')
            # print(model_time)
            plotdata.append((str(model_time.value_in(units.Myr)) + ' Myr', 
                            (BEET.particles.x.value_in(units.kpc), 
                                BEET.particles.y.value_in(units.kpc),
                                BEET.particles.z.value_in(units.kpc)), 
                            (beet_cloud.x.value_in(units.kpc), 
                                beet_cloud.y.value_in(units.kpc),
                                beet_cloud.z.value_in(units.kpc)),
                            (sun.x.value_in(units.kpc),
                                    sun.y.value_in(units.kpc),
                                    sun.z.value_in(units.kpc))
            ))
    
    end_time = datetime.now() - start_time
    print(f'Duration: {end_time}')

    pk.dump(plotdata, open('last_run_plotdata.pk', 'wb'))

    if len(detections) > 0:
        detection_df = pd.DataFrame(detections)
        detection_df.columns = "key", "time", "x", "y", "z", "vx", "vy", "vz"
        print(detection_df)
        pk.dump(detection_df, open('detections.pk', 'wb'))


def make_plots(plotdata=None):
    import matplotlib.patches as patches
    import mpl_toolkits.mplot3d.art3d as art3d
    if not plotdata:
        plotdata = pk.load(open('last_run_plotdata.pk', 'rb'))
    
    fignum = 0
    num_plots = len(plotdata)
    for plot_data in plotdata:
        bar_x_out_of_y(fignum, num_plots)
        tit, beetpos, cloudpos, sunpos = plot_data
        
        colors = ['purple' , 'goldenrod', 'forestgreen', 'steelblue', 'teal']
        colors = colors * (len(cloudpos) // len(colors))

        fig, ax = plt.subplots(1, figsize=(4,4), dpi=200, subplot_kw=dict(projection='3d'))
        ax.set_title(tit)
        # ax.scatter(0, 0, c='maroon')
        ax.scatter(beetpos[0],
                    beetpos[1],
                    beetpos[2],
                    c = 'red', zorder=2,
                    s=2)
        ax.scatter(sunpos[0],
                    sunpos[1],
                    sunpos[2],
                    c = 'gold', zorder=2,
                    s=2)

        # circ = patches.Circle((sunpos[0], sunpos[1]), radius=0.02, transform=ax.transData, fill=False, color='purple', linestyle='--')
        # ax.add_patch(circ)
        
        # for i in ["x","y","z"]:
        #     circle = patches.Circle((sunpos[0], sunpos[1]), radius=0.02, transform=ax.transData, fill=False, color='purple', linestyle='--')
        #     ax.add_patch(circle)
        #     art3d.pathpatch_2d_to_3d(circle, z=sunpos[2], zdir=i)

        ax.scatter(cloudpos[0], 
                    cloudpos[1],
                    cloudpos[2],
                    # c = colors,
                    s = 1
                    )

        # For full galaxy
        # ax.set_xlim(-12, 12)
        # ax.set_ylim(-12, 12)
        # ax.set_zlim(-1, 1)

        # For Sun center
        ax.set_xlim(sunpos[0] - 0.4, sunpos[0] + 0.4)
        ax.set_ylim(sunpos[1] - 0.4, sunpos[1] + 0.4)        
        ax.set_zlim(sunpos[2] - 0.4, sunpos[2] + 0.4)


        # For Beet center
        # ax.set_xlim(beetpos[0] - 0.002, beetpos[0] + 0.002)
        # ax.set_ylim(beetpos[1] - 0.002, beetpos[1] + 0.002)
        # ax.set_zlim(beetpos[2] - 0.002, beetpos[2] + 0.002)

        plt.savefig(f'../figures/fig_{fignum:03d}.png')
        plt.close()
        fignum +=1



def make_movie():
    command = "ffmpeg -y -framerate 25 -i ../figures/fig_%3d.png -c:v libx264 -vb 20M -pix_fmt yuv420p -filter:v 'setpts=2*PTS' -y ../figures/movie.mp4"
    os.system(command)


# %%
if __name__ in '__main__':
    # run_simulation(end_time=(250|units.Myr), 
    #                 timestep=(0.1|units.Myr),
    #                 n_oort_objects=500,
    #                 detection_radius=120|units.pc)
    make_plots()
    make_movie()