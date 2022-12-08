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


class star:
    def __init__(self, 
                particles, 
                timestep, 
                mdot=(0.1|units.MSun / units.yr),
                stellar_evo_code=None):
        self.particles = particles
        self.model_time = 0 | units.yr
        self.timestep = timestep
        self.mdot = mdot

        self.stellar = stellar_evo_code
        if self.stellar:
            self.stellar.particles.add_particles(self.particles)
            self.channel_mass = self.stellar.particles.new_channel_to(self.particles)

    def pre_evolve(self, end_time):
        if not self.stellar:
            print("This star does not have a stellar evo code!")
            return
        self.stellar.evolve_model(end_time)
        self.channel_mass.copy()

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
        if self.stellar:
            self.stellar.evolve_model(self.stellar.model_time + self.timestep)
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
    
    def change_timestep(self, new_timestep):
        self.timestep = new_timestep
        
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
            
    def change_timestep(self, new_timestep):
        self.timestep = new_timestep


def bar_x_out_of_y(x, y, text: str='') -> None:
    maxbars = 20
    nbars = int(x / y * maxbars)
    print('\rProcessing: | ' + '█' * nbars + '-' * (maxbars - nbars) + ' |', end=' ')
    print(f'~ {x}/{y} {text}', end='')
    if x == y:
        print('')


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

def energy_check(kinetic, potential, particles):
    # PRAY that AMUSE works this one out by itsellf.
    k = particles.kinetic_energy()
    p = particles.potential_energy()
    # Append to existing lists
    kinetic.append(k)
    potential.append(p)
    return kinetic, potential

def delete_bound(cloud, beet):
    print("BEET MASS POST SUPERNOVA:", beet.mass)
    relpos = cloud.position - beet.position
    relvel = cloud.velocity - beet.velocity
    distances = relpos.lengths()
    velocities = relvel.lengths()
    bound = []
    j=1 # counter
    for i  in range(len(distances)):
        k = 0.5*velocities[i]**2
        p = constants.G * beet.mass / distances[i]
        if p>k:
            bound.append(i)
            # maybe a-e plot to see who sticks around
            j=j+1
    print('DARKNESS NUMBER: ', j)
    cloud.remove_particle(cloud[bound])
    return bound

def detect_encounters(cloud, sun, model_time, detections, detection_keys, detection_radius):
    relpos = cloud.position - sun.position
    distances = relpos.lengths()

    if distances.min() < detection_radius:
        # print('detection!!')
        cpi = distances.argmin()  # closest particle index
        closest = cloud[cpi]
        relvx, relvy, relvz = closest.velocity - sun.velocity[0]
        relx, rely, relz = relpos[cpi]
        detections.append([closest.key, model_time, relx, rely, relz, relvx, relvy, relvz])
        detection_keys.append(closest.key)
        cloud.remove_particle(cloud[cpi])
        detect_encounters(cloud, sun, model_time, detections, detection_keys, detection_radius)


def run_simulation(end_time=(100|units.Myr), 
                    timestep=(1e-3|units.Myr),
                    timestep_after_sn = (2 | units.Myr),
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
    beet.mass = 21 | units.MSun
    beet.position = (26660, 0, 0) | units.lightyear
    beet.velocity = (0, 250, 0) | units.kms
    
    beet_seba = SeBa()
    beet_seba.set_metallicity(0.024) # From Dolan & Mathews, https://arxiv.org/pdf/1406.3143v2.pdf

    BEET = star(beet, timestep, stellar_evo_code=beet_seba)
    # Pre-evolution:
    BEET.pre_evolve(8.4|units.Myr)
    # An age of 8.4 Myr works pretty well, and that puts the current mass at 18 or so MSun
    # only downside is that rapid mass loss starts at ~8 Myr, which the oort cloud would notice.
    
    # sun
    sun = Particles(1)
    sun.mass = 1 | units.MSun
    sun.position = (26660*1.02, 0, 0) | units.lightyear
    sun.velocity = (0, 250*.98, 0) | units.kms

    # Oort cloud
    beet_cloud = ic.new_isotropic_cloud(number_of_particles=n_oort_objects,
                            m_star = BEET.particles[0].mass,
                            a_min = 60_000 | units.AU,
                            a_max = 300_000 | units.AU,
                            q_min = 20_000 | units.AU,
                            seed=42069)

    beet_cloud.position += beet.position
    beet_cloud.velocity += beet.velocity

    # Make test particles (cloud and sun)
    OORT = test_particles(beet_cloud, timestep)
    SUN = test_particles(sun, timestep)    

    # Make static MWG potential
    MWG = MilkyWay_galaxy()

    # Mix it all together with bridge
    milky_way = bridge.Bridge(use_threading=True)
    milky_way.add_system(OORT, (BEET, MWG))
    milky_way.add_system(BEET, (MWG,))
    milky_way.add_system(SUN, (MWG,))
    milky_way.timestep = timestep

    # Make channels
    channel_out = milky_way.particles.new_channel_to(beet_cloud)
    channel_in = beet_cloud.new_channel_to(milky_way.particles)
    channel_out_sun = milky_way.particles.new_channel_to(sun)
    channel_out_beet = milky_way.particles.new_channel_to(beet)
    
    # Start runnin'
    model_time = 0 | units.yr
    plot_interval = 1 | units.Myr
    last_plot_time = model_time - plot_interval
    plotdata = []
    detections = []
    detection_keys = []
    gone_supernova = False
    start_time = datetime.now()
    while(model_time < end_time):
        model_time += timestep
        
        milky_way.evolve_model(model_time)
            
        # collision detection
        channel_out.copy()
        channel_out_sun.copy()          # I don't know if these two do anything
        channel_out_beet.copy()         # It runs just as well without them two channels
        detect_encounters(beet_cloud, sun, model_time, detections, detection_keys, detection_radius)
        channel_in.copy()

        #  Saving data for future plotting
        if (last_plot_time + plot_interval <= model_time):
            # Progress check
            bar_x_out_of_y(model_time.in_(end_time.unit), end_time, text=f'{len(detections)} detections')
            

            plotdata.append((str(model_time.value_in(units.Myr)) + ' Myr',
                            BEET.particles[0].position.value_in(units.kpc),
                            beet_cloud.position.value_in(units.kpc),
                            sun.position.value_in(units.kpc)[0]))

            last_plot_time = model_time

        # 4 is for Core Helium Burning
        # 14 is black hole
        # Ignore the AMUSE/SeBa documentation

        # print(BEET.stellar.particles[0].stellar_type, BEET.stellar.particles[0].stellar_type.value)

        if BEET.stellar.particles[0].stellar_type.value == 14 and gone_supernova==False:
            print('\n SUPERNOVA IS NOW!!')
            # Flag that ensures timestep changes only once
            gone_supernova = True
            # Delete bounds
            channel_out.copy()
            delete_bound(beet_cloud, beet)
            channel_in.copy()
            # Timestep changes
            OORT.change_timestep(timestep_after_sn)
            BEET.change_timestep(timestep_after_sn)
            SUN.change_timestep(timestep_after_sn)
            milky_way.timestep = timestep_after_sn
            timestep = timestep_after_sn

    duration = datetime.now() - start_time
    print(f'Duration: {duration}')
    
    pk.dump(plotdata, open('last_run_plotdata.pk', 'wb'))

    if len(detections) > 0:
        detection_df = pd.DataFrame(detections)
        detection_df.columns = "key", "time", "x", "y", "z", "vx", "vy", "vz"
        print(detection_df)
        pk.dump(detection_df, open('detections.pk', 'wb'))
    

def make_plots(plotdata=None, focus="beet", zoom=None, mask_outside_detection_radius=False, start_plots_at_num=0, detection_radius=1|units.pc):
    import matplotlib.patches as patches
    import mpl_toolkits.mplot3d.art3d as art3d

    if not plotdata:
        plotdata = pk.load(open('last_run_plotdata.pk', 'rb'))
    
    focus_options = ['beet', 'sun', 'galaxy']
    if not focus in focus_options:
        raise IndexError(f"focus must be in {focus_options}, got: {focus}. aborting...")

    colors = ['purple' , 'goldenrod', 'forestgreen', 'steelblue', 'teal']
    colors = colors * ((len(plotdata[0][2].T[0])) // len(colors)+1)

    plotdata = plotdata[start_plots_at_num:]
    fignum = 0
    num_plots = len(plotdata)
    for plot_data in plotdata:
        tit, beetpos, cloudpos, sunpos = plot_data

        # alphas = np.array([0.1, 1])
        relative_z = np.abs((cloudpos.T[2]-sunpos[2]) / detection_radius.value_in(units.kpc)) # relz in detection radii
        
        alphas = np.ones(relative_z.shape, dtype=float)
        if mask_outside_detection_radius:
            alphas[relative_z >= 1] = 0.1
            alphas[relative_z < 1] = 1

        fig, ax = plt.subplots(1, figsize=(4,4), dpi=200) #, subplot_kw=dict(projection='3d'))
        ax.set_title(tit)
        # ax.scatter(0, 0, c='maroon')
        ax.scatter(beetpos[0],
                    beetpos[1],
                    c = 'red', zorder=2,
                    s=2)
        ax.scatter(sunpos[0],
                    sunpos[1],
                    c = 'gold', zorder=2,
                    s=2)
        ax.scatter(cloudpos.T[0], 
                    cloudpos.T[1],
                    c = colors[:len(cloudpos.T[0])],
                    alpha=alphas,
                    s = 1
                    )
        
        # Circle around the sun to show detection radius:
        circ = patches.Circle((sunpos[0], sunpos[1]), radius=detection_radius.value_in(units.kpc), transform=ax.transData, fill=False, color='purple', linestyle='--')
        ax.add_patch(circ)

        if focus == "galaxy":
            # For full galaxy
            if not zoom:
                zoom = 10
            ax.set_xlim(-zoom, zoom)
            ax.set_ylim(-zoom, zoom)
            # ax.set_zlim(-1, 1)
        
        elif focus == "sun":
            # For Sun center
            if not zoom:
                zoom = 0.4
            ax.set_xlim(sunpos[0] - zoom, sunpos[0] + zoom)
            ax.set_ylim(sunpos[1] - zoom, sunpos[1] + zoom)        
            # ax.set_zlim(sunpos[2] - 0.4, sunpos[2] + 0.4)
        
        elif focus == "beet":
            # For Beet center
            if not zoom:
                zoom = 0.02
            ax.set_xlim(beetpos[0] - zoom, beetpos[0] + zoom)
            ax.set_ylim(beetpos[1] - zoom, beetpos[1] + zoom)
            # ax.set_zlim(beetpos[2] - 0.002, beetpos[2] + 0.002)

        plt.savefig(f'../figures/fig_{fignum:03d}.png')
        plt.close()
        fignum +=1
        bar_x_out_of_y(fignum, num_plots)



def make_movie():
    command = "ffmpeg -y -framerate 25 -i ../figures/fig_%3d.png -c:v libx264 -hide_banner -vb 20M -loglevel panic -pix_fmt yuv420p -filter:v 'setpts=2*PTS' -y ../figures/movie.mp4"
    os.system(command)


if __name__ in '__main__':
    detection_radius = 1|units.pc
    run_simulation(end_time=(240|units.Myr), 
                    timestep=(0.001|units.Myr),
                    timestep_after_sn = (0.1|units.Myr),
                    n_oort_objects=100_000,
                    detection_radius=detection_radius)
    make_plots(focus='sun', zoom=.020, mask_outside_detection_radius=True, start_plots_at_num=220, detection_radius=detection_radius)
    make_movie()
