#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:32:25 2022
@author: Konstantinos, Rahul, Jasper
"""
# Python Modules, gotta have 'em
import numpy as np
import pickle as pk
import pandas as pd
import os

# General AMUSE Stuff
from amuse.units import units, constants
from amuse.lab import Particles
from amuse.community.seba.interface import SeBa
from amuse.couple import bridge
from amuse.ic import isotropic_cloud as ic # Has generate oort cloud

# Field/particle codes:
from custom_classes.star import Star
from custom_classes.test_particles import TestParticles
from custom_classes.galactic_potential import MilkyWay_galaxy

# Helper functions
from run_movie_maker import make_movie, make_plots
from progressbar import bar_x_out_of_y
from datetime import datetime # for timing


def get_beet_posrel_to_sun(phi, d=168):
    '''
    Parameters
    ----------
    phi : float
        angle in radians. Beet's angle to the sun
    d : float, optional
        sun beet distance in parsec. The default is 168.

    Returns
    -------
    Beet_x : float
        beet's x distance to the sun in parsec.
    Beet_y : float
        beet's y distance to the sun in parsec.

    '''
    sign = 1
    if np.abs(phi) < np.pi/2:
        # beet is further away from the GC than the sun
        sign = -1
    Beet_x = sign * np.cos(phi) * d
    Beet_y = np.sin(phi) * d

    return Beet_x, Beet_y


def get_theta(Beet_x, Beet_y, sun_x):
    '''
    Parameters
    ----------
    Beet_x : float
        Beet's X distance from the sun
    Beet_y : Beet's distance from the sun
        Beet's Y distance from the sun.
    sun_x : float
        Sun's distance from the galactic centre.

    Returns
    -------
    theta, Beet's polar angle with the galactic centre

    '''
    tan_theta = (Beet_x + sun_x) / Beet_y
    theta = np.arctan(tan_theta)
    return theta


def delete_bound(cloud: Particles, star: Particles) -> None:
    """
    Checks for - and deletes - gravitationally bound asteroids.

    Parameters
    ----------
    cloud: Particles
        Asteroids to check

    star: Particles 
        Star that the bound asteroids would be orbiting
        
    """
    relpos = cloud.position - star.position
    distances = relpos.lengths()
    
    relvel = cloud.velocity - star.velocity
    velocities = relvel.lengths()
    
    # for-loop, more legible and only has to run once so slow is fine
    bound_indices = []
    for index, (dist, vel) in enumerate(zip(distances, velocities)):
        k = 0.5*vel**2
        p = constants.G * star.mass / dist
        if p > k:   # total energy > 0 means unbound
            bound_indices.append(index)

    # let the good people know how many soldiers we lost:
    print('DARKNESS NUMBER:', len(bound_indices), f'(deleted {len(bound_indices)} bound objects)')
    
    # and remove them from the set
    cloud.remove_particle(cloud[bound_indices])



def detect_encounters(cloud: Particles, star: Particles, model_time, detections: list, detection_radius) -> None:
    """
    Checks for asteroids within the detection radius around a star.
    
    Distances to the sun are computed, and if one of the distances is smaller than the
    detection radius, the particle is removed from the set and its information is stored.

    The function loops recursively until no particles are left within the radius. 

    Parameters
    ----------
    cloud: Particles
        particleset to check
    
    star: Particles
        relevant star
    
    model_time: ScalarQuantity
        current model time. gets stored with detection info for statistics
    
    detections: list
        list containing all past detected particles. new detections will
        get appended to this list. 
        
        Contains: key, model_time, x, y, z, vx, vy, vz (relative to star)
    
    detection_radius: ScalarQuantity
        detection radius around the star
    """
    relpos = cloud.position - star.position
    distances = relpos.lengths()        # array of distances to the Sun

    if distances.min() < detection_radius:  # detection!
        cpi = distances.argmin()        # closest particle index
        relx, rely, relz = relpos[cpi] 
        closest = cloud[cpi]
        relvx, relvy, relvz = closest.velocity - star.velocity[0] 
        detections.append([closest.key, model_time, relx, rely, relz, relvx, relvy, relvz])
        cloud.remove_particle(closest)  # remove it to prevent repeat detections
        # rinse and repeat
        detect_encounters(cloud, star, model_time, detections, detection_radius)



def run_simulation(end_time = 100 | units.Myr, 
                    timestep_pre_sn = 1e-3 | units.Myr,
                    timestep_post_sn = 0.1 | units.Myr,
                    timestep_detection = 0.1 | units.Myr,
                    n_oort_objects: int = 100,
                    phi: float = -20 * np.pi/180,
                    detection_radius = 1. | units.pc,
                    plot_interval = 1. | units.Myr,
                    outdir: str = "../runs/temp",
                    random_seed = None):

    """
    Meat and potatoes function that executes the simulations and saves the results
    """
    # Store run parameters for future reference
    run_params = {'t_end': end_time, 
                    'timesteps': (timestep_pre_sn, timestep_post_sn, timestep_detection),
                    'n_objects': n_oort_objects,
                    'det_rad': detection_radius}

    pk.dump(run_params, open(f'{outdir}/run_params.pk', 'wb'))
    
    # Make static MWG potential
    MWG = MilkyWay_galaxy()
    
    # sun
    sun = Particles(1)
    sun.mass = 1 | units.MSun
    # Values from https://ui.adsabs.harvard.edu/abs/1996AJ....111..794K/abstract
    
    sun_to_GC = -8400 # Suns distance to galactic centre [parsec]
    sun.position = (sun_to_GC, 0, 0) | units.pc

    sun_vel = MWG.vel_circ(sun.x,sun.y)
    sun.velocity = sun_vel * (0., 1., 0.)

    # beet
    beet = Particles(1)
    beet.mass = 21 | units.MSun
    beet_rel_x, beet_rel_y = get_beet_posrel_to_sun(phi, d=135)
    
    # Initialise Beet based on the given phi angle and a circular orbit
    beet.position = ((beet_rel_x, beet_rel_y, 0) | units.pc) + sun.position
    theta = get_theta(beet_rel_x, beet_rel_y, sun_to_GC)
    beet_vel = MWG.vel_circ(beet.x,beet.y)
    beet.velocity = beet_vel * (np.cos(theta), np.abs(np.sin(theta)), 0)

    print("Initialising and pre-evolving SeBa")
    beet_seba = SeBa()
    beet_seba.set_metallicity(0.024) # From Dolan & Mathews, https://arxiv.org/pdf/1406.3143v2.pdf

    # Make Beet
    BEET = Star(beet, timestep_pre_sn, stellar_evo_code=beet_seba)
    
    # Pre-evolve:
    BEET.pre_evolve(8.4|units.Myr)
    # An age of 8.4 Myr works pretty well, and that puts the current mass at ~18 MSun
    # The only downside is that rapid mass loss starts at ~8 Myr, which the oort cloud would notice.
    
    print("Generating Sun and Oort cloud")

    # Use AMUSE function for generating the Oort cloud
    beet_cloud = ic.new_isotropic_cloud(number_of_particles=n_oort_objects,
                            m_star = BEET.particles[0].mass,
                            a_min = 60_000 | units.AU,
                            a_max = 300_000 | units.AU,
                            q_min = 20_000 | units.AU,
                            seed = random_seed)
    
    beet_cloud.position += beet.position
    beet_cloud.velocity += beet.velocity

    # Make test particle objects (cloud and sun)
    OORT = TestParticles(beet_cloud, timestep_pre_sn)
    SUN = TestParticles(sun, timestep_pre_sn)    

    print("Building bridges")
    # Mix it all together with bridge
    milky_way = bridge.Bridge(use_threading=True)
    milky_way.add_system(OORT, (BEET, MWG)) # Oort cloud cares about MW and Beet
    milky_way.add_system(BEET, (MWG,))      # Beet only cares about MW
    milky_way.add_system(SUN, (MWG,))       # Sun only cares about MW
    milky_way.timestep = timestep_pre_sn

    def update_all_timesteps(new_ts):
        # quick way to change all internal timesteps
        OORT.timestep = new_ts
        BEET.timestep = new_ts
        SUN.timestep = new_ts
        milky_way.timestep = new_ts

    # Create all the channels
    channel_out_cloud = milky_way.particles.new_channel_to(beet_cloud)
    channel_in_cloud = beet_cloud.new_channel_to(milky_way.particles)
    channel_out_sun = milky_way.particles.new_channel_to(sun)
    channel_out_beet = milky_way.particles.new_channel_to(beet)
    
    print("starting pre-supernova run")
    # Lists that will be filled during the run:
    detections = []         # contains all the information about detected particles
    plotdata = []           # huge list containing all positions at all plotting timesteps
    
    model_time = 0 | units.yr
    last_plot_time = model_time - plot_interval
    
    # Pre-supernova run:
    start_time = datetime.now()
    while True:       
        model_time += timestep_pre_sn   # Stricter timestep because orbits are tighter
        milky_way.evolve_model(model_time)

        if BEET.stellar.particles[0].stellar_type.value == 14:  # if BEET has gone supernova
            print(f'SUPERNOVA at time {model_time.in_(units.Myr)}!!! mass={beet.mass}')
            break       # exit pre-supernova regime

    print(f'Integration complete! Duration: {datetime.now() - start_time}')
    start_time = datetime.now()         # reset timer
    
    # Delete remaining bound objects
    channel_out_cloud.copy()
    channel_out_beet.copy()
    delete_bound(beet_cloud, beet)
    channel_in_cloud.copy()

    update_all_timesteps(timestep_post_sn)

    # main part of the simulation:
    detecting = False
    while(model_time < end_time):
        model_time += milky_way.timestep
        milky_way.evolve_model(model_time)
        
        # encounter detection part
        if detecting: 
            channel_out_cloud.copy()
            channel_out_sun.copy()
            detect_encounters(beet_cloud, sun, model_time, detections, detection_radius)
            channel_in_cloud.copy()

            # progress bar (including number of detections)
            bar_x_out_of_y(x=model_time.in_(end_time.unit).round(1), 
                            y=end_time, 
                            text=f'{len(detections)} detections')
        
        # plotting part
        if (last_plot_time + plot_interval <= model_time):
            #  save data for plotting:
            channel_out_cloud.copy()
            channel_out_sun.copy()
            channel_out_beet.copy()
            
            plotdata.append((str(model_time.value_in(units.Myr)) + ' Myr',
                            beet.position.value_in(units.kpc)[0],
                            beet_cloud.position.value_in(units.kpc),
                            sun.position.value_in(units.kpc)[0]))

            last_plot_time = model_time

            # less frequent progress bar because nothing interesting happens anyway
            if not detecting:
                bar_x_out_of_y(model_time.in_(end_time.unit).round(1), end_time, text='')
            
            # Check if we should be detecting:
            dist_sun_beet = (sun.position - beet.position).lengths()[0].value_in(units.pc)
            # if sun and beet are close, then start checking for encounters and change timestep
            if dist_sun_beet < 75: # in parsec. 
                if not detecting:
                    detecting = True
                    update_all_timesteps(timestep_detection)

            elif detecting:     # if we are detecting but shouldn't be: stop detecting
                detecting = False
                update_all_timesteps(timestep_post_sn)   

    duration = datetime.now() - start_time
    print(f'Integration complete! Duration: {duration}')

    # write out all of the data:
    pk.dump(plotdata, open(f"{outdir}/plotdata.pk", 'wb'))

    if len(detections) > 0:     # no point storing an empty array
        detection_df = pd.DataFrame(detections)
        detection_df.columns = "key", "time", "x", "y", "z", "vx", "vy", "vz"
        print(detection_df)
        pk.dump(detection_df, open(f'{outdir}/detections.pk', 'wb'))
    

if __name__ in '__main__':
    # quick test run:
    detection_radius = 1 |units.pc
    phi = -34.5 * np.pi/180
    
    run_name="quick_test"
    outdir = f"../runs/{run_name}"
    if not os.path.exists(outdir):
        os.system(f'mkdir {outdir}')
    if not os.path.exists(outdir+"/figs"):
        os.system(f'mkdir {outdir}/figs')
    
    run_simulation(end_time=75.|units.Myr,
                    timestep_pre_sn=0.01|units.Myr,
                    timestep_post_sn=1|units.Myr,
                    timestep_detection=.1|units.Myr,
                    plot_interval=1|units.Myr,
                    phi=phi,
                    n_oort_objects=10_000,
                    detection_radius = detection_radius,
                    outdir=outdir,
                    random_seed=1234)

    make_plots(focus='sun', zoom=0.1, mask_z_distance=True, 
                start_plots_at_num=0, detection_radius=detection_radius,
                outdir=outdir)

    make_movie(outdir + "/")
