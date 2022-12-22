"""
@author: Jasper, Rahul, Konstantinos
"""

import numpy as np
import pandas as pd
from amuse.lab import units, Particles
from custom_classes.star import Star
from custom_classes.test_particles import TestParticles
from amuse.couple import bridge

def run_lonely_planet(detections_df: pd.DataFrame) -> tuple([pd.Series, pd.Series]):
    """
    Runs a lonely planet simulation of the particles described in detections_df.
    Essentially uses the same method as the production runs, only excluding Beet and MWG


    Parameters:
    -----------
    detections_df: pandas.DataFrame
        DataFrame containing detection times and relative positions and velocities


    Returns:
    -----------
    aphelion_time: pandas.Series(float)
        Moment of aphelion passage for each particle in Myr

    aphelion: pandas.Series(float) 
        Aphelion distance for each particle in au
        
    """
    sun = Particles(1)
    sun.mass = 1|units.MSun
    sun.position = (0,0,0) | units.au
    sun.velocity = (0,0,0) | units.kms

    num_detections = len(detections_df)
    asteroids = Particles(num_detections)
    asteroids.mass = 0|units.MSun
    asteroids.position = [(x, y, z) | units.au for x, y, z in 
                            zip(detections_df.x, detections_df.y, detections_df.z)]
    asteroids.velocity = [(vx, vy, vz) | units.kms for vx, vy, vz in 
                            zip(detections_df.vx, detections_df.vy, detections_df.vz)]

    timestep = 10 | units.kyr
    SUN = Star(particles=sun, timestep=1|units.Myr)
    AST = TestParticles(particles=asteroids, timestep=timestep)

    system = bridge.Bridge(use_threading=True)
    system.add_system(AST, (SUN,))
    system.timestep = 1 | units.kyr

    channel_out = system.particles.new_channel_to(asteroids)

    time = 0 | units.yr

    end_time = 1000|units.kyr
    aphelion = np.ones(len(asteroids)) * 4e5
    min_dist_times = np.ones(len(asteroids))
    while time < end_time:
        time += timestep
        system.evolve_model(time)
        channel_out.copy()
        distances = asteroids.position.lengths().value_in(units.au)
        selection = [np.where(distances < aphelion)]
        if np.count_nonzero(selection) > 0:
            aphelion[selection] = distances[selection]
            min_dist_times[selection] = time.value_in(units.Myr)

    start_times= detections_df['time']
    aphelion_time = min_dist_times + start_times.values

    return pd.Series(aphelion_time), pd.Series(aphelion)

