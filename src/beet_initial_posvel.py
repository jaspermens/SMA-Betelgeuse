#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 17:19:56 2022

@author: Rahul, Jasper, Konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from amuse.lab import units

def get_initial_posvel_beet():
    # initial positions and velocities of beet and sun (values from simbad)
    # beet
    beet_ra = "05h55m10.30536s" 
    beet_dec = "+07d24m25.4304s"
    # 27.54
    beet_pm_ra_cosdec = 27.54*u.mas/u.yr
    # 11.3 are normal
    beet_pm_dec = 11.3*u.mas/u.yr
    beet_rad_vel = 21.91*u.km/u.s
    # 168 is normal
    beet_sun_distance = 168*u.pc

    beet_sc = SkyCoord(ra = beet_ra, 
                    dec = beet_dec, 
                    pm_ra_cosdec = beet_pm_ra_cosdec, 
                    pm_dec = beet_pm_dec,
                    radial_velocity = beet_rad_vel)

    # print(beet_sc.galactic.pm_b * beet_sun_distance)
    beet_l = (beet_sc.galactic.l.to(u.rad)).value
    beet_b = (beet_sc.galactic.b.to(u.rad)).value
    pm_l = (beet_sc.galactic.pm_l_cosb.to(u.rad/u.s)).value/beet_b
    pm_b = (beet_sc.galactic.pm_b.to(u.rad/u.s)).value
    beet_pm_l = pm_l/u.s
    beet_pm_b = pm_b/u.s

    beet_x = np.cos(beet_l) * np.cos(beet_b) * beet_sun_distance.to(u.pc)
    beet_y = np.sin(beet_l) * np.cos(beet_b) * beet_sun_distance.to(u.pc)
    beet_z = np.sin(beet_b) * beet_sun_distance.to(u.pc)

    beet_vx = (beet_pm_l * beet_sun_distance.to(u.m) * np.cos(beet_l) * np.cos(beet_b)).to(u.km/u.s)
    beet_vy = (beet_pm_l * beet_sun_distance.to(u.m) * np.sin(beet_l) * np.sin(beet_b)).to(u.km/u.s)
    beet_vz = (beet_pm_b * beet_sun_distance.to(u.m) * np.sin(beet_b)).to(u.km/u.s)
    # print(beet_sc)
    # print(beet_sc.galactic)
    return ((beet_x.value, beet_y.value, beet_z.value), (beet_vx.value, beet_vy.value, beet_vz.value))
    # print(beet_x)
    # print(beet_y)
    # print(beet_z)
    # print(beet_vx)
    # print(beet_vy)
    # print(beet_vz)

if __name__ in '__main__':
    p, v = get_initial_posvel_beet()
    print(p|units.pc, v|units.kms)
