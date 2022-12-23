# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:26:52 2022

@author: Jasper, Rahul, Konstantinos

Detection analytics
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from amuse.lab import units
from read_detections import get_detections_from_runs
from lonely_planet import run_lonely_planet


def plot_relvel(det):
    det.relvel.plot.hist(figsize=[6,4], bins=100)
    plt.title("Relative velocities of detections wrt the Sun")
    plt.xlabel('relative velocity [km/s]')
    plt.show()


def plot_relpos(det):
    det.relpos.plot.hist(figsize=[6,4], bins=100)
    plt.title("Distances between the detection point and the Sun")
    plt.xlabel('distance [10^3 au]')
    plt.show()


def plot_relpos_hex(det):
    det.plot.hexbin(x='x', y='y', gridsize=20, cmap='viridis')
    plt.title("Relative positions of detected particles")
    plt.xlabel("x [$10^3$ au]")
    plt.ylabel("y [$10^3$ au]")
    plt.show()


def plot_relvel_hex(det):
    det.plot.hexbin(x='x', y='y', C='vz', gridsize=20, cmap='viridis')
    plt.title("Relative positions of detected particles")
    plt.xlabel("x [$10^3$ au]")
    plt.ylabel("y [$10^3$ au]")
    plt.show()


def plot_arrival_time(det):
    det.time.plot.hist(figsize=[6,4])
    plt.title("Histogram of arrival time")
    plt.ylabel("Frequency")
    plt.xlabel("Arrival time (Myr)")
    plt.show()


def get_orbital_elements(det):
    from amuse.ext.orbital_elements import orbital_elements_for_rel_posvel_arrays as orbels
    els = []
    for _, d in det.iterrows():
        relpos = [d.x*1e3, d.y*1e3, d.z*1e3] | units.au
        relvel = [d.vx, d.vy, d.vz] | units.kms

        els.append([*orbels(relpos, relvel, 1.|units.MSun)])
    
    return els


def aphelion_cdf_kep(det):    
    """
    Produces a cumulative histogram of the analytical (keplerian) aphelia
    """
    orbels = get_orbital_elements(det)

    smas = [orbels[i][0].value_in(units.au) for i in range(len(orbels))]
    eccs = [orbels[i][1] for i in range(len(orbels))]

    min_seps = pd.Series([float(sma) * (1.-float(ecc)) for sma, ecc in zip(smas, eccs)])
    # print(min_seps)
    # min_seps.plot.kde()
    min_seps.plot.hist(cumulative=True, bins=100, density=True)
    plt.title('CDF of keplerian periapsis')
    plt.axvline((1|units.pc).value_in(units.au), ls='--')
    plt.ylabel('Cum. probability')
    plt.xlabel("Periapsis distance [au]")
    seps = np.arange(0, 206_000, 100)
    plt.plot(seps, 2.37953599e-11*np.power(seps,2), 
                    label='quadratic fit, $P=2.38\cdot10^{-11}a_{p}^2$')
    plt.legend()
    plt.show() 


def aphelion_cdf_lp(detections):
    """
    Produces a cumulative histogram of the aphelion distances using lonely planet data
    """
    _, d = run_lonely_planet(detections_df=detections)
    seps = np.arange(0, 206_000, 100)

    print(len(d[d < 1e3]))
    d.plot.hist(bins=100, cumulative=True, density=True)

    plt.axvline((1|units.pc).value_in(units.au), ls='--', label='1 parsec')
    
    plt.plot(seps, 2.37953599e-11*np.power(seps,2), 
                    label='quadratic fit, $P=2.38\cdot10^{-11}a_{p}^2$')
    
    plt.title('Cumulative histogram of numerical periapsis distance')
    plt.xlabel("Periapsis distance [au]")
    plt.ylabel('Cum. probability')
    
    plt.legend()
    plt.show() 


def arrival_time_histogram_lp(detections):
    """
    Produces a histogram of the aphelion times using lonely planet data
    """
    t, _ = run_lonely_planet(detections_df=detections)

    t.plot.hist(bins=50)
    plt.title("Histogram of arrival time")
    plt.ylabel("Frequency")
    plt.xlabel("Arrival time (Myr)")
    plt.show()
        

if __name__ in '__main__':
    run_name ='milly_1'

    # detections= get_detections_from_run(run_name)
    detections = get_detections_from_runs(*[f'milly_{n}' for n in range(1,12)])
    arrival_time_histogram_lp(detections)
    