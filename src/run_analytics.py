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
import colorcet as cc
fire = cc.fire

#%% Out of report
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

def plot_relpos(det):
    det.relpos.plot.hist(figsize=[6,4], bins=100)
    plt.title("Distances between the detection point and the Sun")
    plt.xlabel('distance [10^3 au]')
    plt.show()

def plot_relvel_hex(det):
    det.plot.hexbin(x='x', y='y', C='vz', gridsize=20, cmap='viridis')
    plt.title("Relative positions of detected particles")
    plt.xlabel("x [$10^3$ au]")
    plt.ylabel("y [$10^3$ au]")
    plt.show()
#%% In the report
def plot_relvel(det):
    det.relvel.plot.hist( bins=100, color='maroon')
    plt.title("Relative velocities of arriving Oort Objects", fontsize=27, pad=10)
    plt.xlabel('Relative Velocity [Km/s]', fontsize=20)
    plt.ylabel('Frequency', fontsize=20)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.xlim(3,3.5)
    plt.show()
    plt.savefig('relvelo.pdf')

def plot_relpos_hex(det):
    # Black background
    xs = [-250, -250, 250, 250]
    ys = [-250, 250, 250, -250]
    plt.fill(xs, ys, color='k')
    # Plot
    plt.hexbin(x='x', y='y', 
               gridsize=50, cmap = 'cet_fire',data=det)
    #det.plot.hexbin(x='x', y='y', gridsize=50, cmap='cet_fire')
    plt.title("Relative positions of arriving Oort Objects", fontsize=27, pad=10)
    plt.ylabel("y [$10^3$ AU]", fontsize=20)
    plt.xlabel("x [$10^3$ AU]", fontsize=20)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    cb = plt.colorbar()
    for t in cb.ax.get_yticklabels():
     t.set_fontsize(18)
    plt.xlim(-220,180)
    plt.ylim(-215,190)
    plt.savefig('relpos.pdf')

def aphelion_cdf_lp(detections):
    """
    Produces a cumulative histogram of the aphelion distances using lonely planet data
    """
    _, d = run_lonely_planet(detections_df=detections)
    seps = np.arange(0, 206_000, 100)
    plt.rcParams['lines.linewidth'] = 2
    print(len(d[d < 1e3]))
    d.plot.hist(bins=100, cumulative=True, density=True, color='maroon')

    plt.axvline((1|units.pc).value_in(units.au), ls='--', label='1 parsec',
                color='yellow')
    
    plt.plot(seps, 2.37953599e-11*np.power(seps,2),  color='cyan',
                    label='quadratic fit, $P=2.38\cdot10^{-11}a_{p}^2$')
    
    plt.title('Cumulative histogram of numerical periapsis distance',
              fontsize=27, pad=10)
    plt.xlabel("Periapsis distance [AU]", fontsize=20)
    plt.ylabel('Cum. probability', fontsize=20)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.legend(loc='best', fontsize=20)
    plt.savefig('cdf.pdf')

def arrival_time_histogram_lp(detections):
    """
    Produces a histogram of the aphelion times using lonely planet data
    """
    t, _ = run_lonely_planet(detections_df=detections)

    t.plot.hist(bins=50, color='maroon')
    plt.title("Histogram of arrival time", fontsize = 27, pad=10)
    plt.ylabel("Frequency", fontsize=20)
    plt.xlabel("Arrival Time (Myr)", fontsize=20)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.savefig('arrival_time.pdf')
        

if __name__ in '__main__':
    detections = get_detections_from_runs(*[f'milly_{n}' for n in range(1,26)])
    # plt.rcParams.update(params)
    plt.rcParams["figure.figsize"] = (10,8)
    plot_relvel(detections)
   # plot_relpos_hex(detections)
   # aphelion_cdf_lp(detections)
    #arrival_time_histogram_lp(detections)
    