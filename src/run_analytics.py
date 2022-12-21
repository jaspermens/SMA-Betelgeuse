# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:26:52 2022

@author: Jasper, Rahul, Konstantinos

Detection analytics
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from amuse.lab import units, constants
import pickle as pk 


def get_detections_from_run(run_name):
    d =  pk.load(open(f'../runs/{run_name}/detections.pk', 'rb'))
    for v in ('vx', 'vy', 'vz'):
        d[v] = [vel.value_in(units.kms) for vel in d[v]]

    for p in ('x', 'y', 'z'):
        d[p] = [pos.value_in(units.au)/1e3 for pos in d[p]]

    d.time = [t.value_in(units.Myr) for t in d.time]

    return d


def plot_det_relvel(det):
    relvel = np.sqrt((det.vx**2 + det.vy**2 + det.vz**2))

    relvel.plot.hist(figsize=[6,4])
    plt.title("Relative velocities of detections wrt the Sun")
    plt.xlabel('relative velocity [km/s]')
    plt.show()


def plot_det_relpos(det):
    relvel = np.sqrt((det.x**2 + det.y**2 + det.z**2))

    relvel.plot.hist(figsize=[6,4])
    plt.title("Distances between the detection point and the Sun")
    plt.xlabel('distance [10^3 au]')
    plt.show()


def plot_relpos_kde(det):
    det.plot.hexbin(x='x', y='y', gridsize=20)
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


def plot_min_sep(det):
    orbels = get_orbital_elements(det)

    smas = [orbels[i][0].value_in(units.au) for i in range(len(orbels))]
    eccs = [orbels[i][1] for i in range(len(orbels))]

    min_seps = pd.Series([float(sma) * (1.-float(ecc)) for sma, ecc in zip(smas, eccs)])
    print(min_seps)
    min_seps.plot.kde()
    plt.axvline((1|units.pc).value_in(units.au), ls='--')
    plt.title("KDE for the keplerian periapsis")
    plt.xlabel("Periapsis distance [au]")
    plt.show()


def combine_detections_from_runs(*runs):
    return pd.concat([get_detections_from_run(run) for run in runs])
        

if __name__ in '__main__':
    run_name ='milly_1'
    # detections= get_detections_from_run(run_name)
    detections = combine_detections_from_runs('milly_1', 'milly_2', 'milly_3')
    plot_relpos_kde(detections)