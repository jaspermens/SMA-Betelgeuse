"""
@author: Jasper, Rahul, Konstantinos

Module for reading the detection data from production runs.

contains: 
get_detections_from_run
get_detections_from_runs
"""

import numpy as np
import pandas as pd
import pickle as pk
from amuse.lab import units

def get_detections_from_run(run_name: str) -> pd.DataFrame:
    """
    Reads and returns the pandas dataframe containing the detection data from
    ../runs/<run_name>/detections.pk
    
    Parameters:
    -----------
    run_name: str
        Directory name containing the detection data.
        ../runs/<run_name> MUST contain detections.pk
    """
    det_df =  pk.load(open(f'../runs/{run_name}/detections.pk', 'rb'))
    for v in ('vx', 'vy', 'vz'):
        det_df[v] = [vel.value_in(units.kms) for vel in det_df[v]]

    for p in ('x', 'y', 'z'):
        det_df[p] = [pos.value_in(units.au)/1e3 for pos in det_df[p]]

    det_df.time = [t.value_in(units.Myr) for t in det_df.time]
    det_df['relvel'] = det_df.apply(lambda x: np.sqrt(x.vx**2+ x.vy**2+ x.vz**2), axis=1)
    det_df['relpos'] = det_df.apply(lambda x: np.sqrt(x.x**2+ x.y**2+ x.z**2), axis=1)
    
    return det_df


def get_detections_from_runs(*runs: list) -> pd.DataFrame:
    """
    Combines detection data from several runs
    """
    return pd.concat([get_detections_from_run(run) for run in runs])

 
if __name__ in '__main__':
    print(get_detections_from_run('temp'))