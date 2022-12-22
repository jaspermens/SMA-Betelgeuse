import os
import numpy as np
from amuse.lab import units

from analytic_solver import run_simulation
from run_movie_maker import make_plots, make_movie

if __name__ in '__main__':
    detection_radius = 1 |units.pc
    # 'true' value is -20, but -34.5 gives detections
    phi = -34.5 * np.pi/180
    
    run_name="temp" 
    outdir = f"../runs/{run_name}"
    if not os.path.exists(outdir):
        os.system(f'mkdir {outdir}')
    if not os.path.exists(outdir+"/figs"):
        os.system(f'mkdir {outdir}/figs')
    
    run_simulation(end_time=75.|units.Myr,
                    timestep_pre_sn=0.001|units.Myr,
                    timestep_post_sn=.1|units.Myr,
                    timestep_detection=.1|units.Myr,
                    plot_interval=1|units.Myr,
                    phi=phi,
                    n_oort_objects=1_000,
                    detection_radius=detection_radius,
                    outdir=outdir,
                    random_seed=None)

    make_plots(focus='sun', zoom=0.1, mask_z_distance=True, 
                start_plots_at_num=0, detection_radius=detection_radius,
                outdir=outdir)

    make_movie(outdir + "/")
