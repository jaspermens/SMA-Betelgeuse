import numpy as np
import matplotlib.pyplot as plt
import pickle as pk
import os
from amuse.lab import units
from progressbar import bar_x_out_of_y 
import matplotlib.patches as patches

def make_plots(focus: str = "sun", 
                zoom: float = None, 
                mask_z_distance: bool = True, 
                start_plots_at_num: int = 0, 
                detection_radius = 1|units.pc,
                outdir: str = "../runs/temp") -> None:
    """
    Creates a scatterplot out of each row of plotdata.
    If plotdata is not provided, it is read from a file. 

    Parameters
    ----------
    focus: str
        Either 'galaxy', 'sun', or 'beet'.
        Dictates where the plots are centered

    zoom: float
        Zoom level. half side length in kpc

    mask_z_distance: bool
        Whether or not to remove asteroids at the wrong z
        i.e. particles that won't get detected anyway

    start_plots_at_num: int
        At what row of plotdata to start

    detection_radius: ScalarQuantity
        Detection radius used in the simulation. 
        Allows for visualising the radius as a circle around the Sun

    outdir: str
        Directory where the figures will be stored. 
        More specifically, figures will end up in <outdir>/figs

        Also, plotdata will be read from <outdir>/plotdata.pk
    """

    print("Generating figures")
    plotdata = pk.load(open(outdir + "/plotdata.pk", 'rb'))
    
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

        cloud_x, cloud_y, cloud_z = cloudpos.T

        if mask_z_distance:
            # relz in detection radii
            relative_z = np.abs((cloud_z-sunpos[2]) / detection_radius.value_in(units.kpc)) 
            if len(cloudpos) > 100_000:
                # Plotting a million particles is way too slow, so cap it at 100k                
                cloud_x = cloud_x[relative_z < 1.]
                cloud_y = cloud_y[relative_z < 1.]
                relative_z = relative_z[relative_z < 1.]

            # These lines technically don't do anything if there's more than 50k particles
            # but otherwise they dim alpha for irrelevant particles
            alphas = np.ones(relative_z.shape, dtype=float)
            alphas[relative_z >= 1] = 0.1
            alphas[relative_z < 1] = 1

        _, ax = plt.subplots(1, figsize=(4,4), dpi=200) 
        ax.set_title(tit)   # Title from plotdata

        if focus == "galaxy":
            # For full galaxy
            if not zoom:
                zoom = 10
            focus_x = 0
            focus_y = 0

        elif focus == "sun":
            # For Sun center
            if not zoom:
                zoom = 0.4
            focus_x = sunpos[0]
            focus_y = sunpos[1]    
        
        elif focus == "beet":
            # For Beet center
            if not zoom:
                zoom = 0.02
            focus_x = beetpos[0]
            focus_y = beetpos[1]

        ax.scatter(beetpos[0] - focus_x,
                    beetpos[1] - focus_y,
                    c = 'red', zorder=2,
                    s=2)
        ax.scatter(sunpos[0] - focus_x,
                    sunpos[1] - focus_y,
                    c = 'gold', zorder=2,
                    s=2)
        ax.scatter(cloud_x - focus_x, 
                    cloud_y - focus_y,
                    c = colors[:len(cloud_x)],
                    alpha=alphas, 
                    s = 1)
        
        # Circle around the sun to show detection radius:
        circ = patches.Circle((sunpos[0] - focus_x, sunpos[1] - focus_y), radius=detection_radius.value_in(units.kpc), transform=ax.transData, fill=False, color='purple', linestyle='--')
        ax.add_patch(circ)

        ax.set_xlim(-zoom, zoom)
        ax.set_ylim(-zoom, zoom)

        plt.savefig(f'{outdir}/figs/fig_{fignum:03d}.png')
        plt.close()
        fignum +=1
        bar_x_out_of_y(fignum, num_plots)


def make_movie(run_dir: str='../runs/temp') -> None:
    """
    Stitches together all of the figures contained within the given directory.

    Assumes the figures are contained within <run_dir>/figs and named fig_000.png,
    fig_001.png, fig_002.png, etc.

    """
    command = f"ffmpeg -y -framerate 25 -i {run_dir}/figs/fig_%3d.png \
                -c:v libx264 -hide_banner -vb 20M -loglevel panic \
                -pix_fmt yuv420p -filter:v 'setpts=2*PTS' -y {run_dir}/movie.mp4"
    os.system(command)


if __name__ in '__main__':
    make_plots(mask_z_distance=True, zoom=0.1)
    make_movie()