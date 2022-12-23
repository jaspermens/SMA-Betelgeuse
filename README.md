# Betelgeuse's expanding Oort Cloud

## Abstract & Contributors

#### Abstract
Betelgeuse (hereafter nicknamed *beet*) is a nearby star in Orion. Beet is going to go supernova. We are interested in the dynamical evolution of its Oort Cloud after the explosion, inside the Galactic Potential. Specifically we want to investigate whether any of its Oort Objects are going to pass through our Solar system

#### Contributors
- Konstantinos Kilmetis (s3745597)
- Jasper Mens (s2015242)
- Rahul Priyadarshan Ravichandran (s3589943)

Together we are known as *Crazy Coincidence*. Our team logo can be found below:\
<img src="https://i.imgur.com/CRZ7iVt.jpg" width="300"> \
 We also want to thank Prof. Alex De Koter for his input.

## Project Specifics
1. Simulate the last few years of Beet's evolution using `SeBa` while utilizing `Hermite` to integrate the orbits of the Oort objects around Beet. Bridge with galactic potential.
2. Resolve the supernova.
3. Integrate the resulting orbits in the Galactic Potential, accounting for the movement of the Solar System (This is going to be our nice movie).

## What we're on the hook for
Listed below you'll find our deliverables. These are notes from the meeting with Simon and Veronica.

#### Minimum
1. Probability of an object arriving to within 1 pc of the Sun.
2. Cumulative distribution of distances and velocities of the Oort objects with respect to the Sun.

There are plenty of interesting ideas that we plan to implement in future.
1. Pruning. Determine from which region of the Oort cloud we're getting objects from, run a simulation on only those objects. Goal is to have better statistics on the distribution.
2. Backpropagation. Try and figure out whether 'Oumuamua could have been an Oort object from some distant supernova. 
3. Orion Arm's most wanted. This simulation should be generalizable to any star in our galactic neighborhood, find out which one poses the greatest threat to Earth.

## Project log
** Week 0:** (Oct 24 - Oct 30) Initial idea, presentation and project proposal\
** Week 1:** (Oct 31 - Nov 6) `beet_test.py` Initializing the problem, scaling the Oort cloud.
** Week 2:** (Nov 7 - Nov 13) `galactic_potential.py` Implementing the galactic potential.
** Week 3:** (Nov 14 - Nov 20) Meeting with Alex de Koter, bridging between galactic potential and beet system. Creation of this very nice `README`.
** Week 4:** (Nov 21 - Nov 27) Stopped using `Hermite`, made a custom class for integrating the orbits around beet analytically. 
** Week 5:** (Nov 28 - Dec 4) Ran some simulations with the custom class and made some movies with simple parameter values.
** Week 6:** (Dec 5 - Dec 11) Modified initial conditions, included inclination angle of Betelgeuse with the Sun, added code for statistics.
** Week 7:** (Dec 12 - Dec 19) Added a robust code for storing and plotting the simulation data. Final presentation happened (and Luigi became everyone's favourite:))
** Week 8:** (Dec 20 - Dec 23) Final fine tuning, preparation of the report, code cleaning.

## Tutorial
The script that runs the main simulation is `analytic_solver.py` in the `/src` directory. There are three functions of interest: `run_simulation()`, `make_plots()` and `make_movie()`.

In the `run_simulation()` function, you can set the values of various parameters governing the simulation. For a simple run, you can set the following:
1. `end_time = 100 | units.Myr`, 
2. `timestep_pre_sn = 1e-3 | units.Myr`
3. `timestep_post_sn = 0.1 | units.Myr` (Speeding up to 1 Myr changes the results very slightly)
4. `timestep_detection = 0.1 | units.Myr` 
(the motivation for choosing these timesteps can be found in our final report) 
5. `n_oort_objects: int = 100`
6. `phi: float = -34.5 * np.pi/180` (this was the value that gave us many detections) 
7. `detection_radius = 1. | units.pc`
8. `plot_interval = 1. | units.Myr` (these plots will later be made into a movie) 

The `run_simulation()` function evolves the system using the above mentioned parameters and dumps the data into a pickle file called `plotdata.pk`. The `make_plot()` function then reads this file and plots the data according to the parameters chosen by the user:
1. The `focus` parameter decides the center of the plots. There are three possible arguments: `sun`, `beet` and `galaxy`. For instance, to see the orbits around the Milky Way center, set  `focus = 'galaxy'`. 
2. The `zoom` parameter is the half-side length in kpc. `zoom = 0.1` sets a view of (0.2 kpc, 0.2 kpc). We suggest 0.1 kpc for beet, 0.01 kpc for sun and 10 kpc for the galaxy

With these plots, `make_movie()` generates a movie based on the frame rate provided in the `run_simulation()` function. The movie is stored in a newly created folder `/figures` created in the working directory. With the above mentioned parameters, the entire process (simulation, plots and movie) should take less than a minute to run.



