# Betelgeuse's expanding Oort Cloud

## Abstract & Contributors

#### Abstract
Betelgeuse (hereafter nicknamed *beet*) is a nearby star in Orion. Beet is going to go supernova. We are interested in the dynamical evolution of its Oort Cloud after the explosion, inside the Galactic Potential. Specifically we want to investigate whether any of its Oort Objects are going to pass through our Solar system

#### Contributors
- Konstantinos Kilmetis (s3745597)
- Jasper Mens (s2015242)
- Rahul Priyadarshan Ravichandran (s3589943)

Together we are known as *Crazy Coincidence*. Our team logo can be found below:
INSERT LOGO HERE
 We also want to thank Prof. Alex De Koter for his input.

## Project Specifics
1. Simulate the last few years of Beet's evolution using `SeBa` while utilizing `Hermite` to integrate the orbits of the Oort objects around Beet. Bridge with Galactic potential.
2. Resolve the Supernova. 
3. Integrate the resulting orbits in the Galactic Potential, accounting for the movement of the Solar system (This is going to be our nice movie).

NOTE: Maybe make a UML diagram here? Would be clearer.

## What we're on the hook for
Listed below you'll find our deliverables. These are notes from the meeting with Simon and Veronica. 

#### Minimum
1. Probability of an object arriving to within 1 parsec of the Sun
2.  Cumulative distribution of distances and velocities of the Oort objects with respect to the sun (Perhaps through a KDE)

There are plenty of interesting ideas, if we're finished on time we intend to implement these three.
#### Extras
1. Pruning. Determine from which region of the Oort Cloud we're getting objects from, run a simulation on only those objects. Goal is to have better statistics on the distributions.
2. Backpropagation. Try and figure out whether ʻOumuamua could have been an Oort object from some distant supernova.
3. Orion Arm's most wanted. This simulation should be generalizable to any star in our galactic neighbourhood, find out which one poses the greatest threat to Earth.
##  Project Log
**Week 0:** (Oct. 24 - Oct. 30) Initial Idea, Presentation and Project Proposal
**Week 1:** (Oct. 31 - Nov. 6) `beet_test.py` Initializing the problem, scaling the Oort Cloud.
**Week 2:** (Nov. 7  - Nov 13) `galactic_potential.py`.
**Week 3:** (Nov. 14 - Nov.  20) Meeting with Alex De Koter, bridging between galactic potential and beet system. Creation of this very nice `README`

```
