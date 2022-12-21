from amuse.lab import units, constants, Particles
import numpy as np
import matplotlib.pyplot as plt
from amuse.community.seba.interface import SeBa

stars1 = Particles(7, mass=[18,18.5,19,19.5,20,20.5,21]|units.MSun)
# stars2 = stars1.copy()
# stars3 = stars1.copy()
print(dir(stars1))
seba_z00 = SeBa()
# seba_z01 = SeBa()
# seba_z02 = SeBa()

seba_z00.set_metallicity(0.02)
# seba_z01.set_metallicity(0.04)
# seba_z02.set_metallicity(0.06)

times = np.arange(0,12,.1) | units.Myr

seba_z00.particles.add_particles(stars1)
# seba_z01.particles.add_particles(stars2)
# seba_z02.particles.add_particles(stars3)

masses_z00 = [] 
# masses_z01 = []
# masses_z02 = [] 
for time in times:
    seba_z00.evolve_model(time)
    # seba_z01.evolve_model(time)
    # seba_z02.evolve_model(time)
    masses_z00.append(seba_z00.particles.mass.value_in(units.MSun))
    
    # masses_z01.append(seba_z01.particles.mass.value_in(units.MSun))
    # masses_z02.append(seba_z02.particles.mass.value_in(units.MSun))

# for m00,m01,m02 in zip(np.array(masses_z00).T, np.array(masses_z01).T, np.array(masses_z02).T):
for m in masses_z00:
    plt.plot(times.value_in(units.Myr), m, c='blue')
    
    # plt.plot(times.value_in(units.Myr), m01, c='orange')
    # plt.plot(times.value_in(units.Myr), m02, c='red')
plt.show()