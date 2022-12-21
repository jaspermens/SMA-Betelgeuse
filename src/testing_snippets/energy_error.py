import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import constants, units
import pickle as pk
from amuse.lab import Particles

dir(Particles)
energy_error_filename = 'last_run_energy_checks.pk'
energy_data = pk.load(open(energy_error_filename, 'rb'))

times = energy_data['time']
total_energies = [ke + pe for ke,pe in zip(energy_data['kinetic'], energy_data['potential'])]
energy_errors = total_energies/total_energies[0]

# print(times[0], energy_data['potential'][-1])
plt.plot(times, energy_data['kinetic'])
# plt.plot(times, energy_data['potential'])
plt.show()