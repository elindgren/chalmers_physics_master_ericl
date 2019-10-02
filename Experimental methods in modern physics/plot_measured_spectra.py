import numpy as np
import matplotlib.pyplot as plt
from load_data import load_file

# Load the data
filename = 'Group_1c/Ag_200nm_10pt'
lamp_spectrum_file = 'Group_1c/CRS_700nm'

# Size of sensor is 1024x256



data, nbr_particles = load_file(filename)
lamp_data, s = load_file(lamp_spectrum_file)


# Plot the data
fig, ax = plt.subplots()
wvl = data["lambda"]
background = data["spectra_0"]
lamp_spectra = lamp_data["spectra_0"]

for i in range(nbr_particles):
    spectra = data[f'spectra_{i}']
    if i == 0:
        c = 'k'
        a = 0.7
        l = "Background"
    elif i == 1 or i==nbr_particles-1:
        # First and last samples
        c = 'r'  # color
        a = 0.6  # alpha
        l = f'Particle {i}'  # Label
    else: 
        # All other samples
        c = 'b'
        a = 0.2
        l = f'Particle {i}'
    corrected_spectra = (spectra - background)/lamp_spectra
    norm_spectra = spectra/spectra.max()
    ax.plot(wvl, corrected_spectra, label=l, color=c, alpha=a)

ax.set_title(f'Corrected spectra for file {filename}')
ax.set_xlabel("Wavelength")
ax.set_ylabel("Counts")
plt.grid()
plt.legend()
plt.savefig(f'{filename}.png')
plt.show()

# TODO Plot shift of peaks as a function of wavelength


# Centroid = center of mass wavelength
# Normalize spectra for a better fit - best around 0 to 1
