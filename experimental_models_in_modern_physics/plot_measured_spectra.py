import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from cycler import cycler
from data_loader import load_file, read_gas_file
from fit_function import lorentzian_fit
import tikzplotlib
import seaborn as sns
sns.set_palette(sns.color_palette("husl", 20))


# Set plot params
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

# Load the data
filename = 'Ag_200nm_10pt'
lamp_spectrum_file = 'CRS_700nm'
background_spectra = 0

# Size of sensor is 1024x256
data, nbr_particles = load_file(filename)
lamp_data, s = load_file(lamp_spectrum_file)

# Plot the data
fig, ax = plt.subplots(figsize=(10,6))
ax = fig.gca(projection='3d')
wvl = data["lambda"]
background = data["spectra_" + str(background_spectra)]
lamp_spectra = lamp_data["spectra_0"]

peak_guesses = [540+i*20 for i in range(nbr_particles)]
peak_guesses.reverse()
peak_positions = []
for i in range(0,nbr_particles):
    spectra = data[f'spectra_{i}']
    if not i == background_spectra:
        print(i)
        # Don't plot background
        corrected_spectra = (spectra - background)/lamp_spectra
        #norm_spectra = spectra/spectra.max()

        peak_pos, fwhm = lorentzian_fit(wvl, corrected_spectra, 780, False)
        peak_positions.append(peak_pos)
        delta_y = np.ones(len(wvl))*(i)
        ax.plot(wvl, delta_y, corrected_spectra, label=f'Particle {i+1}', alpha=0.7)



#ax.set_title(f'Corrected spectra for file {filename}')
ax.tick_params(axis='both', which='major', pad=0)
ax.set_xlabel("Wavelength")
ax.set_ylabel("Measurement")
ax.set_zlabel("Counts")
plt.grid()
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig(f'{filename}.png')
tikzplotlib.save(f'{filename}.tex')

# Plot peak position as well
fig, ax = plt.subplots()
ax.plot(list(range(nbr_particles-1)), peak_positions, 'k--', label="Peak position")
ax.set_xlabel("Particle")
ax.set_ylabel("Peak position")
tikzplotlib.save(f'{filename}_peak_pos.tex')

plt.show()

# Centroid = center of mass wavelength
# Normalize spectra for a better fit - best around 0 to 1
