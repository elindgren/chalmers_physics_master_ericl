import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from data_loader import load_file
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

# Specify files to load
ag_file = 'Ag_200nm_10pt'
au_file = 'Au_200nm_10pt'
pd_file = 'Pd_200nm_10pt'
lamp_file = 'CRS_700nm'

# Load lamp spectra
lamp_data, _ = load_file(lamp_file)
lamp_spectra = lamp_data["spectra_0"]

# Define dict
spectra_dict = {
    "Ag": {
        "file": ag_file,
        "wvl": [],
        "avg_spectra": np.zeros((1024)),
        "peak_pos": 0,
        "color": "slategray"
    },
    "Au": {
        "file": au_file,
        "wvl": [],
        "avg_spectra": np.zeros((1024)),
        "peak_pos": 0,
        "color": "gold"
    },
    "Pd": {
        "file": pd_file,
        "wvl": [],
        "avg_spectra": np.zeros((1024)),
        "peak_pos": 0,
        "color": "dimgray"
    }
}

# Define figure
fig, ax = plt.subplots(figsize=(10,6))

for sample in spectra_dict:
    # Load data from file
    sample_data, _ = load_file(spectra_dict[sample]["file"])
    # Pick out wavelength
    wvl = sample_data["lambda"]
    # Pick out background spectrum
    sample_bground = sample_data["spectra_0"]
    # Calculate average spectra
    measurements = 0
    avg_spectra = spectra_dict[sample]["avg_spectra"]
    for measurement in sample_data:
        if not (measurement == "spectra_0") and not (measurement=="lambda"):
            # Skip background measurement
            avg_spectra += sample_data[measurement]
            measurements += 1
    print(f'Number of measurements: {measurements} for sample: {sample}')
    avg_spectra /= measurements
    # Renormalize, by removing background and then dividing by lamp_spectra
    avg_spectra = (avg_spectra-sample_bground)/lamp_spectra
    # calc peak_pos and fwhm
    peak_pos, fwhm = lorentzian_fit(wvl, avg_spectra, 720, False)
    print(peak_pos)
    # Save results
    spectra_dict[sample]["wvl"] = wvl
    spectra_dict[sample]["avg_spectra"] = avg_spectra
    spectra_dict[sample]["peak_pos"] = peak_pos


# Iterate through and plot all spectras. Also append all peak pos together to a vector and plot that
peak_vector = []
for i, sample in enumerate(spectra_dict):
    wvl = spectra_dict[sample]["wvl"]
    avg_spectra = spectra_dict[sample]["avg_spectra"]
    peak_vector.append(spectra_dict[sample]["peak_pos"])
    ax.plot(wvl, avg_spectra, color=spectra_dict[sample]["color"], label=sample, linewidth=2)
print(peak_vector)

#ax.plot(peak_vector, list(range(len(spectra_dict))), np.zeros(len(spectra_dict)), 'k--', label="Peak position")

#ax.tick_params(axis='both', which='major', pad=0)
ax.set_xlabel("Wavelength")
ax.set_ylabel("Counts")
plt.grid()
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig("ag_au_pd_comparison.png")
tikzplotlib.save('ag_au_pd_comparison.tex')

fig, ax = plt.subplots()
ax.plot(peak_vector, 'k*')
ax.set_ylabel("Peak position")
ax.set_xlabel("Sample (Ag Au Pd)")
plt.savefig("ag_au_pd_peak_pos.png")
tikzplotlib.save('ag_au_pd_peak_pos.tex')

plt.show()

