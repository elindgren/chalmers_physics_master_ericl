import numpy as np
import matplotlib.pyplot as plt
from data_loader import load_file, read_gas_file
import pandas as pd
# Standard libraries
import os

# Third party libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from scipy.optimize import minimize_scalar


# Files
filename = 'Pd_diskrod_200nm_15pt_pulse2'
lamp_file = 'CRS_700nm_glas'
gas_file = 'Lab_TIF295_ArH2_pulses_Pd_2'

# Load lamp spectra
lamp_data, _  = load_file(lamp_file)
lamp_spectra = lamp_data["spectra_0"]

# Load gas data - ndarray [t,h2]
gas_data = read_gas_file(gas_file)
fig_gas, ax_gas = plt.subplots()
ax_gas.plot(gas_data[0], gas_data[1])
# Load data
data = pd.read_table("Group_1c/" + filename + '.asc', header=None)
print(data)
# Drop last column with NaNs
data = data.drop(columns=12)
# Split the data into 120 equal parts
nbr_spectra = 120
first_wavelength = data.iloc[0,0]

# Split 
wvl_rep = data.index[data.iloc[:,0]==first_wavelength]
# Split the repeated measurements
measurements = [data[wvl_rep[i]:wvl_rep[i+1]] for i in range(len(wvl_rep)-1)]  # TODO Add last measurement as well
# Define array to hold max values
peak_vals = {}
# Define figure
fig, ax = plt.subplots()
# Loop over all measurements at different times t
x = np.linspace(450, 900)
for iter, measurement_df in enumerate(measurements):
    l = measurement_df.iloc[:,0]
    b_ground = measurement_df.iloc[:,1]
    for p in measurement_df:
        if p > 1:
            if iter == 0:
                peak_vals[f'particle_{p}'] = []
            p_spectra_t = measurement_df.iloc[:, p]
            corr_spectra_t = (p_spectra_t-b_ground)/l
            peak_pos, peak_val, gamma = lorentzian_fit(l, corr_spectra_t)
            #lor = gen_lorentzian(x, peak_pos, peak_val, gamma)
            #ax.plot(l, corr_spectra_t, label=f'particle_{p}')
            #ax.plot(x, lor, 'r', label=f'fit')
            # Append peak
            peak_vals[f'particle_{p}'].append(peak_pos)
for key, peaks in peak_vals.items():
    ax.plot(peaks, label=key)
plt.grid()
plt.legend()
plt.show()
