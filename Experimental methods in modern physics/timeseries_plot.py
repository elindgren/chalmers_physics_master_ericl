import numpy as np
import matplotlib.pyplot as plt
from load_data import load_file
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


def gen_lorentzian(x, x0, max, gamma):
    return max/(1+((x-x0)/gamma)**2)


def find_max(x, intercept, coeff):
    m_val = intercept
    for i, coeff in enumerate(coeff):
        m_val += coeff * x**i
    return -m_val


def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return idx

def lorentzian_fit(wavelength, data):
    # Fit a polynomial
    order = 20
    wl = np.array(wavelength)
    data = np.array(data)
    # Pick out the index between 600 & 800 nm
    idx_b = wl >= 550
    idx_s = wl <= 750
    idx = idx_b*idx_s
    wl = wl[idx]
    data = data[idx]
    # Reshape
    wl = wl.reshape(-1,1)


    # Construct the design 
    poly = PolynomialFeatures(order, include_bias=False)
    X_poly = poly.fit_transform(wl)

    # Perform the fit
    lin_reg = LinearRegression()
    lin_reg.fit(X_poly, data)
    intercept = lin_reg.intercept_
    coeff = lin_reg.coef_

    # Find maxmimum by optimizing the lin_reg.predict
    x0 = minimize_scalar(find_max, method='bounded', args=(intercept, coeff), bounds=(600,800))
    peak_pos = x0.x
    peak_val = data.max()  # TODO do better
    #x0 = 0
    # Find the FWHM by finding the roots to lin_reg.predict -lin_reg.predict(x0)/2
    nearest_wl_idx = find_nearest_idx(wl, peak_pos)
    gamma_idx = find_nearest_idx(data, data[nearest_wl_idx]/2)
    gamma_wl = wl[gamma_idx]
    gamma = (peak_pos-gamma_wl)[0]
    return peak_pos, peak_val, gamma

# Files
filename = 'Group_1c/Pd_diskrod_200nm_15pt_pulse2'
lamp_file = 'Group_1c/CRS_700nm_glas'

lamp_data, _  = load_file(lamp_file)
lamp_spectra = lamp_data["spectra_0"]

# Load data
data = pd.read_table(filename + '.asc', header=None)
print(data)
# Drop last column with NaNs
data = data.drop(columns=17)
# Split the data into 120 equal parts
nbr_spectra = 120
first_wavelength = data.iloc[0,0]

# Split 
wvl_rep = data.index[data.iloc[:,0]==first_wavelength]
# Split the
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

# Extract peaks for each file for each datapoint

# Plot the peaks for each particle as a function of time - 
# plot curve for each particle in the same graph