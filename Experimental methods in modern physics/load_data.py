import numpy as np


def load_file(filename):
    d = {}
    with open(filename + ".asc") as f:
        for idx, line in enumerate(f):
            val = line.split()
            # If first line, define dict arrays
            if idx == 0:
                nbr_samples = len(val[1:])  # The number of sampels
                d["lambda"] = []  
            # The first row is wavelength
            d["lambda"].append(float(val[0]))
            # Loop through the rest of the lines
            for i, v in enumerate(val[1:]):
                # If first line, define dict arrays
                if idx == 0:
                    d[f'spectra_{i}'] = []
                d[f'spectra_{i}'].append(float(v))
    # Convert all arrays to numpy arrays
    for key, val in d.items():
        d[key] = np.array(val)
    return d, nbr_samples