#!/usr/bin/env python
# coding: utf-8


# Standard libraries
import os

# Third party libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from scipy.optimize import minimize_scalar



# Psuedo

# Load data
# Fit a polynomial
# Calculate the peak position
# Calculate the FWHM



# Fit a polynomial
wavelength = ''
data = ''
order = 20

# Construct the design 
poly = PolynomialFeatures(order, include_bias=False)
X_poly = poly.fit_transform(wavelength)

# Perform the fit
lin_reg = LinearRegression()
lin_reg.fit(wavelength, data)

# Find maxmimum by optimizing the lin_reg.predict
x0 = ''

# Find the FWHM by finding the roots to lin_reg.predict -lin_reg.predict(x0)/2
gamma = '' /2