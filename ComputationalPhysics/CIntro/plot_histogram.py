import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read data
data = pd.read("histogram.txt")

# Plot data
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(data)
plt.show()