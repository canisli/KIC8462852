import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

results = Table.read("./out/" + sys.argv[1] + ".csv", format='csv')
name = sys.argv[2]
filter_type = sys.argv[3]

plt.figure()
plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["font.size"] = "20"
plt.errorbar(
    results["JD"] - np.min(results["JD"]),
    results["Magnitude"],
    yerr=0,
    ls="None",
    marker="o",
    markersize=2,
    color="k",
)
plt.xlabel("Time since first observation (days)")
plt.ylabel(filter_type + " Magnitude of " + name)
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.tight_layout()
plt.show()
