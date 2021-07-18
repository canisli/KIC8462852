"""
@author Canis Li

USAGE:
$ python3 plot_data out_name obs_type
    out_name = name of output file without extension to read starting from ./out
    obs_type = which variation of preset coordinates and values to use
        1 = recent, 2 = archival 2019, 3 = archival 2018
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.time import Time


def mag_to_flux(magnitude):
    return 10 ** (-0.4 * magnitude)


n = len(sys.argv)
if n != 3:
    print("Invalid arguments")
    raise Exception

out_name = sys.argv[1]
obs_type = int(sys.argv[2])

out_path = "./out/" + str(out_name) + ".csv"
comparison_magnitude = 0

if obs_type == 1:  # 7/10/21-7/17/21
    comparison_magnitude = -2.5 * np.log10(
        np.average([mag_to_flux(11.35), mag_to_flux(9.49)])
    )
    # TYC 3162-879-1 (11.35) and HD 191224 (9.49) in V filter
elif obs_type == 2:  # 8/1/19-9/30/19
    comparison_magnitude = -2.5 * np.log10(
        np.average([mag_to_flux(10.880), mag_to_flux(10.5)])
    )
    # TYC 3162-879-1 (10.880) and HD 191224 (10.5) in R filter
elif obs_type == 3:  # 2018
    comparison_magnitude = -2.5 * np.log10(
        np.average([mag_to_flux(10.880), mag_to_flux(10.5)])
    )
    # TYC 3162-879-1 (10.880) and HD 191224 (10.5) in R' filter
else:
    print("Invalid obs_type")
    raise Exception


tabby_results = Table.read("./out/" + str(out_name) + ".csv", format='csv')
tabby_results['time'] = Time(tabby_results['time'])

to_remove = []
print(tabby_results)
threshold = 0.2

# for i in range(len(tabby_results['mag_error'])):
#     # print(tabby_results['mag_error'][i])
#     if tabby_results['mag_error'][i] > threshold:
#        # del tabby_results[i] iterator error
#        to_remove.insert(0, i) #sorted backwards

# for i in to_remove:
#     # print(i,tabby_results['mag_error'][i])
#     del tabby_results[i]

# to_remove = []

# for i in range(len(tabby_results['mag'])):
#     # print(tabby_results['mag_error'][i])
#     if tabby_results['mag'][i] > 13-comparison_magnitude or tabby_results['mag'][i] <11-comparison_magnitude:
#        # del tabby_results[i] iterator error
#        to_remove.insert(0, i) #sorted backwards

# for i in to_remove:
#     # print(i,tabby_results['mag_error'][i])
#     del tabby_results[i]

plt.figure()
plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["font.size"] = "16"
plt.errorbar(
    (tabby_results["time"].jd - np.min(tabby_results["time"].jd)),  # *24,
    tabby_results["mag"] + comparison_magnitude,
    yerr=tabby_results["mag_error"],
    ls="None",
    marker="o",
    markersize=2,
    color="k",
)  # connecting the points is not standard in astronomy
plt.xlabel("Time since first observation (days)")
plt.ylabel("Apparent Magnitude of Tabby's Star")
# plt.title('Light Curve of Tabby\'s Star') put in captions of paper
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.tight_layout()
plt.show()
