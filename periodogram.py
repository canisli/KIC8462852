import sys

import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from astropy.timeseries import TimeSeries

file = sys.argv[1]
ts = TimeSeries.read("./out/" + file + ".csv", format="csv", time_column="time")
frequency, power = LombScargle.from_timeseries(ts, signal_column_name="mag").autopower()

plt.figure()
plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["font.size"] = "20"
plt.xlabel("Frequency (cycles per day)")
plt.ylabel("Lomb-Scargle Power")
plt.plot(frequency, power)
plt.show()