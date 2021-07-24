import sys

import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from astropy.timeseries import TimeSeries

file = sys.argv[1]
data_type = int(sys.argv[2])
ny = int(sys.argv[3])
if data_type == 1:
    time = "time"
    mag = "mag"
else:
    time = "JD"
    mag = "Magnitude"

ts = TimeSeries.read("./out/" + file + ".csv", format="csv", time_column=time, time_format="jd")
frequency, power = LombScargle.from_timeseries(ts, signal_column_name=mag).autopower(nyquist_factor=ny)
ls = LombScargle.from_timeseries(ts, signal_column_name=mag)
print(ls)
print(ls.false_alarm_probability(power.max()))


plt.figure()

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["font.size"] = "20"
# plt.gca().set_ylim(0,1)
plt.xlabel("Frequency (cycles per day)")
plt.ylabel("Lomb-Scargle Power")
plt.plot(frequency, power)
plt.show()
