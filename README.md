# KIC 8462852 (Tabby's Star)
KIC 8462852 A is a variable star featuring irregular light fluctuations, including a dimming of 22% in brightness. The cause is not fully understood, with several hypotheses ranging from an uneven ring of dust orbiting the star to an artificial megastructure orbiting the star. Recent research has also demonstrated that Tabby’s star is extremely likely to be a binary star system.

---
Libraries used: astropy, photutils, numpy, matplotlib

`ap.py` - processes fits files and plots a light curve using aperture photometry, along with writing the processed data to a csv file

`periodogram.py` - generates a Lomb-Scargle periodogram from a csv file containing time and magnitude columns

`plot_data.py` - plots light curve from the csv file created by `ap.py`

`plot_aavso.py` - plots light curve from AAVSO data

`FilterMagnitudes.py` - contains magnitude of the calibration stars, HD 191224 and TYC 3162-879-1, under specific filters

`MoveBadFiles.py` - quick utility to move erroneous files out of the data folder
