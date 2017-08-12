'''
Copyright 2017, Matthew A. Kenworthy

plot_photometry.py - pretty plots for the PDS 110 Monitoring slack channel

'''

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import seaborn as sns   


aavso_file = 'data/aavsodata_598d6894a988c.txt'

ta = ascii.read(aavso_file)

from astropy.time import Time

# read in JD into a Time object for easy conversion

times = Time(ta['JD'], format='jd')

print(ta)

from datetime import datetime
print(datetime.utcnow())


from time import gmtime, strftime
strftime("%Y-%m-%d_%H:%M:%S", gmtime())

# get a list of the unique bandpasses
ta_by_band = ta.group_by('Band')
print(ta_by_band.groups.keys)


# loop over the bands
for band in ta_by_band.groups.keys:
    print(band)
    mask = ta_by_band.groups.keys['Band'] == band

    # loop over each band and estimate out of transit flux
    tb = ta_by_band.groups[mask]

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.errorbar(tb['JD'], tb['Magnitude'], tb['Uncertainty'])
    ax.text(0.1, 0.9, band, ha='center', va='center', transform=ax.transAxes)

    tmax = 2457975.

    t_noecl = (tb['JD'] < tmax)

    # make an out of eclipse average
    t_out = tb[t_noecl]

    mean_mag = np.mean(t_out['Magnitude'])

    
    
