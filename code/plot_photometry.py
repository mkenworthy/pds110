'''
Copyright 2017, Matthew A. Kenworthy

plot_photometry.py - pretty plots for the PDS 110 Monitoring slack channel

'''

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import seaborn as sns   


aavso_file = '../data/aavsodata_598d6894a988c.txt'

ta = ascii.read(aavso_file)

from astropy.time import Time

# read in JD into a Time object for easy conversion

ta['MJD'] = ta['JD'] - 2400000.5

times = Time(ta['JD'], format='jd')

from datetime import datetime
print(datetime.utcnow())

from time import gmtime, strftime
strftime("%Y-%m-%d_%H:%M:%S", gmtime())

# get a list of the unique bandpasses
ta_by_band = ta.group_by('Band')

# the number of different bands
n_bands = ta_by_band.groups.indices.size - 1

fig, axes = plt.subplots(n_bands, 1, figsize=(8, 11), sharex=True)
                         
                         # subplot_kw={'xticks': [], 'yticks': []})

fig.subplots_adjust(hspace=0.05, wspace=0.05)

# loop over the bands
for (ax, band) in zip(axes, ta_by_band.groups.keys):
    print(band)
    mask = ta_by_band.groups.keys['Band'] == band[0]

    # loop over each band and estimate out of transit flux
    tb = ta_by_band.groups[mask]

    ax.errorbar(tb['MJD'], tb['Magnitude'], tb['Uncertainty'], fmt='o')
    ax.text(0.1, 0.9, band[0], ha='center', va='center', transform=ax.transAxes)

    tmax = 57970.

    t_noecl = (tb['MJD'] < tmax)

    # make an out of eclipse average
    t_out = tb[t_noecl]

    mean_mag = np.mean(t_out['Magnitude'])

axes[-1].set_xlabel('Time [MJD]')
fig.suptitle('PDS 110 Photometry', fontsize='large')
 
fout = datetime.today().strftime('pds110_aavso_%Y%m%d_%H%M.png')
plt.savefig(fout)
plt.draw()
plt.show()
