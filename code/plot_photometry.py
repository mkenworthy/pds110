'''
Copyright 2017, Matthew A. Kenworthy

plot_photometry.py - pretty plots for the PDS 110 Monitoring slack channel

'''

import numpy as np
from astropy.io import ascii
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import seaborn as sns   


aavso_file = '../data/aavsodata_598d6894a988c.txt'
aavso_file = '../data/aavsodata_59962839569b7.txt'

ta = ascii.read(aavso_file)

from astropy.time import Time

# read in JD into a Time object for easy conversion

ta['MJD'] = ta['JD'] - 2400000.5

times = Time(ta['JD'], format='jd')

from datetime import datetime

now = Time(datetime.utcnow(),format='datetime')
print('current MJD is {}'.format(now.mjd))

# get a list of the unique bandpasses
ta_by_band = ta.group_by('Band')

# the number of different bands
n_bands = ta_by_band.groups.indices.size - 1

fig, axes = plt.subplots(n_bands, 1, figsize=(8, 11), sharex=True)
                         
                         # subplot_kw={'xticks': [], 'yticks': []})

fig.subplots_adjust(hspace=0.05, wspace=0.05)

conv = {'HKEB':'s',
        'HMB':'o',
        'MGW':'v',
        'DKS':'<',
        'MGAB':'>'}

obscol = {'HKEB':'r',
        'HMB':'g',
        'MGW':'b',
        'DKS':'y',
        'MGAB':'m'}


# loop over the bands
for (ax, band) in zip(axes, ta_by_band.groups.keys):
    print(band)
    mask = ta_by_band.groups.keys['Band'] == band[0]

    # loop over each band and estimate out of transit flux
    tb = ta_by_band.groups[mask]

    tb_by_obs = tb.group_by('Observer Code')
    n_obs = tb_by_obs.groups.indices.size - 1
    print(n_obs)

#    ax.errorbar(tb['MJD'], tb['Magnitude'], tb['Uncertainty'], fmt='o')
#    ax.text(0.1, 0.9, band[0], ha='center', va='center', transform=ax.transAxes)

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


tmax = 57974.

fig, axes = plt.subplots(n_bands, 1, figsize=(8, 11), sharex=True, sharey=True)
                         # subplot_kw={'xticks': [], 'yticks': []})

fig.subplots_adjust(hspace=0.05, wspace=0.05)


for (ax, band) in zip(axes, ta_by_band.groups.keys):
    print(band)
    mask = ta_by_band.groups.keys['Band'] == band[0]

    # loop over each band and estimate out of transit flux
    tb = ta_by_band.groups[mask]
    tb_by_obs = tb.group_by('Observer Code')
    n_obs = tb_by_obs.groups.indices.size - 1
    print(n_obs)
    print('{} observers in filter {}'.format(n_obs,band[0]))
    for nob in tb_by_obs.groups.keys:
        mask2 = tb_by_obs.groups.keys['Observer Code'] == nob[0]        
        tc = tb_by_obs.groups[mask2]
        n_points = tc['JD'].size
        print('In band {} observer {} has {} observations'.format(band[0],nob[0],n_points))



        t_noecl = (tc['MJD'] < tmax)

        # make an out of eclipse average
        t_out = tc[t_noecl]
        print(t_out)
        mean_mag = np.array(t_out['Magnitude']).mean()
        print('mean magnitude is {}'.format(mean_mag))

        # mag to intensity
        tc['dmag'] = tc['Magnitude'] - mean_mag
        tc['I'] = np.power(10, tc['dmag'] / -2.5)

#        ax.errorbar(tc['MJD'], tc['Magnitude'], tc['Uncertainty'], fmt=conv[nob[0]], color=obscol[nob[0]])
        ax.errorbar(tc['MJD'], tc['I'], tc['Uncertainty'], fmt=conv[nob[0]], color=obscol[nob[0]])
        ax.text(0.1, 0.9, band[0], ha='center', va='center', transform=ax.transAxes)
        ax.vlines(now.mjd, 0.9, 1.1,linestyle='dashed')
        ax.hlines(1.0, 0,60000,linestyle='dotted')

ax.set_ylim(0.88,1.12)
ax.set_xlim(now.mjd-20, now.mjd+1)

axes[-1].set_xlabel('Time [MJD]')
fig.suptitle('PDS 110 Photometry', fontsize='large')

fout = datetime.today().strftime('pds110_intens_aavso_%Y%m%d_%H%M.png')
plt.savefig(fout)
plt.draw()
plt.show()
