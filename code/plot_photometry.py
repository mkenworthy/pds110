'''
Copyright 2017, Matthew A. Kenworthy

plot_photometry.py - pretty plots for the PDS 110 Monitoring slack channel

'''

import numpy as np
from astropy.io import ascii
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
#import seaborn as sns   

aavso_file = '../data/aavsodata.txt'

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
n_bands = len(ta_by_band.groups)

conv = {'HKEB':'s',
        'HMB':'o',
        'MGW':'v',
        'DKS':'<',
        'HJW':'.',
        'PCC':'^',
        'LCLC':'*',
        'MGAB':'>'}

obscol = {'HKEB':'r',
        'HMB':'g',
        'LCLC':'darkblue',
        'PCC':'brown',
        'MGW':'b',
        'DKS':'y',
        'HJW':'purple',
        'MGAB':'m'}


tmax = 57998.

mybands = ('B','V','R','I')
n_bands = len(mybands)

fig, axes = plt.subplots(n_bands, 1, figsize=(8, 11), sharex=True, sharey=True)

fig.subplots_adjust(hspace=0.05, wspace=0.05)

# add expected location of eclipse
from astropy.modeling.models import Gaussian1D
#     .. math:: f(x) = A e^{- \frac{\left(x - x_{0}\right)^{2}}{2 \sigma^{2}}}
# 0.5 = exp-x*x / (2*s*s)
# -2 s s ln 0.5 = 1
# s = -1/(2 ln 2) 
# HJD=2458015.5 ± 10 (Sept 9-30 2017
hjd_mid = 2458015.5 - 2400000.5
ecl_fwhm = 7.
ecl_sig = ecl_fwhm / np.sqrt((2*np.log(2)))
f = Gaussian1D(0.26, hjd_mid, ecl_sig)
f_late = Gaussian1D(0.26, hjd_mid+10., ecl_sig)

t_ecl = np.arange(hjd_mid - 50, hjd_mid + 50, 1)


for (ax, band) in zip(axes, mybands):
#for (ax, band) in zip(axes, ta_by_band.groups.keys):
    mask = ta_by_band.groups.keys['Band'] == band[0]

    # loop over each band and estimate out of transit flux
    tb = ta_by_band.groups[mask]
    tb_by_obs = tb.group_by('Observer Code')
    n_obs = len(tb_by_obs.groups)
#    print('{} observers in filter {}'.format(n_obs,band[0]))
    for nob in tb_by_obs.groups.keys:
        mask2 = tb_by_obs.groups.keys['Observer Code'] == nob[0]        
        tc = tb_by_obs.groups[mask2]
        n_points = tc['JD'].size
#        print('In band {} observer {} has {} observations'.format(band[0],nob[0],n_points))

        t_noecl = (tc['MJD'] < tmax)

        # make an out of eclipse average
        t_out = tc[t_noecl]
#        print(t_out)
        mean_mag = np.array(t_out['Magnitude']).mean()
#        print('mean magnitude is {}'.format(mean_mag))

        # mag to intensity
        tc['dmag'] = tc['Magnitude'] - mean_mag
        tc['I'] = np.power(10, tc['dmag'] / -2.5)

        # photometry with errorbars
        ax.errorbar(tc['MJD'], tc['I'], tc['Uncertainty'], fmt=conv.get(nob[0],"*"), color=obscol.get(nob[0],'black'), label=nob[0])

        # photometric band label
        ax.text(0.2, 0.2, band[0], ha='center', va='center', fontsize=24, transform=ax.transAxes)

        # today's date
        ax.vlines(now.mjd, 0.0, 1.1,linestyle='dashed')

        # out of eclipse flux level
        ax.hlines(1.0, 0, 60000, linestyle='dotted')

        # expected location of eclipse
        ax.plot(t_ecl, 1 - f(t_ecl),color='red')

        # expected location of eclipse
        ax.plot(t_ecl, 1 - f_late(t_ecl), color='yellow', alpha=0.5)
        
        ax.legend(loc="lower right", fontsize=8, ncol=3, fancybox=True)
        
ax.set_ylim(0.50,1.08)
ax.set_xlim(now.mjd-20, now.mjd+40)

axes[-1].set_xlabel('Time [MJD]')
fig.suptitle('PDS 110 Photometry', fontsize='large')

fout = datetime.today().strftime('pds110_intens_aavso_%Y%m%d_%H%M.png')
plt.savefig(fout)
plt.draw()
plt.show()
