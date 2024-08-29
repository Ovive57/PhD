# Basic imports
import os
import numpy as np
import random
import time

# Plots
from matplotlib import pyplot as plt
import style
plt.style.use(style.style1)

# Astropy
from astropy.table import Table
from astropy import coordinates as co, units as u
from astropy.io import fits

# Debugging
import ipdb #*ipdb.set_trace()

from dlr_stepbystep import get_sn_angle, get_dDLR, get_dDLR_rotated


seed = int(time.time() * 1000) % (2**32)
np.random.seed(seed)


def dDLR_comparison(RA_sn, Dec_sn, RA_gal, Dec_gal, id_sn, id_gal, ind, radius_deg=2/60):
    sn_ids = []
    galaxy_ids = []
    classic_dDLRs = []
    rotated_dDLRs = []
    times_classic = []
    times_rotated = []

    RA_sn = RA_sn[ind]
    Dec_sn = Dec_sn[ind]
    id_sn = id_sn[ind]
    for i, sn in enumerate(id_sn):
        # Selection of my region:
        #ind = random.randint(0, len(ID_sn)-1)

        center_RA = RA_sn[i]
        center_Dec = Dec_sn[i]
        sn_id = id_sn[i]

        #radius_deg = 2/60 #2 arcmin in deg

        max_RA = center_RA + radius_deg/2 # 1 arcmin in deg
        min_RA = center_RA - radius_deg/2
        max_Dec = center_Dec + radius_deg/2
        min_Dec = center_Dec - radius_deg/2

        # Find the data within this region:

        # We have already the SN
        # The galaxies:
        ind = np.where((RA_gal>min_RA)&(RA_gal<max_RA)&(Dec_gal>min_Dec)&(Dec_gal<max_Dec))
        RA_galaxies_region = RA_gama[ind]
        Dec_galaxies_region = Dec_gama[ind]
        gal_id = id_gal[ind]

        sn_ids.extend([sn_id]*len(gal_id))
        galaxy_ids.extend(gal_id)

        # Preparing the data for the dDLR calculation:
        major_galaxies_region = major_gama[ind]
        minor_galaxies_region = minor_gama[ind]
        position_angle_region = np.radians(position_angle_gama[ind]-90)

        for i,pos in enumerate(RA_galaxies_region):
            classic_dDLR, time_classic = get_dDLR(RA_galaxies_region[i], Dec_galaxies_region[i], position_angle_region[i], major_galaxies_region[i], minor_galaxies_region[i], center_RA, center_Dec, plot=False)
            rotated_dDLR,a,b, time_rotated = get_dDLR_rotated(sn_coords, co.SkyCoord(RA_galaxies_region[i]*u.deg, Dec_galaxies_region[i]*u.deg), major_galaxies_region[i], minor_galaxies_region[i], position_angle_region[i], plot=False)

            classic_dDLRs.append(classic_dDLR)
            rotated_dDLRs.append(rotated_dDLR)
            times_classic.append(time_classic)
            times_rotated.append(time_rotated)


    table = Table([sn_ids, galaxy_ids,classic_dDLRs, rotated_dDLRs, times_classic, times_rotated], names=('sn_id', 'gal_id', 'classic_dDLR', 'rotated_dDLR', 'time_classic', 'time_rotated'))
    table.write(f'files/10_SN_2arcmin_region/med_dec.fits', format='fits', overwrite=True)

    dif_dDLR = np.array(classic_dDLRs) - np.array(rotated_dDLRs)
    mean = np.mean(dif_dDLR)
    median = np.median(dif_dDLR)
    std_dev = np.std(dif_dDLR)
    plt.hist(dif_dDLR, bins=50)
    plt.legend([f'Mean: {mean:9f}\nMedian: {median:9f}\nStd Dev: {std_dev:9f}'])

    plt.xlabel('Difference in dDLR')
    plt.ylabel('Frequency')
    plt.title('Histogram of the difference in dDLR between the classic and rotated methods')
    plt.savefig(f'plots/10_SN/med_dec.png')

    plt.show()

#### The DATA ####
galaxies_table = Table.read('data/gkvScienceCatv02.fits')
#SN_table = Table.read('data/tns_SNIa_20240424_copy.csv')
SN_table = Table.read('crossmatches/specifics/colour_mass_matched465.fits')
## GALAXIES
RA_gama = galaxies_table['RAcen']
Dec_gama = galaxies_table['Deccen']
ID_gama = galaxies_table['CATAID']
position_angle_gama = galaxies_table['ang']
axrat = galaxies_table['axrat']
major_gama = galaxies_table['R50'] # In arcsec
minor_gama = major_gama*axrat

## SN
RA_sn = SN_table['ra']
Dec_sn = SN_table['declination']
sn_coords = co.SkyCoord(RA_sn*u.deg, Dec_sn*u.deg)
ID_sn = SN_table['objid']

print(f'maximum dec SN: {Dec_sn.max()}')
ind_max = np.where(Dec_sn > Dec_sn.max() - 1.2)
print(np.shape(ind_max))
print(f'minimum dec SN: {Dec_sn.min()}')
ind_min = np.where(Dec_sn < Dec_sn.min() + 1.2)
print(np.shape(ind_min))
ind_med = np.where((Dec_sn > Dec_sn.min() + 1.8) & (Dec_sn < Dec_sn.max() - 1.7))
print(np.shape(ind_med))


dDLR_comparison(RA_sn, Dec_sn, RA_gama, Dec_gama, ID_sn, ID_gama, ind_med)