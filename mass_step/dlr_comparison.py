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
    uberIDs = []
    classic_dDLRs = []
    rotated_dDLRs = []
    times_classic = []
    times_rotated = []
    ras_sn = []
    decs_sn = []
    decs_gal = []
    ras_gal = []

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
        if len(gal_id) == 0:
            continue

        sn_ids.extend([sn_id]*len(gal_id))
        ras_sn.extend([center_RA]*len(gal_id))
        decs_sn.extend([center_Dec]*len(gal_id))
        galaxy_ids.extend(gal_id)
        uberIDs.extend(uberID[ind])
        decs_gal.extend(Dec_galaxies_region)
        ras_gal.extend(RA_galaxies_region)

        # Preparing the data for the dDLR calculation:
        major_galaxies_region = major_gama[ind]
        minor_galaxies_region = minor_gama[ind]
        position_angle_region = np.radians(position_angle_gama[ind]-90)

        for i,pos in enumerate(RA_galaxies_region):
            classic_dDLR, sep_classic,time_classic = get_dDLR(RA_galaxies_region[i], Dec_galaxies_region[i], position_angle_region[i], major_galaxies_region[i], minor_galaxies_region[i], center_RA, center_Dec, plot=False)
            rotated_dDLR,sep_rotated, time_rotated = get_dDLR_rotated(RA_galaxies_region[i], Dec_galaxies_region[i], center_RA, center_Dec, major_galaxies_region[i], minor_galaxies_region[i], position_angle_region[i], plot=False)

            classic_dDLRs.append(classic_dDLR)
            rotated_dDLRs.append(rotated_dDLR)
            times_classic.append(time_classic)
            times_rotated.append(time_rotated)


    dif_dDLR = np.array(classic_dDLRs) - np.array(rotated_dDLRs)

    table = Table([sn_ids, ras_sn, decs_sn, galaxy_ids,uberIDs, classic_dDLRs, rotated_dDLRs, times_classic, times_rotated, dif_dDLR, ras_gal, decs_gal], names=('sn_id', 'ra_sn', 'dec_sn','gal_id','uberID', 'classic_dDLR', 'rotated_dDLR', 'time_classic', 'time_rotated', 'dif_dDLR', 'ra_gal', 'dec_gal'))
    table.write(f'files/10_SN_2arcmin_region/test.fits', format='fits', overwrite=True)

    mean = np.mean(dif_dDLR)
    median = np.median(dif_dDLR)
    std_dev = np.std(dif_dDLR)
    plt.hist(dif_dDLR, bins=50)
    plt.legend([f'Mean: {mean:9f}\nMedian: {median:9f}\nStd Dev: {std_dev:9f}'])

    plt.xlabel('Difference in dDLR')
    plt.ylabel('Frequency')
    plt.title('Histogram of the difference in dDLR between the classic and rotated methods')
    plt.savefig(f'plots/10_SN/{seed}.png')

    plt.show()
    return dif_dDLR, sn_ids

#### The DATA ####
galaxies_table = Table.read('data/gkvScienceCatv02.fits')
SN_table = Table.read('data/tns_SNIa_20240424_copy.csv')
#SN_table = Table.read('crossmatches/specifics/colour_mass_matched465.fits')
## GALAXIES
RA_gama = galaxies_table['RAcen']
Dec_gama = galaxies_table['Deccen']
ID_gama = galaxies_table['CATAID']
uberID = galaxies_table['uberID']
position_angle_gama = galaxies_table['ang']
axrat = galaxies_table['axrat']
major_gama = galaxies_table['R50'] # In arcsec
minor_gama = major_gama*axrat

## SN
RA_sn = SN_table['ra']
Dec_sn = SN_table['declination']
sn_coords = co.SkyCoord(RA_sn*u.deg, Dec_sn*u.deg)
ID_sn = SN_table['fullname']

print(f'maximum dec SN: {Dec_sn.max()}')
print(f'maximum dec gama: {Dec_gama.max()}')
ind_max = np.where((Dec_sn < Dec_gama.max())&(Dec_sn > Dec_gama.max() - 1.2))
print(np.shape(ind_max))
print(f'minimum dec SN: {Dec_sn.min()}')
print(f'minimum dec gama: {Dec_gama.min()}')
ind_min = np.where((Dec_sn > Dec_sn.min()+25) & (Dec_sn < Dec_sn.min() + 30))
print(np.shape(ind_min))
ind_med = np.where((Dec_sn > Dec_sn.min() + 25) & (Dec_sn < Dec_gama.max() - 1.7))
print(np.shape(ind_med))


#ind_random = random.randint(0, len(ID_sn)-1)
ind_random = np.random.choice(len(ID_sn), 10)

diff_dlr, sn_ids = dDLR_comparison(RA_sn, Dec_sn, RA_gama, Dec_gama, ID_sn, ID_gama, ind_min)



table_diff = Table.read(f'files/10_SN_2arcmin_region/test.fits')
dec_gama = table_diff['dec_gal']
RA_gama = table_diff['ra_gal']
sn_ids = table_diff['sn_id']
gal_ids = table_diff['gal_id']
uberIDs = table_diff['uberID']
diff_dlr = table_diff['dif_dDLR']
ra_sn = table_diff['ra_sn']
dec_sn = table_diff['dec_sn']
classic_dDLR = table_diff['classic_dDLR']
rotated_dDLR = table_diff['rotated_dDLR']

fig,ax1 = plt.subplots()
ax1.scatter(np.abs(diff_dlr), dec_gama)
ax1.set_xlabel('Difference in dDLR')
ax1.set_ylabel('Dec of the galaxy')
plt.show()

fig,ax2 = plt.subplots()
ax2.scatter(np.abs(diff_dlr), np.abs(RA_gama-ra_sn))
ax2.set_xlabel('Difference in dDLR')
ax2.set_ylabel('RA difference')
plt.show()


print(dec_gama.min(), dec_gama.max())
exit()


percent_diff = (diff_dlr/classic_dDLR)*100
for i, id_sn in enumerate(sn_ids):
    if np.abs(diff_dlr[i]) > 0.005:
        print(f' sn_id: {id_sn}, diff_dlr: {diff_dlr[i]}, %diff_dDLR: {percent_diff[i]}, gal_id: {gal_ids[i]}, uberID: {uberIDs[i]}')#, RA_sn: {ra_sn[i]}, Dec_sn: {dec_sn[i]}')
        #print(id_sn, diff_dlr[i], RA_sn[i], Dec_sn[i])
