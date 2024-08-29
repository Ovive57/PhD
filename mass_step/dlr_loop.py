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

from dlr_stepbystep import get_sn_angle, get_dDLR, get_dDLR_rotated, find_host_galaxy

start = time.time()
seed = int(time.time() * 1000) % (2**32)
np.random.seed(seed)

#### The DATA ####
galaxies_table = Table.read('data/gkvScienceCatv02.fits')
SN_table = Table.read('data/tns_SNIa_20240424_copy.csv')
#SN_table = Table.read('crossmatches/specifics/colour_mass_matched465.fits')
## GALAXIES
RA_gama = galaxies_table['RAcen']
Dec_gama = galaxies_table['Deccen']
ID_gama = galaxies_table['CATAID']
Uber_gama = galaxies_table['uberID']
position_angle_gama = galaxies_table['ang']
axrat = galaxies_table['axrat']
major_gama = galaxies_table['R50'] # In arcsec
minor_gama = major_gama*axrat

## SN
RA_sn = SN_table['ra']
Dec_sn = SN_table['declination']
sn_coords = co.SkyCoord(RA_sn*u.deg, Dec_sn*u.deg)
ID_sn = SN_table['objid']

num_sn = 10

sn_ids = []
galaxy_ids = []
galaxy_ids_uber = []
classic_dDLRs = []
rotated_dDLRs = []
times_classic = []
times_rotated = []

for i, iden_sn in enumerate(ID_sn):
    # Selection of my region:
    #ind = random.randint(0, len(ID_sn)-1)

    center_RA = RA_sn[i]
    center_Dec = Dec_sn[i]
    sn_id = ID_sn[i]

    radius_deg = 1/60 #1 arcmin in deg

    max_RA = center_RA + radius_deg/2 # 0.5 arcmin in deg = 30 arcsec
    min_RA = center_RA - radius_deg/2
    max_Dec = center_Dec + radius_deg/2
    min_Dec = center_Dec - radius_deg/2

    # Find the data within this region:

    # We have already the SN
    # The galaxies:
    ind = np.where((RA_gama>min_RA)&(RA_gama<max_RA)&(Dec_gama>min_Dec)&(Dec_gama<max_Dec))
    RA_galaxies_region = RA_gama[ind]
    Dec_galaxies_region = Dec_gama[ind]
    gal_id = ID_gama[ind]
    gal_id_uber = Uber_gama[ind]

    #sn_ids.extend([sn_id]*len(gal_id))
    #galaxy_ids.extend(gal_id)

    # Preparing the data for the dDLR calculation:
    major_galaxies_region = major_gama[ind]
    minor_galaxies_region = minor_gama[ind]
    position_angle_region = np.radians(position_angle_gama[ind]-90)

    ###### FINDING THE HOST GALAXY ######
    ind_host_galaxy, [dDLR_rotated,time_rotated], [dDLR_classic, time_classic]= find_host_galaxy(RA_galaxies_region, Dec_galaxies_region, major_galaxies_region, minor_galaxies_region, position_angle_region, center_RA, center_Dec)
    if ind_host_galaxy == None:
        continue


    # The host galaxy
    RA_host_galaxy = RA_galaxies_region[ind_host_galaxy]
    Dec_host_galaxy = Dec_galaxies_region[ind_host_galaxy]

    major_host_galaxy = major_galaxies_region[ind_host_galaxy]
    minor_host_galaxy = minor_galaxies_region[ind_host_galaxy]
    position_angle_host = position_angle_region[ind_host_galaxy]
    ID_host = gal_id[ind_host_galaxy]
    ID_host_uber = gal_id_uber[ind_host_galaxy]

    sn_ids.append(sn_id)
    galaxy_ids.append(ID_host)
    galaxy_ids_uber.append(ID_host_uber)
    rotated_dDLRs.append(dDLR_rotated)
    classic_dDLRs.append(dDLR_classic)
    times_rotated.append(time_rotated)
    times_classic.append(time_classic)

print(f'Number of galaxies: {len(galaxy_ids)}')
print(f'Number of uber galaxies: {len(galaxy_ids_uber)}')
print(f'Number of SNe: {len(sn_ids)}')
print(f'Number of classic dDLRs: {len(classic_dDLRs)}')
print(f'Number of rotated dDLRs: {len(rotated_dDLRs)}')
print(f'Number of classic times: {len(times_classic)}')
print(f'Number of rotated times: {len(times_rotated)}')

table = Table([sn_ids, galaxy_ids, galaxy_ids_uber, classic_dDLRs, rotated_dDLRs, times_classic, times_rotated], names=('sn_id', 'CATAID','uberID', 'classic_dDLR', 'rotated_dDLR', 'time_classic', 'time_rotated'))
table.write(f'files/SN_ZTF/7041_SN.fits', format='fits', overwrite=True)
end = time.time()

print(f'Total time: {end-start} s')
print(f'Total time classic: {np.sum(times_classic)} s')
print(f'Total time rotated: {np.sum(times_rotated)} s')