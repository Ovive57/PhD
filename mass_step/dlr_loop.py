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
star_mask = galaxies_table['starmask']
duplicate_mask = galaxies_table['duplicate']
mask = np.where((star_mask == False) & (duplicate_mask == False))
#galaxies_table = galaxies_table[mask]
RA_gama = galaxies_table['RAcen']
Dec_gama = galaxies_table['Deccen']
ID_gama = galaxies_table['CATAID']
Uber_gama = galaxies_table['uberID']
Z_gama = galaxies_table['Z']
position_angle_gama = galaxies_table['ang']
axrat = galaxies_table['axrat']
major_gama = galaxies_table['R100'] # In arcsec
minor_gama = major_gama*axrat

## SN
RA_sn = SN_table['ra']
Dec_sn = SN_table['declination']
sn_coords = co.SkyCoord(RA_sn*u.deg, Dec_sn*u.deg)
ID_sn = SN_table['objid']
sn_name = SN_table['fullname']
z_sn = SN_table['redshift']

num_sn = 10

sn_ids = []
full_sn_names = []
galaxy_ids = []
galaxy_ids_uber = []
classic_dDLRs = []
rotated_dDLRs = []
times_classic = []
times_rotated = []
separations_classic = []
separations_rotated = []

RA_host_galaxies = []
Dec_host_galaxies = []
Z_host_galaxies = []
RA_sn_cross = []
Dec_sn_cross = []
Z_sn_cross = []


for i, iden_sn in enumerate(ID_sn):
    # Selection of my region:
    #ind = random.randint(0, len(ID_sn)-1)

    center_RA = RA_sn[i]
    center_Dec = Dec_sn[i]
    sn_id = ID_sn[i]
    name_sn = sn_name[i]
    center_z = z_sn[i]

    radius_deg = 2/60 #1 arcmin in deg

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
    Z_galaxies_region = Z_gama[ind]
    gal_id = ID_gama[ind]
    gal_id_uber = Uber_gama[ind]

    #sn_ids.extend([sn_id]*len(gal_id))
    #galaxy_ids.extend(gal_id)

    # Preparing the data for the dDLR calculation:
    major_galaxies_region = major_gama[ind]
    minor_galaxies_region = minor_gama[ind]
    position_angle_region = np.radians(position_angle_gama[ind]-90)

    ###### FINDING THE HOST GALAXY ######
    ind_host_galaxy, [dDLR_rotated,time_rotated,sep_rotated], [dDLR_classic, time_classic,sep_classic]= find_host_galaxy(RA_galaxies_region, Dec_galaxies_region, major_galaxies_region, minor_galaxies_region, position_angle_region, center_RA, center_Dec)
    if ind_host_galaxy == None:
        continue


    # The host galaxy
    RA_host_galaxy = RA_galaxies_region[ind_host_galaxy]
    Dec_host_galaxy = Dec_galaxies_region[ind_host_galaxy]
    Z_host_galaxy = Z_galaxies_region[ind_host_galaxy]

    major_host_galaxy = major_galaxies_region[ind_host_galaxy]
    minor_host_galaxy = minor_galaxies_region[ind_host_galaxy]
    position_angle_host = position_angle_region[ind_host_galaxy]
    ID_host = gal_id[ind_host_galaxy]
    ID_host_uber = gal_id_uber[ind_host_galaxy]

    sn_ids.append(sn_id)
    full_sn_names.append(name_sn)
    galaxy_ids.append(ID_host)
    galaxy_ids_uber.append(ID_host_uber)
    rotated_dDLRs.append(dDLR_rotated)
    classic_dDLRs.append(dDLR_classic)
    times_rotated.append(time_rotated)
    times_classic.append(time_classic)
    separations_rotated.append(sep_rotated)
    separations_classic.append(sep_classic)
    RA_host_galaxies.append(RA_host_galaxy)
    Dec_host_galaxies.append(Dec_host_galaxy)
    Z_host_galaxies.append(Z_host_galaxy)
    RA_sn_cross.append(center_RA)
    Dec_sn_cross.append(center_Dec)
    Z_sn_cross.append(center_z)

print(f'Number of galaxies: {len(galaxy_ids)}')
print(f'Number of uber galaxies: {len(galaxy_ids_uber)}')
print(f'Number of SNe: {len(sn_ids)}')
print(f'Number of full SN names: {len(full_sn_names)}')
print(f'Number of classic dDLRs: {len(classic_dDLRs)}')
print(f'Number of rotated dDLRs: {len(rotated_dDLRs)}')
print(f'Number of classic times: {len(times_classic)}')
print(f'Number of rotated times: {len(times_rotated)}')
print(f'Number of classic separations: {len(separations_classic)}')
print(f'Number of rotated separations: {len(separations_rotated)}')

table = Table([sn_ids, full_sn_names, galaxy_ids, galaxy_ids_uber, classic_dDLRs, rotated_dDLRs, times_classic, times_rotated, separations_rotated, separations_classic, RA_host_galaxies, Dec_host_galaxies, Z_host_galaxies, RA_sn_cross, Dec_sn_cross, Z_sn_cross], names=('sn_id','sn_fullname', 'CATAID','uberID', 'classic_dDLR', 'rotated_dDLR', 'time_classic', 'time_rotated', 'sep_rotated', 'sep_classic', 'RA_host', 'Dec_host', 'Z_host', 'RA_sn', 'Dec_sn', 'Z_sn'))
table.write(f'files/SN_ZTF/7041_SN_R100_stars.fits', format='fits', overwrite=True)
end = time.time()

print(f'Total time: {end-start} s')
print(f'Total time classic: {np.sum(times_classic)} s')
print(f'Total time rotated: {np.sum(times_rotated)} s')