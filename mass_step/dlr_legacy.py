# Basic imports
import os
import numpy as np
import random
import time

#! If specific region:
random.seed(0)

# Plots
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.image as mpimg

import style

# Astropy
from astropy.table import Table
from astropy import coordinates as co, units as u
from astropy.io import fits

# Legacy Survey catalogs
import pyvo

# Legacy survey images
import urllib
from PIL import Image
import requests

# Debugging
import ipdb #*ipdb.set_trace()

# Functions
from dlr_stepbystep import get_sn_angle, get_dDLR, get_dDLR_rotated, draw_ellipses, find_host_galaxy, get_jpeg_cutout, get_lsdr10_cat, get_region

'''
def get_region(ra_gal, dec_gal, major_gal, minor_gal, pos_ang_gal, ra_sn, dec_sn, sn_ids, catalog = 'GAMA', radius_deg = 2/60):
    RA_galaxies_region = []
    while len(RA_galaxies_region) == 0:
        # Selection of my region:
        ind = random.randint(0, len(sn_ids)-1)
        center_RA = ra_sn[ind]
        center_Dec = dec_sn[ind]
        sn_id = sn_ids[ind]

        max_RA = center_RA + radius_deg/2 # 1 arcmin in deg
        min_RA = center_RA - radius_deg/2
        max_Dec = center_Dec + radius_deg/2
        min_Dec = center_Dec - radius_deg/2

        # Find the data within this region:
        ind = np.where((ra_gal>min_RA)&(ra_gal<max_RA)&(dec_gal>min_Dec)&(dec_gal<max_Dec))
        RA_galaxies_region = ra_gal[ind]
        Dec_galaxies_region = dec_gal[ind]
        uberID_galaxies_region = uberID[ind]
        cataid_galaxies_region = cataid[ind]
        #print(uberID_galaxies_region, cataid_galaxies_region)
        major_galaxies_region = major_gal[ind]
        minor_galaxies_region = minor_gal[ind]
        #print(position_angle_gama[ind]-90)
        position_angle_region = np.radians(pos_ang_gal[ind]-90)
        #position_angle_region = [round(position_angle_region[i], 3) for i in range(len(position_angle_region))]
    print(f"\nNumber of galaxies in {catalog}: {len(RA_galaxies_region)}")

    if len(RA_galaxies_region) == 0:
        print("No galaxies in the region")
        exit()

    return ind, center_RA, center_Dec, sn_id, max_RA, min_RA, max_Dec, min_Dec
'''
'''
while len(RA_galaxies_region) == 0:
    # Selection of my region:
    ind = random.randint(0, len(sn_ids)-1)



    center_RA = RA_sn[ind]
    center_Dec = Dec_sn[ind]
    sn_id = sn_ids[ind]

    #print(center_RA, center_Dec, RA_gama.min(), Dec_gama.min(), RA_gama.max(), Dec_gama.max())

    radius_deg = 2/60 #2 arcmin in deg


    max_RA = center_RA + radius_deg/2 # 1 arcmin in deg
    min_RA = center_RA - radius_deg/2
    max_Dec = center_Dec + radius_deg/2
    min_Dec = center_Dec - radius_deg/2

    # Find the data within this region:

    # We have already the SN
    # The galaxies:
    ind = np.where((RA_gama>min_RA)&(RA_gama<max_RA)&(Dec_gama>min_Dec)&(Dec_gama<max_Dec))
    RA_galaxies_region = RA_gama[ind]
    Dec_galaxies_region = Dec_gama[ind]
    uberID_galaxies_region = uberID[ind]
    cataid_galaxies_region = cataid[ind]
    print(uberID_galaxies_region, cataid_galaxies_region)
'''

#### The DATA ####
galaxies_table_GAMA = Table.read('data/gkvScienceCatv02.fits')
SN_table = Table.read('data/tns_SNIa_20240424_copy.csv')
#SN_table = Table.read('crossmatches/specifics/colour_mass_matched465.fits')

## GALAXIES
# Take away the stars
star_mask = galaxies_table_GAMA['starmask']
duplicate_mask = galaxies_table_GAMA['duplicate']
mask = np.where((star_mask == False)&(duplicate_mask == False))
galaxies_table_GAMA = galaxies_table_GAMA[mask]

RA_gama = galaxies_table_GAMA['RAcen']
Dec_gama = galaxies_table_GAMA['Deccen']
uberID = galaxies_table_GAMA['uberID']
cataid = galaxies_table_GAMA['CATAID']
position_angle_gama = galaxies_table_GAMA['ang']
axrat = galaxies_table_GAMA['axrat']
major_gama = galaxies_table_GAMA['R50'] # In arcsec
R100 = False
minor_gama = major_gama*axrat

## SN
RA_sn = SN_table['ra']
Dec_sn = SN_table['declination']
sn_coords = co.SkyCoord(RA_sn*u.deg, Dec_sn*u.deg)
sn_ids = SN_table['fullname']

ind_GAMA = np.where((RA_sn>RA_gama.min())&(RA_sn<RA_gama.max())&(Dec_sn>Dec_gama.min())&(Dec_sn<Dec_gama.max()))
RA_sn = RA_sn[ind_GAMA]
Dec_sn = Dec_sn[ind_GAMA]
sn_ids = sn_ids[ind_GAMA]

#ind_outliers = np.where(sn_ids == 'SN 2021ddh')

#RA_sn = RA_sn[ind_outliers]
#Dec_sn = Dec_sn[ind_outliers]
#sn_ids = sn_ids[ind_outliers]

radius_deg = 2/60 #2 arcmin in deg


# Get the region:
#RA_galaxies_region, Dec_galaxies_region, center_RA, center_Dec, sn_id, major_galaxies_region, minor_galaxies_region, position_angle_region, uberID_galaxies_region, cataid_galaxies_region = get_region(RA_gama, Dec_gama, major_gama, minor_gama, position_angle_gama, RA_sn, Dec_sn, sn_ids,catalog = 'GAMA', radius_deg = radius_deg)
#! if random region:
ind_region, center_RA, center_Dec, sn_id, max_RA, min_RA, max_Dec, min_Dec = get_region(RA_gama, Dec_gama, major_gama, minor_gama, position_angle_gama, RA_sn, Dec_sn, sn_ids,catalog = 'GAMA', radius_deg = radius_deg)


'''
sn_id = sn_ids[0]
ind_region = np.where(sn_ids == sn_id)
center_RA = RA_sn[ind_region]
center_Dec = Dec_sn[ind_region]
max_RA = center_RA + radius_deg/2 # 1 arcmin in deg
min_RA = center_RA - radius_deg/2
max_Dec = center_Dec + radius_deg/2
min_Dec = center_Dec - radius_deg/2
'''
uberID_galaxies_GAMA = uberID[ind_region]
cataid_galaxies_GAMA = cataid[ind_region]
RA_galaxies_GAMA = RA_gama[ind_region]
Dec_galaxies_GAMA = Dec_gama[ind_region]
major_galaxies_GAMA = major_gama[ind_region]
minor_galaxies_GAMA = minor_gama[ind_region]
position_angle_GAMA = np.radians(position_angle_gama[ind_region]-90)
table_GAMA = Table([uberID_galaxies_GAMA, cataid_galaxies_GAMA, RA_galaxies_GAMA, Dec_galaxies_GAMA], names=('uberID','CATAID', 'RA_galaxy', 'Dec_galaxy'))
#print('GAMA: ',RA_galaxies_region, major_galaxies_region)


######## LEGACY SURVEY ########

galaxies_table_legacy = get_lsdr10_cat('files/SN_2arcmin_region/legacy_comparison/'+str(sn_id), center_RA, center_Dec, radius_deg/2, overwrite=True, shape_r=None, excess_factor=None, parallax=None)

# We read the table:
shape_r = galaxies_table_legacy['shape_r'] #arcsec (effective radius or the half-light radius)
shape_e1 = galaxies_table_legacy['shape_e1'] #dimensionless
shape_e2 = galaxies_table_legacy['shape_e2'] #dimensionless
excess_factor = galaxies_table_legacy['gaia_phot_bp_rp_excess_factor']
parallax = galaxies_table_legacy['parallax']
#! Test this with Anais, with excess factor and paralax == None in the function
ind = np.where(excess_factor != 0)
print(excess_factor[ind])
#! ind = np.where(parallax != 0)
#! print(parallax[ind])
exit()

id_legacy = galaxies_table_legacy['ls_id']
RA_galaxies_legacy = galaxies_table_legacy['ra']
Dec_galaxies_legacy = galaxies_table_legacy['dec']


############### FIRST PLOT ################
plt.style.use(style.style1)
fig1, ax1 = plt.subplots()
ax1.set_title('SN and galaxies in a 2 arcmin region')
ax1.set(xlabel = ('RA (deg)'), ylabel = ('Dec (deg)'), xlim=(max_RA, min_RA), ylim=(min_Dec, max_Dec))
ax1.set_aspect('equal')
#ax1.invert_xaxis()

### Plot the image of the galaxies in the region
get_jpeg_cutout(2, 'deleteme.jpeg', center_RA, center_Dec)
img = mpimg.imread('deleteme.jpeg')
#img = np.fliplr(img)

plt.imshow(img, extent=[max_RA, min_RA, max_Dec, min_Dec])

### Plot the center of the galaxies over the image for GAMA and legacy

ax1.scatter(RA_galaxies_GAMA, Dec_galaxies_GAMA, label = 'GAMA', color = 'blue', s = 10)
ax1.scatter(RA_galaxies_legacy, Dec_galaxies_legacy, label = 'legacy', color = 'red', s = 10)


############# Parameters calculation #############
ellipticity = np.sqrt(shape_e1**2 + shape_e2**2)
position_angle = (np.arctan2(shape_e2, shape_e1))/2
#print('position angle of the galaxies in deg legacy survey ',np.degrees(position_angle))
position_angle_region = position_angle - np.pi/2

axrat = (1-ellipticity)/(1+ellipticity)

major_galaxies_legacy = shape_r
minor_galaxies_legacy = shape_r*axrat

#print(galaxies_table_legacy)

### We define the valid parameters to plot the ellipses:
ind = np.where(shape_r != 0)

#shape_r_valid = shape_r[ind]
#shape_e1_valid = shape_e1[ind]
#shape_e2_valid = shape_e2[ind]

RA_galaxies_legacy_valid = RA_galaxies_legacy[ind]
Dec_galaxies_legacy_valid = Dec_galaxies_legacy[ind]
id_legacy_valid = id_legacy[ind]
position_angle_region_valid = position_angle_region[ind]
major_galaxies_legacy_valid = major_galaxies_legacy[ind]
minor_galaxies_legacy_valid = minor_galaxies_legacy[ind]

table_legacy = Table([id_legacy_valid, RA_galaxies_legacy_valid, Dec_galaxies_legacy_valid], names=('id_ls', 'RA_galaxy', 'Dec_galaxy'))

# The two concentric ellipses for all galaxies in the region:

#!########################################
#! TO TEST
#!position_angle_region = [np.radians(45) for i in range(len(position_angle_region))]
#!print(np.degrees(position_angle_region))
#!########################################
#print('legacy: ', RA_galaxies_region, major_galaxies_region)

for i,pos in enumerate(RA_galaxies_legacy_valid):
    #print(pos)
    ax1 = draw_ellipses([1, 4], major_galaxies_legacy_valid[i], minor_galaxies_legacy_valid[i], position_angle_region_valid[i], RA_galaxies_legacy_valid[i], Dec_galaxies_legacy_valid[i], ax1, rot_type='anticlockwise', colours = ['red', 'red'])

for i,pos in enumerate(RA_galaxies_GAMA):
    #print(pos)
    ax1 = draw_ellipses([1, 4], major_galaxies_GAMA[i], minor_galaxies_GAMA[i], position_angle_GAMA[i], RA_galaxies_GAMA[i], Dec_galaxies_GAMA[i], ax1, rot_type='anticlockwise', colours=['blue', 'blue'])


angles_from_major = []
classic_dDLRs = []
rotated_dDLRs = []
times_classic = []
times_rotated = []
if len(RA_galaxies_legacy_valid)!=0:
    for i,pos in enumerate(RA_galaxies_legacy_valid):
        angle = get_sn_angle(RA_galaxies_legacy_valid[i], Dec_galaxies_legacy_valid[i], major_galaxies_legacy_valid[i], position_angle_region_valid[i], center_RA, center_Dec,plot=True, ax=ax1)
        classic_dDLR, sep_classic, time_classic = get_dDLR(RA_galaxies_legacy_valid[i], Dec_galaxies_legacy_valid[i], position_angle_region_valid[i], major_galaxies_legacy_valid[i], minor_galaxies_legacy_valid[i], center_RA, center_Dec, plot=False)
        rotated_dDLR, sep_rotated, time_rotated = get_dDLR_rotated(RA_galaxies_legacy_valid[i], Dec_galaxies_legacy_valid[i], center_RA, center_Dec, major_galaxies_legacy_valid[i], minor_galaxies_legacy_valid[i], position_angle_region_valid[i], plot=False)
        angles_from_major.append(angle)
        classic_dDLRs.append(classic_dDLR)
        rotated_dDLRs.append(rotated_dDLR)
        times_classic.append(time_classic)
        times_rotated.append(time_rotated)

    #table['classic_dDLR'] = classic_dDLRs
    #table['rotated_dDLR'] = rotated_dDLRs
    #table['time_classic'] = times_classic
    #table['time_rotated'] = times_rotated

    table_legacy['classic_dDLR'] = classic_dDLRs
    table_legacy['rotated_dDLR'] = rotated_dDLRs
    table_legacy['time_classic'] = times_classic
    table_legacy['time_rotated'] = times_rotated

    if R100:
        #table.write(f'files/SN_2arcmin_region/R100/{sn_id}.fits', overwrite=True)
        table_legacy.write(f'files/SN_2arcmin_region/legacy/{sn_id}.fits', overwrite=True)
    else:
        #table.write(f'files/SN_2arcmin_region/{sn_id}.fits', overwrite=True)
        table_legacy.write(f'files/SN_2arcmin_region/legacy/{sn_id}.fits', overwrite=True)
else:
    print('\nNo legacy galaxies in the region\n')

###### FINDING THE HOST GALAXY ######
#! Maybe add this if in the function later
if len(RA_galaxies_legacy_valid)!=0:
    ind_host_galaxy_legacy, dDLR_rotated_legacy, dDLR_classic_legacy= find_host_galaxy(RA_galaxies_legacy_valid, Dec_galaxies_legacy_valid, major_galaxies_legacy_valid, minor_galaxies_legacy_valid, position_angle_region_valid, center_RA, center_Dec,plot = True, catalogue = 'legacy')

else:
    ind_host_galaxy_legacy = None
    dDLR_rotated_legacy = 999
    dDLR_classic_legacy = 999
ind_host_galaxy_GAMA, dDLR_rotated_GAMA, dDLR_classic_GAMA= find_host_galaxy(RA_galaxies_GAMA, Dec_galaxies_GAMA, major_galaxies_GAMA, minor_galaxies_GAMA, position_angle_GAMA, center_RA, center_Dec,plot = True, catalogue='GAMA')

print('\nThe index of the host galaxies in legacy and GAMA respectively: ')
print(ind_host_galaxy_legacy, ind_host_galaxy_GAMA, '\n')

if ind_host_galaxy_legacy==None:
    RA_host_galaxy_legacy = None
    Dec_host_galaxy_legacy = None
    major_host_galaxy_legacy = None
    minor_host_galaxy_legacy = None
    position_angle_host_legacy = None
else:
    # The host galaxy legacy
    RA_host_galaxy_legacy = RA_galaxies_legacy_valid[ind_host_galaxy_legacy]
    Dec_host_galaxy_legacy = Dec_galaxies_legacy_valid[ind_host_galaxy_legacy]
    major_host_galaxy_legacy = major_galaxies_legacy_valid[ind_host_galaxy_legacy]
    minor_host_galaxy_legacy = minor_galaxies_legacy_valid[ind_host_galaxy_legacy]
    position_angle_host_legacy = position_angle_region_valid[ind_host_galaxy_legacy]

if ind_host_galaxy_GAMA == None:
    RA_host_galaxy_GAMA = None
    Dec_host_galaxy_GAMA = None
    major_host_galaxy_GAMA = None
    minor_host_galaxy_GAMA = None
    position_angle_host_GAMA = None
else:
    RA_host_galaxy_GAMA = RA_galaxies_GAMA[ind_host_galaxy_GAMA]
    Dec_host_galaxy_GAMA = Dec_galaxies_GAMA[ind_host_galaxy_GAMA]
    #cataid_host_galaxy = cataid_galaxies_region[ind_host_galaxy]
    #uberID_host_galaxy = uberID_galaxies_region[ind_host_galaxy]
    major_host_galaxy_GAMA = major_galaxies_GAMA[ind_host_galaxy_GAMA]
    minor_host_galaxy_GAMA = minor_galaxies_GAMA[ind_host_galaxy_GAMA]
    position_angle_host_GAMA = position_angle_GAMA[ind_host_galaxy_GAMA]


print('\nthe major axis of the host galaxy of legacy and GAMA respectively: ')
print(major_host_galaxy_legacy, major_host_galaxy_GAMA)

print('\nthe position angle of the host galaxy of legacy and GAMA respectively: ')
print(position_angle_host_legacy, position_angle_host_GAMA)

# PLOT the host galaxy in the same plot as the rest of the galaxies
ax1.scatter(center_RA, center_Dec, label = f'{sn_id}', color = 'lime')
print("SN coordinates = ",center_RA, center_Dec)

if RA_host_galaxy_GAMA is not None:
    ax1.scatter(RA_host_galaxy_GAMA, Dec_host_galaxy_GAMA, color='royalblue', marker = '^')#, label = f'cataID host {cataid_host_galaxy}')
if RA_host_galaxy_legacy is not None:
    ax1.scatter(RA_host_galaxy_legacy, Dec_host_galaxy_legacy, color='pink', marker = '^')#, label = f'cataID host {cataid_host_galaxy}')

plt.legend()


plt.savefig(f'plots/SN/legacy_comparison/{sn_id}_2arcmin_region_shape_r_0_excess_factor_2_5_parallax_0.png')
plt.savefig(f'plots/SN/legacy_comparison/pdfs/{sn_id}_2arcmin_region_shape_r_0_excess_factor_2_5_parallax_0.pdf')


###### NOW ONLY THE HOST GALAXY AND THE VECTORS ########
### Plot 2 settings ###

fig2, ax2 = plt.subplots()
ax2.set_title('SN and host galaxy')
ax2.set(xlabel = ('RA (deg)'), ylabel = ('Dec (deg)'))
ax2.set_aspect('equal')
ax2.invert_xaxis()


# Calculation of the DLR_host for every angle to do the plot of the ellipse
'''
print('major_host_galaxy_legacy', major_host_galaxy_legacy)
print(type(major_host_galaxy_legacy))
print(np.float32(major_host_galaxy_legacy))
print(type(np.float32(major_host_galaxy_legacy)))
print('major_host_galaxy_GAMA', major_host_galaxy_GAMA)
print(type(major_host_galaxy_GAMA))
print(np.float32(major_host_galaxy_GAMA))
'''
if ind_host_galaxy_legacy is not None:
    #print(ind_host_galaxy_legacy)
    #print('I am here')
    draw_ellipses([1,4], major_host_galaxy_legacy, minor_host_galaxy_legacy, position_angle_host_legacy, RA_host_galaxy_legacy, Dec_host_galaxy_legacy, ax2, colours=['red', 'red'], label = 'legacy')
    #print('works?')
    get_sn_angle(RA_host_galaxy_legacy, Dec_host_galaxy_legacy, major_host_galaxy_legacy, position_angle_host_legacy, center_RA, center_Dec, host=True, ax=ax2)
    get_dDLR_rotated(RA_host_galaxy_legacy, Dec_host_galaxy_legacy, center_RA, center_Dec,major_host_galaxy_legacy, minor_host_galaxy_legacy, position_angle_host_legacy, ax=ax2, plot=True)
    get_dDLR(RA_host_galaxy_legacy, Dec_host_galaxy_legacy, position_angle_host_legacy, major_host_galaxy_legacy, minor_host_galaxy_legacy, center_RA, center_Dec,ax=ax2, plot=True)

if ind_host_galaxy_GAMA is not None:
    #print('else I am here')
    draw_ellipses([1,4], major_host_galaxy_GAMA, minor_host_galaxy_GAMA, position_angle_host_GAMA, RA_host_galaxy_GAMA, Dec_host_galaxy_GAMA, ax2, colours = ['blue','blue'], label = 'GAMA')
    get_sn_angle(RA_host_galaxy_GAMA, Dec_host_galaxy_GAMA, major_host_galaxy_GAMA, position_angle_host_GAMA, center_RA, center_Dec, host=True, ax=ax2)
    get_dDLR_rotated(RA_host_galaxy_GAMA, Dec_host_galaxy_GAMA, center_RA, center_Dec,major_host_galaxy_GAMA, minor_host_galaxy_GAMA, position_angle_host_GAMA, ax=ax2, plot=True)
    get_dDLR(RA_host_galaxy_GAMA, Dec_host_galaxy_GAMA, position_angle_host_GAMA, major_host_galaxy_GAMA, minor_host_galaxy_GAMA, center_RA, center_Dec,ax=ax2, plot=True)

#radec_sn = co.SkyCoord(center_RA*u.deg, center_Dec*u.deg)
#radec_gal = co.SkyCoord(RA_host_galaxy*u.deg, Dec_host_galaxy*u.deg)
# scatter of the galaxy and the SN PLOT
ax2.scatter(center_RA, center_Dec, label = f'{sn_id}', color = 'lime')
#ax2.plot([center_RA,RA_host_galaxy], [center_Dec, Dec_host_galaxy], linestyle='dashed', color = 'green')


plt.savefig(f'plots/SN/legacy_comparison/{sn_id}_host_galaxy_shape_r_0_excess_factor_2_5_parallax_0.png')
plt.savefig(f'plots/SN/legacy_comparison/pdfs/{sn_id}_host_galaxy_shape_r_0_excess_factor_2_5_parallax_0.pdf')

plt.show()
