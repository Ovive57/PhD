# Basic imports
import os
import numpy as np
import random
import time

# Plots
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
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
from dlr_stepbystep import get_sn_angle, get_dDLR, get_dDLR_rotated, draw_ellipses, find_host_galaxy, get_jpeg_cutout
#### The DATA ####
galaxies_table = Table.read('data/gkvScienceCatv02.fits')
#SN_table = Table.read('data/tns_SNIa_20240424_copy.csv')
SN_table = Table.read('crossmatches/specifics/colour_mass_matched465.fits')
## GALAXIES
RA_gama = galaxies_table['RAcen']
Dec_gama = galaxies_table['Deccen']
ID = galaxies_table['CATAID']
position_angle_gama = galaxies_table['ang']
axrat = galaxies_table['axrat']
major_gama = galaxies_table['R50'] # In arcsec
minor_gama = major_gama*axrat

## SN
RA_sn = SN_table['ra']
Dec_sn = SN_table['declination']
sn_coords = co.SkyCoord(RA_sn*u.deg, Dec_sn*u.deg)
sn_ids = SN_table['objid']


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

print("\n number of galaxies: ",len(RA_galaxies_region))

major_galaxies_region = major_gama[ind]
minor_galaxies_region = minor_gama[ind]
#print(position_angle_gama[ind]-90)
position_angle_region = np.radians(position_angle_gama[ind]-90)
#position_angle_region = [round(position_angle_region[i], 3) for i in range(len(position_angle_region))]
ID = ID[ind]

table = Table([ID, RA_galaxies_region, Dec_galaxies_region], names=('CATAID', 'RA_galaxy', 'Dec_galaxy'))

# dDLR: is a distance in DLR units. This metric is used to express the distance of a supernova relative to a galaxy's characteristic length scale (DLR).
# A point with dDLR=1 is one DLR unit away from the galaxy center, while a point with dDLR=4 is four DLR units away.
# To get the ellipse corresponding to a particular dDLR value, I scale the DLR by the dDLR value.

# For dDLR = 1, the radius of the ellipse is the DLR itself, because ang_sep = DLR. It is the size of the galaxy.
# For dDLR = 4, the radius of the ellipse is 4 times the DLR

# If dDLR < 1 : SN inside the galaxy, if dDLR > 1 : SN outside the galaxy.

### Plot 1 settings ###

plt.style.use(style.style1)
fig1, ax1 = plt.subplots()
ax1.set_title('SN and galaxies in a 2 arcmin region')
ax1.set(xlabel = ('RA (deg)'), ylabel = ('Dec (deg)'), xlim=(max_RA, min_RA), ylim=(min_Dec, max_Dec))
ax1.set_aspect('equal')
#ax1.invert_xaxis()


get_jpeg_cutout(2, 'deleteme.jpeg', center_RA, center_Dec)
import matplotlib.image as mpimg
img = mpimg.imread('deleteme.jpeg')
#img = np.fliplr(img)

plt.imshow(img, extent=[max_RA, min_RA, max_Dec, min_Dec])


# The two concentric ellipses for all galaxies in the region:

#!########################################
#! TO TEST
#!position_angle_region = [np.radians(-45) for i in range(len(position_angle_region))]
#!print(np.degrees(position_angle_region))
#!########################################

for i,pos in enumerate(RA_galaxies_region):
    ax1 = draw_ellipses([1, 4], major_galaxies_region[i], minor_galaxies_region[i], position_angle_region[i], RA_galaxies_region[i], Dec_galaxies_region[i], ax1)

angles_from_major = []
classic_dDLRs = []
rotated_dDLRs = []
times_classic = []
times_rotated = []
for i,pos in enumerate(RA_galaxies_region):
    angle = get_sn_angle(RA_galaxies_region[i], Dec_galaxies_region[i], major_galaxies_region[i], position_angle_region[i], center_RA, center_Dec,plot=True, ax=ax1)
    classic_dDLR, time_classic = get_dDLR(RA_galaxies_region[i], Dec_galaxies_region[i], position_angle_region[i], major_galaxies_region[i], minor_galaxies_region[i], center_RA, center_Dec, plot=False)
    rotated_dDLR,a,b, time_rotated = get_dDLR_rotated(sn_coords, co.SkyCoord(RA_galaxies_region[i]*u.deg, Dec_galaxies_region[i]*u.deg), major_galaxies_region[i], minor_galaxies_region[i], position_angle_region[i], plot=False)
    angles_from_major.append(angle)
    classic_dDLRs.append(classic_dDLR)
    rotated_dDLRs.append(rotated_dDLR)
    times_classic.append(time_classic)
    times_rotated.append(time_rotated)

table['classic_dDLR'] = classic_dDLRs
table['rotated_dDLR'] = rotated_dDLRs
table['time_classic'] = times_classic
table['time_rotated'] = times_rotated
table.write(f'files/SN_2arcmin_region/{sn_id}.fits', overwrite=True)

###### FINDING THE HOST GALAXY ######
ind_host_galaxy, dDLR_rotated, dDLR_classic= find_host_galaxy(RA_galaxies_region, Dec_galaxies_region, major_galaxies_region, minor_galaxies_region, position_angle_region, center_RA, center_Dec)

# The host galaxy
RA_host_galaxy = RA_galaxies_region[ind_host_galaxy]
Dec_host_galaxy = Dec_galaxies_region[ind_host_galaxy]

major_host_galaxy = major_galaxies_region[ind_host_galaxy]
minor_host_galaxy = minor_galaxies_region[ind_host_galaxy]
position_angle_host = position_angle_region[ind_host_galaxy]

# PLOT the host galaxy in the same plot as the rest of the galaxies
ax1.scatter(center_RA, center_Dec, label = 'SN', color = 'deeppink')
print("SN coordinates = ",center_RA, center_Dec)
ax1.scatter(RA_host_galaxy, Dec_host_galaxy,color='green', marker = '^', label = 'host galaxy' )
plt.legend()
plt.savefig(f'plots/SN/{sn_id}_2arcmin_region.png')
plt.savefig(f'plots/SN/pdfs/{sn_id}_2arcmin_region.pdf')
#plt.show()

###### NOW ONLY THE HOST GALAXY AND THE VECTORS ########
### Plot 2 settings ###

fig2, ax2 = plt.subplots()
ax2.set_title('SN and host galaxy')
ax2.set(xlabel = ('RA (deg)'), ylabel = ('Dec (deg)'))
ax2.set_aspect('equal')
ax2.invert_xaxis()



# Calculation of the DLR_host for every angle to do the plot of the ellipse
draw_ellipses([1,4], major_host_galaxy, minor_host_galaxy, position_angle_host, RA_host_galaxy, Dec_host_galaxy, ax2)
#* This is just to have the message in the terminal and the vector of the dDLR, not useful:
get_sn_angle(RA_host_galaxy, Dec_host_galaxy, major_host_galaxy, position_angle_host, center_RA, center_Dec, host=True, ax=ax2)
radec_sn = co.SkyCoord(center_RA*u.deg, center_Dec*u.deg)
radec_gal = co.SkyCoord(RA_host_galaxy*u.deg, Dec_host_galaxy*u.deg)
get_dDLR_rotated(radec_sn, radec_gal,major_host_galaxy, minor_host_galaxy, position_angle_host, ax=ax2, plot=True)
get_dDLR(RA_host_galaxy, Dec_host_galaxy, position_angle_host, major_host_galaxy, minor_host_galaxy, center_RA, center_Dec,ax=ax2, plot=True)

# scatter of the galaxy and the SN PLOT
ax2.scatter(center_RA, center_Dec, label = 'SN')
#ax2.plot([center_RA,RA_host_galaxy], [center_Dec, Dec_host_galaxy], linestyle='dashed', color = 'green')

plt.legend()
plt.savefig(f'plots/SN/{sn_id}_host_galaxy.png')
plt.savefig(f'plots/SN/pdfs/{sn_id}_host_galaxy.pdf')
#plt.show()

