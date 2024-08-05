### Script where I will calculate the dDLR step by step ###

# Basic imports
import numpy as np
import random

from matplotlib import pyplot as plt

from astropy.table import Table
from astropy import coordinates as co, units as u

# Constants
seed = random.seed(25)

def dDLR(search_radec, galaxy_radec, galaxy_major, galaxy_minor, galaxy_angle_rad):
    """Calculates the dDLR between a galaxy and a SN without worrying about the angle between the galaxy and the SN

    Args:
        search_radec (_type_): SN or other object coordinates
        galaxy_radec (_type_): Galaxy coordinates
        galaxy_major (_type_): Semi-major axis of the galaxy
        galaxy_minor (_type_): Semi-minor axis of the galaxy
        galaxy_angle_rad (_type_): Position angle of the galaxy. Mathematicians angles convention.

    Returns:
        dDLR (float): dDLR between the object and the galaxy
        da, db (floats): distances to the object in the coordinates of the semi-major and semi-minor axis
    """

    # work out the RA/dec offsets from search pos to galaxy
    alpha, delta = search_radec.spherical_offsets_to(galaxy_radec)

    # choose to work in units of arcsec for convenience
    dx, dy = alpha.arcsec, delta.arcsec

    # now rotate these dx/dy to match the galaxy angle:
    cosphi, sinphi = np.cos(galaxy_angle_rad), np.sin(galaxy_angle_rad)
    # this is just the ordinary rotation matrix
    da, db = cosphi*dx - sinphi*dy, sinphi*dx + cosphi*dy
    #! note this is mathematicians and not astronomers angle
    # now da and db are separations in arcsec in the coord system of the semimaj/minor axes.

    dDLR = np.sqrt((da / galaxy_major)**2. + (db / galaxy_minor)**2. )

    return dDLR, da, db



#### The DATA ####
galaxies_table = Table.read('data/gkvScienceCatv02.fits')
SN_table = Table.read('data/tns_SNIa_20240424_copy.csv')

## GALAXIES
RA_gama = galaxies_table['RAcen']
Dec_gama = galaxies_table['Deccen']
ID = galaxies_table['CATAID']
position_angle_gama = galaxies_table['ang']
axrat = galaxies_table['axrat']
major_gama = galaxies_table['R50']/3600 # From arcsec to deg
minor_gama = major_gama*axrat

## SN
RA_sn = SN_table['ra']
Dec_sn = SN_table['declination']
sn_coords = co.SkyCoord(RA_sn*u.deg, Dec_sn*u.deg)
sn_ids = SN_table['objid']

# Selection of my region:
ind = np.where((RA_sn>=RA_gama.min())&(Dec_sn<=Dec_gama.max()))
sn_ids = sn_ids[ind]
RA_sn = RA_sn[ind]
Dec_sn = Dec_sn[ind]
ind = random.randint(0, len(sn_ids)) #! Later: Take one of the SN with an actual cross-match (host galaxy) with the old method

center_RA = RA_sn[ind]
center_Dec = Dec_sn[ind]

#print(center_RA, center_Dec, RA_gama.min(), Dec_gama.min(), RA_gama.max(), Dec_gama.max())

max_RA = center_RA + 1/60 #1 arcmin in deg
min_RA = center_RA - 1/60
max_Dec = center_Dec + 1/60
min_Dec = center_Dec - 1/60


# Find the data within this region:

# We have already the SN
# The galaxies:
ind = np.where((RA_gama>min_RA)&(RA_gama<max_RA)&(Dec_gama>min_Dec)&(Dec_gama<max_Dec))
RA_galaxies_region = RA_gama[ind]
Dec_galaxies_region = Dec_gama[ind]

major_galaxies_region = major_gama[ind]
minor_galaxies_region = minor_gama[ind]
print(position_angle_gama[ind]-90)
position_angle_region = np.radians(position_angle_gama[ind]-90)

ID = ID[ind]
print(ID)

# Now I am going to plot the zone of influence of each galaxy. The way I am going to do this plotting the zone where dDLR = 1 and dDLR = 4

# dDLR = ang_sep/DLR. If dDLR = 1 -> ang_sep = DLR, if dDLR = 4 -> ang_sep = 4*DLR. For this I would need to calculate the DLR, i.e. take into accound the angle again?
# But I don't need the angle since I want to plot all the ellipse.

ang_sep = 1 # I put it at 1 for simplicity

# I define phi in every angle around a circle:
phi_deg = np.linspace(0,360,1000)
phi_rad = np.radians(phi_deg)

# I calculate the DLR for all this circle:
DLR = [np.sqrt(major_galaxies_region[i] * minor_galaxies_region[i] / np.sqrt((major_galaxies_region[i] * np.sin(phi_rad))**2 + (minor_galaxies_region[i] * np.cos(phi_rad))**2)) for i in range(len(RA_galaxies_region))]
# Adjust the radius for d = 1 and d = 4


#! We need to do a loop, for the moment one by one to do comprobations
#i = 6
plt.figure(figsize=(8, 8))
for i in range(len(major_galaxies_region)):
    r_d1 = ang_sep / 1 * DLR[i]
    r_d4 = ang_sep / 4 * DLR[i]

    x_1 = r_d1 * np.cos(phi_rad)
    y_1 = r_d1 * np.sin(phi_rad)

    x_4 = r_d4 * np.cos(phi_rad)
    y_4 = r_d4 * np.sin(phi_rad)

    # Rotate the ellipse by the position angle of the galaxy
    x1_rotated = x_1 * np.cos(position_angle_region[i]) - y_1 * np.sin(position_angle_region[i]) + RA_galaxies_region[i]
    y1_rotated = x_1 * np.sin(position_angle_region[i]) + y_1 * np.cos(position_angle_region[i]) + Dec_galaxies_region[i]

    x4_rotated = x_4 * np.cos(position_angle_region[i]) - y_4 * np.sin(position_angle_region[i]) + RA_galaxies_region[i]
    y4_rotated = x_4 * np.sin(position_angle_region[i]) + y_4 * np.cos(position_angle_region[i]) + Dec_galaxies_region[i]

    # Plot the results
    plt.plot(x1_rotated, y1_rotated, color='red', linestyle='-')
    plt.plot(x4_rotated, y4_rotated, color='blue', linestyle='--')


### Plot settings ###
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.title('SN and galaxies in a 2 arcmin region')
plt.gca().invert_xaxis()

# PLOT
plt.scatter(center_RA, center_Dec, label = 'SN')
plt.scatter(RA_galaxies_region, Dec_galaxies_region, marker = 'x', label = 'galaxies')
plt.legend()
plt.show()

