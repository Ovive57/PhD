### Script where I will calculate the dDLR step by step ###

# Basic imports
import os
import numpy as np
import random

from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import style

from astropy.table import Table
from astropy import coordinates as co, units as u
from astropy.io import fits

import pyvo
import urllib

import ipdb #*ipdb.set_trace()
# Constants
#seed = random.seed(25)

def draw_ellipses(scale_factor, major_axis, minor_axis, position_angle, ra_gal, dec_gal, ax):
    """Draws an ellipse with the DLR of a galaxy

    Args:
        scale_factor (list shape (2,)): Scale factor of the 2 ellipses in an array, how big they are going to be #! make a default value later
        major_axis (float64): Major axis of the galaxy. In arcsec
        minor_axis (float64): Minor axis of the galaxy. In arcsec
        position_angle (float32): Position angle of the galaxy in radians in mathematicians convention, i.e. starting from the RA axis and going anti-clockwise
        ra_gal (float64): ra of the galaxy
        dec_gal (float64): dec of the galaxy
    Returns:
        ax: the axis of the plot
    """
    # I define xi in every angle around a circle:
    xi_deg = np.linspace(0,360,1000)
    xi_rad = np.radians(xi_deg)

    # I calculate the DLR for all this circle:
    # major axis = major_galaxies_region
    # minor axis = minor_galaxies_region
    # ra_gal = RA_galaxies_region
    # position_angle = position_angle_region
    major_axis = major_axis/3600 #arcsec to deg
    minor_axis = minor_axis/3600 #arcsec to deg

    ### PLOTTING WITH THE DLR DIRECTLY
    DLR_ellipse = major_axis * minor_axis / np.sqrt((major_axis * np.sin(xi_rad))**2 + (minor_axis * np.cos(xi_rad))**2)


    # Adjust the radius for d = 1 and d = 4
    r_1 = scale_factor[0] * DLR_ellipse
    r_2 = scale_factor[1] * DLR_ellipse

    x_1 = r_1 * np.cos(xi_rad)
    y_1 = r_1 * np.sin(xi_rad)

    x_2 = r_2 * np.cos(xi_rad)
    y_2 = r_2 * np.sin(xi_rad)

    # Rotate the ellipse by the position angle of the galaxy
    x1_rotated = x_1 * np.cos(position_angle) + y_1 * np.sin(position_angle) + ra_gal
    y1_rotated = - x_1 * np.sin(position_angle) + y_1 * np.cos(position_angle) + dec_gal

    x2_rotated = x_2 * np.cos(position_angle) + y_2 * np.sin(position_angle) + ra_gal
    y2_rotated = - x_2 * np.sin(position_angle) + y_2 * np.cos(position_angle) + dec_gal

    # Plot the results

    ax.plot(x1_rotated, y1_rotated, color='red', linestyle='-')#, label = 'dDLR = 1/20')
    ax.plot(x2_rotated, y2_rotated, color='blue', linestyle='--')#, label = 'dDLR = 4/20')

    ax.scatter(ra_gal, dec_gal, marker = 'x', color = 'orange')#, label = 'galaxy') #! Can I put the label only once even if I do a for loop? If label is, put it in.
    '''
    # Ellipse with matplotlib is shorter and maybe more efficient
    ellipse = Ellipse(xy=(ra_gal, dec_gal), width=scale_factor[0]*major_axis*2, height=scale_factor[0]*minor_axis*2, angle=np.degrees(position_angle), edgecolor='red',linestyle = '--', fill = False)
    plt.gca().add_patch(ellipse)
    ellipse = Ellipse(xy=(ra_gal, dec_gal), width=scale_factor[1]*major_axis*2, height=scale_factor[1]*minor_axis*2, angle=np.degrees(position_angle), edgecolor='blue',linestyle = '--', fill = False)
    plt.gca().add_patch(ellipse)
    '''
    return ax

def get_sn_angle(ra_gal, dec_gal, major_axis, position_angle, ra_sn, dec_sn, host=False, plot=True, ax = None, plot_angle = False):
    """get the angle between the major axis of the galaxy and the SN.
    if plot=True, plot the major axis of the galaxy and the vector between the galaxy and the SN if we are plotting the host galaxy alone (host=True).


    Args:
        ra_gal (float64): ra of the galaxy
        dec_gal (float64): dec of the galaxy
        major_axis (float64): major axis of the galaxy. In arcsec
        position_angle (float32): position angle of the galaxy in radians in mathematicians convention, i.e. starting from the RA axis and going anti-clockwise
        ra_sn (float64): ra of the SN
        dec_sn (float64): dec of the SN
        host (bool, optional): whether we are plotting the host galaxy or not. Defaults to False.
        plot (bool, optional): whether we are plotting or not. Defaults to True.
        ax (matplotlib axis, optional): axis of the plot. Defaults to None.
    Returns:
        angle_from_major (float): angle from the major axis of the galaxy to the vector between the galaxy and the SN
    """
    major_axis = np.array(major_axis)/3600 #arcsec to deg

    if plot:
        major_x1 = ra_gal - major_axis * np.cos(position_angle)
        major_x2 = ra_gal + major_axis * np.cos(position_angle)

        major_y1 = dec_gal + major_axis * np.sin(position_angle)
        major_y2 = dec_gal - major_axis * np.sin(position_angle)
        ax.plot([major_x1, major_x2], [major_y1, major_y2], '-', color='red', zorder=100)


    # Separation in the 2D
    delta_ra = ra_sn - ra_gal
    delta_dec = dec_sn - dec_gal

    angle_from_ra = np.arctan2(delta_dec, delta_ra) #angle in radians already
    angle_from_major = (position_angle + angle_from_ra) # Because of how angles move, they move in opposite directions so instead of subtracting we add

    x_mid = (ra_gal + ra_sn) / 2
    y_mid = (dec_gal + dec_sn) / 2
    #plt.annotate(f'{np.degrees(angle_from_ra):.1f}°', xy=(x_mid+0.0002, y_mid), color='purple', fontsize=12) #plot the angle from RA
    if plot and host:
        ax.annotate(f'{np.degrees(angle_from_major):.1f}°', xy=(x_mid, y_mid), color='red', fontsize=12) #plot the angle from major axis
        print(f"The angle from the major axis is {np.degrees(angle_from_major):.1f}°, the angle from the RA axis is {np.degrees(angle_from_ra):.1f}° and the position angle of the galaxy is {np.degrees(position_angle):.1f}°")
    #angle_from_major = round(angle_from_major, 2)

    if plot_angle:
        x1_dlr_vector = ra_gal
        y1_dlr_vector = dec_gal

        x2_dlr_vector = ra_gal + (np.cos(angle_from_major))/3600
        y2_dlr_vector = dec_gal + (np.sin(angle_from_major))/3600
        ax.plot([x1_dlr_vector,x2_dlr_vector], [y1_dlr_vector,y2_dlr_vector], linestyle='-', color = 'black')
        ax.annotate(f'{np.degrees(angle_from_major):.1f}°', xy=(x_mid, y_mid), color='black', fontsize=12)

        x1_dlr_vector = ra_gal
        y1_dlr_vector = dec_gal

        x2_dlr_vector = ra_gal + (np.cos(angle_from_ra))/3600
        y2_dlr_vector = dec_gal + (np.sin(angle_from_ra))/3600

        #ax.plot([x1_dlr_vector,x2_dlr_vector], [y1_dlr_vector,y2_dlr_vector], linestyle='-', color = 'deeppink')
        ax.annotate(f'{np.degrees(angle_from_ra):.1f}°', xy=(x_mid+1, y_mid+1), color='deeppink', fontsize=12)
    return angle_from_major

def dDLR_ned(search_radec, galaxy_radec, galaxy_major, galaxy_minor, galaxy_angle_rad):

    """Calculates the dDLR between a galaxy and a SN without worrying about the angle between the galaxy and the SN

    Args:
        search_radec (float64): SN or other object coordinates
        galaxy_radec (float64): Galaxy coordinates
        galaxy_major (float64): Semi-major axis of the galaxy
        galaxy_minor (float64): Semi-minor axis of the galaxy
        galaxy_angle_rad (float32): Position angle of the galaxy. Mathematicians angles convention, i.e. starting from the RA axis and going anti-clockwise

    Returns:
        dDLR (float): dDLR between the object and the galaxy
        da, db (floats): distances to the object in the coordinates of the semi-major and semi-minor axis
    """

    # work out the RA/dec offsets from search pos to galaxy
    alpha, delta = search_radec.spherical_offsets_to(galaxy_radec)

    # choose to work in units of arcsec for convenience
    dx, dy = alpha.arcsec, delta.arcsec

    # now rotate these dx/dy to match the galaxy angle:
    cos_gal_ang, sin_gal_ang = np.cos(galaxy_angle_rad), np.sin(galaxy_angle_rad)
    # this is just the ordinary rotation matrix
    #! I changed again the signs, we are going the other way now
    da, db = cos_gal_ang*dx - sin_gal_ang*dy, +sin_gal_ang*dx + cos_gal_ang*dy
    #! note this is mathematicians and not astronomers angle
    # now da and db are separations in arcsec in the coord system of the semimaj/minor axes.

    dDLR = np.sqrt((da / galaxy_major)**2. + (db / galaxy_minor)**2. )

    return dDLR, da, db

def get_dDLR(ra_gal, dec_gal, position_angle_gal, major_gal, minor_gal, ra_sn, dec_sn, ax=None, plot = True):
    """Get the dDLR with the classic method

    Args:
        ra_gal (float64): RA of the galaxy
        dec_gal (float64): Dec of the galaxy
        position_angle_gal (float32): position angle of the galaxy in radians in mathematicians convention, i.e. starting from the RA axis and going anti-clockwise
        major_gal (float64): semi-major axis of the galaxy
        minor_gal (float64): semi-minor axis of the galaxy
        ra_sn (float64): RA of the SN
        dec_sn (float64): Dec of the SN
        plot (bool, optional): Whether to plot or not. Defaults to True.

    Returns:
        float64: the dDLR of the galaxy with the SN
    """
    gal_radec = co.SkyCoord(ra_gal, dec_gal, unit='deg')
    sn_radec = co.SkyCoord(ra_sn, dec_sn, unit = 'deg')
    sep = gal_radec.separation(sn_radec).arcsec

    phi = get_sn_angle(ra_gal, dec_gal, major_gal, position_angle_gal, ra_sn, dec_sn, host=False, plot=False, plot_angle=False)
    DLR = (major_gal*minor_gal)/(np.sqrt((major_gal*np.sin(phi))**2+(minor_gal*np.cos(phi))**2)) #major and minor are in arcsec
    dDLR_Ned_value = dDLR_ned(sn_radec, gal_radec,major_gal, minor_gal, position_angle_gal)[0]
    DLR_ned = sep/dDLR_Ned_value
    dDLR = sep/DLR
    if plot:  #! HERE is the plot of the vector
        x1_dlr_vector = ra_gal
        y1_dlr_vector = dec_gal
        #x2_dlr = dDLR * np.cos(phi)
        #y2_dlr = dDLR * np.sin(phi)

        #x2_dlr_vector = x2_dlr * np.cos(position_angle_gal) - y2_dlr * np.sin(position_angle_gal) + ra_gal
        #y2_dlr_vector = x2_dlr * np.sin(position_angle_gal) + y2_dlr * np.cos(position_angle_gal) + dec_gal
        unit_dDLR = DLR #! It is 1 DLR
        x2_dlr_vector = ra_gal + (unit_dDLR * np.cos(phi-position_angle_gal))/3600
        y2_dlr_vector = dec_gal + (unit_dDLR * np.sin(phi-position_angle_gal))/3600
        ax.plot([x1_dlr_vector,x2_dlr_vector], [y1_dlr_vector,y2_dlr_vector], linestyle='-', color = 'black')

        unit_dDLR_Ned = DLR_ned #! It is 1 DLR
        x2_dlr_vector = ra_gal + (unit_dDLR_Ned * np.cos(phi - position_angle_gal))/3600
        y2_dlr_vector = dec_gal + (unit_dDLR_Ned * np.sin(phi - position_angle_gal))/3600
        ax.plot([x1_dlr_vector,x2_dlr_vector], [y1_dlr_vector,y2_dlr_vector], linestyle='-', color = 'deeppink')

    return dDLR

def find_host_galaxy(ra_gal, dec_gal, major_gal, minor_gal, position_angle_gal, ra_sn, dec_sn):
    """Finding the index of host galaxy with the minimum dDLR

    Args:
        ra_gal (astropy.table.column.Column): ra of the galaxies
        dec_gal (astropy.table.column.Column): dec of the galaxies
        major_gal (astropy.table.column.Column): semi-major axis of the galaxies
        minor_gal (astropy.table.column.Column): semi-minor axis of the galaxies
        position_angle_gal (astropy.table.column.Column): position angle of the galaxies. Mathematicians angles convention, i.e. starting from the RA axis and going anti-clockwise
        ra_sn (float64): ra of the SN
        dec_sn (float64): dec of the SN

    Returns:
        numpy.int64: index of the host galaxy
    """

    dDLRs_ned = []
    dDLRs = []

    for i,ra in enumerate(ra_gal):
        gal_radec = co.SkyCoord(ra_gal[i], dec_gal[i], unit='deg')
        sn_radec = co.SkyCoord(ra_sn, dec_sn, unit = 'deg')
        sep = gal_radec.separation(sn_radec).arcsec
        # I calculate the dDLR if we are within 15"
        if sep < 30:
            dDLR_value_ned, da, db = dDLR_ned(sn_radec, gal_radec, major_gal[i], minor_gal[i], position_angle_gal[i])
            dDLR_value = get_dDLR(ra_gal[i], dec_gal[i], position_angle_gal[i], major_gal[i], minor_gal[i], ra_sn, dec_sn, plot=False)
            dDLRs_ned.append(dDLR_value_ned)
            dDLRs.append(dDLR_value)
            #separation = np.sqrt((da)**2+(db)**2)
        # else, I put something big so it is not going to be the minimum.
        else:
            dDLRs_ned.append(np.inf)
            dDLRs.append(np.inf)

    dDLRs_ned = np.array(dDLRs_ned)
    dDLRs = np.array(dDLRs)

    ind = np.argmin(dDLRs)
    print('neds dDLR = ',dDLRs_ned[ind], 'dDLR = ', dDLRs[ind])

    return ind


def get_lsdr10_cat( identifier, ra, dec, search_radius_deg, i_limit=18., overwrite=False ):
    cat_filename = f'{identifier}.ls_dr10.cat.fits'

    if os.path.exists( cat_filename ) and not overwrite:
        return Table.read( cat_filename )

    print( f'querying legacy with a search radius of {search_radius_deg} deg.' )

    tap_service = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
    ex_query = f"""
            SELECT tr.ls_id, tr.ra, tr.dec, tr.ref_cat, tr.ref_id, tr.gaia_phot_bp_rp_excess_factor, 
                tr.mag_g, tr.mag_r, tr.mag_i, tr.mag_z, tr.mag_w1, tr.mag_w2, 
                tr.sersic, tr.sersic_ivar,  tr.shape_e1, tr.shape_e1_ivar, tr.shape_e2, tr.shape_e2_ivar, 
                tr.shape_r,  tr.shape_r_ivar, zp.z_spec, zp.survey, zp.z_phot_mean, zp.z_phot_mean_i, 
                zp.z_phot_std, zp.z_phot_l68, zp.z_phot_u68, zp.z_phot_l95, zp.z_phot_u95
            FROM ls_dr10.tractor as tr JOIN ls_dr10.photo_z as zp ON tr.ls_id = zp.ls_id
            WHERE 't' = Q3C_RADIAL_QUERY(tr.ra,tr.dec,{ra:.6f},{dec:+.6f},{search_radius_deg:.6f}) AND tr.mag_i < {i_limit}
            """

    # NB. uses gaia 'excess factor' to exclude gaia point sources (i.e. stars)
    # Can effectively exclude stars using
    # AND (gaia_phot_bp_rp_excess_factor > 3.5 or gaia_phot_bp_rp_excess_factor <= 0)
    #         AND shape_e1 > 0 AND shape_e2 > 0  --- to exclude saturated stars
    print( ex_query )

    result = tap_service.search(ex_query)
    tblResult = result.to_table()
    print(f'query result has {len(tblResult)} rows.' )

    print('writing query result to:', cat_filename )
    tblResult.write(cat_filename, overwrite=overwrite )

    return tblResult
from PIL import Image
import requests
def get_jpeg_cutout(size_arcmin, filename, ra, dec, pixscale = 0.262):
    #! pixscale = 0.262 # this is the LS DR10 standard; i.e. no resampling
    num_pixels = int( (size_arcmin*60.) / pixscale )+1
    imagequery = (
        f'https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&size={num_pixels}&pixscale={pixscale}&layer=ls-dr10')  

    # and here is actually issuing the request and opening image
    f = urllib.request.urlopen(imagequery)
    im = Image.open(requests.get(imagequery, stream=True).raw)
    #mirrored_img = im.transpose(Image.FLIP_LEFT_RIGHT)
    # NB. this does not actually save the image to disk
    # but this does!
    im.save(filename)


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
#ind = np.where((RA_sn>=RA_gama.min())&(Dec_sn<=Dec_gama.max()))
#sn_ids = sn_ids[ind]
#RA_sn = RA_sn[ind]
#Dec_sn = Dec_sn[ind]
ind = random.randint(0, len(sn_ids)-1)

center_RA = RA_sn[ind]
center_Dec = Dec_sn[ind]


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
print(ID)
print(Dec_galaxies_region)
#! I need to think more about this:
# dDLR = ang_sep/DLR. If dDLR = 1 -> ang_sep = DLR, if dDLR = 4 -> ang_sep = 4*DLR. For this I would need to calculate the DLR, i.e. take into accound the angle again?
# But I don't need the angle since I want to plot all the ellipse.

# dDLR: is a distance in DLR units. This metric is used to express the distance of a supernova relative to a galaxy's characteristic length scale (DLR).
# A point with dDLR=1 is one DLR unit away from the galaxy center, while a point with dDLR=4 is four DLR units away.
# To get the ellipse corresponding to a particular dDLR value, I scale the DLR by the dDLR value.

# For dDLR = 1, the radius of the ellipse is the DLR itself, because ang_sep = DLR
# For dDLR = 4, the radius of the ellipse is 4 times the DLR

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

#legacy_survey_table = get_lsdr10_cat('deleteme', center_RA, center_Dec, radius_deg, overwrite=True)
#hdu = fits.open('deleteme.ls_dr10.cat.fits')
#print(hdu.info())
#image_data = hdu.data
#print(image_data)
#ipdb.set_trace()
#ax3.imshow(hdu, origin='lower', extent=[min_RA, max_RA, min_Dec, max_Dec])
#exit()
#ax1.imshow(image_data, origin='lower', extent=[min_RA, max_RA, min_Dec, max_Dec])

# The two concentric ellipses for all galaxies in the region:
#########################################
#! TO TEST
#!position_angle_region = [np.radians(-45) for i in range(len(position_angle_region))]
#!print(np.degrees(position_angle_region))
#########################################

for i,pos in enumerate(RA_galaxies_region):
    ax1 = draw_ellipses([1, 4], major_galaxies_region[i], minor_galaxies_region[i], position_angle_region[i], RA_galaxies_region[i], Dec_galaxies_region[i], ax1)

angles_from_major = []
for i,pos in enumerate(RA_galaxies_region):
    angle = get_sn_angle(RA_galaxies_region[i], Dec_galaxies_region[i], major_galaxies_region[i], position_angle_region[i], center_RA, center_Dec,plot=True, ax=ax1)
    angles_from_major.append(angle)

###### FINDING THE HOST GALAXY ######
# Another loop again, the loop before was only for some tests plottings, this one is going to be the real calculation (?)
#position_angle_region = [round(position_angle_region[i], 3) for i in range(len(position_angle_region))]
ind_host_galaxy= find_host_galaxy(RA_galaxies_region, Dec_galaxies_region, major_galaxies_region, minor_galaxies_region, position_angle_region, center_RA, center_Dec)

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
plt.show()

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
get_dDLR(RA_host_galaxy, Dec_host_galaxy, position_angle_host, major_host_galaxy, minor_host_galaxy, center_RA, center_Dec,ax=ax2, plot=True)

# scatter of the galaxy and the SN PLOT
ax2.scatter(center_RA, center_Dec, label = 'SN')
#ax2.plot([center_RA,RA_host_galaxy], [center_Dec, Dec_host_galaxy], linestyle='dashed', color = 'green')

plt.legend()
plt.show()

