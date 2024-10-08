### Script where I will calculate the dDLR step by step ###

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

# Constants
#seed = random.seed(25)

def rotation(x, y, angle, rot_type = 'anticlockwise'):
    if rot_type == 'anticlockwise':
        x_rotated = x * np.cos(angle) + y * np.sin(angle)
        y_rotated = - x * np.sin(angle) + y * np.cos(angle)
    elif rot_type == 'clockwise':
        x_rotated = x * np.cos(angle) - y * np.sin(angle)
        y_rotated = x * np.sin(angle) + y * np.cos(angle)
    return x_rotated, y_rotated

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

def draw_ellipses(scale_factor, major_axis, minor_axis, position_angle, ra_gal, dec_gal, ax, rot_type = 'anticlockwise', colours = ['red', 'blue'], label = None):
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
    '''
    print("\nTHE MAJOR AXIS",major_axis, "\n")
    if major_axis == None:
        ax.scatter(ra_gal, dec_gal, marker = 'x', color = 'orange')
        return ax
    '''
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
    #x1_rotated = x_1 * np.cos(position_angle) + y_1 * np.sin(position_angle) + ra_gal
    #y1_rotated = - x_1 * np.sin(position_angle) + y_1 * np.cos(position_angle) + dec_gal

    x1_rotated, y1_rotated = rotation(x_1, y_1, position_angle, rot_type = rot_type)
    x1_rotated = x1_rotated + ra_gal
    y1_rotated = y1_rotated + dec_gal

    #x2_rotated = x_2 * np.cos(position_angle) + y_2 * np.sin(position_angle) + ra_gal
    #y2_rotated = - x_2 * np.sin(position_angle) + y_2 * np.cos(position_angle) + dec_gal

    x2_rotated, y2_rotated = rotation(x_2, y_2, position_angle, rot_type = rot_type)
    x2_rotated = x2_rotated + ra_gal
    y2_rotated = y2_rotated + dec_gal

    # Plot the results
    if label:
        ax.plot(x1_rotated, y1_rotated, color=colours[0], linestyle='-', label = label)
        ax.plot(x2_rotated, y2_rotated, color=colours[1], linestyle='--')
    else:
        ax.plot(x1_rotated, y1_rotated, color=colours[0], linestyle='-')#, label = 'dDLR = 1/20')
        ax.plot(x2_rotated, y2_rotated, color=colours[1], linestyle='--')#, label = 'dDLR = 4/20')

    #ax.scatter(ra_gal, dec_gal, marker = 'x', color = 'orange')#, label = 'galaxy') #! Can I put the label only once even if I do a for loop? If label is, put it in.

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
    delta_ra = (ra_sn - ra_gal)*np.cos(np.radians(dec_sn)) # This is an approximation for small angles
    delta_dec = dec_sn - dec_gal

    angle_from_ra = np.arctan2(delta_dec, delta_ra) #angle in radians already
    angle_from_major = (position_angle + angle_from_ra) # Because of how angles move, they move in opposite directions so instead of subtracting we add

    x_mid = (ra_gal + ra_sn) / 2
    y_mid = (dec_gal + dec_sn) / 2
    #plt.annotate(f'{np.degrees(angle_from_ra):.1f}°', xy=(x_mid+0.0002, y_mid), color='purple', fontsize=12) #plot the angle from RA
    if plot and host:
        ax.annotate(f'{np.degrees(angle_from_major):.1f}°', xy=(x_mid, y_mid), color='red', fontsize=12) #plot the angle from major axis
        #print(f"The angle from the major axis is {np.degrees(angle_from_major):.1f}°, the angle from the RA axis is {np.degrees(angle_from_ra):.1f}° and the position angle of the galaxy is {np.degrees(position_angle):.1f}°")


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

def get_dDLR_rotated(ra_gal, dec_gal, ra_sn, dec_sn, galaxy_major, galaxy_minor, galaxy_angle_rad, ax=None, plot = True):

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
    start = time.time()
    gal_radec = co.SkyCoord(ra_gal, dec_gal, unit='deg')
    sn_radec = co.SkyCoord(ra_sn, dec_sn, unit = 'deg')
    # work out the RA/dec offsets from search pos to galaxy
    alpha, delta = sn_radec.spherical_offsets_to(gal_radec)

    # choose to work in units of arcsec for convenience
    dx, dy = alpha.arcsec, delta.arcsec

    # now rotate these dx/dy to match the galaxy angle:
    cos_gal_ang, sin_gal_ang = np.cos(galaxy_angle_rad), np.sin(galaxy_angle_rad)
    # this is just the ordinary rotation matrix
    #da, db = cos_gal_ang*dx - sin_gal_ang*dy, + sin_gal_ang*dx + cos_gal_ang*dy
    da = rotation(dx, dy, galaxy_angle_rad, rot_type = 'clockwise')[0]
    db = rotation(dx, dy, galaxy_angle_rad, rot_type = 'clockwise')[1]
    #! note this is mathematicians and not astronomers angle
    # now da and db are separations in arcsec in the coord system of the semimaj/minor axes.

    dDLR = np.sqrt((da / galaxy_major)**2. + (db / galaxy_minor)**2. )
    end = time.time()
    #! Change this when I discover why:
    #print(f'shape dDLR = {np.shape(dDLR)}')
    dDLR = dDLR.min()

    sep = np.sqrt(da**2 + db**2)
    DLR = sep/dDLR
    phi = np.arctan2(db,da)
    phi2 = np.arcsin(db/sep) # This put the angle always in the first and forth quadrant (-90° to 90°)
    phi3 = np.arccos(da/sep) # This put the angle always in the first and second quadrant (0° to 180°)

    #print(f' Rotated axis angle = {np.degrees(phi):.1f}°, other angle = {np.degrees(phi2):.1f}°, other angle = {np.degrees(phi3):.1f}°')
    #ra_gal = galaxy_radec.ra.deg
    #dec_gal = galaxy_radec.dec.deg
    if plot:
        unit_dDLR = DLR #! It is 1 DLR
        x1_dlr_vector = ra_gal
        y1_dlr_vector = dec_gal
        x2_dlr_vector = ra_gal + (unit_dDLR * np.cos(phi - galaxy_angle_rad + np.pi))/3600 # I added 180° to have the angle in the good quadrant
        y2_dlr_vector = dec_gal + (unit_dDLR * np.sin(phi - galaxy_angle_rad + np.pi))/3600
        ax.plot([x2_dlr_vector,x1_dlr_vector], [y2_dlr_vector,y1_dlr_vector], linestyle='-', color = 'deeppink')

        print(f'rotated axis dDLR calculation took {end-start} seconds')
    return dDLR, sep, end-start

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
    start = time.time()
    gal_radec = co.SkyCoord(ra_gal, dec_gal, unit='deg')
    sn_radec = co.SkyCoord(ra_sn, dec_sn, unit = 'deg')
    sep = gal_radec.separation(sn_radec).arcsec

    phi = get_sn_angle(ra_gal, dec_gal, major_gal, position_angle_gal, ra_sn, dec_sn, host=False, plot=False, plot_angle=False)
    DLR = (major_gal*minor_gal)/(np.sqrt((major_gal*np.sin(phi))**2+(minor_gal*np.cos(phi))**2)) #major and minor are in arcsec

    dDLR = sep/DLR
    end = time.time()
    if plot:  #! HERE is the plot of the vector
        x1_dlr_vector = ra_gal
        y1_dlr_vector = dec_gal
        unit_dDLR = DLR #! It is 1 DLR
        x2_dlr_vector = ra_gal + (unit_dDLR * np.cos(phi-position_angle_gal))/3600
        y2_dlr_vector = dec_gal + (unit_dDLR * np.sin(phi-position_angle_gal))/3600
        ax.plot([x1_dlr_vector,x2_dlr_vector], [y1_dlr_vector,y2_dlr_vector], linestyle='-', color = 'black')
        print(f'classic dDLR calculation took {end-start} seconds')
    return dDLR, sep, end-start

def find_host_galaxy(ra_gal, dec_gal, major_gal, minor_gal, position_angle_gal, ra_sn, dec_sn, plot=False, catalogue = 'GAMA'):
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

    dDLRs_rotated = []
    dDLRs = []
    times_rotated = []
    times_classic = []
    separation_classic = []
    separation_rotated = []

    for i,ra in enumerate(ra_gal):
        # I calculate the dDLR for all the galaxies because I asume that a previous selection was done, i.e. the function is called after a selection of a region near the SN
        dDLR_value_rotated, sep_rotated, time_rotated = get_dDLR_rotated(ra_gal[i], dec_gal[i], ra_sn, dec_sn, major_gal[i], minor_gal[i], position_angle_gal[i], plot = False)
        dDLR_value, sep_classic, time_classic = get_dDLR(ra_gal[i], dec_gal[i], position_angle_gal[i], major_gal[i], minor_gal[i], ra_sn, dec_sn, plot=False)
        if sep_classic < 47:
            dDLRs_rotated.append(dDLR_value_rotated)
            dDLRs.append(dDLR_value)
            times_rotated.append(time_rotated)
            times_classic.append(time_classic)
            separation_classic.append(sep_classic)
            separation_rotated.append(sep_rotated)
            # else, I put something big so it is not going to be the minimum.
        else:
            dDLRs_rotated.append(np.inf)
            dDLRs.append(np.inf)
            times_rotated.append(None)
            times_classic.append(None)
            separation_classic.append(np.inf)
            separation_rotated.append(np.inf)

    dDLRs_rotated = np.array(dDLRs_rotated)
    dDLRs = np.array(dDLRs)
    times_rotated = np.array(times_rotated)
    times_classic = np.array(times_classic)
    separation_classic = np.array(separation_classic)
    separation_rotated = np.array(separation_rotated)

    if len(dDLRs) > 0:
        ind = np.argmin(dDLRs)
        if dDLRs[ind] > 4:
            ind = None
            if plot:
                print(f'\nNo host galaxy found for {catalogue}\n')
            return ind, [None, None, None], [None, None, None]
        else:
            return ind, [dDLRs_rotated[ind], times_rotated[ind], separation_rotated[ind]], [dDLRs[ind], times_classic[ind], separation_classic[ind]]


    else:
        ind = None
        return ind, [None, None, None], [None, None, None]
    #print('rotated axis dDLR = ',dDLRs_rotated[ind], ', classic dDLR = ', dDLRs[ind])




def get_lsdr10_cat(identifier, ra, dec, search_radius_deg, overwrite=False, shape_r = 0, excess_factor = 2.5, parallax = 5):


    print( f'querying legacy with a search radius of {search_radius_deg} deg.' )

    tap_service = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")

    if shape_r!= None:
        ex_query = f"""
                SELECT tr.ls_id, tr.ra, tr.dec, tr.ref_cat, tr.ref_id, tr.gaia_phot_bp_rp_excess_factor, tr.parallax,
                    tr.mag_g, tr.mag_r, tr.mag_i, tr.mag_z, tr.mag_w1, tr.mag_w2, 
                    tr.sersic, tr.sersic_ivar,  tr.shape_e1, tr.shape_e1_ivar, tr.shape_e2, tr.shape_e2_ivar, 
                    tr.shape_r,  tr.shape_r_ivar, tr.gaia_duplicated_source, zp.z_spec, zp.survey, zp.z_phot_mean, zp.z_phot_mean_i, 
                    zp.z_phot_std, zp.z_phot_l68, zp.z_phot_u68, zp.z_phot_l95, zp.z_phot_u95
                FROM ls_dr10.tractor as tr JOIN ls_dr10.photo_z as zp ON tr.ls_id = zp.ls_id
                WHERE 't' = Q3C_RADIAL_QUERY(tr.ra,tr.dec,{ra:.6f},{dec:+.6f},{search_radius_deg:.6f}) AND tr.shape_r > {shape_r}
                """
        cat_filename = f'{identifier}_shape_{shape_r}.ls_dr10.cat.fits'
        if os.path.exists( cat_filename ) and not overwrite:
            return Table.read( cat_filename )
    if shape_r != None and excess_factor != None:
        ex_query = f"""
                SELECT tr.ls_id, tr.ra, tr.dec, tr.ref_cat, tr.ref_id, tr.gaia_phot_bp_rp_excess_factor, tr.parallax,
                    tr.mag_g, tr.mag_r, tr.mag_i, tr.mag_z, tr.mag_w1, tr.mag_w2, 
                    tr.sersic, tr.sersic_ivar,  tr.shape_e1, tr.shape_e1_ivar, tr.shape_e2, tr.shape_e2_ivar, 
                    tr.shape_r,  tr.shape_r_ivar, tr.gaia_duplicated_source, zp.z_spec, zp.survey, zp.z_phot_mean, zp.z_phot_mean_i, 
                    zp.z_phot_std, zp.z_phot_l68, zp.z_phot_u68, zp.z_phot_l95, zp.z_phot_u95
                FROM ls_dr10.tractor as tr JOIN ls_dr10.photo_z as zp ON tr.ls_id = zp.ls_id
                WHERE 't' = Q3C_RADIAL_QUERY(tr.ra,tr.dec,{ra:.6f},{dec:+.6f},{search_radius_deg:.6f}) AND tr.shape_r > {shape_r} AND tr.gaia_phot_bp_rp_excess_factor > {excess_factor}
                """
        cat_filename = f'{identifier}_shape_{shape_r}_excessfactor_{excess_factor}.ls_dr10.cat.fits'
        if os.path.exists( cat_filename ) and not overwrite:
            return Table.read( cat_filename )
    if shape_r != None and excess_factor != None and parallax != None:
        ex_query = f"""
                SELECT tr.ls_id, tr.ra, tr.dec, tr.ref_cat, tr.ref_id, tr.gaia_phot_bp_rp_excess_factor, tr.parallax,
                    tr.mag_g, tr.mag_r, tr.mag_i, tr.mag_z, tr.mag_w1, tr.mag_w2, 
                    tr.sersic, tr.sersic_ivar,  tr.shape_e1, tr.shape_e1_ivar, tr.shape_e2, tr.shape_e2_ivar, 
                    tr.shape_r,  tr.shape_r_ivar, tr.gaia_duplicated_source, zp.z_spec, zp.survey, zp.z_phot_mean, zp.z_phot_mean_i, 
                    zp.z_phot_std, zp.z_phot_l68, zp.z_phot_u68, zp.z_phot_l95, zp.z_phot_u95
                FROM ls_dr10.tractor as tr JOIN ls_dr10.photo_z as zp ON tr.ls_id = zp.ls_id
                WHERE 't' = Q3C_RADIAL_QUERY(tr.ra,tr.dec,{ra:.6f},{dec:+.6f},{search_radius_deg:.6f}) AND tr.shape_r > {shape_r} AND tr.gaia_phot_bp_rp_excess_factor > {excess_factor} AND tr.parallax < {parallax}
                """
        cat_filename = f'{identifier}_shape_{shape_r}_excessfactor_{excess_factor}_parallax_{parallax}.ls_dr10.cat.fits'
        if os.path.exists( cat_filename ) and not overwrite:
            return Table.read( cat_filename )
    if shape_r == None and excess_factor == None and parallax == None:
        ex_query = f"""
                SELECT tr.ls_id, tr.ra, tr.dec, tr.ref_cat, tr.ref_id, tr.gaia_phot_bp_rp_excess_factor, tr.parallax, tr.parallax_ivar, 
                    tr.mag_g, tr.mag_r, tr.mag_i, tr.mag_z, tr.mag_w1, tr.mag_w2, 
                    tr.sersic, tr.sersic_ivar,  tr.shape_e1, tr.shape_e1_ivar, tr.shape_e2, tr.shape_e2_ivar, 
                    tr.shape_r,  tr.shape_r_ivar, tr.gaia_duplicated_source, zp.z_spec, zp.survey, zp.z_phot_mean, zp.z_phot_mean_i, 
                    zp.z_phot_std, zp.z_phot_l68, zp.z_phot_u68, zp.z_phot_l95, zp.z_phot_u95
                FROM ls_dr10.tractor as tr JOIN ls_dr10.photo_z as zp ON tr.ls_id = zp.ls_id
                WHERE 't' = Q3C_RADIAL_QUERY(tr.ra,tr.dec,{ra:.6f},{dec:+.6f},{search_radius_deg:.6f})
                """
        cat_filename = f'{identifier}_without_restrictions.ls_dr10.cat.fits'
        if os.path.exists( cat_filename ) and not overwrite:
            return Table.read( cat_filename )
    #    AND tr.parallax < 5
    # NB. uses gaia 'excess factor' to exclude gaia point sources (i.e. stars)
    # Can effectively exclude stars using
    # AND (gaia_phot_bp_rp_excess_factor > 3.5 or gaia_phot_bp_rp_excess_factor <= 0)
    #         AND shape_e1 > 0 AND shape_e2 > 0  --- to exclude saturated stars
    #print( ex_query )

    result = tap_service.search(ex_query)
    tblResult = result.to_table()
    print(f'query result has {len(tblResult)} rows.' )

    print('writing query result to:', cat_filename )
    tblResult.write(cat_filename, overwrite=overwrite )

    return tblResult

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

