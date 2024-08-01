from astropy.table import Table, vstack, unique,join
from astropy import coordinates as co
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import style

plt.style.use(style.style1)

np.random.seed(23)
#from plots import plot_color_mass

# Load the combined table
#result_table = Table.read('crossmatches/specifics/colour_SFR_matched465.fits')
result_table = Table.read('data/gkvScienceCatv02.fits')
cross_matched_table = Table(names=result_table.colnames + ['dDLR'], dtype=[result_table[col].dtype for col in result_table.colnames] + ['f8'])

# Group the table by SN_ID
#sn_ids = unique(result_table, keys='objid')['objid']
#print('works')
#uber_ids = []
#dDLRs = []

RA_gama = result_table['RAcen']
Dec_gama = result_table['Deccen']
ID = result_table['CATAID']
axrat = result_table['axrat']
#ind_ra = np.where((RA_gama>180)&(RA_gama<220))
#ind_dec = np.where((Dec_gama>-1.5)&(Dec_gama<1.5))

# Selection of galaxies:
#print(RA_gama.min(), RA_gama.max())
#print(Dec_gama.min(),Dec_gama.max())
ra_min = 176#179.3 #176.5 #179
ra_max = 184 #179.82 #177.5 #186.5
dec_min = 1.625 #-3.05 #-1
dec_max = 1.635 #1.820 #-3.02 #1

ind = np.where((RA_gama>ra_min)&(RA_gama<ra_max)&(Dec_gama>dec_min)&(Dec_gama<dec_max)&(ID!=0)&(axrat>0.0)&(axrat<0.5))

#ind = np.where((ID==288679)|(ID == 550086)|(ID==219946))
#print(np.shape(ind))

# Scatter of the selection of galaxies:

#ax1.scatter(RA_gama[ind], Dec_gama[ind], s=0.01)

#plt.show()

major_axis = result_table['R50'][ind]/10 #! I would need to do: arcsec -> degree, but too small
axrat = result_table['axrat'][ind]
minor_axis = axrat*major_axis # axial ratio = minor/major -> minor = axial ratio * major
#plt.errorbar(x = RA_gama[ind], y=Dec_gama[ind], xerr=major_axis, yerr=minor_axis, fmt='.')

position_angle = result_table['ang'][ind]

def plot_galaxies(coord1, coord2, semi_major_axes, semi_minor_axes, angles):
    fig1, ax1 = plt.subplots()
    ax1.set_title('Scatter of the galaxies')
    ax1.set(xlabel = ('RA'), ylabel = ('Dec'))#, xlim = (ra_min-0.1, ra_max+0.1))#, ylim = (dec_min-0.05,dec_max+0.05))#, yscale = ('log'))#, xlim = (0,6))
    ax1.set_aspect('equal')
    ax1.invert_xaxis()

    ######## The random SN ########

    # Generate random RA and Dec values within the specified ranges
    SN_ra = np.random.uniform(ra_min, ra_max,2)
    SN_dec = np.random.uniform(dec_min, dec_max,2)
    #ax1.plot(SN_ra, SN_dec, 'o', color = 'blue')

    # Plot each galaxy
    for i in range(len(coord1)):
        x = coord1[i]
        y = coord2[i]
        semi_major = semi_major_axes[i]
        semi_minor = semi_minor_axes[i]
        angle = angles[i]
        print(angle,x,y)

        # Plot the center of the galaxy
        ax1.plot(x, y, 'o', color='black')

        # Calculate the end points of the semi-major and semi-minor axes
        theta = np.radians(angle)
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        major_x1 = x + semi_major * cos_theta
        major_y1 = y + semi_major * sin_theta
        major_x2 = x - semi_major * cos_theta
        major_y2 = y - semi_major * sin_theta

        minor_x1 = x + semi_minor * (-sin_theta)
        minor_y1 = y + semi_minor * cos_theta
        minor_x2 = x - semi_minor * (-sin_theta)
        minor_y2 = y - semi_minor * cos_theta

        ############ THE GOOD ONE ############################
        major_x1 = x + semi_major * sin_theta
        major_y1 = y + semi_major * cos_theta
        major_x2 = x - semi_major * sin_theta
        major_y2 = y - semi_major * cos_theta

        minor_x1 = x + semi_minor * (-cos_theta)
        minor_y1 = y + semi_minor * sin_theta
        minor_x2 = x - semi_minor * (-cos_theta)
        minor_y2 = y - semi_minor * sin_theta
        ########################################################

        #ax1.plot([major_x1, major_x2], [major_y1, major_y2], '-', color='red', zorder=100)
        #ax1.plot([minor_x1, minor_x2], [minor_y1, minor_y2], '-', color='black',zorder=100)

        '''

        major_x1 = x + semi_major * cos_theta
        major_y1 = y + semi_major * sin_theta
        major_x2 = x - semi_major * cos_theta
        major_y2 = y - semi_major * sin_theta

        minor_x1 = x + semi_minor * (-sin_theta)
        minor_y1 = y + semi_minor * cos_theta
        minor_x2 = x - semi_minor * (-sin_theta)
        minor_y2 = y - semi_minor * cos_theta
        '''
        '''
        major_x1 = x - semi_major * sin_theta
        major_y1 = y + semi_major * cos_theta
        major_x2 = x + semi_major * sin_theta
        major_y2 = y - semi_major * cos_theta

        minor_x1 = x - semi_minor * (cos_theta)
        minor_y1 = y + semi_minor * -sin_theta
        minor_x2 = x + semi_minor * (cos_theta)
        minor_y2 = y - semi_minor * -sin_theta
        '''
        # Plot the semi-major and semi-minor axes
        ax1.plot([major_x1, major_x2], [major_y1, major_y2], '-', color='red', label = 'major')
        ax1.plot([minor_x1, minor_x2], [minor_y1, minor_y2], '-', color='black', label = 'minor')

    #exit()

    radec1 = co.SkyCoord(coord1, coord2, unit='deg')
    radec2 = co.SkyCoord(SN_ra, SN_dec, unit='deg')
    for i, SN_point in enumerate(radec2):
        distances = SN_point.separation(radec1)
        ind_min = distances.argmin()
        print(f"\nDistances from SN {i+1} to Galaxies (in arcseconds):")
        print(distances.arcsec)

        delta_ra = [SN_point.ra.deg - radec1.ra.deg[j] for j, d in enumerate(distances)]
        delta_dec = [SN_point.dec.deg - radec1.dec.deg[j] for j, d in enumerate(distances)]

        angle_from_ra = np.degrees(np.arctan2(delta_dec, delta_ra))
        galaxy_angle = angles[ind_min]

        angle_from_major = 180 + (angle_from_ra[ind_min] - galaxy_angle) # The angle with the SN at the smallest distance

        angles_from_major = 180 + (angle_from_ra - angles)

        DLR = (semi_major_axes*semi_minor_axes)/(np.sqrt((semi_major_axes*np.sin(angles_from_major))**2+(semi_minor_axes*np.cos(angles_from_major))**2))
        dDLR = distances.arcsec/DLR
        ind_min_dlr = dDLR.argmin()

        #for j, d in enumerate(distances):
            #ax1.plot([SN_point.ra.deg,radec1.ra.deg[j]], [SN_point.dec.deg,radec1.dec.deg[j]], '-', color = 'pink')
            #ax1.plot([SN_point.ra.deg,radec1.ra.deg[ind_min_dlr]], [SN_point.dec.deg,radec1.dec.deg[ind_min_dlr]], linestyle='dashed', color = 'green')

        #ax1.plot([radec1.ra.deg[ind_min] - 5, radec1.ra.deg[ind_min] + 5], [radec1.dec.deg[ind_min], radec1.dec.deg[ind_min]], color='black', linestyle=':')

        x_mid = (radec1.ra.deg[ind_min] + SN_point.ra.deg) / 2
        y_mid = (radec1.dec.deg[ind_min] + SN_point.dec.deg) / 2
        #ax1.annotate(f'{angle_from_ra:.1f}°', xy=(x_mid, y_mid), color='purple', fontsize=12) #plot the angle from RA
        #ax1.annotate(f'{angle_from_major:.1f}°', xy=(x_mid, y_mid), color='red', fontsize=12) #plot the angle from major axis

    #plt.legend()


plot_galaxies(RA_gama[ind],Dec_gama[ind], major_axis, minor_axis, position_angle)
plot_galaxies([RA_gama[ind][0]],[Dec_gama[ind][0]], major_axis, minor_axis, [130])
print(ID[ind], RA_gama[ind], Dec_gama[ind], position_angle)
plt.show()
exit()
################################################################### NED's CODE ####################################################################################
import os, sys, time
import pyvo
from astropy.table import Table
from astropy import coordinates as co, units as u

from matplotlib import pyplot as plt
import numpy as np

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

#cat = get_lsdr10_cat( 'deleteme', 179.7174, 1.728273, 0.1 )
#plt.clf()
#plt.scatter( cat['ra'], cat['dec'], 4 )
#plt.show()
#exit()


############################################################################

def dDLR(search_radec, galaxy_radec, galaxy_major, galaxy_minor, galaxy_angle_rad):
    # # work out the separation between search and galaxy coordinates
    # separation = search_radec.separation(galaxy_radec).arcsec
    # # work out the angle from search pos to galaxy
    # sep_angle = 90 - search_radec.position_angle(galaxy_radec).deg 

    # work out the RA/dec offsets from search pos to galaxy
    alpha, delta = search_radec.spherical_offsets_to( galaxy_radec )
    # choose to work in units of arcsec for convenience
    dx, dy = alpha.arcsec, delta.arcsec

    # now rotate these dx/dy to match the galaxy angle:
    # todo check sign!
    cosphi, sinphi = np.cos(galaxy_angle_rad), np.sin(galaxy_angle_rad)
    # this is just the ordinary rotation matrix
    # todo: check sign!!
    da, db = cosphi*dx - sinphi*dy, sinphi*dx + cosphi*dy
    # note this is mathematicians and not astronomers angle

    # now da and db are separations in arcsec 
    # in the coord system of the semimaj/minor axes.
    
    dDLR = np.sqrt((da / galaxy_major)**2. + (db / galaxy_minor)**2. )

    return dDLR, da, db


ra_galaxies = RA_gama[ind]#[180.,182]
dec_galaxies = Dec_gama[ind]#[0.,1.]
#sizes = [4.,5]
#boveras = [0.8,0.7]
#galaxy_major, galaxy_minor = major_axis, minor_axis #size, bovera*size
position_angles = position_angle
#galaxy_angle = np.radians(position_angle - 90) #np.radians( 0. ) # should be in radians
# note this is currently mathematicians not astronomers angle

#galaxy_radec = co.SkyCoord(ra=ra, dec=dec, unit='deg')#co.SkyCoord(ra=RA_gama[ind], dec=Dec_gama[ind], unit='deg')

#dx, dy = np.meshgrid( np.arange( -20, 20.1, 1 ), 
#                      np.arange( -20, 20.1, 1 ), ) #
ra_grid, dec_grid = np.meshgrid(np.linspace(ra_min, ra_max, 100), np.linspace(dec_min, dec_max, 100))
#dx,dy = np.meshgrid(np.linspace(176, 184, 41 ), np.linspace(1.625, 1.635, 41 ), ) #np.meshgrid( np.arange( -20, 20.1, 1 ), np.arange( -20, 20.1, 1 ), )
#print(np.shape(dx), np.shape(dy))
grid_coords = co.SkyCoord(ra=ra_grid, dec=dec_grid, unit='deg')


closest_galaxy = np.zeros(ra_grid.shape, dtype=int)

for i, (ra, dec) in enumerate(zip(ra_galaxies, dec_galaxies)):
    #size = sizes[i]
    #bovera = boveras[i]
    position_angle = position_angles[i]

    galaxy_major, galaxy_minor = major_axis[i], minor_axis[i]
    galaxy_angle = np.radians(position_angle - 90)
    galaxy_radec = co.SkyCoord(ra=ra, dec=dec, unit='deg')

    dDLR_values,da,db = dDLR(grid_coords, galaxy_radec, galaxy_major, galaxy_minor, galaxy_angle)
    print(np.shape(dDLR_values))
    if i == 0:
        min_dDLR = dDLR_values
    else:
        mask = dDLR_values < min_dDLR
        closest_galaxy[mask] = i
        min_dDLR[mask] = dDLR_values[mask]

# Plot the results
plt.figure()
plt.gca().invert_xaxis()

# Plot the grid points colored by the closest galaxy
plt.scatter(ra_grid, dec_grid, c=closest_galaxy, cmap='tab10', s=10, alpha=0.5)

# Plot the galaxies
for i in range(len(ra_galaxies)):
    plt.scatter(ra_galaxies[i], dec_galaxies[i], s=50)

plt.colorbar(label='Closest Galaxy')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.legend()
plt.title('Closest Galaxy by dDLR')
plt.show()

exit()

for i in range(len(ra_galaxies)):
    ra = ra_galaxies[i]
    dec = dec_galaxies[i]
    size = sizes[i]
    bovera = boveras[i]
    position_angle = position_angles[i]

    galaxy_major, galaxy_minor = size, bovera * size
    galaxy_angle = np.radians(position_angle - 90)

    galaxy_radec = co.SkyCoord(ra=ra, dec=dec, unit='deg')
    search_radec = galaxy_radec.spherical_offsets_by(-dx * u.arcsec, dy * u.arcsec)

    result, da, db = dDLR(search_radec, galaxy_radec, galaxy_major, galaxy_minor, galaxy_angle)

    plt.scatter(search_radec.ra.deg, search_radec.dec.deg, 10, result, label=f'Galaxy {i+1}')
    plt.scatter(galaxy_radec.ra.deg, galaxy_radec.dec.deg, 50, label=f'Galaxy {i+1} center')

plt.colorbar(label='dDLR')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.legend()
plt.title('DLR Contours for Multiple Galaxies')
plt.show()

exit()

search_radec = galaxy_radec.spherical_offsets_by( -dx*u.arcsec, dy*u.arcsec) 






result, da, db = dDLR( search_radec, galaxy_radec, galaxy_major, galaxy_minor, galaxy_angle )

plt.scatter( search_radec.ra.deg, search_radec.dec.deg, 10, result )
plt.colorbar()
plt.scatter( galaxy_radec.ra.deg, galaxy_radec.dec.deg, 50, label='galaxy' )

plt.gca().invert_xaxis()
plt.legend()
plt.show()
#######################################################################################################################################################################

exit()
sn_ids = result_table['objid']
# Print the separations for each supernova
for i, sn_id in enumerate(sn_ids):
    print(i, sn_id)
    continue
    exit()
    sn_mask = result_table['objid'] == sn_id
    uberID = result_table['uberID'][sn_mask]
    separations = result_table['Separation_arcsec'][sn_mask]

    major_axis = result_table['R50'][sn_mask] #arcsec
    axrat = result_table['axrat'][sn_mask]
    minor_axis = axrat*major_axis # axial ratio = minor/major -> minor = axial ratio * major

    position_angle = result_table['ang'][sn_mask]

    DLR = (major_axis*minor_axis)/(np.sqrt((major_axis*np.sin(position_angle))**2+(minor_axis*np.cos(position_angle))**2 ))

    dDLR = separations/DLR
    print(dDLR)
    result_table2 = result_table[sn_mask]
    #result_table2['dDLR'][i] = dDLR
    result_table3 = result_table[sn_mask]
    #result_table3['dDLR'][i] = dDLR.min()

    ind = np.argmin(dDLR)
    uberID = uberID[ind]
    uber_ids.append(uberID)
    dDLR = dDLR[ind]
    dDLRs.append(dDLR)


#print(len(uber_ids))
#print(len(dDLRs))
#print(len(result_table3))
exit()
#mask = [uber_id in uber_ids for uber_id in result_table['uberID']]
#reduced_table = result_table[mask]
#print(len(reduced_table))

for uber_id, dDLR in zip(uber_ids, dDLRs):
    cross_matched_table.add_row(result_table2[result_table2['uberID'] == uber_id][0])
    cross_matched_table[-1]['dDLR'] = dDLR
print(cross_matched_table)
exit()
reduced_table['dDLR'] = dDLRs
exit()
# For the colours and masses:
gama_masses_file = "data/StellarMassesGKVv24.fits"
mass_cat = Table.read(gama_masses_file)
selected_mass = mass_cat[['uberID', 'Z', 'logmstar', 'dellogmstar', 'mstar', 'delmstar', 'gminusi', 'delgminusi', 'gminusi_stars', 'delgminusi_stars', 'uminusr', 'deluminusr']]
selected_matched = reduced_table[['uberID', 'ra', 'declination', 'redshift','Separation_arcsec', 'dDLR']]

colour_mass_cat = selected_mass
matched_color_mass_cat = join(selected_mass, selected_matched, keys='uberID')
matched_color_mass_cat.write('crossmatches/dlr_colour_mass.fits', overwrite = True)
