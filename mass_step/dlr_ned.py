################################################################### NED's CODE ####################################################################################
import os, sys, time
import pyvo
from astropy.table import Table
from astropy import coordinates as co, units as u

from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
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
    #separation = search_radec.separation(galaxy_radec).arcsec
    #print(np.shape(separation))
    #ind = np.where(separation<=15)
    #print(separation.min())

    #print(ind)

    # # work out the angle from search pos to galaxy
    # sep_angle = 90 - search_radec.position_angle(galaxy_radec).deg

    # work out the RA/dec offsets from search pos to galaxy
    alpha, delta = search_radec.spherical_offsets_to(galaxy_radec)
    #print(np.shape(search_radec))
    #print(type(search_radec))
    #print(type(galaxy_radec))
    # choose to work in units of arcsec for convenience
    dx, dy = alpha.arcsec, delta.arcsec

    # now rotate these dx/dy to match the galaxy angle:
    cosphi, sinphi = np.cos(galaxy_angle_rad), np.sin(galaxy_angle_rad)
    # this is just the ordinary rotation matrix
    da, db = cosphi*dx - sinphi*dy, sinphi*dx + cosphi*dy
    # note this is mathematicians and not astronomers angle
    #print(np.shape(da))
    # now da and db are separations in arcsec
    # in the coord system of the semimaj/minor axes.

    dDLR = np.sqrt((da / galaxy_major)**2. + (db / galaxy_minor)**2. )

    return dDLR, da, db

#### DATA ####

result_table = Table.read('data/gkvScienceCatv02.fits')

RA_gama = result_table['RAcen']
Dec_gama = result_table['Deccen']
ID = result_table['CATAID']
axrat = result_table['axrat']
position_angle = result_table['ang']
major = result_table['R50']

#### DATA SELECTION ####
ra_min = 176
ra_max = 184
dec_min = 1#1.625
dec_max = 2#1.635

ind = np.where((RA_gama>ra_min)&(RA_gama<ra_max)&(Dec_gama>1.625)&(Dec_gama<1.635)&(ID!=0)&(axrat>0.0)&(axrat<0.5))

ra_galaxies = RA_gama[ind]#[180.,182]
dec_galaxies = Dec_gama[ind]#[0.,1.]

position_angles = position_angle[ind] - 90 #! Corrected from GAMA to DLR
major_axis = major[ind]
axrat = axrat[ind]
minor_axis = axrat*major_axis # axial ratio = minor/major -> minor = axial ratio * major

#### Random grid SN ####

ra_grid, dec_grid = np.meshgrid(np.linspace(ra_min, ra_max, 200), np.linspace(dec_min, dec_max, 200))

#print(np.shape(ra_grid))


grid_coords = co.SkyCoord(ra=ra_grid, dec=dec_grid, unit='deg')

closest_galaxy = np.full(ra_grid.shape, -1, dtype=int)  # -1 indicates "outside everything"
min_dDLR = np.full(ra_grid.shape, np.inf)
#closest_galaxy = np.zeros(ra_grid.shape, dtype=int)

for i, (ra, dec) in enumerate(zip(ra_galaxies, dec_galaxies)):
    #print(i)
    position_angle = position_angles[i]

    galaxy_major, galaxy_minor = major_axis[i], minor_axis[i]
    galaxy_angle = np.radians(position_angle)
    galaxy_radec = co.SkyCoord(ra=ra, dec=dec, unit='deg')

    dDLR_values,da,db = dDLR(grid_coords, galaxy_radec, galaxy_major, galaxy_minor, galaxy_angle)
    #print(np.shape(da))
    #print(np.shape(dDLR))
    angular_separation = np.sqrt((da**2)+(db**2))

    mask = angular_separation<15
    valid_mask = dDLR_values<min_dDLR
    combined_mask = mask & valid_mask
    closest_galaxy[combined_mask] = i
    min_dDLR[combined_mask] = dDLR_values[combined_mask]


# Plot the results
#plt.figure(1) ; plt.clf()
#plt.gca().set_aspect('equal')
#colors = plt.cm.tab20(np.linspace(0, 1, len(ra_galaxies)))

# Convert closest_galaxy indices to colors for imshow
image = np.zeros(closest_galaxy.shape + (3,))
colors = plt.cm.tab20(np.linspace(0, 1, len(ra_galaxies)))

# Plot grid points not associated with any galaxy in grey
#mask = closest_galaxy == -1
#plt.scatter(ra_grid[mask], dec_grid[mask], color='grey', s=10)

# Plot grid points with associated galaxies
#for i in range(len(ra_galaxies)):
#    mask = closest_galaxy == i
#    plt.scatter(ra_grid[mask], dec_grid[mask], color=colors[i], s=10)

for i in range(len(ra_galaxies)):
    mask = closest_galaxy == i
    image[mask] = colors[i][:3]

# Grey color for unassociated SN
mask = closest_galaxy == -1
image[mask] = [0.5, 0.5, 0.5]

# Plotting with imshow
plt.figure(figsize=(10, 8))
plt.imshow(image, origin='lower', extent=[ra_min, ra_max, dec_min, dec_max])


# Plot galaxies as ellipses
for i, (ra, dec) in enumerate(zip(ra_galaxies, dec_galaxies)):
    galaxy_major = major_axis[i]
    galaxy_minor = minor_axis[i]

    angle = position_angles[i]

    # Adjust the angle to counteract the x-axis inversion
    adjusted_angle = 180 - angle
    # Create an ellipse patch
    ellipse = Ellipse(xy=(ra, dec), width=galaxy_major*2/3600, height=galaxy_minor*2/3600, angle=adjusted_angle, edgecolor='black', facecolor=colors[i], alpha=0.5)
    #print(angle)

    plt.gca().add_patch(ellipse)



plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
#plt.legend()
plt.title('Galaxies and Associated Supernovae')

plt.gca().invert_xaxis()
plt.show()

exit()



# Plot the SN points colored by the closest galaxy
#plt.imshow( closest_galaxy, extent=(ra_max,ra_min,dec_min,dec_max), cmap='tab10', aspect='auto' )
plt.scatter(ra_grid, dec_grid, c=closest_galaxy, cmap='tab10', s=1, alpha=0.5)


# Plot the galaxies
for i in range(len(ra_galaxies)):
    plt.scatter(ra_galaxies[i], dec_galaxies[i], s=50, edgecolors='black')

#plt.colorbar(label='Closest Galaxy')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.xlim(ra_max, ra_min)
plt.ylim(dec_min, dec_max)
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