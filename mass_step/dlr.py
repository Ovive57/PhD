from astropy.table import Table, vstack, unique,join
from astropy import coordinates as co
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import style

plt.style.use(style.style1)

np.random.seed(5)
#from plots import plot_color_mass

# Load the combined table
result_table = Table.read('crossmatches/specifics/colour_SFR_matched465.fits')

cross_matched_table = Table(names=result_table.colnames + ['dDLR'], dtype=[result_table[col].dtype for col in result_table.colnames] + ['f8'])

# Group the table by SN_ID
#sn_ids = unique(result_table, keys='objid')['objid']
#print('works')
#uber_ids = []
#dDLRs = []

RA_gama = result_table['RAcen']
Dec_gama = result_table['Deccen']
ID = result_table['CATAID']
ind_ra = np.where((RA_gama>180)&(RA_gama<220))
ind_dec = np.where((Dec_gama>-1.5)&(Dec_gama<1.5))

# Selection of galaxies:
ind = np.where((RA_gama>179)&(RA_gama<203)&(Dec_gama>-1)&(Dec_gama<1))

# Scatter of the selection of galaxies:

#ax1.scatter(RA_gama[ind], Dec_gama[ind], s=0.01)

#plt.show()

major_axis = result_table['R50'][ind] #! I would need to do: arcsec -> degree, but too small
axrat = result_table['axrat'][ind]
minor_axis = axrat*major_axis # axial ratio = minor/major -> minor = axial ratio * major
#plt.errorbar(x = RA_gama[ind], y=Dec_gama[ind], xerr=major_axis, yerr=minor_axis, fmt='.')

position_angle = result_table['ang'][ind]

def plot_galaxies(coord1, coord2, semi_major_axes, semi_minor_axes, angles):
    fig1, ax1 = plt.subplots()
    ax1.set_title('Scatter of the galaxies')
    ax1.set(xlabel = ('RA'), ylabel = ('Dec'))#, yscale = ('log'))#, xlim = (0,6))
    ax1.set_aspect('equal')
    ax1.invert_xaxis()

    ######## The random SN ########
    ra_min = 179
    ra_max = 186.5
    dec_min = -1
    dec_max = 1

    # Generate random RA and Dec values within the specified ranges
    SN_ra = np.random.uniform(ra_min, ra_max,2)
    SN_dec = np.random.uniform(dec_min, dec_max,2)
    ax1.plot(SN_ra, SN_dec, 'o', color = 'blue')

    # Plot each galaxy
    for i in range(len(coord1)):
        x = coord1[i]
        y = coord2[i]
        semi_major = semi_major_axes[i]
        semi_minor = semi_minor_axes[i]
        angle = angles[i]

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

        #ax1.plot([major_x1, major_x2], [major_y1, major_y2], '-', color='red', zorder=100)
        #ax1.plot([minor_x1, minor_x2], [minor_y1, minor_y2], '-', color='black',zorder=100)

        cos_theta, sin_theta = np.cos(theta+90), np.sin(theta+90)

        major_x1 = x + semi_major * cos_theta
        major_y1 = y + semi_major * sin_theta
        major_x2 = x - semi_major * cos_theta
        major_y2 = y - semi_major * sin_theta

        minor_x1 = x + semi_minor * (-sin_theta)
        minor_y1 = y + semi_minor * cos_theta
        minor_x2 = x - semi_minor * (-sin_theta)
        minor_y2 = y - semi_minor * cos_theta
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


    radec1 = co.SkyCoord(coord1, coord2, unit='deg')
    radec2 = co.SkyCoord(SN_ra, SN_dec, unit='deg')
    for i, SN_point in enumerate(radec2):
        distances = SN_point.separation(radec1)
        ind_min = distances.argmin()
        #print(f"\nDistances from SN {i+1} to Galaxies (in arcseconds):")
        #print(distances.arcsec)

        delta_ra = SN_point.ra.deg - radec1.ra.deg[ind_min]
        delta_dec = SN_point.dec.deg - radec1.dec.deg[ind_min]

        angle_from_ra = np.degrees(np.arctan2(delta_dec, delta_ra))

        galaxy_angle = angles[ind_min]
        angle_from_major = (angle_from_ra - galaxy_angle)

        ax1.plot([SN_point.ra.deg,radec1.ra.deg[ind_min]], [SN_point.dec.deg,radec1.dec.deg[ind_min]], '-', color = 'green')
        ax1.plot([radec1.ra.deg[ind_min] - 0.5, radec1.ra.deg[ind_min] + 0.5], [radec1.dec.deg[ind_min], radec1.dec.deg[ind_min]], color='black', linestyle='dashed')

        x_mid = (radec1.ra.deg[ind_min] + SN_point.ra.deg) / 2
        y_mid = (radec1.dec.deg[ind_min] + SN_point.dec.deg) / 2
        #ax1.annotate(f'{angle_from_ra:.1f}°', xy=(x_mid, y_mid), color='purple', fontsize=12)
        ax1.annotate(f'{angle_from_major:.1f}°', xy=(x_mid, y_mid), color='red', fontsize=12)

    #plt.legend()


plot_galaxies(RA_gama[ind],Dec_gama[ind], major_axis, minor_axis, position_angle)
#plot_galaxies([RA_gama[ind][0]],[Dec_gama[ind][0]], major_axis, minor_axis, [130])
print(ID[ind])
plt.show()
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
