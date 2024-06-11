from astropy.table import Table, vstack, unique,join
from astropy import coordinates as co
import astropy.units as u
import numpy as np

#from plots import plot_color_mass

# Load the combined table
result_table = Table.read('crossmatches/specifics/colour_SFR_matched465.fits')

cross_matched_table = Table(names=result_table.colnames + ['dDLR'], dtype=[result_table[col].dtype for col in result_table.colnames] + ['f8'])

# Group the table by SN_ID
#sn_ids = unique(result_table, keys='objid')['objid']
#print('works')
#uber_ids = []
#dDLRs = []

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
