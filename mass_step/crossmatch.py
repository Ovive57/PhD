import numpy as np
import matplotlib.pyplot as plt
import style
from astropy.table import Table, hstack, vstack, join
from astropy import coordinates as co, units as u
from astropy.cosmology import Planck18, FlatLambdaCDM
import time


plt.style.use(style.style1)

gama_file = "data/gkvScienceCatv02.fits"
ztf_file = "data/tns_SNIa_20240424_copy.csv"

'''
ztf_ra = np.char.strip(np.loadtxt(ztf_file, dtype=str, skiprows=1, usecols=(2), unpack=True, delimiter=","),'"') #hourangle
ztf_dec = np.char.strip(np.loadtxt(ztf_file, dtype=str, skiprows=1, usecols=(3), unpack=True, delimiter=","),'"') #hourangle
ztf_ra1 = ztf_ra[0]
ztf_dec1 = ztf_dec[0]
'''

#https://docs.astropy.org/en/stable/coordinates/index.html

def separations(cat_file1, cat_file2, max_sep, frame = 'icrs', test = False):
    start = time.time()
    cat1 = Table.read(cat_file1)
    cat2 = Table.read(cat_file2)
    radec1 = co.SkyCoord(cat1['RAcen'], cat1['Deccen'], frame = frame, unit='deg') # GAMA (ICRS): (ra, dec) in deg
    radec2 = co.SkyCoord(cat2['ra']*u.deg, cat2['declination']*u.deg, frame=frame) # ZTF (ICRS): (ra, dec) in deg

    max_sep = max_sep*u.arcsec

    result_table = Table(names=cat1.colnames + ['objid', 'Separation_arcsec'], dtype=[cat1[col].dtype for col in cat1.colnames] + ['i8', 'f8'])

    for i, c2 in enumerate(radec2):
        sep = radec1.separation(c2)

        mask = sep < max_sep
        nearby_galaxies = cat1[mask]
        if len(nearby_galaxies) > 0:
            sn_id = cat2[i]['objid']
            separations = sep[mask].arcsec
            nearby_galaxies['Separation_arcsec'] = separations
            nearby_galaxies['objid'] = sn_id  # Add the SN_ID to the nearby galaxies table
            result_table = vstack([result_table, nearby_galaxies])

        if test:
            if i>300:
                break
    end = time.time()
    print(end-start)
    result_table_final = join(result_table, cat2, keys='objid')
    result_table_final.write('nearby_galaxies_all_supernovae_with.fits', overwrite=True)

#separations(gama_file, ztf_file, max_sep=465, test = False)
#exit()

def positional_crossmatch(cat_file1, cat_file2, dim, max_sep, max_sep3d=None, cosmology = Planck18,frame = 'icrs', verbose = True):

    cat1 = Table.read(cat_file1)
    cat2 = Table.read(cat_file2)


    # Masks:
    if dim == 2:
        radec1 = co.SkyCoord(cat1['RAcen'], cat1['Deccen'], frame = frame, unit='deg') # GAMA (ICRS): (ra, dec) in deg
        radec2 = co.SkyCoord(cat2['ra']*u.deg, cat2['declination']*u.deg, frame=frame) # ZTF (ICRS): (ra, dec) in deg
        a2b, sep, _ = radec2.match_to_catalog_sky(radec1) # d2d : on-sky distances
        sep_constraint = sep.arcsec < max_sep
        a2b, b2a, sep = a2b[sep_constraint], np.where(sep_constraint)[0], sep[sep_constraint]
        #* a2b are the index of the matches in the first catalogue
        #* b2a are the index of the matches in the second catalogue
        #* sep is the separation between the matches
        if verbose:
            print(f"There are {len(a2b)} matches")
            print("Galaxies: ", radec1[a2b])
            print("SN: ", radec2[b2a])
            print("Separation: in arcsec \n", sep.arcsec)
        match = hstack((cat1[a2b], cat2[b2a]))
        match['sep_arcsec'] = sep.arcsec

        output_path = f"crossmatches/crossmatch_GAMA_ZTF_maxsep_{max_sep}arcsec_dim{dim}.fits"
        match.write(output_path, overwrite = True)
        return match
    if dim ==3:
        z1 = cat1['Z']
        ind1 = np.where(z1>=0)
        comoving_distance1 = cosmology.comoving_distance(z1[ind1]).value
        comoving_distance_kpc1 = comoving_distance1*1e3 #Mpc to kpc

        z2 = cat2['redshift']
        ind2 = np.where(z2>=0)
        comoving_distance2 = cosmology.comoving_distance(z2[ind2]).value
        comoving_distance_kpc2 = comoving_distance2*1e3 #Mpc to kpc

        radecz1 = co.SkyCoord(ra=cat1['RAcen'][ind1], dec=cat1['Deccen'][ind1], distance=comoving_distance_kpc1*u.kpc) # GAMA (ICRS): (ra, dec) in deg, distance in kpc
        #radec2 = co.SkyCoord(ra=cat1['RAcen'][ind], dec=cat1['Deccen'][ind], distance = comoving_distance_kpc,unit = (u.deg, u.deg, u.kpc))
        radecz2 = co.SkyCoord(ra=cat2['ra']*u.deg, dec=cat2['declination']*u.deg, distance=comoving_distance_kpc2*u.kpc) # ZTF (ICRS): (ra, dec) in deg
        #radec2 = co.SkyCoord(cat2['ra'], cat2['declination'],unit = (u.deg, u.deg), frame=frame) # ZTF (ICRS): (ra, dec) in deg

        a2b, sep, sep3d = radecz2.match_to_catalog_3d(radecz1)

        #### ON SKY ####
        sep_constraint = sep.arcsec < max_sep

        a2b, b2a, sep = a2b[sep_constraint], np.where(sep_constraint)[0], sep[sep_constraint]
        #* a2b are the index of the matches in the first catalogue
        #* b2a are the index of the matches in the second catalogue
        #* sep is the separation between the matches
        if verbose:
            print(f"There are {len(a2b)} matches on sky with a max separation of {max_sep} arcsec")
            print("Galaxies: ", radecz1[a2b])
            print("SN: ", radecz2[b2a])
            print("Separation: in arcsec \n", sep.arcsec)

        match1 = hstack((cat1[a2b], cat2[b2a]))
        match1['sep_arcsec'] = sep.arcsec

        output_path = f"crossmatches/crossmatch_GAMA_ZTF_maxsep_{max_sep}arcsec_dim{dim}_onsky.fits"
        match1.write(output_path, overwrite = True)

        #### 3D ####
        a2b, sep, sep3d = radecz2.match_to_catalog_3d(radecz1)

        sep_constraint = sep3d.value < max_sep3d
        a2b, b2a, sep = a2b[sep_constraint], np.where(sep_constraint)[0], sep[sep_constraint]
        #* a2b are the index of the matches in the first catalogue
        #* b2a are the index of the matches in the second catalogue
        #* sep is the separation between the matches
        if verbose:
            print(f"There are {len(a2b)} matches in 3D with a max separation of {max_sep3d} kpc")
            print("Galaxies: ", radecz1[a2b])
            print("SN: ", radecz2[b2a])
            print("Separation: in kpc \n", sep3d.value)

        match2 = hstack((cat1[a2b], cat2[b2a]))
        match2['sep_arcsec'] = sep.arcsec

        output_path = f"crossmatches/crossmatch_GAMA_ZTF_maxsep{max_sep}arcsec_maxsep3d_{max_sep3d}kpc_dim{dim}_3D.fits"
        match2.write(output_path, overwrite = True)

        return match1, match2



#positional_crossmatch(gama_file, ztf_file, dim=3, max_sep = 4665.8898, max_sep3d = 100)
#positional_crossmatch(gama_file, ztf_file, dim=2, max_sep = 465)
#exit()

#! I'm ging to do a rough crossmatch between all catalogues, to have a big catalogue with all the information. I propose to create a function with dictionaries
#! to create this kind of catalogue quikly

# For the colours and masses:
gama_masses_file = "data/StellarMassesGKVv24.fits" # 'uberID', 'Z', 'logmstar', 'dellogmstar', 'mstar', 'delmstar', 'gminusi', 'delgminusi', 'uminusr', 'deluminusr'

# For the sersic:
gama_sersic_file_VIKING = "data/SersicCatVIKINGv09.fits" # 'CATAID', 'RA', 'DEC', 'GALINDEX_i', 'GALINDEX_g', 'GALINDEX_u', 'GALINDEX_r'
gama_sersic_file_SDSS = "data/SersicCatSDSSv09.fits" # 'CATAID', 'RA', 'DEC', 'GALINDEX_i', 'GALINDEX_g', 'GALINDEX_u', 'GALINDEX_r'

# For the ID and RA, Dec, Z:
gama_science_file = "data/gkvScienceCatv02.fits" # 'uberID', 'RAcen', 'Deccen', 'CATAID', 'Z'

# For the sSFR:
gama_SFR_file = "data/MagPhysv06.fits" # 'CATAID', 'Z', 'sSFR_0_1Gyr_best_fit', 'SFR_0_1Gyr_best_fit'


#separations = [1.5]#,5,47,465]
#matched_file = [f"crossmatches/crossmatch_GAMA_ZTF_maxsep_{sep}arcsec_dim2.fits" for sep in separations] # 'uberID', 'RAcen', 'Deccen', 'CATAID', 'Z', 'ra', 'declination', 'redshift'
matched_file = "crossmatches/crossmatch_GAMA_ZTF_maxsep_465arcsec_dim2.fits" # 'uberID', 'RAcen', 'Deccen', 'CATAID', 'Z', 'ra', 'declination', 'redshift'

#! I need 3 catalogues:
    #! plot 1: (colour - mass) + host(colour - mass) EASY, DONE
    #! plot 2: (colour - sersic) + host(colour - sersic) -> gamma_masses_file + gama_sersic_file : but I don't have the same ID, so I add gama_science_file
    #! plot 3: (colour - SFR) + host(colour - SFR) -> gamma_masses_file + gama_SFR_file : but I don't have the same ID, so I add gama_science_file

science_cat = Table.read(gama_science_file)
mass_cat = Table.read(gama_masses_file)
sersic_cat_SDSS = Table.read(gama_sersic_file_SDSS)
sersic_cat_VIKING = Table.read(gama_sersic_file_VIKING)
SFR_cat = Table.read(gama_SFR_file)
matched_cat = Table.read(matched_file)

selected_science = science_cat[['uberID', 'RAcen', 'Deccen', 'CATAID', 'flux_it', 'flux_il', 'flux_err_it', 'flux_gt', 'flux_gl', 'flux_err_gt', 'flux_Zt', 'flux_Zl', 'flux_err_Zt']]
selected_mass = mass_cat[['uberID', 'Z', 'logmstar', 'dellogmstar', 'mstar', 'delmstar', 'gminusi', 'delgminusi', 'gminusi_stars', 'delgminusi_stars', 'uminusr', 'deluminusr']]
selected_sersic_SDSS = sersic_cat_SDSS[['CATAID', 'RA', 'DEC', 'GALINDEX_i', 'GALINDEX_g', 'GALINDEX_u', 'GALINDEX_r', 'GALMAG_i', 'GALMAG_g', 'GALMAG_u', 'GALMAG_r']]
selected_sersic_VIKING = sersic_cat_VIKING[['CATAID', 'RA', 'DEC', 'GALINDEX_Z', 'GALMAG_Z']]

selected_SFR = SFR_cat[['CATAID', 'sSFR_0_1Gyr_best_fit', 'SFR_0_1Gyr_best_fit']]
selected_matched = matched_cat[['uberID', 'ra', 'declination', 'redshift']]


# plot 1:
'''
colour_mass_cat = selected_mass
colour_mass_cat.write('crossmatches/join_colour_mass_cat_all.fits', overwrite = True)
matched_color_mass_cat = join(selected_mass, selected_matched, keys='uberID')
matched_color_mass_cat.write('crossmatches/join_colour_mass_cat_matched465.fits', overwrite = True)

'''
# plot 2:

colour_sersic_cat_viking1 = join(selected_science, selected_mass, keys='uberID')
colour_sersic_cat_viking = join(colour_sersic_cat_viking1, selected_sersic_VIKING, keys='CATAID')
colour_sersic_cat_viking.write('crossmatches/join_colour_sersic_cat_all_VIKING.fits', overwrite = True)

matched_colour_sersic_cat_viking = join(colour_sersic_cat_viking, selected_matched, keys='uberID')
matched_colour_sersic_cat_viking.write('crossmatches/join_colour_sersic_cat_matched465_VIKING.fits', overwrite = True)

colour_sersic_cat_sdss1 = join(selected_science, selected_mass, keys='uberID')
colour_sersic_cat_sdss = join(colour_sersic_cat_sdss1, selected_sersic_SDSS, keys='CATAID')
colour_sersic_cat_sdss.write('crossmatches/join_colour_sersic_cat_all_SDSS.fits', overwrite = True)

matched_colour_sersic_cat_sdss = join(colour_sersic_cat_sdss, selected_matched, keys='uberID')
matched_colour_sersic_cat_sdss.write('crossmatches/join_colour_sersic_cat_matched465_SDSS.fits', overwrite = True)

exit()
# plot 3:

colour_SFR_cat1 = join(selected_science, selected_mass, keys='uberID')
colour_SFR_cat = join(colour_SFR_cat1, selected_SFR, keys='CATAID')
colour_SFR_cat.write('crossmatches/join_colour_SFR_cat_all.fits', overwrite = True)

matched_colour_SFR_cat = join(colour_SFR_cat, selected_matched, keys='uberID')
matched_colour_SFR_cat.write('crossmatches/join_colour_SFR_cat_matched465.fits', overwrite = True)
