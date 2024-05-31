import numpy as np
from astropy.cosmology import Planck18, FlatLambdaCDM

cosmo = Planck18
cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)

def physical2onsky(cosmology, z, physical_dist):
    comoving_distance = cosmology.comoving_distance(z).value
    comoving_distance_kpc = comoving_distance*1e3 #Mpc to kpc

    angular_distance_rad = physical_dist/comoving_distance_kpc
    angular_distance_arcsec = np.rad2deg(angular_distance_rad)*3600 #Radians to arcsec
    return angular_distance_arcsec

print(physical2onsky(cosmo, 0.01, 10))