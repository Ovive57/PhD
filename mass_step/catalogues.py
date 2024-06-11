
from astropy.table import Table, hstack, vstack, join
import numpy as np

gama_science_file = "data/gkvScienceCatv02.fits" # 'uberID', 'RAcen', 'Deccen', 'CATAID', 'Z'

science_cat = Table.read(gama_science_file)

ang = science_cat['ang']
ind = np.where((ang>44.9999)&(ang<45.001))

cat2 = science_cat[ind]
cat2.write('crossmatches/angles/45.fits', overwrite = True)
print(science_cat['CATAID'][ind])
ind = np.where((ang>89.9999)&(ang<90.001))
print(science_cat['CATAID'][ind])

cat2 = science_cat[ind]
cat2.write('crossmatches/angles/90.fits', overwrite = True)

ind = np.where((ang>134.9999)&(ang<135.001))
print(science_cat['CATAID'][ind])

cat2 = science_cat[ind]
cat2.write('crossmatches/angles/135.fits', overwrite = True)

ind = np.where((ang>59.9999)&(ang<60.001))
print(science_cat['CATAID'][ind])

cat2 = science_cat[ind]
cat2.write('crossmatches/angles/60.fits', overwrite = True)

ind = np.where((ang>29.9999)&(ang<30.001))
print(science_cat['CATAID'][ind])

cat2 = science_cat[ind]
cat2.write('crossmatches/angles/30.fits', overwrite = True)
