
# coding: utf-8

# In[1]:


import os
import numpy as np
from astropy.io import fits, ascii as asc
from astropy.table import Table
from utilities_az import supernova, spectroscopy as spec


OPTICAL_DIR = '../data/spectra/EFOSC/'
optical_fname = 'tASASSN-15oz_20151003_Gr13_Free_slit1.0_57720_1_e.fits'
filter_dir = '/Users/bostroem/Dropbox/DLT40_all/script/scalespectra/filters/'

optical_ofile = fits.open(os.path.join(OPTICAL_DIR, optical_fname))
flux = optical_ofile[0].data[0,0,:]
pix = np.arange(len(flux))
wl = spec.calc_wavelength(optical_ofile[0].header, pix)
optical_spec = spec.spectrum1d(wl, flux)

sn15oz = supernova.LightCurve2('asassn-15oz')
redshift = 0.0069

cur_dir = os.getcwd()
#deredshift
rest_wl = spec.apply_redshift(optical_spec.wave, redshift=redshift)
#deredden
obs_spec = spec.spectrum1d(rest_wl, optical_spec.flux)
obs_spec_dered = spec.correct_for_galactic_extinction(obs_spec, sn15oz.ebv_mw, R_V=3.1)
obs_spec_dered = spec.correct_for_galactic_extinction(obs_spec_dered, sn15oz.ebv_host, R_V=3.1)
tbdata_dered = Table([obs_spec_dered.wave, obs_spec_dered.flux], names=['wave', 'flux'])
tbdata_dered.write(os.path.join(OPTICAL_DIR, '{}_rest_dustcorr.dat'.format(optical_fname.strip('.fits')).replace('-', '_')), 
                    format='ascii.no_header', overwrite=True)
optical_date = optical_ofile[0].header['date-obs']
os.chdir(OPTICAL_DIR)
supernova.scale_spectra_quba('asassn-15oz', 
                             '{}_rest_dustcorr.dat'.format(optical_fname.strip('.fits')), 
                             filter_dir=filter_dir, 
                             date_kw=optical_date, 
                             max_sep=2, 
                             header_date=False)

os.chdir(cur_dir)
