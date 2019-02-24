
# coding: utf-8

# In[1]:


import os
import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table
from utilities_az import supernova, spectroscopy as spec


# In[2]:


UV_DIR = '../data/swiftuvot/reduced_default/'
OPTICAL_DIR = '../data/spectra/lco/'


# # UV Data

# In[3]:


optical_fname = 'asassn15oz_20150904_redblu_122216.314.ascii'
UV_fname = 'combine_epoch1.csv'
filter_dir = '/Users/bostroem/Dropbox/DLT40_all/script/scalespectra/filters/'


# In[4]:


uv_tbdata = asc.read(os.path.join(UV_DIR, UV_fname))
optical_tbdata = asc.read(os.path.join(OPTICAL_DIR, optical_fname), names=['wave', 'flux'])


# In[5]:


sn15oz = supernova.LightCurve2('asassn-15oz')
redshift = 0.0069


# In[ ]:

cur_dir = os.getcwd()
pixels = np.arange(len(uv_tbdata['flux']))+1
#deredshift
uv_rest_wl = spec.apply_redshift(uv_tbdata['wave'], redshift=redshift)
#deredden
uv_obs_spec = spec.spectrum1d(uv_rest_wl, uv_tbdata['flux'])
uv_obs_spec_dered = spec.correct_for_galactic_extinction(uv_obs_spec, sn15oz.ebv_mw, R_V=3.1)
uv_obs_spec_dered = spec.correct_for_galactic_extinction(uv_obs_spec_dered, sn15oz.ebv_host, R_V=3.1)
uv_tbdata_dered = Table([uv_obs_spec_dered.wave, uv_obs_spec_dered.flux], names=['wave', 'flux'])
uv_tbdata_dered.write(os.path.join(UV_DIR, 'combine_epoch1_rest_dustcorr.dat'), format='ascii.no_header', overwrite=True)
uv_date = '2015-09-05'
os.chdir(UV_DIR)
supernova.scale_spectra_quba('asassn-15oz', 
                             os.path.join('combine_epoch1_rest_dustcorr.dat'), 
                             filter_dir=filter_dir, 
                             date_kw=uv_date, 
                             max_sep=2, 
                             header_date=False)

os.chdir(cur_dir)


pixels = np.arange(len(optical_tbdata['flux']))+1
#deredshift
optical_rest_wl = spec.apply_redshift(optical_tbdata['wave'], redshift=redshift)
#deredden
optical_obs_spec = spec.spectrum1d(optical_rest_wl, optical_tbdata['flux'])
optical_obs_spec_dered = spec.correct_for_galactic_extinction(optical_obs_spec, sn15oz.ebv_mw, R_V=3.1)
optical_obs_spec_dered = spec.correct_for_galactic_extinction(optical_obs_spec_dered, sn15oz.ebv_host, R_V=3.1)
optical_tbdata_dered = Table([optical_obs_spec_dered.wave, optical_obs_spec_dered.flux], names=['wave', 'flux'])
optical_tbdata_dered.write(os.path.join(OPTICAL_DIR, 'asassn15oz_20150904_redblu_rest_dustcorr.dat'), format='ascii.no_header', overwrite=True)
optical_date = '2015-09-04'
os.chdir(OPTICAL_DIR)
supernova.scale_spectra_quba('asassn-15oz', 
                             os.path.join('asassn15oz_20150904_redblu_rest_dustcorr.dat'), 
                             filter_dir=filter_dir, 
                             date_kw=optical_date, 
                             max_sep=2, 
                             header_date=False)
                             
os.chdir(cur_dir)