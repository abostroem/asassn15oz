 
# coding: utf-8

# Creates:
#     * spec_obs_tab.tex

 
# In[32]:


import os
import glob
from astropy.io import ascii as asc
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from utilities_az import connect_to_sndavis, supernova


 
 
# In[33]:


sn15oz = supernova.LightCurve2('asassn-15oz')


 
 
# In[34]:


IRTF_DIR = '../../data/spectra/IRTF/'
SOFI_DIR = '../../data/spectra/SOFI'

DATA_DIR_LCO = '../../data/spectra/lco/'
DATA_DIR_EFOSC = '../../data/spectra/EFOSC/'
DATA_DIR_XSHOOT = '../../data/spectra/xshooter/'
DATA_DIR_GEM = '/Users/bostroem/Desktop/research/asassn15oz/data/spectra/gmos/'


 
# # IR

 
# In[35]:


sofi1 = os.path.join(SOFI_DIR, 'asassn15oz_20150905_2457270.58657_1.fits')
sofi2= os.path.join(SOFI_DIR, 'asassn15oz_20151005_2457300.50252_1.fits')
irtf = os.path.join(IRTF_DIR, 'A15oz_merge.txt')
#IRTF date from commented header
ir_spec_dates = [Time(fits.getval(sofi1, 'date-obs', 0), out_subfmt='date'), Time(fits.getval(sofi2, 'date-obs', 0), out_subfmt='date'), Time(57305.1963045, format='mjd', out_subfmt='date')]
ir_spec_telescopes = ['NTT', 'NTT', 'IRTF']
ir_spec_instruments = ['SOFI', 'SOFI', 'SpeX']


 
# # Optical

 
# In[36]:


spectra_files = [
         ('asassn15oz_20150904_redblu_122216.314.fits', DATA_DIR_LCO),
         ('asassn-15oz_20150906_redblu_105042.698a.fits', DATA_DIR_LCO),
         ('asassn15oz_20150907_redblu_123835.277.fits', DATA_DIR_LCO),
         ('asassn15oz_20150911_redblu_105336.349.fits', DATA_DIR_LCO),
         ('asassn15oz_20150916_redblu_120911.274.fits', DATA_DIR_LCO),
         ('asassn-15oz_20150920_redblu_135034.512.fits',DATA_DIR_LCO),
         ('ASASSN15oz_VLT_20150921.txt', DATA_DIR_XSHOOT),
         ('asassn-15oz_20150924_redblu_123847.580.fits',DATA_DIR_LCO),
         ('asassn-15oz_20150930_redblu_122858.217.fits',DATA_DIR_LCO),
         ('tASASSN-15oz_20151003_Gr13_Free_slit1.0_57720_1_e.fits',DATA_DIR_EFOSC),
         ('asassn15oz_20151006_redblu_101906.800.fits', DATA_DIR_LCO),
         ('asassn15oz_20151014_redblu_112918.305.fits', DATA_DIR_LCO),
         ('asassn-15oz_20151025_redblu_102221.833.fits', DATA_DIR_LCO),
         ('asassn-15oz_20151107_redblu_101210.833.fits', DATA_DIR_LCO),
         ('tASAS-SN_15oz_20151107_Gr13_Free_slit1.5_57723_1_e.fits', DATA_DIR_EFOSC),
         ('tASAS-SN_15oz_20151118_Gr13_Free_slit1.0_57723_1_e.fits', DATA_DIR_EFOSC),
         ('tASASSN-15oz_20160410_Gr13_Free_slit1.5_57723_1_e.fits', DATA_DIR_EFOSC),
         ('comb20160610_R400.fits', DATA_DIR_GEM),
         ('tASASSN-15oz_20160802_Gr13_Free_slit1.0_57723_1_e.fits', DATA_DIR_EFOSC),
         ('tASASSN-15oz_20160918_Gr13_Free_slit1.5_57723_2_e.fits', DATA_DIR_EFOSC),
                ]

optical_spec_dates = []
optical_spec_telescopes = []
optical_spec_instruments = []

for ifile, idir in spectra_files:
    filename = os.path.join(idir, ifile)
    if idir == DATA_DIR_LCO:
        optical_spec_dates.append(Time(fits.getval(filename, 'date-obs', 0),out_subfmt='date'))
        optical_spec_telescopes.append('LCO')
        optical_spec_instruments.append('FLOYDS')
    elif idir == DATA_DIR_EFOSC:
        optical_spec_dates.append(Time(fits.getval(filename, 'date-obs', 0),out_subfmt='date'))
        optical_spec_telescopes.append('NTT')
        optical_spec_instruments.append('EFOSC')
    elif idir == DATA_DIR_XSHOOT:
        #date from ESO XSHOOTER archive
        optical_spec_dates.append(Time(57286.981832762, format='mjd', out_subfmt='date'))
        optical_spec_telescopes.append('VLT')
        optical_spec_instruments.append('X-SHOOTER')
    elif idir == DATA_DIR_GEM:
        optical_spec_dates.append(Time(fits.getval(filename, 'mjd-obs', 0), format='mjd', out_subfmt='date'))
        optical_spec_telescopes.append('Gemini-S')
        optical_spec_instruments.append('GMOS')
        


 
# # UV

 
# In[37]:


SWIFT_DIR = '../../data/swiftuvot/reduced_default/'


 
 
# In[38]:


uv_spec_dates = []
uv_spec_telescopes = []
uv_spec_instruments = []
obsid_list1 = ['00034040001', '00034040002', '00034040005', '00034040007', '00034040009']

for obsid in obsid_list1:
    flist1 = glob.glob(os.path.join(SWIFT_DIR, obsid, 'uvot', 'image', '*.pha'))
    for ifile in flist1:
        hdr = fits.getheader(ifile, 1)
        uv_spec_dates.append(Time(fits.getval(ifile, 'date-obs', 1), out_subfmt='date'))
        uv_spec_telescopes.append('Swift')
        uv_spec_instruments.append('UVOTA')


 
 
# In[39]:


dateobs = []
jd = []
phase = []
instrument = []
telescope = []

for itel, iinst, idate in zip(uv_spec_telescopes, uv_spec_instruments, uv_spec_dates):
    dateobs.append('{}'.format(idate.iso))
    jd.append('{:7.1f}'.format(idate.jd))
    phase.append('{:3.1f}'.format(idate.jd - sn15oz.jdexpl))
    instrument.append(iinst)
    telescope.append(itel)
for itel, iinst, idate in zip(optical_spec_telescopes, optical_spec_instruments, optical_spec_dates):
    dateobs.append('{}'.format(idate.iso))
    jd.append('{:7.1f}'.format(idate.jd))
    phase.append('{:3.1f}'.format(idate.jd - sn15oz.jdexpl))
    instrument.append(iinst)
    telescope.append(itel)    
for itel, iinst, idate in zip(ir_spec_telescopes, ir_spec_instruments, ir_spec_dates):
    dateobs.append('{}'.format(idate.iso))
    jd.append('{:7.1f}'.format(idate.jd))
    phase.append('{:3.1f}'.format(idate.jd - sn15oz.jdexpl))
    instrument.append(iinst)
    telescope.append(itel)

tbdata = Table([dateobs, jd, phase, telescope, instrument], names=['Date', 'JD', 'Phase', 'Observatory', 'Instrument']) 
tbdata.sort(['Date'])
tbdata.write('../spec_obs_tab.tex', format='ascii.latex', overwrite=True,
            latexdict={'preamble':r'\centering',
                       'caption':r'Spectroscopic Observations of ASASSN-15oz',
                       'data_start':r'\hline',
                      'label':r'tab:SpecObs'},)


 