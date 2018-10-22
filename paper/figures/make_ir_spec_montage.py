 
# coding: utf-8

# Creates:
#     * ir_spec_montage_log.pdf

 
# In[1]:


import os

import numpy as np
from astropy.io import ascii as asc
from astropy.io import fits
import astropy.constants as c
import astropy.units as u
from astropy.time import Time
from astropy.modeling import fitting, models

from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', '')

#https://github.com/abostroem/utilities
import supernova
import spectroscopy as spec
import visualization


 
 
# In[2]:


IRTF_DIR = '../../data/spectra/IRTF/'
SOFI_DIR = '../../data/spectra/SOFI'
LCO_DIR = '../../data/spectra/lco/'
EFOSC_DIR = '../../data/spectra/EFOSC/'

FIG_DIR = '.'
VEL_DATA_DIR = '../../data/line_info'


 
 
# In[3]:


z = 0.0069 #15oz redshift


 
 
# In[4]:


#Date-obs: 2015-10-10
tbdata_irtf = asc.read(os.path.join(IRTF_DIR, 'A15oz_merge.txt'), names=['wave', 'flux', 'err', 'junk'])


 
 
# In[5]:


tbdata_sofi1_blue = fits.getdata(os.path.join(SOFI_DIR, 'asassn15oz_20150905_2457270.58657_1.fits'), 1) #wave, flux, err, skyback
tbdata_sofi1_red = fits.getdata(os.path.join(SOFI_DIR, 'asassn15oz_20150905_2457270.62701_1.fits'), 1)
tbdata_sofi2_blue = fits.getdata(os.path.join(SOFI_DIR, 'asassn15oz_20151005_2457300.50252_1.fits'), 1)
tbdata_sofi2_red = fits.getdata(os.path.join(SOFI_DIR, 'asassn15oz_20151005_2457300.54299_1.fits'), 1)


 
 
# In[6]:


import matplotlib as mpl
font = mpl.font_manager.FontProperties()
font.set_family('cursive')
font.set_name(mpl.rcParams['font.monospace'][6:])
font.set_size('5')
print(font.get_name())
print(mpl.rcParams['font.monospace'][6:])


 
 
# In[9]:


spec_sofi1_blue = spec.spectrum1d(spec.apply_redshift(tbdata_sofi1_blue['WAVE'][0],z), 
                                  tbdata_sofi1_blue['FLUX'][0], 
                                  tbdata_sofi1_blue['ERR'][0])
spec_sofi1_red = spec.spectrum1d( spec.apply_redshift(tbdata_sofi1_red['WAVE'][0], z), 
                                  tbdata_sofi1_red[ 'FLUX'][0], 
                                  tbdata_sofi1_red[ 'ERR'][0])
spec_sofi2_blue = spec.spectrum1d(spec.apply_redshift(tbdata_sofi2_blue['WAVE'][0], z), 
                                  tbdata_sofi2_blue['FLUX'][0], 
                                  tbdata_sofi2_blue['ERR'][0])
spec_sofi2_red = spec.spectrum1d( spec.apply_redshift(tbdata_sofi2_red['WAVE'][0], z), 
                                 tbdata_sofi2_red[ 'FLUX'][0], 
                                 tbdata_sofi2_red[ 'ERR'][0])
spec_irtf = spec.spectrum1d(spec.apply_redshift(tbdata_irtf['wave']*10**4, z), 
                            tbdata_irtf['flux'], 
                            tbdata_irtf['err'])


 
 
# In[10]:


scale_sofi1_blue = spec.scale_spectra(spec_sofi1_blue, spec_sofi1_blue, wlmin=15000, wlmax=16000)
scale_sofi1_red = spec.scale_spectra(spec_sofi1_red, spec_sofi1_blue,   wlmin=15000, wlmax=16000)
scale_sofi2_blue = spec.scale_spectra(spec_sofi2_blue, spec_sofi1_blue, wlmin=15000, wlmax=16000)
scale_sofi2_red = spec.scale_spectra(spec_sofi2_red, spec_sofi1_blue,   wlmin=15000, wlmax=16000)
scale_irtf = spec.scale_spectra(spec_irtf, spec_sofi1_blue,             wlmin=15000, wlmax=16000)


 
 
# In[14]:


plt.style.use(['seaborn-paper', 'az-paper-twocol'])
fig = plt.figure()
ax_ir = fig.add_subplot(1,1,1)

telluric_blue = (13300, 14500)
telluric_red = (17800, 19250)


segment1_indx = scale_irtf.wave<telluric_blue[0]
segment2_indx = (scale_irtf.wave>telluric_blue[1]) & (scale_irtf.wave<telluric_red[0])
segment3_indx = scale_irtf.wave>telluric_red[1]

ax_ir.axvspan(telluric_blue[0], telluric_blue[1], color='LightGray', alpha=0.5)
ax_ir.axvspan(telluric_red[0], telluric_red[1], color='LightGray', alpha=0.5)
#Plot Second SOFI spec deredshifted
l2, = ax_ir.plot(scale_sofi1_blue.wave, np.log10(scale_sofi1_blue.flux))
l3, = ax_ir.plot(scale_sofi2_blue.wave, np.log10(scale_sofi2_blue.flux)-0.5)
l1, = ax_ir.plot(scale_irtf.wave[segment1_indx], np.log10(scale_irtf.flux[segment1_indx])-1)
ax_ir.plot(scale_irtf.wave[segment2_indx], np.log10(scale_irtf.flux[segment2_indx])-1, color=l1.get_color())


ax_ir.plot(scale_sofi1_red.wave, np.log10(scale_sofi1_red.flux), color=l2.get_color())
ax_ir.plot(scale_sofi2_red.wave, np.log10(scale_sofi2_red.flux)-0.5, color=l3.get_color())
ax_ir.plot(scale_irtf.wave[segment3_indx], np.log10(scale_irtf.flux[segment3_indx])-1, color=l1.get_color())

ax_ir.set_ylim(-17.75,-14.5)
ax_ir.set_xlim(9000, 24900)

annotation_dict_blue = {'CaII 9906': (9686,       -16.15),
                   'MgII 10092': (9857,           -16.15),
                   r'P-$\delta$':(10232,          -16.15),
                   'CI 10691/HeI 10830': (10500, -16.4),
                   r'P-$\gamma$': (10675,         -16.5),
                   'OI 11290': (11075,            -16.5),
                   'CaII 11836': (11572,          -16.25),
                   'CaII 11947': (11685,          -16.9),
                   r'P-$\beta$':(12503,           -16.5)}

for key in annotation_dict_blue.keys():
    x,y = annotation_dict_blue[key]
    yoffset=0.1
    textoffset = 0.05
    ax_ir.axvline(x, color='k', alpha=0.05, zorder=0)
    ax_ir.plot([x,x], [y, y-yoffset], color='k', lw=0.5)
    ax_ir.text(x, y-yoffset-textoffset, key, rotation='vertical', ha='center', va='top',
        fontproperties=font)
    
    
annotation_dict_red = {'SiI 15893': (15632,   -15.25),
                       'CI 16600': (16351,    -15.25), 
                       'CI 16890': (16655,    -15.25),
                       r'Br-$\gamma$':(21480, -15.75),
                       'HeI 20587': (20041,   -15.5)}
for key in annotation_dict_red.keys():
    x,y = annotation_dict_red[key]
    yoffset=0.1
    textoffset = 0.05
    ax_ir.axvline(x, color='k', alpha=0.05, zorder=0)
    ax_ir.plot([x,x], [y, y+yoffset], color='k', lw=0.5)
    ax_ir.text(x, y+yoffset+textoffset, key, rotation='vertical', ha='center', va='bottom', 
        fontproperties=font)
ax_ir.set_ylabel(r'log(Flux) ($\rm erg/cm^2/s/\AA$)+ offset')
ax_ir.set_xlabel(r'Wavelength ($\rm \AA$)')
ax_ir.set_yticks([])

ax_ir.text(25000, -16.0, '8d' , ha='left')
ax_ir.text(25000, -16.65, '38d', ha='left')
ax_ir.text(25000, -17.0, '43d', ha='left')
plt.savefig(os.path.join(FIG_DIR, 'ir_spec_montage_log.pdf'))


 