 
# coding: utf-8

# Creates:
#     * cachito_evolution.pdf

 
# In[17]:


import os

from astropy.io import fits
from astropy.io import ascii as asc
from astropy.time import Time
import astropy.constants as c
import astropy.units as u
from astropy.modeling import models, fitting
import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from cycler import cycler

from utilities_az import spectroscopy as spec, visualization, supernova


 
 
# In[18]:


plt.style.use(['seaborn-paper','az-paper-onecol'])


 
 
# In[19]:


cm_rainbow = visualization.make_rainbow_cm()


 
# # Plot HA, Hb evolution Series

 
# In[20]:


sn15oz = supernova.LightCurve2('asassn-15oz')


 
 
# In[21]:


def read_iraf_spectrum(filename):
    ofile = fits.open(filename)
    flux = ofile[0].data[0,0,:]
    wave = spec.calc_wavelength(ofile[0].header, np.arange(len(flux))+1)
    rest_wave = spec.apply_redshift(wave, 0.0069)
    return(spec.spectrum1d(rest_wave, flux))


 
 
# In[22]:


def calc_velocity(obs_wl, rest_wl):
    velocity = c.c*(obs_wl/rest_wl - 1)
    return velocity


 
 
# In[23]:


def calc_obs_wave(velocity, rest_wl):
    obs_wl = (velocity/c.c+1.0)*rest_wl
    return obs_wl


 
 
# In[24]:


VEL_DATA_DIR =   '../../data/line_info/'
DATA_DIR_LCO =   '../../data/spectra/lco'
DATA_DIR_EFOSC = '../../data/spectra/EFOSC'
DATA_DIR_XSHOOT = '../../data/spectra/xshooter/'
IRTF_DIR =       '../../data/spectra/IRTF/'
SOFI_DIR =       '../../data/spectra/SOFI'
TEST_FILE_DIR = '../../data/line_info/testing/'
HBETA_DIR = '../../data/line_info'


 
 
# In[25]:


z = 0.0069 #15oz redshift
Ha = 6563.0
Hb = 4850.0
HeI = 10830.0
CI = 10691.0
paschen_g = 10940.0
paschen_d = 10050.0
FeII = 5169.0
IR_dates = Time(['2015-09-05','2015-10-05', '2015-10-10'])


 
 
# In[26]:


spectra_files = [
         ('asassn15oz_20150904_redblu_122216.314.fits', DATA_DIR_LCO),
         ('asassn-15oz_20150906_redblu_105042.698a.fits', DATA_DIR_LCO),
         ('asassn15oz_20150907_redblu_123835.277.fits', DATA_DIR_LCO),
         ('asassn15oz_20150911_redblu_105336.349.fits', DATA_DIR_LCO),
         ('asassn15oz_20150916_redblu_120911.274.fits', DATA_DIR_LCO),
         ('asassn-15oz_20150920_redblu_135034.512.fits',DATA_DIR_LCO),
        #('ASASSN15oz_VLT_20150921.txt', DATA_DIR_XSHOOT),
         ('asassn-15oz_20150924_redblu_123847.580.fits',DATA_DIR_LCO),
         #('asassn-15oz_20150930_redblu_122858.217.fits',DATA_DIR_LCO),
         ('tASASSN-15oz_20151003_Gr13_Free_slit1.0_57720_1_e.fits',DATA_DIR_EFOSC),
         #('asassn15oz_20151006_redblu_101906.800.fits', DATA_DIR_LCO),
         ('asassn15oz_20151014_redblu_112918.305.fits', DATA_DIR_LCO),
         ('asassn-15oz_20151025_redblu_102221.833.fits', DATA_DIR_LCO),
         #('asassn-15oz_20151107_redblu_101210.833.fits', DATA_DIR_LCO),
         ('tASAS-SN_15oz_20151107_Gr13_Free_slit1.5_57723_1_e.fits', DATA_DIR_EFOSC),
         ('tASAS-SN_15oz_20151118_Gr13_Free_slit1.0_57723_1_e.fits', DATA_DIR_EFOSC)
                ]


 
 
# In[27]:


tbdata_sofi1_blue = fits.getdata(os.path.join(SOFI_DIR, 'asassn15oz_20150905_2457270.58657_1.fits'), 1) #wave, flux, err, skyback
tbdata_sofi2_blue = fits.getdata(os.path.join(SOFI_DIR, 'asassn15oz_20151005_2457300.50252_1.fits'), 1)
#Date-obs: 2015-10-10
tbdata_irtf = asc.read(os.path.join(IRTF_DIR, 'A15oz_merge.txt'), names=['wave', 'flux', 'err', 'junk'])


 
 
# In[28]:


spec_sofi1 = spec.spectrum1d(spec.apply_redshift(tbdata_sofi1_blue['WAVE'][0], z), tbdata_sofi1_blue['flux'][0], tbdata_sofi1_blue['err'])
spec_sofi2 = spec.spectrum1d(spec.apply_redshift(tbdata_sofi2_blue['WAVE'][0], z), tbdata_sofi2_blue['flux'][0], tbdata_sofi2_blue['err'])
spec_irtf = spec.spectrum1d(spec.apply_redshift(tbdata_irtf['wave']*10**4, z), tbdata_irtf['flux'], tbdata_irtf['err'])


 
 
# In[29]:


new_fit_cachito = asc.read(os.path.join(TEST_FILE_DIR, 'cachito.tab'))
phase_cachito = (Time(new_fit_cachito['date'])-Time(sn15oz.jdexpl, format='jd')).value
velocity_cachito = -1*calc_velocity(new_fit_cachito['vel0'], Ha).to(u.km/u.s).value
fitter_power = fitting.LevMarLSQFitter()
fitter_linear = fitting.LinearLSQFitter()
power_model = models.PowerLaw1D()
power_fit_cachito = fitter_power(power_model, phase_cachito, velocity_cachito)


 
 
# In[30]:


new_fit_ha = asc.read(os.path.join(TEST_FILE_DIR, 'Ha.tab'))
phase_ha = (Time(new_fit_ha['date'])-Time(sn15oz.jdexpl, format='jd')).value
velocity_ha = -1*calc_velocity(new_fit_ha['vel0'], Ha).to(u.km/u.s).value
power_fit_ha = fitter_power(power_model, phase_ha, velocity_ha)


 
 
# In[31]:


new_fit_hb = asc.read(os.path.join(HBETA_DIR, 'Hbeta.tab'))
phase_hb = (Time(new_fit_hb['date'])-Time(sn15oz.jdexpl, format='jd')).value
velocity_hb = -1*calc_velocity(new_fit_hb['velocity'], Hb).to(u.km/u.s).value
power_fit_hb = fitter_power(power_fit_ha, phase_hb, velocity_hb)


 
 
# In[33]:


plt.close()
fig = plt.figure()
fig.set_figheight(5.0)

width = 0.21
left = 0.2
height = 0.86
bottom = 0.09
pad = 0.025

ax1 = plt.axes([left, bottom, width, height])
ax2 = plt.axes([left+width+pad, bottom, width, height], sharex=ax1)
ax3 = plt.axes([left+2*width+2*pad, bottom, width, height], sharex=ax1)
offset = np.arange(len(spectra_files))*0.4E-14
ref_spec = read_iraf_spectrum(os.path.join(spectra_files[0][1], spectra_files[0][0]))
ha_velocity = calc_velocity(ref_spec.wave, Ha).to(u.km/u.s)
ref_norm_indx = np.where(ha_velocity<(calc_velocity(6000, Ha).to(u.km/u.s)))[0]
#Consider making this all one color
colors = visualization.make_color_wheel(spectra_files+[1], cmap=cm_rainbow)

#Scale phase by this when offsetting
ha_scale = 0.1
hb_scale = 0.6
he_scale = 0.1

#Plot the dates for the HV features
HV_H_limits = np.array((-40,-80))
HV_He_limits = np.array((-20, -60))
ax1.axhspan(*(HV_H_limits*ha_scale), color='#777777', alpha=0.2)
ax2.axhspan(*(HV_H_limits*hb_scale), color='#777777', alpha=0.2)
ax3.axhspan(*(HV_He_limits*he_scale),color='#777777', alpha=0.2)

#Find the velocity range plotted in Ha for use in Hb and HeI
ha_vel_xmin = calc_velocity(6000, Ha)
hb_obs_wl_min = calc_obs_wave(ha_vel_xmin, Hb)

#Plot Ha and Hb
for indx, ifile in enumerate(spectra_files):
    filename, idir = ifile
    ispec = read_iraf_spectrum(os.path.join(idir, filename))
    ha_velocity = calc_velocity(ispec.wave, Ha).to(u.km/u.s)
    hb_velocity = calc_velocity(ispec.wave, Hb).to(u.km/u.s)    
    scale_ispec = spec.scale_spectra(ispec, ref_spec)
    scale_factor_ha = np.interp(6000, scale_ispec.wave, scale_ispec.flux)
    scale_factor_hb = np.interp(hb_obs_wl_min, scale_ispec.wave, scale_ispec.flux)
    date = Time(fits.getval(os.path.join(idir, filename), 'date-obs', 0), out_subfmt='date')
    phase = int((date-Time(sn15oz.jdexpl, format='jd')).value)
    
    if int(phase) == 8:
        ax1.text(-26, -ha_scale*phase*0.85, '{}d'.format(int(phase)), 
                  ha='right', va='bottom',  color='k')
        ax1.plot(ha_velocity/1000, scale_ispec.flux/scale_factor_ha-1.0-ha_scale*phase*0.8, color=colors[indx])
        ax2.plot(hb_velocity/1000, scale_ispec.flux/scale_factor_hb-1.0-hb_scale*phase*0.8, color=colors[indx])
    elif int(phase) == 9:
        ax1.text(-26, -ha_scale*phase*0.98, '{}d'.format(int(phase)), 
                  ha='right', va='bottom',  color='k')
        ax1.plot(ha_velocity/1000, scale_ispec.flux/scale_factor_ha-1.0-ha_scale*phase, color=colors[indx])
        ax2.plot(hb_velocity/1000, scale_ispec.flux/scale_factor_hb-1.0-hb_scale*phase, color=colors[indx])
    else:
        ax1.plot(ha_velocity/1000, scale_ispec.flux/scale_factor_ha-1.0-ha_scale*phase, color=colors[indx])
        ax1.text(-26, -ha_scale*phase, '{}d'.format(int(phase)), 
                  ha='right', va='bottom',  color='k')
        ax2.plot(hb_velocity/1000, scale_ispec.flux/scale_factor_hb-1.0-hb_scale*phase, color=colors[indx])

    
#Plot Cachito
dates_cachito = np.arange(4, 95)
model_cachito_vel = power_fit_cachito(dates_cachito)
ax1.plot(-model_cachito_vel/1000, -dates_cachito*ha_scale, c='#AE76A3', ls='--')
ax2.plot(-model_cachito_vel/1000, -dates_cachito*hb_scale, c='#AE76A3', ls='--')

model_ha_vel = power_fit_ha(dates_cachito)
ax1.plot(-model_ha_vel/1000, -dates_cachito*ha_scale, color='#A5170E', ls=':')

model_hb_vel = power_fit_hb(dates_cachito)
ax2.plot(-model_hb_vel/1000, -dates_cachito*hb_scale, color='#A5170E', ls=':')

#Plot the IR
heI_obs_wl_min = calc_obs_wave(ha_vel_xmin, HeI)
scale_sofi1 = spec.scale_spectra(spec_sofi1, spec_sofi2) #SOFI1
scale_irtf = spec.scale_spectra(spec_irtf, spec_sofi2) #IRTF
heI_velocity_sofi2 = calc_velocity(spec_sofi2.wave, HeI).to(u.km/u.s)
heI_velocity_sofi1 = calc_velocity(scale_sofi1.wave, HeI).to(u.km/u.s)
heI_velocity_irtf = calc_velocity(scale_irtf.wave, HeI).to(u.km/u.s)

scale_factor_HeI_sofi1 = np.interp(heI_obs_wl_min, scale_sofi1.wave, scale_sofi1.flux)
phase_sofi1 = (Time('2015-09-05')-Time(sn15oz.jdexpl, format='jd')).value
ax3.plot(heI_velocity_sofi1/1000, scale_sofi1.flux/scale_factor_HeI_sofi1-1.0-he_scale*phase_sofi1, color=colors[0])
ax3.text(3, -he_scale*phase_sofi1, '{}d'.format(int(phase_sofi1)), 
              ha='left', va='bottom',  color='k')

scale_factor_HeI_sofi2 = np.interp(heI_obs_wl_min, spec_sofi2.wave, spec_sofi2.flux)
phase_sofi2 = (Time('2015-10-05')-Time(sn15oz.jdexpl, format='jd')).value
scale_factor_HeI_sofi22 = np.interp(heI_obs_wl_min, spec_sofi2.wave, spec_sofi2.flux)
ax3.plot(heI_velocity_sofi2/1000, spec_sofi2.flux/scale_factor_HeI_sofi2-1.0-he_scale*phase_sofi2, color=colors[-6])
ax3.text(3, -he_scale*phase_sofi2, '{}d'.format(int(phase_sofi2)), 
              ha='left', va='bottom',  color='k')

scale_factor_HeI_irtf = np.interp(heI_obs_wl_min, scale_irtf.wave, scale_irtf.flux)
phase_irtf = (Time('2015-10-10')-Time(sn15oz.jdexpl, format='jd')).value
ax3.plot(heI_velocity_irtf/1000, scale_irtf.flux/scale_factor_HeI_irtf-1.0-he_scale*phase_irtf, color=colors[-5])
ax3.text(3, -he_scale*phase_irtf, '{}d'.format(int(phase_irtf)), 
              ha='left', va='bottom', color='k')

ylimits = np.array((-90, -4))
ax1.set_xticks([-20, 0])
ax1.set_yticks([])
ax1.set_xlim(calc_velocity(6000, Ha).to(u.km/u.s).value/1000, calc_velocity(6600, Ha).to(u.km/u.s).value/1000)
ax1.set_ylim(ha_scale*ylimits)
ax1.set_ylabel(r'Flux + offset', labelpad=25.0)
ax1.set_title(r'H $\rm \alpha$')

ax2.set_ylim(hb_scale*ylimits)
ax2.set_yticks([])
ax2.set_xlabel(r'Velocity (1000 km $\rm s^{-1}$)')
ax2.set_title(r'H $\rm \beta$')

ax3.set_yticks([])
ax3.set_ylim(he_scale*ylimits)
ax3.set_title(r'He I', ha='center')

plt.subplots_adjust(wspace=0.1)
plt.savefig('cachito_evolution.pdf')



 