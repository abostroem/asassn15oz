 
# coding: utf-8

# Creates:
#     * nebular_spectra_OI.pdf

 
# In[8]:


import os
import copy

from astropy.io import ascii as asc
from astropy.io import fits
from astropy.time import Time
import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', '')
import matplotlib as mpl
from utilities_az import spectroscopy as spec, supernova, define_filters


 
 
# In[9]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
# Confirmed that scaled spectra files have not been de-redshifted

 
# In[10]:


FIG_DIR = '.'
neb_repo = '/Users/bostroem/Desktop/research/asassn15oz/data/spectra/EFOSC/scaled'
GEMINI_DIR = '/Users/bostroem/Desktop/research/asassn15oz/data/spectra/gmos/'
model_repo = '../../../nebular_spectra_OI/models/'


 
# ## Read in spectra and de-redshift

 
# In[11]:


specfile_212 = os.path.join(neb_repo, 'tASASSN_15oz_20160410_Gr13_Free_slit1.5_57723_1_esca.asci') 
tbdata_212 = asc.read(specfile_212, names = ['wavelength', 'flux'])
wl_212 = spec.apply_redshift(tbdata_212['wavelength'], redshift=0.0069)
spectrum_212 = tbdata_212['flux']


 
 
# In[12]:


specfile_340 = os.path.join(neb_repo, 'tASASSN_15oz_20160802_Gr13_Free_slit1.0_57723_1_esca.asci') 
#specfile = os.path.join(neb_repo, 'asassn15oz_20160919a_387dayssca.fits')#low S/N
tbdata_340 = asc.read(specfile_340, names = ['wavelength', 'flux'])
wl_340 = spec.apply_redshift(tbdata_340['wavelength'], redshift=0.0069)
spectrum_340 = tbdata_340['flux']


 
 
# In[103]:


specfile_306 = os.path.join(GEMINI_DIR, 'gmos_merge_rest_dustcorrsca.dat')
#specfile_306 = os.path.join(GEMINI_DIR, 'gmos_merge_rest_dustcorr.dat')
tbdata_306 = asc.read(specfile_306, names=['wave', 'flux'])
wl_306 = tbdata_306['wave']
spectrum_306 = tbdata_306['flux']


 
# ## Read in Models

 
# In[105]:


modelfile12_212 = os.path.join(model_repo, 'mzams12_212d.dat')
modelfile15_212 = os.path.join(model_repo, 'mzams15_212d.dat')
modelfile19_212 = os.path.join(model_repo, 'mzams19_212d.dat')
modelfile25_212 = os.path.join(model_repo, 'mzams25_212d.dat')


 
 
# In[106]:


modelfile12_340 = os.path.join(model_repo, 'mzams12_306d.dat')
modelfile15_340 = os.path.join(model_repo, 'mzams15_350d.dat')
modelfile19_340 = os.path.join(model_repo, 'mzams19_369d.dat')
modelfile25_340 = os.path.join(model_repo, 'mzams25_369d.dat')


 
 
# In[107]:


modelfile12_306 = os.path.join(model_repo, 'mzams12_306d.dat')
modelfile15_306 = os.path.join(model_repo, 'mzams15_306d.dat')
modelfile19_306 = os.path.join(model_repo, 'mzams19_250d.dat')
modelfile25_306 = os.path.join(model_repo, 'mzams25_306d.dat')


 
 
# In[108]:


mod12_212 = asc.read(modelfile12_212, names = ['wavelength', 'flux'])
mod15_212 = asc.read(modelfile15_212, names = ['wavelength', 'flux'])
mod19_212 = asc.read(modelfile19_212, names = ['wavelength', 'flux'])
mod25_212 = asc.read(modelfile25_212, names = ['wavelength', 'flux'])


 
 
# In[109]:


mod12_340 = asc.read(modelfile12_340, names = ['wavelength', 'flux'])
mod15_340 = asc.read(modelfile15_340, names = ['wavelength', 'flux'])
mod19_340 = asc.read(modelfile19_340, names = ['wavelength', 'flux'])
mod25_340 = asc.read(modelfile25_340, names = ['wavelength', 'flux'])


 
 
# In[110]:


mod12_306 = asc.read(modelfile12_306, names = ['wavelength', 'flux'])
mod15_306 = asc.read(modelfile15_306, names = ['wavelength', 'flux'])
mod19_306 = asc.read(modelfile19_306, names = ['wavelength', 'flux'])
mod25_306 = asc.read(modelfile25_306, names = ['wavelength', 'flux'])


 
# ## Create Spectrum Objects

 
# In[111]:


mod_spec12_212 = spec.spectrum1d(mod12_212['wavelength'], mod12_212['flux'])
mod_spec15_212 = spec.spectrum1d(mod15_212['wavelength'], mod15_212['flux'])
mod_spec19_212 = spec.spectrum1d(mod19_212['wavelength'], mod19_212['flux'])
mod_spec25_212 = spec.spectrum1d(mod25_212['wavelength'], mod25_212['flux'])
data_spec_212 = spec.spectrum1d(wl_212-10, spectrum_212)


 
 
# In[112]:


mod_spec12_340 = spec.spectrum1d(mod12_340['wavelength'], mod12_340['flux'])
mod_spec15_340 = spec.spectrum1d(mod15_340['wavelength'], mod15_340['flux'])
mod_spec19_340 = spec.spectrum1d(mod19_340['wavelength'], mod19_340['flux'])
mod_spec25_340 = spec.spectrum1d(mod25_340['wavelength'], mod25_340['flux'])
data_spec_340 = spec.spectrum1d(wl_340-20, spectrum_340)


 
 
# In[113]:


mod_spec12_306 = spec.spectrum1d(mod12_306['wavelength'], mod12_306['flux'])
mod_spec15_306 = spec.spectrum1d(mod15_306['wavelength'], mod15_306['flux'])
mod_spec19_306 = spec.spectrum1d(mod19_306['wavelength'], mod19_306['flux'])
mod_spec25_306 = spec.spectrum1d(mod25_306['wavelength'], mod25_306['flux'])
data_spec_306 = spec.spectrum1d(wl_306-20, spectrum_306)


 
# # Scale Models to spectrum

 
# In[114]:


Ni_mass_mod = 0.062 #Msun
d_mod = 5.5 #Mpc
d_15oz = 28.83 #Mpc, NED Hubble + Virgo Infall
Co_halflife = 111.4

t_obs_306 = 288.0
t_mod_12_306 = 306.0
t_mod_15_306 = 306.0
t_mod_19_306 = 250.0
t_mod_25_306 = 306.0

t_obs_340 = 340.0
t_mod_12_340 = 306.0
t_mod_15_340 = 350.0
t_mod_19_340 = 369.0
t_mod_25_340 = 369.0

t_obs_212 = 228.0
t_mod_12_212 = 212.0
t_mod_15_212 = 212.0
t_mod_19_212 = 212.0
t_mod_25_212 = 212.0


 
# ##### Scale by time difference

 
# In[115]:


#Create new object
scale_time_mod_spec12_212 = copy.deepcopy(mod_spec12_212)
scale_time_mod_spec15_212 = copy.deepcopy(mod_spec15_212)
scale_time_mod_spec19_212 = copy.deepcopy(mod_spec19_212)
scale_time_mod_spec25_212 = copy.deepcopy(mod_spec25_212)
scale_time_mod_spec12_212.flux = scale_time_mod_spec12_212.flux*np.exp((t_mod_12_212-t_obs_212)/Co_halflife)
scale_time_mod_spec15_212.flux = scale_time_mod_spec15_212.flux*np.exp((t_mod_15_212-t_obs_212)/Co_halflife)
scale_time_mod_spec19_212.flux = scale_time_mod_spec19_212.flux*np.exp((t_mod_19_212-t_obs_212)/Co_halflife)
scale_time_mod_spec25_212.flux = scale_time_mod_spec25_212.flux*np.exp((t_mod_25_212-t_obs_212)/Co_halflife)


 
 
# In[116]:


#Create new object
scale_time_mod_spec12_340 = copy.deepcopy(mod_spec12_340)
scale_time_mod_spec15_340 = copy.deepcopy(mod_spec15_340)
scale_time_mod_spec19_340 = copy.deepcopy(mod_spec19_340)
scale_time_mod_spec25_340 = copy.deepcopy(mod_spec25_340)
scale_time_mod_spec12_340.flux = scale_time_mod_spec12_340.flux*np.exp((t_mod_12_340-t_obs_340)/Co_halflife)
scale_time_mod_spec15_340.flux = scale_time_mod_spec15_340.flux*np.exp((t_mod_15_340-t_obs_340)/Co_halflife)
scale_time_mod_spec19_340.flux = scale_time_mod_spec19_340.flux*np.exp((t_mod_19_340-t_obs_340)/Co_halflife)
scale_time_mod_spec25_340.flux = scale_time_mod_spec25_340.flux*np.exp((t_mod_25_340-t_obs_340)/Co_halflife)


 
 
# In[117]:


#Create new object
scale_time_mod_spec12_306 = copy.deepcopy(mod_spec12_306)
scale_time_mod_spec15_306 = copy.deepcopy(mod_spec15_306)
scale_time_mod_spec19_306 = copy.deepcopy(mod_spec19_306)
scale_time_mod_spec25_306 = copy.deepcopy(mod_spec25_306)
scale_time_mod_spec12_306.flux = scale_time_mod_spec12_306.flux*np.exp((t_mod_12_306-t_obs_306)/Co_halflife)
scale_time_mod_spec15_306.flux = scale_time_mod_spec15_306.flux*np.exp((t_mod_15_306-t_obs_306)/Co_halflife)
scale_time_mod_spec19_306.flux = scale_time_mod_spec19_306.flux*np.exp((t_mod_19_306-t_obs_306)/Co_halflife)
scale_time_mod_spec25_306.flux = scale_time_mod_spec25_306.flux*np.exp((t_mod_25_306-t_obs_306)/Co_halflife)


 
# #### Scale model to observed spectrum by distance and Ni mass (empirically)

 
# In[118]:


scale_mod12_212, scale_factor_12_212 =  spec.scale_spectra(scale_time_mod_spec12_212, data_spec_212, scale_factor=True)
scale_mod15_212, scale_factor_15_212 =  spec.scale_spectra(scale_time_mod_spec15_212, data_spec_212, scale_factor=True)
scale_mod19_212, scale_factor_19_212 =  spec.scale_spectra(scale_time_mod_spec19_212, data_spec_212, scale_factor=True)
scale_mod25_212, scale_factor_25_212 =  spec.scale_spectra(scale_time_mod_spec25_212, data_spec_212, scale_factor=True)


 
 
# In[119]:


scale_mod12_340, scale_factor_12_340 =  spec.scale_spectra(scale_time_mod_spec12_340, data_spec_340, scale_factor=True)
scale_mod15_340, scale_factor_15_340 =  spec.scale_spectra(scale_time_mod_spec15_340, data_spec_340, scale_factor=True)
scale_mod19_340, scale_factor_19_340 =  spec.scale_spectra(scale_time_mod_spec19_340, data_spec_340, scale_factor=True)
scale_mod25_340, scale_factor_25_340 =  spec.scale_spectra(scale_time_mod_spec25_340, data_spec_340, scale_factor=True)


 
 
# In[120]:


scale_mod12_306, scale_factor_12_306 =  spec.scale_spectra(scale_time_mod_spec12_306, data_spec_306, scale_factor=True)
scale_mod15_306, scale_factor_15_306 =  spec.scale_spectra(scale_time_mod_spec15_306, data_spec_306, scale_factor=True)
scale_mod19_306, scale_factor_19_306 =  spec.scale_spectra(scale_time_mod_spec19_306, data_spec_306, scale_factor=True)
scale_mod25_306, scale_factor_25_306 =  spec.scale_spectra(scale_time_mod_spec25_306, data_spec_306, scale_factor=True)


 
 
# In[121]:


fig = plt.figure()
fig.set_figheight(6)
fig.subplotpars.update(left=0.15, top=0.97, right=0.97)
fig.subplotpars.update(bottom=0.12)
ax1 = fig.add_subplot(3,1,1)
l4, = ax1.plot(data_spec_212.wave, data_spec_212.flux/10**-16,  lw=1.0, label='Day 228 Spectrum')
l0, = ax1.plot(scale_mod12_212.wave, scale_mod12_212.flux/10**-16, lw=0.5, alpha=0.8, label='12M$\odot$')
l1, = ax1.plot(scale_mod15_212.wave, scale_mod15_212.flux/10**-16, lw=0.5, alpha=0.8, label='15M$\odot$')
l2, = ax1.plot(scale_mod19_212.wave, scale_mod19_212.flux/10**-16, lw=0.5, alpha=0.8, label='19M$\odot$')
l3, = ax1.plot(scale_mod25_212.wave, scale_mod25_212.flux/10**-16, lw=0.5, alpha=0.8, label='25M$\odot$')


ax1.set_xlim(3500, 9200)
ax1.set_ylim(0, 9E-16/10**-16)
ax1.legend(bbox_to_anchor=[0.54,0.64, 0.39, 0.39], framealpha=1.0)

#Inset
ax1_inset = plt.axes([.25, .85, .2, .1])
ax1_inset.plot(data_spec_212.wave-10, data_spec_212.flux/10**-16, color=l4.get_color())
ax1_inset.plot(scale_mod12_212.wave, scale_mod12_212.flux/10**-16, alpha=0.8, color=l0.get_color())
ax1_inset.plot(scale_mod15_212.wave, scale_mod15_212.flux/10**-16, alpha=0.8, color=l1.get_color())
ax1_inset.plot(scale_mod19_212.wave, scale_mod19_212.flux/10**-16, alpha=0.8, color=l2.get_color())
ax1_inset.plot(scale_mod25_212.wave, scale_mod25_212.flux/10**-16, alpha=0.8, color=l3.get_color())


ax1_inset.set_xlim(6200,6500)
ax1_inset.set_ylim(0, 2.5E-16/10**-16)
ax1_inset.set_xticks(())
ax1_inset.set_yticks(())


ax2 = fig.add_subplot(3,1,2)
ax2.plot(data_spec_306.wave, data_spec_306.flux/10**-16, lw=1.0, label='Day 306 Spectrum', color=l4.get_color())
ax2.plot(scale_mod12_306.wave, scale_mod12_306.flux/10**-16, label='12M$\odot$',   lw=0.75, alpha=0.8,color=l0.get_color())
ax2.plot(scale_mod15_306.wave, scale_mod15_306.flux/10**-16, label='15M$\odot$',   lw=0.75, alpha=0.8,color=l1.get_color())
ax2.plot(scale_mod19_306.wave, scale_mod19_306.flux/10**-16, label='19M$\odot$',   lw=0.75, alpha=0.8,color=l2.get_color())
ax2.plot(scale_mod25_306.wave, scale_mod25_306.flux/10**-16, label='25M$\odot$',   lw=0.75, alpha=0.8,color=l3.get_color())



#plt.xlim(4000, 6500)
ax2.set_xlim(3500, 9200)
ax2.set_ylim(0, 7E-16/10**-16)
ax2.set_xlabel(r'Wavlength ($\rm \AA$)')
ax2.set_ylabel('\t\t\t\t\tFlux ($x10^{-16}$'+r'$\rm erg/cm^2/s/\AA$)')
ax2.legend(bbox_to_anchor=[0.54,0.64, 0.39, 0.39], framealpha=1.0)
#Inset
ax2_inset = plt.axes([.25, .55, .2, .1])
ax2_inset.plot(data_spec_306.wave, data_spec_306.flux/10**-16, color=l4.get_color())
ax2_inset.plot(scale_mod12_306.wave, scale_mod12_306.flux/10**-16,alpha=0.8, color=l0.get_color())
ax2_inset.plot(scale_mod15_306.wave, scale_mod15_306.flux/10**-16,alpha=0.8, color=l1.get_color())
ax2_inset.plot(scale_mod19_306.wave, scale_mod19_306.flux/10**-16,alpha=0.8, color=l2.get_color())
ax2_inset.plot(scale_mod25_306.wave, scale_mod25_306.flux/10**-16,alpha=0.8, color=l3.get_color())


ax2_inset.set_xlim(6200,6500)
ax2_inset.set_ylim(0, 3E-16/10**-16)
#xticks = ax2_inset.get_xticklabels()
#xticks = [tick.set_fontsize(6) for tick in xticks]
#yticks = ax2_inset.get_yticklabels()
#yticks = [tick.set_fontsize(6) for tick in yticks]
ax2_inset.set_xticks(())
ax2_inset.set_yticks(())


ax3 = fig.add_subplot(3,1,3)
ax3.plot(data_spec_340.wave, data_spec_340.flux/10**-16, lw=1.0, label='Day 340 Spectrum', color=l4.get_color())
ax3.plot(scale_mod12_340.wave, scale_mod12_340.flux/10**-16, label='12M$\odot$',   lw=0.75, alpha=0.8,color=l0.get_color())
ax3.plot(scale_mod15_340.wave, scale_mod15_340.flux/10**-16, label='15M$\odot$',   lw=0.75, alpha=0.8,color=l1.get_color())
ax3.plot(scale_mod19_340.wave, scale_mod19_340.flux/10**-16, label='19M$\odot$',   lw=0.75, alpha=0.8,color=l2.get_color())
ax3.plot(scale_mod25_340.wave, scale_mod25_340.flux/10**-16, label='25M$\odot$',   lw=0.75, alpha=0.8,color=l3.get_color())



#plt.xlim(4000, 6500)
ax3.set_xlim(3500, 9200)
ax3.set_ylim(0, 3E-16/10**-16)
ax3.set_xlabel(r'Wavlength ($\rm \AA$)')
ax3.set_ylabel('\t\t\t\t\tFlux ($x10^{-16}$'+r'$\rm erg/cm^2/s/\AA$)')
ax3.legend(bbox_to_anchor=[0.54,0.64, 0.39, 0.39], framealpha=1.0)
#Inset
ax3_inset = plt.axes([.25, .25, .2, .1])
ax3_inset.plot(data_spec_340.wave, data_spec_340.flux/10**-16, color=l4.get_color())
ax3_inset.plot(scale_mod12_340.wave, scale_mod12_340.flux/10**-16,alpha=0.8, color=l0.get_color())
ax3_inset.plot(scale_mod15_340.wave, scale_mod15_340.flux/10**-16,alpha=0.8, color=l1.get_color())
ax3_inset.plot(scale_mod19_340.wave, scale_mod19_340.flux/10**-16,alpha=0.8, color=l2.get_color())
ax3_inset.plot(scale_mod25_340.wave, scale_mod25_340.flux/10**-16,alpha=0.8, color=l3.get_color())


ax3_inset.set_xlim(6200,6500)
ax3_inset.set_ylim(0, 1.75E-16/10**-16)
#xticks = ax2_inset.get_xticklabels()
#xticks = [tick.set_fontsize(6) for tick in xticks]
#yticks = ax2_inset.get_yticklabels()
#yticks = [tick.set_fontsize(6) for tick in yticks]
ax3_inset.set_xticks(())
ax3_inset.set_yticks(())

#ax2.set_title('August 03, 2016 - Day 340')


fontsize=6

ax1_inset.annotate('[OI]', xy=(6415,0.95E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            horizontalalignment='center', fontsize=fontsize)

ax1.annotate('Mg I]', xy=(4571,1E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax1.annotate('Na I', xy=(5896,1.5E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax1.annotate('[OI]', xy=(6300,2.5E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax1.annotate(r'H$\alpha$', xy=(6400,7E-16/10**-16), xycoords='data',
            xytext=(-20, 0), textcoords='offset points', ha='left',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center',verticalalignment='center', fontsize=fontsize)
ax1.annotate('[Ca II]', xy=(7303,2E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax1.annotate('Ca II', xy=(8662,3E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)


ax2_inset.annotate('[OI]', xy=(6415,0.95E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            horizontalalignment='center', fontsize=fontsize)

ax2.annotate('Mg I]', xy=(4571,1E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax2.annotate('Na I', xy=(5896,1.5E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax2.annotate('[OI]', xy=(6300,2.75E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax2.annotate(r'H$\alpha$', xy=(6400,5.5E-16/10**-16), xycoords='data',
            xytext=(-20, 0), textcoords='offset points', ha='left',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center',verticalalignment='center', fontsize=fontsize)
ax2.annotate('[Ca II]', xy=(7303,2E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax2.annotate('Ca II', xy=(8662,2E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)

ax3_inset.annotate('[OI]', xy=(6415,0.7E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            horizontalalignment='center', fontsize=fontsize)

ax3.annotate('Mg I]', xy=(4571,0.5E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax3.annotate('Na I', xy=(5896,0.5E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax3.annotate('[OI]', xy=(6300,1.75E-16/10**-16), xycoords='data',
            xytext=(0, 20), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax3.annotate(r'H$\alpha$', xy=(6400,2.75E-16/10**-16), xycoords='data',
            xytext=(-20, 0), textcoords='offset points', ha='left',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center',verticalalignment='center', fontsize=fontsize)
ax3.annotate('[Ca II]', xy=(7303,0.9E-16/10**-16), xycoords='data',
            xytext=(0, 15), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
ax3.annotate('Ca II', xy=(8662,0.45E-16/10**-16), xycoords='data',
            xytext=(0, 15), textcoords='offset points',
            arrowprops=dict(facecolor='black', headlength=2, headwidth=5, width=1),
            horizontalalignment='center', fontsize=fontsize)
plt.savefig(os.path.join(FIG_DIR,'nebular_spectra_OI.pdf'))


 
# ## Calculate Ni mass referenced in paper 

 
# In[70]:


Ni_mass_12_212 = scale_factor_12_212*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 12 Msun @ 212 days = {}'.format(Ni_mass_12_212))
Ni_mass_15_212 = scale_factor_15_212*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 15 Msun @ 212 days = {}'.format(Ni_mass_15_212))
Ni_mass_19_212 = scale_factor_19_212*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 19 Msun @ 212 days = {}'.format(Ni_mass_19_212))
Ni_mass_25_212 = scale_factor_25_212*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 25 Msun @ 212 days = {}'.format(Ni_mass_25_212))
print('Average Ni mass for day 212 = {}'.format(np.mean(np.array([Ni_mass_12_212, Ni_mass_15_212, Ni_mass_19_212, Ni_mass_25_212]))))


 
 
# In[95]:


Ni_mass_12_306 = scale_factor_12_306*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 12 Msun @ 212 days = {}'.format(Ni_mass_12_306))
Ni_mass_15_306 = scale_factor_15_306*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 15 Msun @ 212 days = {}'.format(Ni_mass_15_306))
Ni_mass_19_306 = scale_factor_19_306*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 19 Msun @ 212 days = {}'.format(Ni_mass_19_306))
Ni_mass_25_306 = scale_factor_25_306*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 25 Msun @ 212 days = {}'.format(Ni_mass_25_306))
print('Average Ni mass for day 212 = {}'.format(np.mean(np.array([Ni_mass_12_306, Ni_mass_15_306, Ni_mass_19_306, Ni_mass_25_306]))))


 
 
# In[71]:


Ni_mass_12_340 = scale_factor_12_340*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 12 Msun @ 340 days = {}'.format(Ni_mass_12_340))
Ni_mass_15_340 = scale_factor_15_340*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 15 Msun @ 340 days = {}'.format(Ni_mass_15_340))
Ni_mass_19_340 = scale_factor_19_340*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 19 Msun @ 340 days = {}'.format(Ni_mass_19_340))
Ni_mass_25_340 = scale_factor_25_340*Ni_mass_mod * (d_15oz/d_mod)**2
print('Ni mass for 25 Msun @ 340 days = {}'.format(Ni_mass_25_340))
print('Average Ni mass for day 340 = {}'.format(np.mean(np.array([Ni_mass_12_340, Ni_mass_15_340, Ni_mass_19_340, Ni_mass_25_340]))))


 