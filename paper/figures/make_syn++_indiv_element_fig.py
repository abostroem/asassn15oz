 
# coding: utf-8

# Creates:
#     * syn++indiv_elements.pdf

 
# In[1]:


import os
import glob

import numpy as np
from astropy.io import ascii as asc
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
from cycler import cycler

from utilities_az import spectroscopy as spec


 
 
# In[2]:


plt.style.use(['seaborn-paper', 'az-paper-twocol'])


 
 
# In[3]:


SYN_DIR = '../../data/syn++/'
DATA_DIR = '../../data/spectra'
FIG_DIR = '.'


 
 
# In[4]:


flist = glob.glob(os.path.join(SYN_DIR, 'syn_spec*only.txt'))
spectrum = asc.read(os.path.join(DATA_DIR, 'lco','asassn15oz_20151006_redblu_101906.800_dust_corr.asc'))
#spectrum = asc.read(os.path.join(DATA_DIR, 'EFOSC','ASASSN-15oz_20151003_Gr13_Free_slit1.0_57720_1_e_rest_dustcorrsca.dat'), names=['rest_wavelength', 'dust_corr_flux'])
scale_indx = spectrum['rest_wavelength']>5000


 
 
# In[5]:


tbdata_all = asc.read(os.path.join(SYN_DIR, 'syn_spec_Hb_all.txt'), names=['wave', 'flux', '?'])
all_scale_indx = tbdata_all['wave']>5000
scale_all = np.trapz(spectrum['dust_corr_flux'][scale_indx],spectrum['rest_wavelength'][scale_indx])/np.trapz(tbdata_all['flux'][all_scale_indx],tbdata_all['wave'][all_scale_indx])


 
 
# In[6]:


fig, axlist = plt.subplots(nrows=11)
fig.set_figheight(6)
fig.subplotpars.update(left=0.08, bottom=0.1, top=0.98, right=0.98)
elements = [r'H-$\alpha$', r'H-$\beta$', 'Na I','O I' ,'Ca II','Sc II', 'Ti II'  , 'Fe I' , 'Fe II', 'Ba I']
filename = {r'H-$\alpha$': 'H' ,
            r'H-$\beta$': 'HB_',
            'Na I': 'Na_',
            'O I': 'O_' ,
            'Ca II': 'Ca_',
            'Sc II': 'Sc_', 
            'Ti II': 'Ti_', 
            'Fe I': 'FeI_' , 
            'Fe II': 'FeII_', 
            'Ba I': 'Ba_'}
indx=1
default_cycler = plt.rcParams['axes.prop_cycle']
ax_cycler = default_cycler[indx:indx+1]
indx+=1
axlist[0].set_prop_cycle(ax_cycler)
l1 = axlist[0].plot(spectrum['rest_wavelength'], spectrum['dust_corr_flux']*10**15, label='observed spectrum', color='k')
l2 = axlist[0].plot(tbdata_all['wave'], tbdata_all['flux']*scale_all*10**15, label='syn++ all ions')
#first_legend = axlist[0].legend(borderaxespad=0., bbox_to_anchor=(0.7, 0.85, 0.25, .102))
axlist[0].legend(bbox_to_anchor=(0.025, 0.7, 0.95, .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0., framealpha=0)

#axlist[0].legend(handles=l2, borderaxespad=0.,bbox_to_anchor=(0.75, 0.75, 0.25, .102))
axlist[0].set_xticks([])

for ion, ax in zip(elements, axlist[1:]):
    ax_cycler = default_cycler[indx:indx+1]
    ax.set_prop_cycle(ax_cycler)
    ifile = os.path.join(SYN_DIR, 'syn_spec_{}only.txt'.format(filename[ion]))
    ax.plot(spectrum['rest_wavelength'], spectrum['dust_corr_flux']*10**15, color='k')
    tbdata = asc.read(ifile, names=['wave', 'flux', '?'])
    iscale_indx = tbdata['wave']>5000
    scale = np.trapz(spectrum['dust_corr_flux'][scale_indx],spectrum['rest_wavelength'][scale_indx])/np.trapz(tbdata['flux'][iscale_indx],tbdata['wave'][iscale_indx])
    line = ax.plot(tbdata['wave'], tbdata['flux']*scale*10**15)
    ax.legend([line[0]], [ion], framealpha=0, bbox_to_anchor=(0.84, 0.95, 0.175, .102), mode='expand')
    ax.set_xticks([])
    indx+=1
    
plt.subplots_adjust(hspace=0.0)
axlist[5].set_ylabel(r'Flux ($x10^{-15}$ $\rm erg/cm^2/s/\AA$)')
axlist[10].set_xlabel(r'Wavelength ($\rm \AA$)')
axlist[10].set_xticks(np.arange(3000, 11000, 1000))
plt.savefig(os.path.join(FIG_DIR, 'syn++indiv_elements.pdf'))


 
 
# In[15]:


spectrum_old = asc.read(os.path.join(DATA_DIR, 'lco','asassn15oz_20151006_redblu_101906.800_dust_corr.asc'))
spectrum_old.colnames
plt.figure()
plt.plot(spectrum_old['rest_wavelength'], spectrum_old['dust_corr_flux'])
plt.plot(spectrum['rest_wavelength'], spectrum['dust_corr_flux'])


 