 
# coding: utf-8

# Creates:
#     * lightcurve_snec.pdf

 
# In[5]:


from astropy.io import ascii as asc
import astropy.units as u
from matplotlib import pyplot as plt
import supernova
import glob
import os
import sys
import numpy as np
import visualization
#get_ipython().run_line_magic('matplotlib', 'inline')

import matplotlib as mpl


 
 
# In[6]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[7]:


sys.path.append('/Users/bostroem/Desktop/research/not_my_code/SNEC-1.01/')
import chisq_analysis


 
 
# In[8]:


DARK_DIR = '/Users/bostroem/dark'


 
# # Get Photometry

 
# In[9]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz.get_photometry()
sn15oz.get_abs_mag()


 
# # Get Model

 
# In[10]:


ni_mass = [0.083, 0.0965, 0.11]
energies = [0.5, 0.8, 1.1, 1.4, 1.7, 2.0]
masses = [11, 13, 15, 16, 17, 18, 21]
ni_mixing = [5.0]
time_offsets = np.arange(-4, 4, 1)
Kvalues = [10, 20, 30, 35, 40, 50, 60]
radii = [1500, 1800, 2100, 2400, 2700, 3000, 3300]

snec_models = os.path.join(DARK_DIR,'SNEC/snec_models/')
snname = 'asassn15oz'
S2_start = 50
S2_end = 88  #Vary this parameter
snec_15oz = chisq_analysis.SnecAnalysis(snname, snec_models, S2_start, S2_end, 
                 ni_mass, ni_mixing, masses, energies, time_offsets, 
                 Kvalues, radii, fig_dir='../figures')


 
 
# In[11]:


best_model_csm_dir = os.path.join(snec_models, 
                                 'Ni_mass_{:1.4f}'.format(0.083),
                                 'Ni_mixing_{:1.1f}'.format(5.0),
                                 'M{:2.1f}'.format(18.0),
                                 'E_{:1.3f}'.format(1.400),
                                 'K_{:2.1f}'.format(10.0), 
                                 'R_{}'.format(2400),
                                 'Data')
print(best_model_csm_dir)
best_model_csm_data = snec_15oz.prepare_model_data(best_model_csm_dir)


 
 
# In[12]:


best_model_bare_dir = os.path.join(DARK_DIR,'bostroem/research/not_my_code/SNEC/asassn15oz/mixing_5.0/M18/E_1.4/Data')
best_model_bare_data = snec_15oz.prepare_model_data(best_model_bare_dir)


 
 
# In[14]:


fig = plt.figure()
fig.set_figheight(6)
fig.subplotpars.update(bottom=.09)

ax = fig.add_subplot(111)
toffset = -4
ax.axhspan(-4, -14.3, hatch='/', facecolor='none', edgecolor='k')
ax.axhspan(-4, -14.3, color='white', alpha = 0.7)


ax.set_title('Light Curve and SNEC Model')

li, = ax.plot(sn15oz.phase['i']+toffset, sn15oz.abs_mag['i']-4, 'o', label='i-4')
ax.plot(best_model_csm_data['time'], best_model_csm_data['i']-4,color=li.get_color())
ax.plot(best_model_bare_data['time'], best_model_bare_data['i']-4, ls='--', color=li.get_color())

lr, = ax.plot(sn15oz.phase['r']+toffset, sn15oz.abs_mag['r']-1.2, 'o', label='r-1.2')
ax.plot(best_model_csm_data['time'], best_model_csm_data['r']-1.2, color=lr.get_color())
ax.plot(best_model_bare_data['time'], best_model_bare_data['r']-1.2, ls='--', color=lr.get_color())

lv, = ax.plot(sn15oz.phase['V']+toffset, sn15oz.abs_mag['V'], 'o', label='V+0')
ax.plot(best_model_csm_data['time'], best_model_csm_data['V'],color=lv.get_color() )
ax.plot(best_model_bare_data['time'], best_model_bare_data['V'], ls='--', color=lv.get_color())

lg, = ax.plot(sn15oz.phase['g']+toffset, sn15oz.abs_mag['g']+1, 'o', label='g+1')
ax.plot(best_model_csm_data['time'], best_model_csm_data['g']+1, color=lg.get_color())
ax.plot(best_model_bare_data['time'], best_model_bare_data['g']+1, ls='--', color=lg.get_color())

lb, = ax.plot(sn15oz.phase['B']+toffset, sn15oz.abs_mag['B']+4, 'o', alpha=0.75, label='b+4')
ax.plot(best_model_csm_data['time'], best_model_csm_data['B']+4, color=lb.get_color(), alpha=0.75)
ax.plot(best_model_bare_data['time'], best_model_bare_data['B']+4, ls='--', color=lb.get_color(), alpha=0.75)

lu, = ax.plot(sn15oz.phase['U']+toffset, sn15oz.abs_mag['U']+5.8, 'o', alpha=0.75, label='U+5.8')
ax.plot(best_model_csm_data['time'], best_model_csm_data['U']+5.8, color=lu.get_color(),  alpha=0.75)
ax.plot(best_model_bare_data['time'], best_model_bare_data['U']+5.8, ls='--', color=lu.get_color(), alpha=0.75)

ax.barh(-9.25, 97.5, 1.,6.6, facecolor='w', edgecolor='LightGrey', alpha=0.75)
ax.annotate(s='With CSM (wind)', xy=(12, -9.5), xytext=(32, -9.4), xycoords='data', 
            arrowprops={'arrowstyle':'-',  'linestyle':'--','linewidth':1.5}, fontsize=6.0, 
            backgroundcolor='none')
ax.annotate(s='Without CSM (wind)', xy=(12, -9), xytext=(32, -8.9), xycoords='data', fontsize=6.0, 
            arrowprops={'arrowstyle':'-', 'linewidth':1.5}, 
            backgroundcolor='none')


ax.set_ylim(-7.5, -23)
ax.set_xlim(0, 140)
ax.set_xlabel('Phase (day)')
ax.set_ylabel('Absolute magnitude + offset')
#plt.grid()
handles, labels = ax.get_legend_handles_labels()
leg = plt.legend(handles[::3], labels[::3], ncol=3, loc='lower left')
plt.savefig('lightcurve_snec.pdf')


 