 
# coding: utf-8

# Creates:
#     * ni_mass_lc.pdf

 
# In[2]:


import os
from astropy.io import ascii as asc
from astropy.time import Time
from astropy.table import Table, vstack
from matplotlib import pyplot as plt
import numpy as np
#get_ipython().run_line_magic('matplotlib', '')
import astropy.units as u
from astropy.modeling import fitting, models

from utilities_az import supernova, connect_to_sndavis


 
 
# In[3]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[4]:


FIG_DIR = '.'


 
 
# In[5]:


tbdata = asc.read('../../data/asassn15oz_bolo_UBgVrRiI.txt', names=['phase', 'logL', 'err'])
tbdata_87A = asc.read('../../data/bol_lum_1987A_extrap.txt', names=['phase', 'logL', 'err'])


 
 
# In[6]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz.get_photometry('V')
sn15oz.get_slope('s2', band='V')
sn15oz.get_slope('s50', band='V')
sn15oz.get_slope('tail')
Time(sn15oz.jdexpl, format='jd').iso


 
 
# In[7]:


lin_model = models.Linear1D()
fitter = fitting.LinearLSQFitter()

s2_indx = (tbdata['phase']>50) & (tbdata['phase']<100)
s2_fit = fitter(lin_model, tbdata['phase'][s2_indx], tbdata['logL'][s2_indx], weights=1/tbdata['err'][s2_indx])

tail_indx = (tbdata['phase']>100) & (tbdata['phase']<300)
tail_fit = fitter(lin_model, tbdata['phase'][tail_indx], tbdata['logL'][tail_indx], weights=1/tbdata['err'][tail_indx])


 
 
# In[8]:


fall_from_plateau_length = 20 #days
start_fall_phase_upper = tbdata['phase'][s2_indx][-1]
end_fall_phase_upper = start_fall_phase_upper + fall_from_plateau_length

conservative_tpt = 125 #days
end_fall_phase_lower = conservative_tpt+fall_from_plateau_length/2.
start_fall_phase_lower = conservative_tpt-fall_from_plateau_length/2.


 
 
# In[9]:


sn87A_indx = tbdata_87A['phase']>75
sn87A_fit = fitter(lin_model, tbdata_87A['phase'][sn87A_indx], tbdata_87A['logL'][sn87A_indx])#, weights=1/tbdata_87A['err'][sn87A_indx])
sn87A_phase = tbdata_87A['phase'][sn87A_indx]
sn87A_logL = sn87A_fit(sn87A_phase)
sn87A_upper_scale_factor = tail_fit(end_fall_phase_upper)/np.interp(end_fall_phase_upper, sn87A_phase, sn87A_logL)
sn87A_upper_scale_flux = sn87A_logL*sn87A_upper_scale_factor

sn87A_lower_scale_factor = tail_fit(end_fall_phase_lower)/np.interp(end_fall_phase_lower, sn87A_phase, sn87A_logL)
sn87A_lower_scale_flux = sn87A_logL*sn87A_lower_scale_factor


 
 
# In[10]:


fig = plt.figure()
fig.set_figheight(4)
fig.subplotpars.update(bottom=0.12)
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)


tail_phase = np.arange(end_fall_phase_upper, 400)
ax1.errorbar(tbdata['phase'], tbdata['logL'], tbdata['err'], fmt='.', capsize=1, linestyle='none', label='Observation')
#ax1.plot(tbdata['phase'], tbdata['logL'], '.', color='#BBBBBB', label='Observations')
ax1.plot(tbdata['phase'][s2_indx], s2_fit(tbdata['phase'][s2_indx]),  label='S2 fit')
ax1.plot([start_fall_phase_upper, end_fall_phase_upper], [s2_fit(start_fall_phase_upper), tail_fit(end_fall_phase_upper)],  label='Artificial fall')
ax1.plot(tail_phase, tail_fit(tail_phase),  label='Tail fit')
ax1.plot(sn87A_phase, sn87A_upper_scale_flux, ls=':', label='Complete trapping')
ax1.plot(end_fall_phase_upper, tail_fit(end_fall_phase_upper), '*', markersize=8,  label='Ni Scale Luminosity', mec='k')

ax1.legend(loc='lower left', framealpha=0)
ax1.set_xlim(-5, 250)
ax1.set_ylim(40.25, 42.75)
ax1.set_xlabel('Phase (days)')
#ax1.set_ylabel('Log(L$_{bol}$) ($ergs/cm^2/s$)')
#ax1.set_title('Scale luminosity to determine upper limit on Ni mass')


tail_phase = np.arange(end_fall_phase_lower, 300)
s2_phase = np.arange(tbdata['phase'][s2_indx][0], start_fall_phase_lower)

ax2.errorbar(tbdata['phase'], tbdata['logL'], tbdata['err'], fmt='.', capsize=1, linestyle='none', label='Observation')
#ax2.plot(tbdata['phase'], tbdata['logL'], '.',  label='Observations')
ax2.plot(s2_phase, s2_fit(s2_phase), label='S2 fit')
ax2.plot([s2_phase[-1], tail_phase[0]], [s2_fit(s2_phase[-1]), tail_fit(tail_phase[0])], label='Artificial fall')
ax2.plot(tail_phase, tail_fit(tail_phase), label='Tail fit')
ax2.plot(sn87A_phase, sn87A_lower_scale_flux, ls=':',  label='Complete trapping')
ax2.plot(tail_phase[0], tail_fit(end_fall_phase_lower), '*', markersize=8, label='Ni Scale Luminosity', mec='k')
ax2.legend(loc='lower left', framealpha=0)
ax2.set_xlim(-5, 250)
ax2.set_ylim(40.25, 42.75)
ax2.set_xlabel('Phase (day)')
ax2.set_ylabel('\t\t\t\t\tLog(L$_{bol}$) (erg/s)')
#ax2.set_title('Scale luminosity to determine lower limit on Ni mass')
plt.savefig(os.path.join(FIG_DIR, 'ni_mass_lc.pdf'))


 
# Check whether $L_{complete}/L_{obs} = M_{Ni,mod}/M_{Ni, obs}$

 
# In[11]:


phase = 228
Lcomplete = np.interp(phase, sn87A_phase, sn87A_lower_scale_flux)
Lobs = tail_fit(phase)
print('Lcomplete/Lobs = {}'.format(Lcomplete/Lobs))


 
 
# In[12]:


phase = 340
Lcomplete = 10**np.interp(phase, sn87A_phase, sn87A_lower_scale_flux)
Lobs = 10**tail_fit(phase)
print('Lcomplete/Lobs = {}'.format(Lcomplete/Lobs))


 
# # What would the luminosity discrepancy have to be to get the Ni mass ratio of 4 (or 3.75) ?

 
# In[13]:


fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

tail_phase = np.arange(end_fall_phase_lower, 300)
s2_phase = np.arange(tbdata['phase'][s2_indx][0], start_fall_phase_lower)

ax1.errorbar(tbdata['phase'], tbdata['logL'], tbdata['err'], fmt='.', capsize=1, linestyle='none', label='Observation')
ax1.plot(s2_phase, s2_fit(s2_phase), label='S2 fit')
ax1.plot([s2_phase[-1], tail_phase[0]], [s2_fit(s2_phase[-1]), tail_fit(tail_phase[0])], label='Artificial fall')
ax1.plot(tail_phase, tail_fit(tail_phase), label='Tail fit')
ax1.plot(sn87A_phase, sn87A_lower_scale_flux, ls=':',  label='Complete trapping')
ax1.plot(tail_phase[0], tail_fit(end_fall_phase_lower), '*', markersize=8, label='Ni Scale Luminosity', mec='k')
ax1.plot(228, tail_fit(228)+np.log10(0.08/0.03), 's')
ax1.plot(340, tail_fit(340)+np.log10(0.08/0.02), 's')
ax1.legend(loc='lower left', framealpha=0)
ax1.set_xlim(-5, 350)
ax1.set_ylim(40, 42.75)
ax1.set_xlabel('Phase (day)')
ax1.set_ylabel('Log(L$_{bol}$) (erg/s)')
ax1.grid()


 