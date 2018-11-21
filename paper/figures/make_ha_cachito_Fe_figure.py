 
# coding: utf-8

# Creates:
# * cachito_fe_vel_comp.pdf

 
# In[1]:


import os

import numpy as np
import yaml
from astropy.io import ascii as asc
from astropy.time import Time
import astropy.units as u
import astropy.constants as c
from astropy.modeling import models, fitting

from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

from utilities_az import supernova


 
 
# In[2]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[3]:


TEST_FILE_DIR = '../../data/line_info/testing/'
FIG_DIR = './'
DATA_DIR = '../../data/line_info'


 
 
# In[4]:


HA = 6563.0
SiII = 6355.0
FeII = 5169.0
IR_dates = Time(['2015-09-05','2015-10-05', '2015-10-10'])


 
 
# In[5]:


sn15oz = supernova.LightCurve2('asassn-15oz')
texpl = Time(sn15oz.jdexpl, format='jd')


 
 
# In[6]:


new_fit_cachito = asc.read(os.path.join(TEST_FILE_DIR, 'cachito.tab'))


 
 
# In[7]:


def calc_velocity(obs_wl, rest_wl):
    velocity = c.c*(obs_wl/rest_wl - 1)
    return velocity


 
 
# In[8]:


phase_cachito = (Time(new_fit_cachito['date'])-texpl).value
velocity_cachito = -1*calc_velocity(new_fit_cachito['vel0'], HA).to(u.km/u.s).value


 
 
# In[9]:


#tbdata_feII = asc.read(os.path.join(DATA_DIR, 'FeII_multi.tab'))
#tbdata_feII.remove_columns(['vel1', 'vel_err_left_1', 'vel_err_right_1', 'vel_pew_1', 'vel_pew_err1'])
tbdata_feII = asc.read(os.path.join(DATA_DIR, 'FeII_5169.tab'))
tbdata_feII.rename_column('vel0', 'velocity')
tbdata_feII.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_feII.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_feII.rename_column('vel_pew_0', 'pew')
tbdata_feII.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[10]:


phase_feII = (Time(tbdata_feII['date'])-texpl).value
velocity_feII = -1*calc_velocity(tbdata_feII['velocity'], FeII).to(u.km/u.s)


 
 
# In[11]:


fig = plt.figure()
fig.subplotpars.update(left=.17)
ax_Fe = fig.add_subplot(1,1,1)

ax_Fe.plot((Time(new_fit_cachito['date'])-texpl).value, -1*calc_velocity(new_fit_cachito['vel0'], SiII).to(u.km/u.s)/1000, '^', label='Cachito (as SiII 6533)')
ax_Fe.plot(phase_feII, velocity_feII/1000, 'o', label='FeII (5169)') 

ax_Fe.set_xticks(np.arange(0, 90, 10))
ax_Fe.legend()
ax_Fe.set_ylim(5, 11)
ax_Fe.set_xlim(0, 40)
ax_Fe.set_xlabel('Phase (day)')
ax_Fe.set_ylabel('Velocity (1000 km/s)')
plt.savefig(os.path.join(FIG_DIR, 'cachito_fe_vel_comp.pdf'))


 