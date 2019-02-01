 
# coding: utf-8

# Creates:
#     * lightcurve_snec.pdf

 
# In[1]:


from astropy.io import ascii as asc
import astropy.units as u
from matplotlib import pyplot as plt
import matplotlib as mpl
import glob
import os
import sys
import numpy as np
from utilities_az import visualization, supernova
#get_ipython().run_line_magic('matplotlib', 'inline')

import matplotlib as mpl


 
 
# In[2]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[3]:


DARK_DIR = '/Users/bostroem/dark'


 
 
# In[4]:


def get_breakout_time(model_dir):
    ofile = open(os.path.join(model_dir, 'info.dat'), 'r')
    all_lines = ofile.readlines()
    if len(all_lines)>6:
        time_breakout = float((all_lines[5].split('=')[1]).strip('seconds\n'))
    else: #SN never got to breakout
        time_breakout = None
    return time_breakout


 
 
# In[5]:


def prepare_model_data(model_dir, type='mag'):
    if type=='mag':
        model_tbdata = asc.read(os.path.join(model_dir,'magnitudes.dat'),
                                names=['time', 'temp', 'weird_mag', 
                                       'u', 'g', 'r', 'i', 'z', 'U', 
                                       'B', 'V', 'R', 'I'])
    elif type=='lum':
        model_tbdata = asc.read(os.path.join(model_dir, 'lum_observed.dat'), names=['time', 'luminosity'])
    else:
        print('Type must be "lum" or "mag", {} not understood')
        sys.exit()
    #Observers call t_explosion the time when the explosion is first visible, therefore
    #the models should use this same conventions
    time_breakout = get_breakout_time(model_dir)
    if time_breakout is not None:
        model_tbdata['time'] = ((model_tbdata['time']-time_breakout)*u.second).to(u.day).value
    else:
        model_tbdata=None
    return model_tbdata


 
# # Get Photometry

 
# In[6]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz.get_photometry()
sn15oz.get_abs_mag()


 
# # Get Model

 
# In[7]:


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


 
# **best model parameters for color:**  
# bostroem@dark:/dark/bostroem/research/ASASSN15oz/snec_models  
#  mass: 17  
#  energy: 1.4   
#  K: 40   
#  R: 1800   
#  t_offset: 3   
#  Ni mass: 0.083  
# M_CSM = 1.5 solMass  

# **best model parameters for bare color**  
#  mass: 18  
#  energy: 2.0   
#  K: 0   
#  R: 0   
#  t_offset: 3   
#  Ni mass: 0.083  

# **best model parameters for black body Bolometric Luminosity:**  
# bostroem@dark:/dark/bostroem/research/ASASSN15oz/snec_models  
#  mass: 18  
#  energy: 1.4   
#  K: 60   
#  R: 1500   
#  t_offset: 2   
#  Ni mass: 0.083  
#  M_CSM = 1.4 sol  
# 

# **best model parameters for bare BB bolometric luminosity**  
#  mass: 17  
#  energy: 2.0   
#  K: 0   
#  R: 0   
#  t_offset: 1   
#  Ni mass: 0.083  

 
# In[8]:


best_model_csm_dir_color = os.path.join(snec_models, 
                                 'Ni_mass_{:1.4f}'.format(0.083),
                                 'Ni_mixing_{:1.1f}'.format(5.0),
                                 'M{:2.1f}'.format(17.0),
                                 'E_{:1.3f}'.format(1.400),
                                 'K_{:2.1f}'.format(40.0), 
                                 'R_{}'.format(1800),
                                 'Data')
best_model_csm_toff_color = 3


 
 
# In[9]:


best_model_bare_dir_color = os.path.join(snec_models, 
                                 'Ni_mass_{:1.4f}'.format(0.083),
                                 'Ni_mixing_{:1.1f}'.format(5.0),
                                 'M{:2.1f}'.format(18.0),
                                 'E_{:1.3f}'.format(2.00),
                                 'K_{:2.1f}'.format(0.0), 
                                 'R_{}'.format(0),
                                 'Data')
best_model_bare_toff_color = 3


 
 
# In[11]:


tbdata_model_color_csm = prepare_model_data(best_model_csm_dir_color)
tbdata_model_color_bare = prepare_model_data(best_model_bare_dir_color)


 
# # Make Bolometric luminosity plot

 
# In[12]:


best_model_csm_dir_bolo = os.path.join(snec_models, 
                                 'Ni_mass_{:1.4f}'.format(0.083),
                                 'Ni_mixing_{:1.1f}'.format(5.0),
                                 'M{:2.1f}'.format(18.0),
                                 'E_{:1.3f}'.format(1.400),
                                 'K_{:2.1f}'.format(60.0), 
                                 'R_{}'.format(1500),
                                 'Data')
best_model_csm_toff_bolo = 2


 
 
# In[13]:


best_model_bare_dir_bolo = os.path.join(snec_models, 
                                 'Ni_mass_{:1.4f}'.format(0.083),
                                 'Ni_mixing_{:1.1f}'.format(5.0),
                                 'M{:2.1f}'.format(17.0),
                                 'E_{:1.3f}'.format(2.000),
                                 'K_{:2.1f}'.format(0.0), 
                                 'R_{}'.format(0),
                                 'Data')
best_model_bare_toff_bolo = 1


 
 
# In[14]:


tbdata_model_bolo_csm = prepare_model_data(best_model_csm_dir_bolo, type='lum')
tbdata_model_bolo_bare = prepare_model_data(best_model_bare_dir_bolo, type='lum')


 
 
# In[15]:


sn15oz_tbdata = asc.read('../../data/asassn-15oz_bolo_BB.txt', names=['phase', 'logL'], delimiter=' ')


 
 
# In[16]:


def plot_lightcurve(ax, data, model_csm, model_bare, toff_csm, toff_bare, mag_off, ifilter, alpha=1):
    if mag_off < 0:
        label = '{}+{}'.format(ifilter, abs(mag_off))
    else:
        label = '{}-{}'.format(ifilter, mag_off)
    l, = ax.plot(data.phase[ifilter], data.abs_mag[ifilter]-mag_off, 'o', label=label, alpha=alpha)
    ax.plot(model_csm['time']-toff_csm, 
            model_csm[ifilter]-mag_off,
            color=l.get_color(), alpha=alpha)
    ax.plot(model_bare['time']-toff_bare, 
            model_bare[ifilter]-mag_off, 
            ls='--', color=l.get_color(), alpha=alpha)    
    return l, ax


 
 
# In[17]:


fig = plt.figure()
fig.set_figheight(5.5)
fig.subplotpars.update(bottom=.09, left=0.22)

#ax = fig.add_subplot(1,1,1)
ax = fig.add_axes([0.2, 0.25, 0.75, 0.65])
ax_bolo = fig.add_axes([0.2, 0.1, 0.75, 0.15])

toffset = -4
ax.axhspan(-4, -14.3, hatch='/', facecolor='none', edgecolor='k')
ax.axhspan(-4, -14.3, color='white', alpha = 0.7)


ax.set_title('Light Curve and SNEC Model')


#ax_bolo = ax.twinx()
lbolo, = ax_bolo.plot(sn15oz_tbdata['phase'], sn15oz_tbdata['logL'], 'o',  label='Bolometric Luminosity')
ax_bolo.plot(tbdata_model_bolo_csm['time']-best_model_csm_toff_bolo, 
             np.log10(tbdata_model_bolo_csm['luminosity']), label='With CSM (wind)', color='gray')
ax_bolo.plot(tbdata_model_bolo_bare['time']-best_model_bare_toff_bolo, 
             np.log10(tbdata_model_bolo_bare['luminosity']), ls='--', label='Without CSM (wind)', color='gray')

ax.plot([], [])
li, ax = plot_lightcurve(ax, sn15oz, tbdata_model_color_csm, tbdata_model_color_bare,
                        best_model_csm_toff_color, best_model_bare_toff_color,
                        3, 'i')


lr, ax = plot_lightcurve(ax, sn15oz, tbdata_model_color_csm, tbdata_model_color_bare,
                        best_model_csm_toff_color, best_model_bare_toff_color,
                        1.2, 'r')

lv, ax = plot_lightcurve(ax, sn15oz, tbdata_model_color_csm, tbdata_model_color_bare,
                        best_model_csm_toff_color, best_model_bare_toff_color,
                        0, 'V')

lg, ax = plot_lightcurve(ax, sn15oz, tbdata_model_color_csm, tbdata_model_color_bare,
                        best_model_csm_toff_color, best_model_bare_toff_color,
                        -1, 'g')

lb, ax = plot_lightcurve(ax, sn15oz, tbdata_model_color_csm, tbdata_model_color_bare,
                        best_model_csm_toff_color, best_model_bare_toff_color,
                        -4, 'B', alpha=0.75)

lu, ax = plot_lightcurve(ax, sn15oz, tbdata_model_color_csm, tbdata_model_color_bare,
                        best_model_csm_toff_color, best_model_bare_toff_color,
                        -5.8, 'U', alpha=0.75)
#lbolo, = ax.plot([], [])

#ax.barh(-9.80, 85, 1.,5.5, facecolor='w', edgecolor='LightGrey', alpha=0.75) #y, width, height, left
#ax.annotate(s='Without CSM (wind)', xy=(12, -10), xytext=(32, -9.9), xycoords='data', 
#            arrowprops={'arrowstyle':'-',  'linestyle':'--','linewidth':1.5}, fontsize=6.0, 
#            backgroundcolor='none')
#ax.annotate(s='With CSM (wind)', xy=(12, -9.5), xytext=(32, -9.4), xycoords='data', fontsize=6.0, 
#            arrowprops={'arrowstyle':'-', 'linewidth':1.5}, 
#            backgroundcolor='none')
top = 43.2
ax_bolo.barh(top-0.25, 85, 0.7,60, facecolor='w', edgecolor='LightGrey', alpha=1) #y, width, height, left
ax_bolo.annotate(s='Without CSM (wind)', xy=(62, top-0.26), xytext=(80, top-0.26), xycoords='data', 
            arrowprops={'arrowstyle':'-',  'linestyle':'--','linewidth':1.5}, fontsize=6.0, 
            backgroundcolor='w')
ax_bolo.annotate(s='With CSM (wind)', xy=(62, top-0.5), xytext=(80, top-0.5), xycoords='data', fontsize=6.0, 
            arrowprops={'arrowstyle':'-', 'linewidth':1.5}, 
            backgroundcolor='none')




ax.set_ylim(-7.5, -21.75)
ax.set_xlim(0, 140)
ax.set_xticklabels([])
ax.set_ylabel('Absolute magnitude + offset')
handles, labels = ax.get_legend_handles_labels()

ax_bolo.set_xlim(0, 140)
ax_bolo.set_ylim(41.5, 43.2)


leg = ax.legend(handles[::3]+[lbolo], labels[::3]+[lbolo.get_label()], ncol=3, bbox_to_anchor=[0.3, -.02, 0.75, 0.2])
#ax_bolo.legend([lbolo], [lbolo.get_label()], loc='upper left')
ax_bolo.set_ylabel(r'$\rm Log_{10}(L_{bolo})$ (erg $\rm s^{-1})$')
ax_bolo.set_xlabel('Phase (day)')
plt.savefig('lightcurve_snec.pdf')


 