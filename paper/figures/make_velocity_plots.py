 
# coding: utf-8

# Creates:
# * line_velocity.pdf
# * velocity_compare_guitierrez.pdf
# * snec_velocity_comp.pdf
# * phot_velocity_table.csv

 
# In[1]:


import os
import glob
from collections import namedtuple

from astropy.io import ascii as asc
from astropy.time import Time
from astropy import table
import astropy.units as u
import astropy.constants as c

import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

from utilities_az import spectroscopy as spec, supernova


 
 
# In[2]:


plt.style.use(['seaborn-paper','az-paper-onecol'])


 
 
# In[3]:


FIG_DIR = '.'
DATA_DIR = '../../data/line_info/'
TEST_FILE_DIR = '../../data/line_info/testing/'


 
# ### Wavelengths are already corrected for redshift!

 
# In[4]:


#sn = namedtuple('sn', ('jdexpl', 'redshift'))
#jd_expl = Time(2457262., format='jd')
#redshift = 0.006931 #Ref: NED object: HIPASS J1919-33
#sn15oz = sn(jdexpl=jd_expl, redshift=redshift)


 
 
# In[5]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz.jdexpl = Time(sn15oz.jdexpl, format='jd')


 
 
# In[6]:


flist = [#'hbeta.tab', 'hbeta_2.tab', 'hbeta_3.tab', 
         #'NaI5898.tab', 'NaI5898_2.tab', 
         'ScII5239.tab', 'ScII5239_2.tab', 
         #'ScII5526.tab', 
         'ScII5526_2.tab', 'ScII5526_3.tab', 
         'ScII5662.tab', 'ScII5662_2.tab', 
         'ScII6262.tab', 'ScII6262_2.tab', 'ScII6262_3.tab']

tbdata = asc.read(os.path.join(DATA_DIR,flist[0]))
for ifile in flist[1:]:
    tbdata = table.vstack((tbdata, asc.read(os.path.join(DATA_DIR,ifile))))


 
 
# In[7]:


for iline in tbdata:
    iline['line'] = iline['line'].split('_')[0]


 
 
# In[8]:


rest_wavelengths = {'NaI': 5898,
                    'HA-cachito': 6561,
                    'SiII':6355.0,
                    'HB': 4861,
                    #'hbeta': 4861, #Better fit with Gaussian than Silverman method
                    #'NaI5898':5893,  #Better fit with Gaussian than Silverman method 
                    'ScII5662':5662, 
                    'ScII6262': 6262,
                    'ScII5526': 5526,
                    #'ScII5239': 5239, #Better fit with Gaussian than Silverman method
                    'FeII_multi': 5169, 
                    'OI': 7774,
                    'CaII': 8498.
                   }
colors = {'HA-cachito': '#4477AA',
          'HB': '#332288', 
          'SiII': 'b',
          'hbeta': 'c',
          'NaI5898':'gold', 
          'NaI': '#DDCC77',
          'ScII5662':'#44AA99', 
          'ScII6262': '#117733',
          'ScII5526': '#999933',
          'ScII5239': 'r',
          'FeII_multi':'#882255', 
         'OI': '#AA4466', 
         'CaII': '#CC6677'}
markers = {'HA-cachito': 'o',
           'HB': '.',
           'SiII':'s',
          'hbeta': '.',
           'NaI': 's',
          'NaI5898':'s', 
          'ScII5662':'^', 
          'ScII6262': '>',
          'ScII5526': '<',
          'ScII5239': 'v',
          'FeII_multi':'d', 
          'OI': 'p',
          'CaII': '*'}


 
 
# In[9]:


tbdata_Ha_cach = asc.read(os.path.join(TEST_FILE_DIR, 'HA-cachito.tab'))
#tbdata_feII = asc.read(os.path.join(DATA_DIR, 'FeII_multi.tab'))
tbdata_feII = asc.read(os.path.join(DATA_DIR, 'FeII_5169.tab'))
tbdata_HB = asc.read(os.path.join(DATA_DIR, 'HB.tab'))
tbdata_NaI = asc.read(os.path.join(DATA_DIR, 'NaI.tab'))
tbdata_OI = asc.read(os.path.join(DATA_DIR, 'OI.tab'))
tbdata_CaII = asc.read(os.path.join(DATA_DIR, 'CaII.tab'))


 
 
# In[10]:


tbdata_halpha = tbdata_Ha_cach.copy()
tbdata_halpha.remove_columns(['vel0', 'vel_err_left_0', 'vel_err_right_0', 'vel_pew_0', 'vel_pew_err0'])
tbdata_halpha.remove_row(3) #row three is anomolously high

tbdata_halpha.rename_column('vel1', 'velocity')
tbdata_halpha.rename_column('vel_err_left_1', 'vel_err_left')
tbdata_halpha.rename_column('vel_err_right_1', 'vel_err_right')
tbdata_halpha.rename_column('vel_pew_1', 'pew')
tbdata_halpha.rename_column('vel_pew_err1', 'pew_err')


 
 
# In[11]:


tbdata_cachito = tbdata_Ha_cach.copy()
tbdata_cachito.remove_rows([3, 8, 9, 10, 11, 12, 13]) #row 3 is anomolously high and cachito isn't reliable after row 7
tbdata_cachito.remove_columns(['vel1', 'vel_err_left_1', 'vel_err_right_1', 'vel_pew_1', 'vel_pew_err1'])
tbdata_cachito.rename_column('vel0', 'velocity')
tbdata_cachito.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_cachito.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_cachito.rename_column('vel_pew_0', 'pew')
tbdata_cachito.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[12]:


#tbdata_feII.remove_columns(['vel1', 'vel_err_left_1', 'vel_err_right_1', 'vel_pew_1', 'vel_pew_err1'])
tbdata_feII.rename_column('vel0', 'velocity')
tbdata_feII.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_feII.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_feII.rename_column('vel_pew_0', 'pew')
tbdata_feII.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[13]:


tbdata_HB.rename_column('vel0', 'velocity')
tbdata_HB.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_HB.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_HB.rename_column('vel_pew_0', 'pew')
tbdata_HB.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[14]:


tbdata_NaI.rename_column('vel0', 'velocity')
tbdata_NaI.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_NaI.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_NaI.rename_column('vel_pew_0', 'pew')
tbdata_NaI.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[15]:


tbdata_OI.rename_column('vel0', 'velocity')
tbdata_OI.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_OI.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_OI.rename_column('vel_pew_0', 'pew')
tbdata_OI.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[16]:


tbdata_CaII.remove_columns(['vel1', 'vel_err_left_1', 'vel_err_right_1', 'vel_pew_1', 'vel_pew_err1',
                           'vel2', 'vel_err_left_2', 'vel_err_right_2', 'vel_pew_2', 'vel_pew_err2'])
tbdata_CaII.rename_column('vel0', 'velocity')
tbdata_CaII.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_CaII.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_CaII.rename_column('vel_pew_0', 'pew')
tbdata_CaII.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[17]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])
fig = plt.figure()
plt.ion()
#plt.figure(figsize=[5,5])
plot_list = ['HA-cachito', 'HB','CaII','NaI','OI','ScII5662','ScII6262','ScII5526','FeII_multi','SiII']
for i, feature in enumerate(plot_list):
    if feature == 'HA-cachito':
        feature_tbdata = tbdata_halpha
        label = r'H$\alpha$'
    elif feature == 'SiII':
        feature_tbdata = tbdata_cachito
        label='SiII 6355'
    elif feature == 'FeII_multi':
        feature_tbdata = tbdata_feII
        label = 'FeII 5169'
    elif feature =='HB':
        feature_tbdata = tbdata_HB
        label = r'H$\beta$'
    elif feature == 'NaI':
        feature_tbdata = tbdata_NaI
        label = 'NaI 5898'
    elif feature == 'OI':
        feature_tbdata = tbdata_OI
        label = 'OI 7774'
    elif feature == 'CaII':
        feature_tbdata = tbdata_CaII
        label = 'CaII 8498'
    else:
        feature_tbdata = tbdata[tbdata['line']==feature]
        label = '{} {}'.format(feature[:-4], feature[-4:])
    if len(feature_tbdata) != 0:
        dates = Time(feature_tbdata['date'])
        sort_indx = np.argsort(dates)
        dates = dates[sort_indx]
    else:
        print(feature)
    v = (((feature_tbdata['velocity'] - rest_wavelengths[feature])/rest_wavelengths[feature])*c.c.to(u.km/u.s))
    v = v[sort_indx]
    #plt.plot((dates-sn15oz.jdexpl).value+i*0.1, np.abs(v.value), '.', color=colors[feature], label=feature)
    v_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths[feature])*feature_tbdata['vel_err_left'])**2)
    v_err_l = v_err_l[sort_indx]
    v_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths[feature])*feature_tbdata['vel_err_right'])**2)
    v_err_r = v_err_r[sort_indx]
    l = plt.errorbar((dates-sn15oz.jdexpl).value+i*0.1, np.abs(v.value)/1000, yerr=(v_err_l.value/1000, v_err_r.value/1000), 
                  linestyle='none', capsize=0, alpha=0.25, lw=1.0)
    if feature != 'HB':
        plt.plot((dates-sn15oz.jdexpl).value+i*0.1, np.abs(v.value)/1000, 
                       marker=markers[feature], label=label,
                      color=l[0].get_color(), markeredgecolor='k', markersize=3, alpha=1.0, lw=0.5)
    else:
        plt.plot((dates-sn15oz.jdexpl).value+i*0.1, np.abs(v.value)/1000, 
                       marker=markers[feature], label=label,
                      color=l[0].get_color(), markeredgecolor='k', markersize=3, alpha=1.0, zorder=100, lw=0.5)
    
    
    

plt.legend(loc='best', ncol=2, bbox_to_anchor=[0.52, 0.55, 0.5, 0.5])
plt.xlabel('Phase (day)')
plt.ylabel(r'Velocity (1000 $\rm km\,s^{-1}$)')
plt.savefig(os.path.join(FIG_DIR,'line_velocity.pdf'))


 
# # Compare to Guitierrez 2017

 
# In[18]:


tbdata_halpha = asc.read(os.path.join(DATA_DIR, 'HA.tab'))
tbdata_halpha.rename_column('vel0', 'velocity')
tbdata_halpha.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_halpha.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_halpha.rename_column('vel_pew_0', 'pew')
tbdata_halpha.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[19]:


tbdata_HB = asc.read(os.path.join(DATA_DIR, 'HB.tab'))
tbdata_HB.rename_column('vel0', 'velocity')
tbdata_HB.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_HB.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_HB.rename_column('vel_pew_0', 'pew')
tbdata_HB.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[20]:


#tbdata_feII = asc.read(os.path.join(DATA_DIR, 'FeII_multi.tab'))
tbdata_feII = asc.read(os.path.join(DATA_DIR, 'FeII_5169.tab'))
#tbdata_feII.remove_columns(['vel1', 'vel_err_left_1', 'vel_err_right_1', 'vel_pew_1', 'vel_pew_err1'])
tbdata_feII.rename_column('vel0', 'velocity')
tbdata_feII.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_feII.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_feII.rename_column('vel_pew_0', 'pew')
tbdata_feII.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[21]:


v_halpha = (((tbdata_halpha['velocity'] - rest_wavelengths['HA-cachito'])/rest_wavelengths['HA-cachito'])*c.c.to(u.km/u.s))
v_feII = (((tbdata_feII['velocity'] - rest_wavelengths['FeII_multi'])/rest_wavelengths['FeII_multi'])*c.c.to(u.km/u.s))
v_hbeta = (((tbdata_HB['velocity'] - rest_wavelengths['HB'])/rest_wavelengths['HB'])*c.c.to(u.km/u.s))

v_feII_calc = 0.805*v_hbeta
v_feII_full = np.empty(len(tbdata_HB))
early_indx = Time(tbdata_HB['date']) < np.min(Time(tbdata_feII['date']))
late_indx = Time(tbdata_HB['date']) >= np.min(Time(tbdata_feII['date'])) 
v_feII_full[early_indx] = v_feII_calc[early_indx]
v_feII_full[late_indx] = v_feII
dates_full = Time(tbdata_HB['date'])
v50_phot = np.abs(np.interp(50, (dates_full-sn15oz.jdexpl).value, v_feII_full))


 
# ## Velocity spread error

 
# In[22]:


v_halpha_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HA-cachito'])*tbdata_halpha['vel_err_left'])**2)
v_halpha_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HA-cachito'])*tbdata_halpha['vel_err_right'])**2)

v_feII_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*tbdata_feII['vel_err_left'])**2)
v_feII_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*tbdata_feII['vel_err_right'])**2)

v_hbeta_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*tbdata_HB['vel_err_left'])**2)
v_hbeta_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*tbdata_HB['vel_err_right'])**2)


 
# ## 2 $\AA$ error

 
# In[23]:


delta_wave = np.ones(len(v_halpha))*2.0
v_halpha_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HA-cachito'])*delta_wave)**2)
v_halpha_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HA-cachito'])*delta_wave)**2)
delta_wave = np.ones(len(v_feII))*2.0
v_feII_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*delta_wave)**2)
v_feII_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*delta_wave)**2)
delta_wave = np.ones(len(v_hbeta))*2.0
v_hbeta_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*delta_wave)**2)
v_hbeta_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*delta_wave)**2)


 
# ## Guitierrez Empirical

 
# In[24]:


t_gutierrez = np.array([4,8.6,12.8,18.1,23.1,27.7,33.1,38.1,42.8,47.8,53.1,58.6,63.3,68,72.8,78.2,83.5,87.5,93.3,98.2,103,108.2,115.7])
vhalpha_gutierrez = np.array([12845,10702 ,9468 ,8987 ,7798 ,8369 ,7745 ,7478 ,6551 ,7145 ,6535 ,6615 ,6619 ,6613 ,6720 ,6061 ,6490 ,6197 ,5666 ,5788 ,4422 ,5625 ,5805])
vhalpha_err_gutierrez = np.array([ 950 , 138,1787 ,1430 ,1966 ,1825 ,1194 ,1548 ,1745 ,1772 ,1755 ,1900 ,1565 ,1038 ,1564 ,1660 ,1883 ,1930 ,2016 ,2073 ,1353 ,1226 ,1176])

vhbeta_gutierrez = np.array([11379 ,9605 ,8748 ,8364 ,7083 ,7426 ,6668 ,6297 ,5267 ,5541 ,5004 ,5025 ,4836 ,4909 ,4725 ,4460 ,4386 ,4448 ,4261 ,4372 ,3014 ,4128 ,4536 ])
vhbeta_err_gutierrez = np.array([1800 ,1574 ,1653 ,1478 ,1690 ,1527 ,1260 ,1582 ,1600 ,1688 ,1785 ,1774 ,1548 ,1146 ,1665 ,1553 ,1492 ,1514 ,1493 ,1505 , 993 , 885 ,1025])

vfeII_gutierrez = np.array([10447,6871 ,6300 ,5274 ,5440 ,4942,4428 ,3760 ,3938,3537,3631,3401,3397,3374,3078,3074,3253,2946,2476,2119,2625,2451])
vfeII_err_gutierrez = np.array([ 1050 ,2234 ,1174 ,1254 ,1098 , 892 ,1065 ,1045 , 990 , 851 , 973 , 758 , 639 , 828 , 820 , 919 , 785 , 718 , 631 , 525 , 457 , 679])

v50_phot_gutierrez = np.interp(50, t_gutierrez[1:], vfeII_gutierrez)


 
# ## Plot in velocity (without normalizing by V50)

 
# In[25]:


plt.style.use(['seaborn-paper', 'az-paper-twocol'])
fig=plt.figure()
width = 0.27
left = 0.08
height = 0.75
bottom = 0.16
pad = 0.05

ax1 = plt.axes([left, bottom, width, height])
ax2 = plt.axes([left+width+pad, bottom, width, height], sharex=ax1, sharey=ax1)
ax3 = plt.axes([left+2*width+2*pad, bottom, width, height], sharex=ax1, sharey=ax1)
dates_halpha = Time(tbdata_halpha['date'])



l2, = ax1.plot(t_gutierrez, vhalpha_gutierrez/1000)
ax1.fill_between(t_gutierrez, 
                 (vhalpha_gutierrez+vhalpha_err_gutierrez)/1000,
                 (vhalpha_gutierrez-vhalpha_err_gutierrez)/1000, 
                 alpha=0.25, color=l2.get_color())
l1, = ax1.plot((dates_halpha-sn15oz.jdexpl).value, np.abs(v_halpha.value)/1000, 
                  linestyle='none', 
                  marker=markers['HA-cachito'], label=r'H-$\alpha$', zorder=10)
ax1.set_xlim(0, 120)
ax1.set_title(r'H $\alpha$')
ax1.set_ylabel(r'Velocity (1000 $\rm km s^{-1}$)')
#ax1.set_ylim(0, 2.5)

dates_hbeta = Time(tbdata_HB['date'])
ax2.plot((dates_hbeta-sn15oz.jdexpl).value, np.abs(v_hbeta.value)/1000,
                 linestyle='none', marker=markers['HB'], 
                 label=r'H-$\beta$')

#ax2.errorbar((dates_hbeta-sn15oz.jdexpl).value, np.abs(v_hbeta.value),
#                 yerr = ((v_hbeta_err_l).value, (v_hbeta_err_r).value),
#                 color='white', linestyle='none', zorder=9)
#

ax2.plot(t_gutierrez,  vhbeta_gutierrez/1000, color=l2.get_color())
ax2.fill_between(t_gutierrez, 
                 (vhbeta_gutierrez+vhbeta_err_gutierrez)/1000,
                 (vhbeta_gutierrez-vhbeta_err_gutierrez)/1000, 
                 alpha=0.25, color=l2.get_color() )
ax2.plot((dates_hbeta-sn15oz.jdexpl).value, np.abs(v_hbeta.value)/1000,
                 linestyle='none', marker=markers['HB'], 
                 label=r'H-$\beta$', zorder=10, color=l1.get_color())

ax2.set_title(r'H $\beta$')
ax2.set_xlabel('Phase (day)')

dates_feII = Time(tbdata_feII['date'])


ax3.set_title('Fe II')
ax3.plot(t_gutierrez[1:],  vfeII_gutierrez/1000, color=l2.get_color())
ax3.fill_between(t_gutierrez[1:], 
                 (vfeII_gutierrez+vfeII_err_gutierrez)/1000,
                 (vfeII_gutierrez-vfeII_err_gutierrez)/1000, 
                 alpha=0.25, color=l2.get_color() )
ax3.plot((dates_feII-sn15oz.jdexpl).value, np.abs(v_feII.value)/1000, 
                  linestyle='none', 
                 marker=markers['FeII_multi'], label=r'FeII 5159', zorder=10, color=l1.get_color())

plt.savefig(os.path.join(FIG_DIR, 'velocity_compare_guitierrez.pdf'))


 
# # Compare to SNEC

 
# In[26]:


#Multi gaussian fit
tbdata_feII_multi = asc.read(os.path.join(DATA_DIR, 'FeII_multi.tab'))
tbdata_feII_multi.remove_columns(['vel1', 'vel_err_left_1', 'vel_err_right_1', 'vel_pew_1', 'vel_pew_err1'])
tbdata_feII_multi.rename_column('vel0', 'velocity')
tbdata_feII_multi.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_feII_multi.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_feII_multi.rename_column('vel_pew_0', 'pew')
tbdata_feII_multi.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[27]:


#Fit centered on 5169 line
tbdata_feII = asc.read(os.path.join(DATA_DIR, 'FeII_5169.tab'))
tbdata_feII.rename_column('vel0', 'velocity')
tbdata_feII.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_feII.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_feII.rename_column('vel_pew_0', 'pew')
tbdata_feII.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[28]:


tbdata_HB = asc.read(os.path.join(DATA_DIR, 'HB.tab'))
tbdata_HB.rename_column('vel0', 'velocity')
tbdata_HB.rename_column('vel_err_left_0', 'vel_err_left')
tbdata_HB.rename_column('vel_err_right_0', 'vel_err_right')
tbdata_HB.rename_column('vel_pew_0', 'pew')
tbdata_HB.rename_column('vel_pew_err0', 'pew_err')


 
 
# In[29]:


snec_dir = '../../../snec_models/Ni_mass_0.0830/Ni_mixing_5.0/M18.0/E_1.400/K_10.0/R_2400/'
tbdata_vel = asc.read(os.path.join(snec_dir, 'Data/vel_photo.dat'), names=['time', 'velocity'])
ofile_info = open(os.path.join(snec_dir, 'Data/info.dat'))

model_offset = -4 #days - best fit phase offset


 
 
# In[30]:


for iline in ofile_info:
    sline = iline.split('=')
    if "Time of breakout" in sline[0]:
        breakout_time = float(sline[1].strip().split(' ')[0])*u.second


 
 
# In[31]:


tbdata_tau_sob = asc.read(os.path.join(DATA_DIR, 'vel_tau_sob_1.0.dat'), names=['time', 'velocity'])
phase_tau_sob = (tbdata_tau_sob['time']*u.second - breakout_time).to(u.day)
velocity_tau_sob = (tbdata_tau_sob['velocity']*u.km/u.second)/1000
start_phase_tau_sob = Time('2015-09-13') - Time(sn15oz.jdexpl, format='jd')
indx_tau_sob = phase_tau_sob > start_phase_tau_sob


 
# ### Infer FeII from H-beta following Faran 2014

 
# In[32]:


v_feII = (((tbdata_feII['velocity'] - rest_wavelengths['FeII_multi'])/rest_wavelengths['FeII_multi'])*c.c.to(u.km/u.s))
v_feII_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*tbdata_feII['vel_err_left'])**2)
v_feII_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*tbdata_feII['vel_err_right'])**2)
dates_feII = Time(tbdata_feII['date'])

v_hbeta = (((tbdata_HB['velocity'] - rest_wavelengths['HB'])/rest_wavelengths['HB'])*c.c.to(u.km/u.s))
v_hbeta_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*tbdata_HB['vel_err_left'])**2)
v_hbeta_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*tbdata_HB['vel_err_right'])**2)

v_feII_calc = 0.805*v_hbeta
v_feII_err_l_calc = 0.805*v_hbeta_err_l
v_feII_err_r_calc = 0.805*v_hbeta_err_r
v_feII_full = np.empty(len(tbdata_HB))
v_feII_err_l_full = np.empty(len(tbdata_HB))
v_feII_err_r_full = np.empty(len(tbdata_HB))
early_indx = Time(tbdata_HB['date']) < np.min(Time(tbdata_feII['date']))
late_indx = Time(tbdata_HB['date']) >= np.min(Time(tbdata_feII['date'])) 
v_feII_full[early_indx] = v_feII_calc[early_indx]
v_feII_full[late_indx] = v_feII
v_feII_err_l_full[early_indx] = v_feII_err_l_calc[early_indx]
v_feII_err_r_full[early_indx] = v_feII_err_r_calc[early_indx]
v_feII_err_l_full[late_indx] = v_feII_err_l
v_feII_err_r_full[late_indx] = v_feII_err_r
dates_full = Time(tbdata_HB['date'])
v50_phot = np.abs(np.interp(50, (dates_full-sn15oz.jdexpl).value, v_feII_full))


 
 
# In[33]:


v_feII_multi = (((tbdata_feII_multi['velocity'] - rest_wavelengths['FeII_multi'])/rest_wavelengths['FeII_multi'])*c.c.to(u.km/u.s))
v_feII_err_l_multi = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*tbdata_feII_multi['vel_err_left'])**2)
v_feII_err_r_multi = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['FeII_multi'])*tbdata_feII_multi['vel_err_right'])**2)

v_hbeta = (((tbdata_HB['velocity'] - rest_wavelengths['HB'])/rest_wavelengths['HB'])*c.c.to(u.km/u.s))
v_hbeta_err_l = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*tbdata_HB['vel_err_left'])**2)
v_hbeta_err_r = np.sqrt(((c.c.to(u.km/u.s)/rest_wavelengths['HB'])*tbdata_HB['vel_err_right'])**2)

v_feII_calc_multi = 0.805*v_hbeta
v_feII_calc_err_l_multi = 0.805*v_hbeta
v_feII_full_multi = np.empty(len(tbdata_HB))
early_indx_multi = Time(tbdata_HB['date']) < np.min(Time(tbdata_feII_multi['date']))
late_indx_multi = Time(tbdata_HB['date']) >= np.min(Time(tbdata_feII_multi['date'])) 
v_feII_full_multi[early_indx_multi] = v_feII_calc_multi[early_indx_multi]
v_feII_full_multi[late_indx_multi] = v_feII_multi
dates_full_multi = Time(tbdata_HB['date'])
v50_phot_multi = np.abs(np.interp(50, (dates_full_multi-sn15oz.jdexpl).value, v_feII_full_multi))


 
# # Per Referee's request, use Si II rather than Hb for FeII velocity at early times

 
# In[34]:


new_fit_cachito = asc.read(os.path.join(TEST_FILE_DIR, 'cachito.tab'))
phase_cachito = (Time(new_fit_cachito['date'])-sn15oz.jdexpl).value
velocity_siII = -1*((new_fit_cachito['vel0']-rest_wavelengths['SiII'])/rest_wavelengths['SiII'])*c.c.to(u.km/u.s)
visible_indx = phase_cachito < 40 #Keep only the fits that we believe
len(phase_cachito[visible_indx])


 
 
# In[35]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])
fig = plt.figure()
fig.subplotpars.update(left=0.22)
ax = fig.add_subplot(1,1,1)


l1, = ax.plot(phase_cachito[visible_indx]+4, (velocity_siII[visible_indx])/1000, 's', label='Si II Velocity (Photospheric)')
ax.plot((dates_feII-sn15oz.jdexpl).value+4, (v_feII*-1)/1000, 'o', label='Fe II Velocity (Photospheric)')
ax.plot((tbdata_vel['time']*u.second - breakout_time).to(u.day), (tbdata_vel['velocity']*u.cm/u.second).to(u.km/u.second)/1000, label = 'SNEC Photospheric Velocity')
ax.plot(phase_tau_sob[indx_tau_sob], velocity_tau_sob[indx_tau_sob], label = 'SNEC Velocity for $\\tau_{sob}=1.0$')
ax.set_xlim(0, 100)
#ax.plot((dates_full-sn15oz.jdexpl).value+4, (v_feII_full*-1)/1000, 'o', label='Fe II Velocity (Photospheric)')


ax.plot(())
ax.set_xlabel('Phase (day)')
ax.set_ylabel(r'Velocity (1000 $\rm km\,s^{-1}$)')
ax.legend(bbox_to_anchor=[0.15, 0.2, 0.4, 0.2])
fig.savefig('snec_velocity_comp.pdf')


 
 
# In[36]:


from astropy.table import Table
tbdata = Table([(dates_full -sn15oz.jdexpl).value, v_feII_full*-1], names = ['phase', 'v_phot'])
tbdata.write('phot_velocity_table.csv', overwrite=True)


 