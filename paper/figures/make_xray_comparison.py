 
# coding: utf-8

# Creates:
# * xray_comp.pdf

 
# In[1]:


from astropy.table import Table
from astropy.io import ascii as asc
import astropy.units as u
from astropy.time import Time
from matplotlib import pyplot as plt
import numpy as np
#get_ipython().run_line_magic('matplotlib', '')
import matplotlib as mpl
from utilities_az import supernova


 
 
# In[2]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
# # Read in SNaX data

 
# In[3]:


#        0       1           2                  3             4              5          6        7           8                                  9              10         11                           12           13           14          15                       16                              17                          18        19          
names=["name","type","Date Exploded (JD)","Coordinates","Distance (Mpc)","galaxy","redshift","redshiftErr","Observation Start Date (JD)","instrument","Age (days)","Flux (10^-13 erg cm^-2 s^-1)","isUpperBound","fluxErrL","fluxErrH","Energy Lower Bound (KeV)","Energy Upper Bound (KeV)","Luminosity (10^39 erg s^-1)","lumErrL","lumErrH","model","dateExplodedRef","redshiftRef","dateExplodedRef","coordsRef","distRef","dateObservedRef","fluxRef", "junk"]
name = []
sntype = []
age = []# (days)
isUpperBound = []
luminosity = [] # (10^39 erg s^-1) 
lumErrL =[]
lumErrH =[]

ofile = open('../../data/xray/SNaX.TSV', 'r')
for iline in ofile:
    if iline.startswith('SN'):
        sline = iline.split('\t')
        name.append(sline[0])
        sntype.append(sline[1])
        age.append(float(sline[10]))
        isUpperBound.append(bool(int(sline[12])))
        luminosity.append(float(sline[17]))
        lumErrL.append(float(sline[18]))
        lumErrH.append(float(sline[19]))
tbdata = Table([name, sntype, age, isUpperBound, luminosity, lumErrL, lumErrH], 
              names=['name', 'sntype', 'age', 'isUpperBound', 'luminosity', 'lumErrL', 'lumErrH'])


 
# # Read in data from Stefano

 
# In[4]:


ofile_ximage = asc.read('../../data/xray/upper_lim_table.txt', data_start=1,
                       names=['visit','date-obs','exptime','obsid','cts/s','Flux','unabsorbed-flux','NH'],
                       format='fixed_width')
distance = (28.83*u.Mpc).to(u.cm)
ofile_ximage['luminosity'] = (4*np.pi*distance**2)*ofile_ximage['unabsorbed-flux']
ofile_ximage['luminosity'][ofile_ximage['luminosity']==0]=np.nan
ofile_ximage['MJD'] = Time(ofile_ximage['date-obs']).mjd


 
 
# In[5]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz_phase = Time(['2015-09-05', '2015-11-18']) - Time(sn15oz.jdexpl, format='jd')


 
 
# In[11]:


fig = plt.figure()
fig.subplotpars.update(left=0.25)
ax = fig.add_subplot(1,1,1)

leg_lines, leg_labels = ax.get_legend_handles_labels()

snnames = set(tbdata['name'])
for snname in snnames:
    indx = (tbdata['name']==snname) & (tbdata['age']<100)
    if (indx==True).any():
        sntype = tbdata[indx]['sntype'][0]

        if 'L' in sntype:
            iline = ax.errorbar(tbdata[indx]['age'], tbdata[indx]['luminosity'], fmt='s',
                         yerr=[tbdata[indx]['lumErrH'], tbdata[indx]['lumErrL']], ls='--')
            leg_lines.append(iline[0])
            leg_labels.append('{}, {}'.format(snname, tbdata[indx]['sntype'][0]))
            
            ax.errorbar(tbdata[indx]['age'][tbdata[indx]['isUpperBound']], tbdata[indx]['luminosity'][tbdata[indx]['isUpperBound']], fmt='s',
                         yerr=[0.25]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]),
                         uplims = [True]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]), ecolor=iline[0].get_color())
            
            
            
        elif 'n' not in sntype:
            iline = ax.errorbar(tbdata[indx]['age'], tbdata[indx]['luminosity'], fmt='d',
                         yerr=[tbdata[indx]['lumErrH'], tbdata[indx]['lumErrL']], ls=':')
            leg_lines.append(iline[0])
            leg_labels.append('{}, {}'.format(snname, tbdata[indx]['sntype'][0]))
            ax.errorbar(tbdata[indx]['age'][tbdata[indx]['isUpperBound']], tbdata[indx]['luminosity'][tbdata[indx]['isUpperBound']], fmt='d',
                         yerr=[0.25]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]),
                         uplims = [True]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]), ecolor=iline[0].get_color())
ax.set_yscale("log")
iline = ax.errorbar(ofile_ximage['MJD']-Time(sn15oz.jdexpl, format='jd').mjd, np.array(ofile_ximage['luminosity'])/10**39, np.array(ofile_ximage['luminosity'])/10**39*0.2, 
                    uplims=True, fmt='.')
leg_lines.append(iline[0])
leg_labels.append('ASASSN-15oz, IIL')
ax.set_xlim(0,160)
#ax.set_ylim(-10, 20)
ax.set_xlabel('Phase (day)')
ax.set_ylabel(r'Luminosity (x10$^{39}$ erg $\rm s^{-1}$)', position=(1,0.38))
leg = ax.legend(leg_lines, leg_labels, bbox_to_anchor=[0.55, 0.65, 0.4, 0.4], framealpha=1.0, frameon=True)
plt.savefig('xray_comp.pdf')


 