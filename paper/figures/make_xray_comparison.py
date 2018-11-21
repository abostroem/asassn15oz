 
# coding: utf-8

# Creates:
# * xray_comp.pdf

 
# In[7]:


from astropy.table import Table
from astropy.io import ascii as asc
from astropy.time import Time
from matplotlib import pyplot as plt
import numpy as np
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib as mpl
from utilities_az import supernova


 
 
# In[8]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[9]:


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


 
 
# In[10]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz_phase = Time(['2015-09-05', '2015-11-18']) - Time(sn15oz.jdexpl, format='jd')
#sn15oz_phase = Time(['2015-09-05', '2015-11-18']) - Time('2015-08-27')


 
 
# In[11]:


fig = plt.figure()
fig.subplotpars.update(left=0.15)
snnames = set(tbdata['name'])
leg_lines = []
leg_labels = []
for snname in snnames:
    indx = (tbdata['name']==snname) & (tbdata['age']<100)
    if (indx==True).any():
        sntype = tbdata[indx]['sntype'][0]

        if 'L' in sntype:
            iline = plt.errorbar(tbdata[indx]['age'], tbdata[indx]['luminosity'], fmt='.',
                         yerr=[tbdata[indx]['lumErrH'], tbdata[indx]['lumErrL']], ls='--')
            leg_lines.append(iline[0])
            leg_labels.append('{}, {}'.format(snname, tbdata[indx]['sntype'][0]))
            plt.errorbar(tbdata[indx]['age'][tbdata[indx]['isUpperBound']], tbdata[indx]['luminosity'][tbdata[indx]['isUpperBound']], fmt='.',
                         yerr=[0.25]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]),
                         uplims = [True]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]), ecolor=iline[0].get_color())
        elif 'n' not in sntype:
            iline = plt.errorbar(tbdata[indx]['age'], tbdata[indx]['luminosity'], fmt='d',
                         yerr=[tbdata[indx]['lumErrH'], tbdata[indx]['lumErrL']], ls=':')
            leg_lines.append(iline[0])
            leg_labels.append('{}, {}'.format(snname, tbdata[indx]['sntype'][0]))
            plt.errorbar(tbdata[indx]['age'][tbdata[indx]['isUpperBound']], tbdata[indx]['luminosity'][tbdata[indx]['isUpperBound']], fmt='.',
                         yerr=[0.25]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]),
                         uplims = [True]*len(tbdata[indx]['age'][tbdata[indx]['isUpperBound']]), ecolor=iline[0].get_color())
iline = plt.errorbar(np.array([np.mean(sn15oz_phase).value]), np.array([3.72]), fmt='s',
            yerr=np.array([0.25]), uplims=np.array([True]))
leg_lines.append(iline[0])
leg_labels.append('ASASSN-15oz, IIL')
plt.xlim(0,100)
plt.xlabel('Phase (day)')
plt.ylabel('Luminosity (x10$^{39}$ erg/s)')
plt.legend(leg_lines, leg_labels)
plt.savefig('xray_comp.pdf')


 