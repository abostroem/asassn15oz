 
# coding: utf-8

# creates: 
#    * full_lc.pdf

 
# In[1]:


import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl
from utilities_az import define_filters, supernova, visualization

#get_ipython().run_line_magic('matplotlib', 'inline')


 
 
# In[2]:


plt.style.use(['seaborn-paper','az-paper-onecol'])


 
 
# In[3]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz.get_photometry()


 
 
# In[4]:


bands = np.array([iband for iband in sn15oz.jd.keys()])
filters = define_filters.define_filters()
cenwave = []
for iband in bands:
    cenwave.append(filters[iband][2])
cenwave = np.array(cenwave)
sort_indx = np.argsort(cenwave)
sort_bands = bands[sort_indx]


 
 
# In[5]:


offsets = [8, 6, 4.5, 3.5, 2.75, 1.5, 0.5, 0, -0.5, -1, -1.5, -2, -3, -3.25, -4.0, -4.25, -4.5]


 
 
# In[6]:


rainbow_cm = visualization.make_rainbow_cm()
colors = visualization.make_color_wheel(offsets, cmap=rainbow_cm)


 
 
# In[10]:


fig = plt.figure()
fig.set_figheight(5.5)
fig.subplotpars.update(bottom=0.08, top=0.98, right=0.77, left=0.22)
ax = fig.add_subplot(1,1,1)
for ioffset, icolor, iband in zip(offsets, colors, sort_bands):
    if iband in ['J', 'H', 'K']:
        symbol = 'd'
    elif iband.isupper():
        symbol = 's'
    elif iband in ['uw2', 'um2', 'uw1', 'us', 'bs', 'vs']:
        symbol = '^'
    else:
        symbol = 'o'
    ax.errorbar(sn15oz.phase[iband], sn15oz.apparent_mag[iband]+ioffset, sn15oz.apparent_mag_err[iband],
                ecolor=icolor, capsize=2, linestyle='none')
    if ioffset < 0:
        ax.plot(sn15oz.phase[iband], sn15oz.apparent_mag[iband]+ioffset, symbol, color=icolor, mec='k', label='{}-{:1.2f}'.format(iband, abs(ioffset)))
    else:
        ax.plot(sn15oz.phase[iband], sn15oz.apparent_mag[iband]+ioffset, symbol, color=icolor, mec='k', label='{}+{:1.2f}'.format(iband, ioffset))
ax.invert_yaxis()
leg = ax.legend(loc='best', ncol=1, bbox_to_anchor=[0.535, 0.61, 0.2, 0.3], framealpha=1)
ax.set_xlim(0, 425)
ax.set_ylim(30, 7.6)
ax.set_xlabel('Phase (day)')
ax.set_ylabel('Apparent magnitude + offset')
yticks = ax.get_yticks()
abs_mag, abs_mag_err = supernova.calc_abs_mag(yticks, 
                                              sn15oz.dist_mod, 
                                              sn15oz.A_mw['V'],
                                              A_host=sn15oz.A_host['V'],
                                              )
ax_abs_mag = ax.twinx()
ax_abs_mag.set_ylim(30, 7.6)
ax_abs_mag.set_yticklabels(['{:2.1f}'.format(i) for i in abs_mag])
ax_abs_mag.set_ylabel('Absolute magnitdue + offset')
leg.set_zorder(1000)
plt.savefig('full_lc.pdf')


 