 
# coding: utf-8

# creates: 
#    * full_lc.pdf

 
# In[1]:


import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl
import define_filters
import supernova
import visualization

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


offsets = [8, 6, 4.5, 3.5, 2.75, 1.5, 0.5, 0, -0.5, -1, -1.5, -2, -3, -3]


 
 
# In[6]:


rainbow_cm = visualization.make_rainbow_cm()
colors = visualization.make_color_wheel(offsets, cmap=rainbow_cm)


 
 
# In[11]:


fig = plt.figure()
fig.set_figheight(6)
fig.subplotpars.update(bottom=0.08, top=0.98)
ax = fig.add_subplot(1,1,1)
for ioffset, icolor, iband in zip(offsets, colors, sort_bands):
    ax.errorbar(sn15oz.phase[iband], sn15oz.apparent_mag[iband]+ioffset, sn15oz.apparent_mag_err[iband],
                ecolor=icolor, capsize=2, linestyle='none')
    ax.plot(sn15oz.phase[iband], sn15oz.apparent_mag[iband]+ioffset, 'o', color=icolor, label='{}+{:1.2f}'.format(iband, ioffset))
ax.invert_yaxis()
ax.legend(loc='best', ncol=2)
ax.set_xlim(0, 425)
ax.set_ylim(30, 10)
ax.set_xlabel('Phase (day)')
ax.set_ylabel('Apparent Magnitude + offset')
plt.savefig('full_lc.pdf')


 