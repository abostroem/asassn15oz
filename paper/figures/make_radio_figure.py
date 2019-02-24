 
# coding: utf-8

 
# In[1]:


import os

from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

from astropy.io import ascii as asc


 
 
# In[2]:


plt.style.use('az-paper-onecol')


 
 
# In[3]:


RADIO_DIR = '../../data/radio/'


 
# ##### From Assaf Horesh:

 
# In[4]:


high_freq_22 = (28, 0.1)
low_freq_4_8 =((28, 42), (0.122, 0.21))
high_freq_15 = (42, 0.08)


 
 
# In[5]:


radio_tbdata = asc.read(os.path.join(RADIO_DIR, '15oz_fit_lines.txt'),names=['x1', 'y1', 'x2', 'y2', 'x3', 'y3'], data_start=1)


 
 
# In[21]:


fig = plt.figure()
fig.subplotpars.update(left=0.23)
ax = fig.add_subplot(1,1,1)
l1,a,b = ax.errorbar(high_freq_22[0], high_freq_22[1], high_freq_22[1]*0.2, uplims=True, label= '22 GHz limit')
l2, = ax.loglog(high_freq_15[0], high_freq_15[1], 's', label= '15 GHz detection')
l3, = ax.loglog(low_freq_4_8[0], low_freq_4_8[1], 'o', label='4.8 GHz detection')

l4, = ax.loglog(radio_tbdata['x1'], radio_tbdata['y1'], ':', label= '22 GHz fit', color=l1.get_color())
l5, = ax.loglog(radio_tbdata['x2'], radio_tbdata['y2'], ':', label= '15 GHz fit', color=l2.get_color())
l6, = ax.loglog(radio_tbdata['x3'], radio_tbdata['y3'], ':', label='4.8 GHz fit', color=l3.get_color())
ax.set_ylim(10**-3, 10**0)
ax.set_xlim(0.8, 100)
ax.set_xlabel('Time (day)')
ax.set_ylabel(r'$\rm F_{\nu}$ (mJy)')
#plt.legend(ncol=2, loc='lower right', framealpha=1)
plt.legend([l1,l2,l3,l4,l5,l6], ['22 GHz limit',l2.get_label(),l3.get_label(),l4.get_label(),l5.get_label(),l6.get_label()], 
           bbox_to_anchor=[0.04, 0.48, 0.5, 0.5], framealpha=0)
plt.savefig('radio_fig.pdf')


 
 
# In[20]:





 