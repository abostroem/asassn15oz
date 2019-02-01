 
# coding: utf-8

 
# In[3]:


import os
import glob

from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

from astropy.io import ascii as asc


 
 
# In[4]:


plt.style.use(['seaborn-paper','az-paper-onecol'])


 
 
# In[5]:


DATA_DIR = '../../data/rho_snapshots/'
FIG_DIR = '.'


 
 
# In[8]:


flist = [ '../../data/rho_snapshots/rho_14d.csv', 
          '../../data/rho_snapshots/rho_43d.csv', 
          '../../data/rho_snapshots/rho_73d.csv',
         '../../data/rho_snapshots/rho_102d.csv']


 
 
# In[10]:


fig = plt.figure()
fig.subplotpars.update(left=0.23, right=0.9)
ax = fig.add_subplot(1,1,1)
for ifile in flist:
    tbdata = asc.read(ifile, names=['velocity', 'logrho'])
    label = os.path.basename(ifile).split('_')[1].split('d')[0]
    ax.plot(tbdata['velocity'][:-5], tbdata['logrho'][:-5], label='{} days'.format(label))
ax.set_xlim(1500, 10000)
ax.set_ylim(-17.25, -9.75)
ax.set_xlabel(r'v (km$\rm s^{-1}$)')
ax.set_ylabel(r'log$\rm _{10}(\rho)$ (g $\rm cm^{-3}$)')
ax.legend()
plt.savefig(os.path.join(FIG_DIR, 'log_rho.pdf'))


 