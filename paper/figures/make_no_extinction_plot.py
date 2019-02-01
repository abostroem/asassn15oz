 
# coding: utf-8

# Creates:
#     * extinction.pdf

 
# In[1]:


from astropy.io import ascii as asc
import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
from utilities_az import spectroscopy as spec
import matplotlib as mpl


 
 
# In[2]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[3]:


filename = '../../data/spectra/xshooter/ASASSN15oz_VLT_20150921.txt'
z = 0.006929
tbdata = asc.read(filename, names=['wavelength', 'flux', 'error'])
spectrum = spec.spectrum1d(tbdata['wavelength'], tbdata['flux'], error=tbdata['error'])


 
 
# In[4]:


sodium_D_lines = [5889.950, 5895.924]
redshifted_sodium = spec.apply_redshift(np.array(sodium_D_lines), z, redden=True )


 
 
# In[9]:


fig = plt.figure()
fig.subplotpars.update(left=0.21)
ax = fig.add_subplot(1,1,1)
ax.plot(spectrum.wave-0.5, spectrum.flux/10**-15, label='Spectrum')
ylim = ax.get_ylim()
ax.vlines(redshifted_sodium,*ylim, linestyle=':', label='Galactic extinction', alpha=0.25)
ax.vlines(sodium_D_lines, *ylim, label='Host extinction', alpha=0.25)
ax.set_xlim(5880, 5950)
ax.set_ylim(1, 3.25)
ax.set_ylabel(r'Flux ($\rm x10^{-15} erg\,cm^{-2}\,s^{-1}\,\AA^{-1}$)', position=(1,0.38))
ax.set_xlabel(r'Wavelength ($\rm \AA$)')
ax.legend(loc=(0.24, 0.1))
plt.savefig('extinction.pdf')


 