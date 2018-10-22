 
# coding: utf-8

# Creates:
#     * extinction.pdf

 
# In[2]:


from astropy.io import ascii as asc
import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import spectroscopy as spec
import matplotlib as mpl


 
 
# In[3]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[4]:


filename = '../../data/spectra/xshooter/ASASSN15oz_VLT_20150921.txt'
z = 0.006929
tbdata = asc.read(filename, names=['wavelength', 'flux', 'error'])
spectrum = spec.spectrum1d(tbdata['wavelength'], tbdata['flux'], error=tbdata['error'])


 
 
# In[5]:


sodium_D_lines = [5889.950, 5895.924]
redshifted_sodium = spec.apply_redshift(np.array(sodium_D_lines), z, redden=True )


 
 
# In[6]:


plt.plot(spectrum.wave-0.5, spectrum.flux/10**-15, label='Spectrum')
ylim = plt.ylim()
plt.vlines(redshifted_sodium,*ylim, linestyle=':', label='Galactic extinction', alpha=0.25)
plt.vlines(sodium_D_lines, *ylim, label='Host extinction', alpha=0.25)
plt.xlim(5880, 5950)
plt.ylim(1, 3.25)
plt.ylabel(r'Flux ($\rm x10^{-15} erg/cm^2/s/\AA$)')
plt.xlabel(r'Wavelength ($\rm \AA$)')
plt.legend(loc=(0.24, 0.1))
plt.savefig('extinction.pdf')


 