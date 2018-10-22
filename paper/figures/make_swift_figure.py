 
# coding: utf-8

# Creates:
#     * swift_spectra.pdf

 
# In[2]:


import os
import glob
import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from astropy.time import Time
import astropy.constants as c
import astropy.units as u


 
 
# In[3]:


plt.style.use(['seaborn-paper', 'az-paper-twocol'])


 
# ### Roll Angle 1

 
# In[4]:


SWIFT_DIR = '../../data/swiftuvot/reduced_default/'
obsid_list1 = ['00034040001']
data_list1 = []
exptime_list1 = []
for obsid in obsid_list1:
    flist1 = glob.glob(os.path.join(SWIFT_DIR, obsid, 'uvot', 'image', '*.pha'))
    print(os.path.join(SWIFT_DIR, obsid, 'uvot', 'image', '*.pha'))
    for ifile in flist1:
        
        tbdata = fits.getdata(ifile, 2)
        data_list1.append(tbdata)
        exptime_list1.append(fits.getval(ifile, 'exposure',2))


 
 
# In[5]:


wl1 = data_list1[0]['LAMBDA']
flux_list1 = []
wht_list1 = []
for exptime, tbdata in zip(exptime_list1, data_list1):
    # interpolate to first wavelength range
    # take average weighted by flux error
    wht1 = 1/tbdata['FLUXERR']
    flux_list1.append(np.interp(wl1, tbdata['LAMBDA'], tbdata['FLUX']))
    wht_list1.append(np.interp(wl1, tbdata['LAMBDA'], wht1)*exptime)
flux_arr1 = np.array(flux_list1)
wht_arr1 = np.array(wht_list1)
final_spec1 = np.sum(flux_arr1*wht_arr1, axis=0)/np.sum(wht_arr1, axis=0)
smoothed_signal1 = convolve(final_spec1, Box1DKernel(11))


 
# ### Roll Angle 2

 
# In[6]:


SWIFT_DIR = '../../data/swiftuvot/reduced_default/'
obsid_list2 = ['00034040002']
data_list2 = []
exptime_list2=[]
for obsid in obsid_list2:
    flist2 = glob.glob(os.path.join(SWIFT_DIR, obsid, 'uvot', 'image', '*.pha'))
    print(os.path.join(SWIFT_DIR, obsid, 'uvot', 'image', '*.pha'))
    for ifile in flist2:
        
        tbdata = fits.getdata(ifile, 2)
        data_list2.append(tbdata)
        exptime_list2.append(fits.getval(ifile, 'exposure', 2))


 
 
# In[7]:


wl2 = data_list2[0]['LAMBDA']
flux_list2 = []
wht_list2 = []
for exptime, tbdata in zip(exptime_list2, data_list2):
    # interpolate to first wavelength range
    # take average weighted by flux error
    wht2 = 1/tbdata['FLUXERR']
    flux_list2.append(np.interp(wl2, tbdata['LAMBDA'], tbdata['FLUX']))
    wht_list2.append(np.interp(wl2, tbdata['LAMBDA'], wht2)*exptime)
flux_arr2 = np.array(flux_list2)
wht_arr2 = np.array(wht_list2)
final_spec2 = np.sum(flux_arr2*wht_arr2, axis=0)/np.sum(wht_arr2, axis=0)
smoothed_signal2 = convolve(final_spec2, Box1DKernel(11))


 
# ### Combine Roll Angles

 
# In[8]:


SWIFT_DIR = '../../data/swiftuvot/reduced_default/'
obsid_list = ['00034040001', '00034040002']
data_list = []
exptime_list = []
for obsid in obsid_list:
    flist = glob.glob(os.path.join(SWIFT_DIR, obsid, 'uvot', 'image', '*.pha'))
    print(os.path.join(SWIFT_DIR, obsid, 'uvot', 'image', '*.pha'))
    for ifile in flist:
        tbdata = fits.getdata(ifile, 2)
        data_list.append(tbdata)
        exptime_list.append(fits.getval(ifile, 'exposure', 2))


 
 
# In[9]:


wl = data_list[0]['LAMBDA']
flux_list = []
wht_list = []
for exptime, tbdata in zip(exptime_list,data_list):
    # interpolate to first wavelength range
    # take average weighted by flux error
    # weight by exposure time too
    wht = 1/tbdata['FLUXERR']
    flux_list.append(np.interp(wl, tbdata['LAMBDA'], tbdata['FLUX']))
    wht_list.append(np.interp(wl, tbdata['LAMBDA'], wht)*exptime)
flux_arr = np.array(flux_list)
wht_arr = np.array(wht_list)
total_exptime = np.sum(np.array(exptime_list))
final_spec = np.sum(flux_arr*wht_arr, axis=0)/np.sum(wht_arr, axis=0)
smoothed_signal = convolve(final_spec, Box1DKernel(11))


 
 
# In[12]:


fig = plt.figure()
fig.subplotpars.update(left=0.075, right=0.98)
ax1 = fig.add_subplot(3,1,1)
ax2 = fig.add_subplot(3,1,2, sharey=ax1)
ax3 = fig.add_subplot(3,1,3, sharey=ax1)
#ax.grid()
default_cycler = plt.rcParams['axes.prop_cycle']

ax3.set_xlabel('Wavelength ($\AA$)')
ax2.set_ylabel('Flux (weighted mean)')
#ax1.set_title('Roll Angle 248')
l1 = ax1.plot(wl1,final_spec1/10**-14,   label='Roll Angle 260')
ax1.set_ylim(-0.5,2)
ax1.set_yticks([0, 1, 2])
ax1.legend(loc='upper right')
ax1.set_xticks([])

ax2_cycler = default_cycler[1:2]

#ax2.set_title('Roll Angle 260')
ax2.set_prop_cycle(ax2_cycler)
l2 = ax2.plot(wl2,final_spec2/10**-14,  label='Roll Angle 260')
#ax2.set_ylim(-1E-14,2)
ax2.legend(loc='upper right')
ax2.set_xticks([])

#ax3.set_title('Combined Spectrum')
ax3_cycler = default_cycler[2:3]
ax3.set_prop_cycle(ax3_cycler)
l3 = ax3.plot(wl, final_spec/10**-14 ,  label='Combined Spectrum')
#ax3.set_ylim(-1E-14,2)
ax3.legend(loc='upper right')

#plt.tight_layout()
ax3.set_xlabel(r'Wavelength ($\rm \AA$)')
ax2.set_ylabel(r'Flux ($\rm x10^{-14}$ $\rm erg/cm^2/s/\AA$)')

for ax in [ax1, ax2, ax3]:
    ax.set_xlim(2000, 6700)

plt.savefig('swift_spectra.pdf')


 
# 
# 
