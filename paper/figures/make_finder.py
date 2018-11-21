 
# coding: utf-8

 
# In[1]:


import os

import numpy as np

from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
#get_ipython().run_line_magic('matplotlib', 'inline')
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

from utilities_az import visualization


 
 
# In[2]:


DATA_DIR = '../../data/images'
filename_lco = 'lsc1m005-kb78-20150903-0153-e90.fits'
filename_dss = 'dss_search.fits'


 
 
# In[4]:


ofile = fits.open(os.path.join(DATA_DIR, filename_lco))
hdr = ofile[0].header
img = ofile[0].data
ra = hdr['RA']
dec = hdr['DEC']
coord = SkyCoord(ra, dec, unit=(u.hourangle, u.degree))
hdr_wcs = WCS(hdr)
x, y = hdr_wcs.all_world2pix(coord.ra, coord.dec, 0)


 
 
# In[5]:


ofile_dss = fits.open(os.path.join(DATA_DIR, 'dss_search.fits'))

hdr_dss = ofile_dss[0].header
img_dss = ofile_dss[0].data

hdr_wcs_dss = WCS(hdr_dss)
x_dss, y_dss = hdr_wcs_dss.all_world2pix(coord.ra, coord.dec, 0)


 
 
# In[6]:


plt.style.use(['seaborn-paper','az-paper-onecol'])
fig = plt.figure()
fig.set_figheight(5)
fig.subplotpars.update(left=0.27, bottom=0.1, top=0.97)
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)

#Make cutout in LCO pixel coordinates
cutout_size_pix = 100
left = x-cutout_size_pix
right = x+ cutout_size_pix
bottom = y-cutout_size_pix
top = y+cutout_size_pix

#Make cross in LCO pixel coordinates
cross_size = 10
cross_offset = 10
left_cross =   x- cross_offset - cross_size
right_cross =  x- cross_offset
bottom_cross = y- cross_offset - cross_size
top_cross =    y - cross_offset

#Convert cutout from LCO pixel coordinates to DSS pixel coordinates
edge_ra, edge_dec = hdr_wcs.all_pix2world((left, right), (top, bottom), 0)
edge_x_dss, edge_y_dss = hdr_wcs_dss.all_world2pix(edge_ra, edge_dec, 0)
edge_x_dss.sort()
edge_y_dss.sort()

#Convert cross from LCO pixel coordinates to DSS pixel coordinates
cross_ra, cross_dec = hdr_wcs.all_pix2world((left_cross, right_cross), (top_cross, bottom_cross), 0)
cross_x_dss, cross_y_dss = hdr_wcs_dss.all_world2pix(cross_ra, cross_dec, 0)

vmin, vmax = visualization.zscale(img[int(bottom):int(top), int(left):int(right)])
vmin_dss, vmax_dss = visualization.zscale(img_dss[int(edge_y_dss[0]):int(edge_y_dss[1]), int(edge_x_dss[0]):int(edge_x_dss[1])])

ax1.imshow(img_dss,vmin=vmin_dss, vmax=vmax_dss, cmap='bone_r')
ax1.plot(cross_x_dss, (y_dss, y_dss), c='#CC6677')
ax1.plot((x_dss, x_dss), cross_y_dss, c='#CC6677')
ax1.set_xlim(*edge_x_dss)
ax1.set_ylim(*edge_y_dss)

ax2.plot((left_cross, right_cross), (y, y), c='#CC6677')
ax2.plot((x, x), (top_cross, bottom_cross), c='#CC6677')
ax2.imshow(img, vmin=vmin, vmax=vmax, interpolation='nearest', cmap='bone_r')
ax2.set_xlim(left, right)
ax2.set_ylim(bottom, top)

#Convert ticklabels from pixels to degrees
xticks_pix = ax2.get_xticks()[1:-1:2]
yticks_pix = ax2.get_yticks()[1:-1:2]
xticks_world, yticks_world = hdr_wcs.all_pix2world(xticks_pix, yticks_pix, 0)
xticklabels_world = ['{:3.2f}'.format(i) for i in xticks_world]
ax2.set_xticklabels(xticklabels_world)
ax2.set_xticks(xticks_pix)
yticklabels_world = ['{:3.2f}'.format(i) for i in yticks_world]
ax2.set_yticklabels(yticklabels_world)
ax2.set_yticks(yticks_pix)

xticks_pix_dss, yticks_pix_dss = hdr_wcs_dss.all_world2pix(xticks_world, yticks_world, 0)
ax1.set_xticks(xticks_pix_dss)
ax1.set_yticks(yticks_pix_dss[::-1])
ax1.invert_yaxis()
ax1.set_xticklabels(xticklabels_world)
ax1.set_yticklabels(yticklabels_world[::-1])

ax2.set_xlabel('RA (degree)')
ax2.set_ylabel('Dec (degree)')
ax1.set_ylabel('Dec (degree)')
plt.savefig('finder_chart.pdf')


 
 
# In[7]:


vmin, vmin_dss
vmax, vmax_dss


 
 
# In[8]:


ra_test, dec_test = hdr_wcs.all_pix2world((x-100, x, x+100), (y-100, y, y+100), 0)
print(((ra_test[0] - ra_test[-1])*u.degree).to(u.arcsec))
print(((dec_test[0] - dec_test[-1])*u.degree).to(u.arcsec))
c1 = SkyCoord(ra_test[0]*u.degree, dec_test[0]*u.degree, 0)
c2 = SkyCoord(ra_test[-1]*u.degree, dec_test[0]*u.degree, 0)
c3 = SkyCoord(ra_test[0]*u.degree, dec_test[-1]*u.degree, 0)
print(c1.separation(c2).to(u.arcsec))
print(c1.separation(c2).to(u.arcsec))


 
 
# In[9]:


Time(hdr['date-obs'])- Time(2457261.5, format='jd')


 