 
# coding: utf-8

# Creates:
#     * spectra_montage.pdf

 
# In[10]:


import os

from astropy.io import fits
from astropy.io import ascii as asc
from astropy.time import Time
import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib

from utilities_az import spectroscopy as spec, visualization, supernova

#plt.ioff()


 
 
# In[2]:


plt.style.use(['seaborn-paper', 'az-paper-twocol'])


 
# This is the only way I could get the bars on my capital I's, I'm not sure why it works, but I tried a million other things that didn't work

 
# In[3]:


import matplotlib as mpl
font = mpl.font_manager.FontProperties()
font.set_family('cursive')
font.set_name(mpl.rcParams['font.monospace'][6:])
font.set_size('5')


 
 
# In[4]:


cols = [(0,0,0)]
for x in np.linspace(0,1, 254):
    rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
    gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
    bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
    cols.append((rcol, gcol, bcol))

cols.append((1,1,1))
cm_rainbow = matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)


 
# # Make Spectra Montage

 
# In[5]:


def read_iraf_spectrum(filename):
    ofile = fits.open(filename)
    flux = ofile[0].data[0,0,:]
    wave = spec.calc_wavelength(ofile[0].header, np.arange(len(flux))+1)
    rest_wave = spec.apply_redshift(wave, 0.0069)
    return(spec.spectrum1d(rest_wave, flux))


 
 
# In[6]:


DATA_DIR_LCO = '../../data/spectra/lco'
DATA_DIR_EFOSC = '../../data/spectra/EFOSC'
DATA_DIR_XSHOOT = '../../data/spectra/xshooter/'


 
 
# In[7]:


spectra_files = [
         ('asassn15oz_20150904_redblu_122216.314.fits', DATA_DIR_LCO),
         ('asassn-15oz_20150906_redblu_105042.698a.fits', DATA_DIR_LCO),
         ('asassn15oz_20150907_redblu_123835.277.fits', DATA_DIR_LCO),
         ('asassn15oz_20150911_redblu_105336.349.fits', DATA_DIR_LCO),
         ('asassn15oz_20150916_redblu_120911.274.fits', DATA_DIR_LCO),
         #('asassn-15oz_20150920_redblu_135034.512.fits',DATA_DIR_LCO),
         ('ASASSN15oz_VLT_20150921.txt', DATA_DIR_XSHOOT),
         ('asassn-15oz_20150924_redblu_123847.580.fits',DATA_DIR_LCO),
         ('asassn-15oz_20150930_redblu_122858.217.fits',DATA_DIR_LCO),
         ('tASASSN-15oz_20151003_Gr13_Free_slit1.0_57720_1_e.fits',DATA_DIR_EFOSC),
         #('asassn15oz_20151006_redblu_101906.800.fits', DATA_DIR_LCO),
         #('asassn15oz_20151014_redblu_112918.305.fits', DATA_DIR_LCO),
         ('asassn-15oz_20151025_redblu_102221.833.fits', DATA_DIR_LCO),
         #('asassn-15oz_20151107_redblu_101210.833.fits', DATA_DIR_LCO),
         ('tASAS-SN_15oz_20151107_Gr13_Free_slit1.5_57723_1_e.fits', DATA_DIR_EFOSC),
         ('tASAS-SN_15oz_20151118_Gr13_Free_slit1.0_57723_1_e.fits', DATA_DIR_EFOSC)
                ]


 
 
# In[84]:


line_list_bottom = [
            ('OI',        9100, -30),
            ('OI',        7615, -30), 
            ('CaII',      8300, -30), 
            ('CaII',      8475, -30),
            ('ScII',      6100, -30), 
            ('NaI',       5770, -30), 
            ('ScII',      5560, -30), 
            ('ScII',      5415, -30),
            ('FeI',       5303, -30), 
            ('ScII/FeII', 5170, -35),
            ('FeII/TiII', 5075, -35), 
            ('BaI/FeII',  4455, -20), 
]


 
 
# In[154]:


line_list_top = [
                    (r'H-$\alpha$', 6350, 30),
                    (r'H-$\beta$',  4690, 30),
                    ('cachito',     6145, 30)
                   ]


 
 
# In[155]:


sn15oz = supernova.LightCurve2('asassn-15oz')
jdexpl = Time(sn15oz.jdexpl, format='jd', out_subfmt='date')


 
 
# In[156]:


fig = plt.figure()
fig.set_figheight(6)
fig.subplotpars.update(left=0.06, bottom=0.1, top=0.97)
ax = fig.add_subplot(1,1,1)
#offset = np.append(np.arange(0,5)*0.09, np.arange(5,len(spectra_files))*0.15)[::-1]
offset = (np.arange(0,len(spectra_files))*0.9)[::-1]
#offset[5] = offset[5]*1.1
ref_spec = read_iraf_spectrum(os.path.join(spectra_files[0][1], spectra_files[0][0]))
#Consider making this all one color
colors = visualization.make_color_wheel(spectra_files+[1], cmap=cm_rainbow)
for indx, ifile in enumerate(spectra_files):
    filename, idir = ifile
    if idir == DATA_DIR_XSHOOT:
        tbdata = asc.read(os.path.join(idir, filename), names=['wave', 'flux', 'err'])
        ispec = spec.spectrum1d(tbdata['wave'], tbdata['flux'], error=tbdata['err'])
        date = Time('2015-09-21', out_subfmt='date')
    else:
        ispec = read_iraf_spectrum(os.path.join(idir, filename))
        date = Time(fits.getval(os.path.join(idir, filename), 'date-obs', 0), out_subfmt='date')
    scale_ispec = spec.scale_spectra(ispec, ref_spec)
    plt_indx = (scale_ispec.wave<9200) & (scale_ispec.wave >3300)
    plt.plot(scale_ispec.wave[plt_indx], np.log10(scale_ispec.flux[plt_indx]/10**-13)+offset[indx], color=colors[indx])
    phase_str = (date-jdexpl).value
    phase_yloc = np.log10(scale_ispec.flux[plt_indx][-2]/10**-13)+offset[indx]
    if np.isfinite(phase_yloc):
        ax.text(9300, phase_yloc, '{:3.2f}d'.format(phase_str))
    else:
        print(date)

ax.set_xlim(3200, 10200)
ax.set_ylim(-4,10.5)
ax.set_yticks([])
ax.set_xlabel(r'Rest Wavelength ($\rm \AA$)')
ax.set_ylabel(r'log(Flux) + offset')


#-----------------------
for ilabel, iwave, text_offset in line_list_top:
    closest_wl_indx = np.argmin(np.abs(ref_spec.wave-iwave))
    line_flux = np.log10(ref_spec.flux[closest_wl_indx]/10**-13)+offset[0]+0.1
    plt.annotate(
        ilabel,
        xy=(iwave, line_flux), xytext=(0, text_offset), 
        arrowprops=dict(arrowstyle='-'),
        rotation='vertical', 
        xycoords='data', 
        textcoords='offset points', 
        ha='center', 
        fontproperties=font)
#-----------------------    
for ilabel, iwave, text_offset in line_list_bottom:
    closest_wl_indx = np.argmin(np.abs(scale_ispec.wave-iwave))
    line_flux = np.log10(scale_ispec.flux[closest_wl_indx]/10**-13)-0.1
    plt.annotate(
        ilabel,
        xy=(iwave, line_flux), xytext=(0, text_offset), 
        arrowprops=dict(arrowstyle='-'),
        rotation='vertical', 
        xycoords='data', 
        textcoords='offset points', 
        ha='center', 
        fontproperties=font)
plt.savefig('spectra_montage.pdf')


 