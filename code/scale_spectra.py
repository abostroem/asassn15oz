import os
import numpy as np

from astropy.io import fits
from astropy.table import Table

from utilities_az import supernova
from utilities_az import spectroscopy as spec

GMOS_DIR = '../data/spectra/gmos'
EFOSC_DIR = '../data/spectra/EFOSC'

flist = ['tASASSN-15oz_20160410_Gr13_Free_slit1.5_57723_1_e.fits', 
         'tASASSN-15oz_20160802_Gr13_Free_slit1.0_57723_1_e.fits']
dates = [None, None]
path_flist = [EFOSC_DIR, EFOSC_DIR]

ra = 289.889800 #from SNEX
dec = -33.767000 #from SNEX
redshift = 0.0069

cur_dir = os.getcwd()


sn = supernova.LightCurve2('asassn-15oz')
#for date, idir, ifile in zip(dates,path_flist, flist):
#    os.chdir(idir)
#    ofile = fits.open(ifile)
#    flux = ofile[0].data
#    if len(flux.shape) == 3:  #EFOSC file
#        flux = flux[0,0,:]
#    else: #GMOS
#        flux = flux*10**-15
#    pixels = np.arange(len(flux))+1
#    wl = spec.calc_wavelength(ofile[0].header, pixels)
#    #deredshift
#    rest_wl = spec.apply_redshift(wl, redshift=redshift)
#    #deredden
#    obs_spec = spec.spectrum1d(rest_wl, flux)
#    obs_spec_dered = spec.correct_for_galactic_extinction(obs_spec, sn.ebv_mw, R_V=3.1)
#    obs_spec_dered = spec.correct_for_galactic_extinction(obs_spec_dered, sn.ebv_host, R_V=3.1)
#    #scale to photometry
#    filter_dir = '/Users/bostroem/Dropbox/DLT40_all/script/scalespectra/filters/'
#    tbdata = Table([obs_spec_dered.wave, obs_spec_dered.flux], names=['wave', 'flux'])
#    output_filename = ifile.replace('-', '_').split('.')[0] #iraf doesn't like -
#    tbdata.write('{}_rest_dustcorr.dat'.format(output_filename), format='ascii.no_header', overwrite=True)
#    if not date:
#        date = fits.getval(ifile, 'date-obs', 0)
#    supernova.scale_spectra_quba('asassn-15oz', '{}_rest_dustcorr.dat'.format(output_filename), filter_dir=filter_dir, date_kw=date, max_sep=7, header_date=False)
#
#    os.chdir(cur_dir)
    
files = ['comb20160609_R400.fits', 'comb20160610_R400.fits', 'comb20160612_B600.fits']
os.chdir(GMOS_DIR)
ofile1 = fits.open(files[0])
ofile2 = fits.open(files[1])
ofile3 = fits.open(files[2])
flux1 = ofile1[0].data*10**-15
flux2 = ofile2[0].data*10**-15
flux3 = ofile3[0].data*10**-15
pixels1 = np.arange(len(flux1))+1
pixels2 = np.arange(len(flux2))+1
pixels3 = np.arange(len(flux3))+1
wl1 = spec.calc_wavelength(ofile1[0].header, pixels1)
wl2 = spec.calc_wavelength(ofile2[0].header, pixels2)
wl3 = spec.calc_wavelength(ofile3[0].header, pixels3)
#deredshift
rest_wl1 = spec.apply_redshift(wl1, redshift=redshift)
rest_wl2 = spec.apply_redshift(wl2, redshift=redshift)
rest_wl3 = spec.apply_redshift(wl3, redshift=redshift)
#scale and merge
obs_spec1 = spec.spectrum1d(rest_wl1, flux1)
obs_spec2 = spec.spectrum1d(rest_wl2, flux2)
obs_spec3 = spec.spectrum1d(rest_wl3, flux3)
obs_spec2_interp = spec.spectrum1d(obs_spec1.wave, np.interp(obs_spec1.wave, obs_spec2.wave, obs_spec2.flux))
scale_obs_spec3 = spec.spectrum1d(obs_spec3.wave, obs_spec3.flux-0.049*10**-15)
scale_obs_spec3_interp = spec.spectrum1d(obs_spec1.wave, np.interp(obs_spec1.wave, scale_obs_spec3.wave, scale_obs_spec3.flux))
overlap_min = 6000
overlap_max = scale_obs_spec3.wave[-5]
seg1_indx = scale_obs_spec3.wave<overlap_min
seg1_wl = scale_obs_spec3.wave[seg1_indx]
seg1_flux = scale_obs_spec3.flux[seg1_indx] 
seg2_indx = (obs_spec2.wave>=overlap_min)&(obs_spec2.wave<=overlap_max)
seg2_wl = obs_spec1.wave[seg2_indx]
seg2_flux = np.mean(np.vstack((obs_spec1.flux[seg2_indx], obs_spec2_interp.flux[seg2_indx], scale_obs_spec3_interp.flux[seg2_indx])), axis=0)
seg3_indx = obs_spec1.wave>overlap_max
seg3_wl = obs_spec1.wave[seg3_indx]
seg3_flux = np.mean(np.vstack((obs_spec1.flux[seg3_indx], obs_spec2_interp.flux[seg3_indx])), axis=0)
merge_spec = spec.spectrum1d(np.hstack((seg1_wl, seg2_wl, seg3_wl)), np.hstack((seg1_flux, seg2_flux, seg3_flux)))
tbdata = Table([merge_spec.wave, merge_spec.flux], names=['wave', 'flux'])
tbdata.write('gmos_merge.csv', overwrite=True)
#deredden
obs_spec_dered1 = spec.correct_for_galactic_extinction(merge_spec, sn.ebv_mw, R_V=3.1)
obs_spec_dered1 = spec.correct_for_galactic_extinction(obs_spec_dered1, sn.ebv_host, R_V=3.1)
#scale to photometry
filter_dir = '/Users/bostroem/Dropbox/DLT40_all/script/scalespectra/filters/'
tbdata1 = Table([obs_spec_dered1.wave, obs_spec_dered1.flux], names=['wave', 'flux'])
output_filename1 = 'gmos_merge' #iraf doesn't like -
tbdata1.write('{}_rest_dustcorr.dat'.format(output_filename1), format='ascii.no_header', overwrite=True)
date1 = '2016-06-11' 
supernova.scale_spectra_quba('asassn-15oz', '{}_rest_dustcorr.dat'.format(output_filename1), filter_dir=filter_dir, date_kw=date1, max_sep=7, header_date=False)
