import os
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import glob
import shutil

def calc_wavelength(header, pixels):
    assert header['CTYPE1'] == 'LINEAR'
    CRVAL1 = header['CRVAL1']
    CRPIX1 = header['CRPIX1']
    CD1_1 = header['CD1_1']
    
    wavelength = CRVAL1 + CD1_1*(pixels - CRPIX1)
    return wavelength

def efosc_file_info(directory, search_str = 'EFOSC*.fits'):
	flist = glob.glob(os.path.join(directory, search_str))
	tbdata = Table(data=None, names = ['filename', 'date-obs', 'target', 'obs-type', 'filter', 'slit', 'grism'], dtype = ['S35', Time, 'S20', 'S20', 'S20', 'S20', 'S20'])
	for ifile in flist:
		hdr = fits.getheader(ifile)
		date_obs = Time(hdr['date-obs'], format='isot')
		tbdata.add_row([os.path.basename(ifile), date_obs, hdr['object'], hdr['HIERARCH ESO DPR TECH'], 
			hdr['HIERARCH ESO INS FILT1 NAME'], hdr['HIERARCH ESO INS SLIT1 NAME'], hdr['HIERARCH ESO INS GRIS1 NAME']])
	tbdata.sort('date-obs')
	return tbdata

def append_to_filename(search_str, new_str):
	flist = glob.glob(search_str)
	for ifile in flist:
		shutil.move(ifile, ifile.replace('.fits', '_{}.fits'.format(new_str)))
