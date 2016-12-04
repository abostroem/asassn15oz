import os
import glob
import shutil
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import spectroscopy
import numpy as np



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

def make_fits_bintable(input_filename, output_filename='multi', clobber=False):
    array_HDU = fits.open(input_filename)
    array_table = array_HDU[0].data
    pix = np.arange(max(array_table.shape))+1
    wl = spectroscopy.calc_wavelength(array_HDU[0].header, pix)
    primary_hdu = fits.PrimaryHDU(header=array_HDU[0].header)
    
    columns = fits.ColDefs([fits.Column(wl, name='WAVELENGTH', format='E'), 
                            fits.Column(array_table[0,0,:], name='OPTIMAL_FLUX', format='E'),
                            fits.Column(array_table[1,0,:] , name='FLUX', format='E'),
                            fits.Column(array_table[2,0,:], name='SKY', format='E'),
                            fits.Column(array_table[3,0,:], name='ERROR', format='E')])
    table_hdu = fits.BinTableHDU.from_columns(columns)
    hdu_list = fits.HDUList([primary_hdu, table_hdu])
    
    if output_filename == 'multi':
        output_filename = input_filename.split('.fits')[0]+'multi.fits'
    hdu_list.writeto(output_filename, clobber=clobber)  
    