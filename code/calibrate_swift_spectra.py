import os
import glob

from astropy.io import fits

from uvotpy.uvotgetspec import getSpec

ra = 289.889800 
dec = -33.767000

calfile='/Users/bostroem/Desktop/research/not_my_code/uvotpy/calfiles/swugu0160wcal20041120v002.fits'

obsid_base = '000340400'
observations = range(1, 12)
os.chdir('../data/swiftuvot/reduced')

original_dir = os.getcwd()
for obs in observations:
    if obs != 3:
        obsid = '{}{:2}'.format(obsid_base, obs).replace(' ', '0')
        print('=======================\n{}\n=======================\n'.format(obsid))
        os.chdir(os.path.join(obsid, 'uvot', 'image'))
        flist = glob.glob('*dt.img*')
        lenticular_flist = glob.glob('sw*uw*.img*')
        if len(flist) != 0:
            grism_ext = len(fits.open(flist[0]))-1
            for ifile in lenticular_flist:
                next = len(fits.open(ifile))-1
                if 'uw1' in ifile:
                    uw1_next = next
                else:
                    uw1_next = None
                if 'uw2' in ifile:
                    uw2_next = next
                else:
                    uw2_next = next

            getSpec(ra, dec, obsid, grism_ext, 
                    fit_second=False, 
                    clobber=True,
                    optimal_extraction=False,
                    interactive=True,
                    background_lower=[20, 40], 
                    background_upper=[10, 30], 
                    calfile=calfile,
                    lfilt1_ext=uw1_next,
                    lfilt2_ext=uw2_next)

        os.chdir(original_dir)