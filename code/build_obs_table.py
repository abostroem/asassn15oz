import os
import glob
from astropy.io import fits, ascii
from astropy.table import Table, Column

obsid_base = '000340400'
observations = range(1, 12)
os.chdir('../data/swiftuvot/reduced')
original_dir = os.getcwd()
tbdata = Table(names=['visit', 'date-obs', 'exptime','expnum', 'PA', 'uw1', 'uw2', 'obsid'],
              dtype=('S2', 'S19', 'f8', 'i4', 'i4', 'bool', 'bool', 'S11'))
for obs in observations:
    if obs != 3:
        obsid = '{}{:2}'.format(obsid_base, obs).replace(' ', '0')
        print('=======================\n{}\n=======================\n'.format(obsid))
        os.chdir(os.path.join(obsid, 'uvot', 'image'))
        flist = glob.glob('*pha*')
        if len(flist)> 0:
            if (os.path.exists('sw{}uw1_ex.img'.format(obsid)) or 
                os.path.exists('sw{}uw1_ex.img.gz'.format(obsid))):
                uw1 = True
            else:
                uw1 = False
            if (os.path.exists('sw{}uw2_ex.img'.format(obsid)) or 
                os.path.exists('sw{}uw2_ex.img.gz'.format(obsid))):
                uw2 = True
            else:
                uw2 = False
            for ifile in flist:
                expnum = int(ifile.split('_')[-2])
                hdr = fits.getheader(ifile, 1)
                exptime = hdr['EXPOSURE']
                dateobs = hdr['DATE-OBS']
                visit = obsid[-2:]
                try:
                    roll_angle = fits.getval('sw{}ugu_dt.img'.format(obsid), 'PA_PNT', 0)
                except IOError:
                    roll_angle = fits.getval('sw{}ugu_dt.img.gz'.format(obsid), 'PA_PNT', 0)
                tbdata.add_row([visit, dateobs, exptime, expnum, int(round(roll_angle,0)), uw1, uw2, obsid])
        os.chdir(original_dir)
tbdata.write('../obs_table.txt', format = 'ascii.fixed_width', overwrite=True)
    

