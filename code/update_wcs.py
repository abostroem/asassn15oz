import os
import glob
from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from astropy.io import fits
from astropy import wcs

from visualization import zscale

FILTERS = ['uw1', 'uw2', 'ugu', 'um2']
DATA_DIR = '../data/swiftuvot/wcs_fix'
FIG_DIR = '../figures'


def fix_wcs(obsid, obj_ra, obj_dec):
    style.use('seaborn-pastel')
    with PdfPages(os.path.join(FIG_DIR, 'swift_wcs_changes_{}.pdf'.format(obsid))) as pp:
        for ifilter in FILTERS:
            filename1 = os.path.join(DATA_DIR, obsid, 'uvot', 'image', 'sw{}{}_sk.img.gz'.format(obsid, ifilter))
            filename2 = os.path.join(DATA_DIR, obsid, 'uvot', 'image', 'sw{}{}_sk.img'.format(obsid, ifilter))
            if os.path.exists(filename1):
                filename = filename1
            elif os.path.exists(filename2):
                filename = filename2
            else: 
                filename = 'junkxx'
            if os.path.exists(filename):
                # read in data
                ofile = fits.open(filename, mode='update')
                next = len(ofile)
                print(filename)
                print(np.arange(1, next))
                print(ofile.info())
                for iext in np.arange(1,next):
                    # plot it
                    img = ofile[iext].data
                    hdr = ofile[iext].header
                    vmin, vmax = zscale(img)
                    plt.figure()
                    plt.imshow(img, vmin=vmin, vmax=vmax)
                    plt.title('{}, ext={}'.format(os.path.basename(filename), iext))
                    # plot location of SN from wcs
                    img_wcs = wcs.WCS(hdr)
                    x, y = img_wcs.all_world2pix(obj_ra, obj_dec, 0)
                    plt.plot(x, y, 'o')
                    # ask user to zoom
                    adjust_wcs = input('Zoom in on object. Would you like to adjust the WCS? (y), n ')
                    if adjust_wcs !='n':
                        print('Click on actual object location')
                        # Ask user to click on actual SN location
                        x_new,y_new = plt.ginput(1, 0, show_clicks=True)[0]
                        # Calculate difference
                        x_shift = x_new - x
                        y_shift = y_new - y
                        # Update header history
                        hdr.add_history('Updated CRPIX1 from {} to {}'.format(img_wcs.wcs.crpix[0], img_wcs.wcs.crpix[0]+x_shift))
                        hdr.add_history('Updated CRPIX2 from {} to {}'.format(img_wcs.wcs.crpix[1], img_wcs.wcs.crpix[1]+y_shift))
                        img_wcs.wcs.crpix[0] = img_wcs.wcs.crpix[0]+x_shift
                        img_wcs.wcs.crpix[1] = img_wcs.wcs.crpix[1]+y_shift
                        x, y = img_wcs.all_world2pix(obj_ra, obj_dec, 1)
                        plt.plot(x, y, 'o')
                        # Update CRPIX1 and CRPIX2 with offset
                        hdr['CRPIX1'] = img_wcs.wcs.crpix[0]
                        hdr['CRPIX2'] = img_wcs.wcs.crpix[1]
                        update_hdr = input('Update header? (y), n ')
                        if update_hdr != 'n':
                            ofile[iext].header = hdr
                            ofile.flush()
                    pp.savefig()
                ofile.close()
                
if __name__ == "__main__":
    obsid_list = ['00034040001', '00034040002', '00034040005', '00034040007', '00034040009', '00034040011']
    #obsid_list = ['00034040013', '00034040015']
    #obsid_list = ['00034040001']
    RA = 289.889800 
    DEC = -33.767000
    for obsid in obsid_list:
        fix_wcs(obsid, RA, DEC)




