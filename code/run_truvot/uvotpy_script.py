#!!!! YOU MUST LAUNCH HEASOFT PRIOR TO EXECUTING THIS SCRIPT!!!!
#  ('heainit' is the default method of launching)

#import some python packages
import os
import glob
from uvotpy import uvotgrism, uvotgetspec
from astropy.io import fits
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

import sys
sys.path.append('/Users/bostroem/Desktop/research/asassn15oz/code')
import calibrate_swift_spectra

DATA_DIR = '/Users/bostroem/Desktop/research/asassn15oz/data/swiftuvot/test_truvot_extraction/truvot'
FIG_DIR = DATA_DIR


uvotgetspec.give_new_result = True

def extract_spectra(ra, dec, observation, obsid, ext, uw1_next, uw2_next, yoffset=None):
    CUR_DIR = os.getcwd()
    os.chdir(os.path.join(DATA_DIR, obsid, 'uvot', 'image'))
    #READ IN THE TEMPLATE FOR THIS SNAPSHOT AS AN ARRAY
    e=fits.getdata('sw{}ugu_dt_template_{}.fits'.format(obsid, int(ext+1))) 

    Y0, Y1, Y2, Y3, Y4 = calibrate_swift_spectra.run_getspec(ra, dec, observation, obsid, ext, interactive=True,
                                        yoffset=yoffset,
                                        bkg_templ=e, 
                                        uw1_next=uw1_next,
                                        uw2_next=uw2_next)
    #CONVERT THE UVOTPY OUTPUT SPECTRA INTO A COLUMNATED FORMAT
    thr = fits.open('sw{}ugu_1ord_{}_f.pha'.format(obsid, ext+1))
    data = thr[2].data
    np.savetxt('sw{}ugu_1ord_{}.dat'.format(obsid, ext+1), data)
    os.chdir(CUR_DIR)
    raw_input('Zoom before saving')
    fig1 = plt.figure(1)
    fig2 = plt.figure(2)
    pp = PdfPages(os.path.join(FIG_DIR, '{}_ext{}_extraction.pdf'.format(obsid, ext+1)))
    pp.savefig(fig1)
    pp.savefig(fig2)
    pp.close()
    calibrate_swift_spectra.pickle_output(Y0, Y1, Y2, Y3, Y4, obsid, ext+1)
    
if __name__ == "__main__":
    ra = 289.889800 
    dec = -33.767000  #SET VARIABLES REQUIRED FOR UVOTPY (TARGET LOCATION ON SKY)
    obsid_base = '000340400'
    observations=[1, 2] #First visit is the only one with signal
    original_dir = os.getcwd()
    calibrate_swift_spectra.check_environment()
    #yoffset={1:[[91,1], [91,1], [91,1]], 2:[[80,1], [93,0], [89, 0]]}

    for obs in observations:
        obsid = '{}{:2}'.format(obsid_base, obs).replace(' ', '0')
        print('=======================\n{}\n=======================\n'.format(obsid))
        flist = glob.glob(os.path.join(DATA_DIR, obsid, 'uvot', 'image', '*dt.img*'))
        with fits.open(flist[0]) as ofile:
            next = len(ofile)-1
        uw1_next, uw2_next = calibrate_swift_spectra.get_lenticular_images(obsid, DATA_DIR)
        for ext in range(next):
            extract_spectra(ra, dec, obs, obsid, ext, uw1_next, uw2_next, yoffset=None)

    #yout = extract_spectra('00034040001', 1, DATA_DIR, offset = [98, 1])
    #yout = extract_spectra('00034040001', 2, DATA_DIR, offset=[98,1])
    #yout = extract_spectra('00034040001', 3, DATA_DIR, offset=[98,1])

    #yout = extract_spectra('00034040002', 1, DATA_DIR, offset=[87, 1])
    #yout = extract_spectra('00034040002', 2, DATA_DIR, offset=[100, 1])
    #yout = extract_spectra('00034040002', 3, DATA_DIR, offset=[96, 1])