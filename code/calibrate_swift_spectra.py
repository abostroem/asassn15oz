import os
import glob
from argparse import ArgumentParser
import pickle

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.io.fits import Column, ColDefs

from uvotpy import uvotgetspec
from uvotpy.uvotgetspec import getSpec
#from uvotpy.uvotspec import sum_PHAspectra





#matplotlib.rcParams['image.cmap'] = 'bone'
uvotgetspec.give_result = True #Return things from getSpec


def check_environment():
    #Check that environment variables are set:
    environ = os.environ
    assert environ.has_key('CALDB'), 'CALDB path must be set before running'
    assert environ.has_key('HEADAS'), 'HEADAS path must be set before running'
    assert environ.has_key('UVOTPY'), 'UVOTPY path must be set before running'
    assert environ.has_key('FTOOLSOUTPUT'), 'please run heainit before running'

def save_background_data(Y0, Y1, obsid, ext):
    '''
    Save background
    '''
    bkg_img = Y0[5]
    bkg_upper = Y1[1][1]
    bkg_lower = Y1[1][2]
    bkg_mean = Y1[1][0]
    col1 = Column(name='bkg_lower', array=bkg_lower, format='f8')
    col2 = Column(name='bkg_upper', array=bkg_upper, format='f8')
    col3 = Column(name='bkg_mean', array=bkg_mean, format='f8')
    cols = ColDefs([col1, col2, col3])
    hdu1 = fits.PrimaryHDU([bkg_img])
    hdu2 = fits.BinTableHDU.from_columns(cols)
    hdulist = fits.HDUList([hdu1, hdu2])
    hdulist.writeto('bkg_{}_{}.fits'.format(obsid, ext), overwrite=True)

def pickle_output(Y0, Y1, Y2, Y3, Y4, obsid, ext):
    '''
    (specfile, lfilt1_, lfilt1_ext_, lfilt2_, lfilt2_ext_, attfile), (method), \
        (Xphi, Yphi, date1), (dist12, ankerimg, ZOpos), expmap, bgimg, bg_limits_used, \
         bgextra = Y0
    (dis,spnet,angle,anker,anker2,anker_field,ank_c), (bg,bg1,bg2,extimg,spimg,spnetimg,offset), 
        (C_1,C_2,img),  hdr,m1,m2,aa,wav1 ) = Y1 
         
    fit,(coef0,coef1,coef2,coef3),(bg_zeroth,bg_first,bg_second,bg_third),\
        (borderup,borderdown),apercorr,expospec=Y2 
    counts, variance, borderup, borderdown, (fractions,cnts,vars,newsigmas) = Y3
    wav2p, dis2p, flux2p, qual2p, dist12p = Y4[0]
    '''
    
    output_dict = {'Y0': Y0, 'Y1': Y1, 'Y2': Y2, 'Y3': Y3, 'Y4': Y4}
    
    ofile = open('getspec_output_{}_{}.pkl'.format(obsid, ext), 'wb')
    pickle.dump(output_dict, ofile)

def run_getspec(ra, dec, observation, obsid, ext, interactive=False, yoffset=None, bkg_lower=[None, None],
                bkg_upper=[None, None], bkg_templ=None, uw1_next=None, uw2_next=None):
    calfile='/Users/bostroem/Desktop/research/not_my_code/uvotpy/calfiles/swugu0160wcal20041120v002.fits'
    if yoffset is not None:
        if observation in yoffset.keys():
            print('using Y offset')
            (Y0, Y1, Y2, Y3, Y4) = getSpec(ra, dec, obsid, ext+1, 
                fit_second=False, 
                clobber=True,
                optimal_extraction=False,
                interactive=True,
                background_lower=bkg_lower, 
                background_upper=bkg_upper,
                background_template=bkg_templ, 
                calfile=calfile,
                lfilt1_ext=uw1_next,
                lfilt2_ext=uw2_next,
                offsetlimit=yoffset[observation][ext],
                chatter=3)

        else:
            (Y0, Y1, Y2, Y3, Y4) = getSpec(ra, dec, obsid, ext+1, 
                                            fit_second=False, 
                                            clobber=True,
                                            optimal_extraction=False,
                                            interactive=interactive,
                                            background_lower=bkg_lower, 
                                            background_upper=bkg_upper,
                                            background_template=bkg_templ,  
                                            calfile=calfile,
                                            lfilt1_ext=uw1_next,
                                            lfilt2_ext=uw2_next,
                                            chatter=3)
    else:
        (Y0, Y1, Y2, Y3, Y4) = getSpec(ra, dec, obsid, ext+1, 
                                            fit_second=False, 
                                            clobber=True,
                                            optimal_extraction=False,
                                            interactive=interactive,
                                            background_lower=bkg_lower, 
                                            background_upper=bkg_upper, 
                                            background_template=bkg_templ,
                                            calfile=calfile,
                                            lfilt1_ext=uw1_next,
                                            lfilt2_ext=uw2_next,
                                            chatter=3)
    return Y0, Y1, Y2, Y3, Y4
    

def get_lenticular_images(obsid, data_dir):
    '''
    Get number of extensions in lenticular images
    '''
    lenticular_flist = glob.glob(os.path.join(data_dir, obsid, 'uvot', 'image', 'sw*uw*.img*'))
    #Loop over lenticular images
    if len(lenticular_flist) == 0:
        print('WARNING: No lenticular files found in {}'.format(os.path.join(data_dir, obsid, 'uvot', 'image')))
        return None, None
    for ifile in lenticular_flist:
        next = len(fits.open(ifile))-1
        if 'uw1' in ifile:
            uw1_next = next
        else:
            uw1_next = None
        if 'uw2' in ifile:
            uw2_next = next
        else:
            uw2_next = None
    return uw1_next, uw2_next

def calibrate_spectra(interactive=False, yoffset=None, obs_end='', bkg_lower=[None, None], 
                      bkg_upper=[None, None], bkg_templ=None):
    #Loop over visits
    for obs in observations:
        obsid = '{}{:2}'.format(obsid_base, obs).replace(' ', '0')
        uw1_next, uw2_next = get_lenticular_images(obsid, DATA_DIR)
        print('=======================\n{}\n=======================\n'.format(obsid))
        os.chdir(os.path.join(obsid+obs_end, 'uvot', 'image'))
        flist = glob.glob('*dt.img*')
        #If grism files exist
        if len(flist) != 0:
            grism_next = len(fits.open(flist[0]))-1
            print('EXT = {}'.format(grism_next))
            #Loop over lenticular images
            
            pp = PdfPages(os.path.join(FIG_DIR,'{}_plots.pdf'.format(obsid)))
            #Loop over grism extensions
            for ext in range(grism_next):
                Y0, Y1, Y2, Y3, Y4 = run_getspec(ra, dec, obs, obsid, ext, 
                                                 interactive=interactive, 
                                                 yoffset=yoffset, 
                                                 bkg_lower=[None, None],
                                                 bkg_upper=[None, None], 
                                                 bkg_templ=None,
                                                 uw1_next=uw1_next, 
                                                 uw2_next = uw2_next)
                save_background_data(Y0, Y1, obsid, ext+1)
                pickle_output(Y0, Y1, Y2, Y3, Y4, obsid, ext+1)
                if interactive is True:
                    raw_input('Press enter to continue')
                for i in range(3):
                    fig = plt.figure(i+1)
                    fig.suptitle('{}'.format(obsid))
                    pp.savefig(fig)
                    plt.close(fig)
            pp.close()
        os.chdir(original_dir)
        
def combine_spectra(obs_end=''):
    for obs in observations:
        if obs != 3:
            obsid = '{}{:2}'.format(obsid_base, obs).replace(' ', '0')
            print('=======================\n{}\n=======================\n'.format(obsid))
            os.chdir(os.path.join(obsid+obs_end, 'uvot', 'image'))
            flist = glob.glob('*pha*')
            if len(flist) > 1:
                sum_PHAspectra(flist, outfile='{}_x1dsum.fits'.format(obsid))
            os.chdir(original_dir)
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--show", dest='show', action='store_true',
                        help='--show \t [%(default)s], if set interactive plots are displayed')
    args = parser.parse_args()
        
    if args.show:
        plt.ion()
    ra = 289.889800 
    dec = -33.767000

    obsid_base = '000340400'
    #observations = range(1, 12)
    #observations=[1, 5, 9]
    observations=[1, 2] #First visit, only one with extractable signal
    #os.chdir('../data/swiftuvot/reduced_default')
    DATA_DIR = './'
    #os.chdir('../data/swiftuvot/reduced_center_zeroth')

    original_dir = os.getcwd()
    check_environment()
    
    #--------------------------------
    #Default extraction parameters
    #Template matches first extraction
    #--------------------------------
#     FIG_DIR = '../../../../../../figures'
#     FIG_DIR = './'
#     yoffset={1:[[91,1], [91,1], [91,1]], 2:[[80,1], [93,0], [89, 0]]}
#     calibrate_spectra(interactive=args.show, yoffset=yoffset)
    #combine_spectra()
    #--------------------------------
    #Test TRUVOT against original extraction
    #Test default extraction
    # No need for y-offsets, default extraction works well
    #--------------------------------
    os.chdir('../data/swiftuvot/test_truvot_extraction/original') 
    FIG_DIR = './'
    calibrate_spectra(interactive=args.show, yoffset=None)
    #combine_spectra()
    
    
    #--------------------------------
    #Offsets to exclude star just above spectrum at PA=248
    #--------------------------------
    #FIG_DIR = '../../../../../../figures'
    #yoffset={1:[102, 0],5:[99, 0], 9:[102, 0]}
    #uvotgetspec.trackwidth = 1.7
    #calibrate_spectra(interactive=args.show, yoffset=yoffset)
    #combine_spectra()
    #--------------------------------
    #Default extraction  
    #-------------------------------- 
    #FIG_DIR = './' 
    #calibrate_spectra(interactive=args.show, obs_end='_default')
    #combine_spectra(obs_end='_default')
    
    #--------------------------------
    #Offsets to extract start just above spectrum at PA=248
    #--------------------------------
    #FIG_DIR = './' 
    #yoffset_star={1:[118, 0],5:[111, 0], 9:[113, 0]}
    #uvotgetspec.trackwidth = 1.7
    #calibrate_spectra(interactive=args.show, yoffset=yoffset_star, obs_end='_star')
    #combine_spectra(obs_end='_star')
    
    #--------------------------------
    #Try different background regions
    #--------------------------------
    #FIG_DIR = './'
    #yoffset={1:[102, 0],5:[99, 0], 9:[102, 0]}
    #uvotgetspec.trackwidth = 1.7
    #calibrate_spectra(interactive=args.show, yoffset=yoffset, bkg_lower=[25, 75], bkg_upper=[25, 75], obs_end='_diffbkg')
    #combine_spectra()