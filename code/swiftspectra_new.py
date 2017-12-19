#!/usr/bin/env python

import os,string,re,sys
from astropy.io import fits
from matplotlib import pyplot as plt
from numpy import interp as ninterp
import numpy as np
import glob
from argparse import ArgumentParser
from matplotlib.backends.backend_pdf import PdfPages

from uvotpy import uvotgetspec
from uvotpy.uvotgetspec import getSpec

description="> reduce swift spectra  " 
usage= "%prog  image ra dec  \t [listfile]"

def xpa(arg):
    import subprocess
    subproc = subprocess.Popen('xpaset -p ds9 '+arg,shell=True)
    subproc.communicate()

def rebin_spec(wave, specin, wavnew):
    from pysynphot import observation
    from pysynphot import spectrum
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
    return obs.binflux

def updateheader(image,dimension,headerdict):
    print 'update header'
    try:
        imm=fits.open(image,mode='update')
        _header=imm[dimension].header
        for i in headerdict.keys():
           _header.update(i,headerdict[i][0],headerdict[i][1])
        imm.flush()
        imm.close()
    except:
#        from floyds.util import correctcard
#        print '\nwarning: problem to update header, try to correct header format ....'
#        correctcard(image)
        try:
            imm=fits.open(image,mode='update')
            _header=imm[dimension].header
            for i in headerdict.keys():
               _header.update(i,headerdict[i][0],headerdict[i][1])
            _header.update(_headername,_value,commento)
            imm.flush()
            imm.close()
        except:
           print '\n# still problems to update header'
#           import sys
#            sys.exit('error: not possible update header')
#################################################################################################
def main(obsid, ra, dec, output=None, directory = './', angle = None,  slit='', show=False, 
         calfile=None, chatter=0, lfilt1_ext=None):
    '''
    This function runs UVOTPY's getSpec function to extract UV Grism spectra. Given a directory
    all extensions of the sw*dt* files will be extracted. These files are merged into a single spectrum
    using IRAFs scombine.

    Parameters
    ----------
    obsid : str
        the observation ID associated with your observation e.g. '00010130004'
    ra : float
        right ascension of observed target in decimal degrees
    dec : float
        declination of observed target in decimal degrees
    directory : str , optional
        directory where grism and lenticular files are located relative to current directory; default: './'
    angle : float, optional
        angle of extracted trace; passed to getSpec as fixed_angle; default: None; use getSpec default
    slit : str, optional
        extraction box height used to set uvotgetspec.trackwidth; default: None; use getSpec default
    show : bool, optional
        set to True to show interactive plots. default: False
    calfile : str, optional
        path to swugu0160wcal20041120v002.fits calibration file (including filename). default: None; look for file in CALDB directory
    chatter : int, optional
        set the level of output messages; default: 0

    '''
    #Check that environment variables are set:
    environ = os.environ
    assert environ.has_key('CALDB'), 'CALDB path must be set before running'
    assert environ.has_key('HEADAS'), 'HEADAS path must be set before running'
    assert environ.has_key('UVOTPY'), 'UVOTPY path must be set before running'
    assert environ.has_key('FTOOLSOUTPUT'), 'please run heainit before running'
        
    if not output: 
        output = re.sub('.img','',obsid)+'_out'


    current=os.getcwd()    
    
    if directory:
        print 'change directory'
        os.chdir(directory)

    try:
        _file=glob.glob('sw'+obsid+'*dt*')[0]
        f=fits.open(_file)
        ext=[]
        for i in range(2,len(f)): 
            ext.append(i)
    except:
        ext=[2]        
        
    if show:
        plt.ion()


    fileoutput=[]
    os.system('rm '+output+'*.pha')
    pp = PdfPages('{}_plots.pdf'.format(obsid))
    for jj in ext:
        _output1=output+'_'+str(jj)

        if slit:     
            uvotgetspec.trackwidth=float(slit)

        print obsid
        uvotgetspec.getSpec(ra, dec, obsid, jj, fit_second=False, clobber=False,outfile=_output1,optimal_extraction=False,interactive=True,\
                            background_lower=[20, 40], background_upper=[10, 30], calfile=calfile, chatter=chatter, fixed_angle=angle, lfilt1_ext=lfilt1_ext)
        #,spextwidth=13)
        _output1=re.sub('.pha','',glob.glob(_output1+'*pha')[0])
        f = fits.open(_output1+'.pha')
        d = f[2].data
        h = f[2].header
        w = d.field('LAMBDA') 
        net=d.field('NETRATE') 
        flux=d.field('FLUX')
        fluxerr=d.field('FLUXERR')
        wave=np.arange(w[0],w[-1],1)
        #Rebin to 1A scale?
        fluxnew=rebin_spec(w, flux, wave)
        #Write to _1d.fits file?
        os.system('rm '+_output1+'_1d.fits')
        out_fits = fits.PrimaryHDU(header=h,data=fluxnew)
        out_fits.writeto(_output1+'_1d.fits', clobber=True, output_verify='fix')
        #Add wavelength solution in fits header standard?
        dict={'CTYPE1':['LINEAR',''], 'CTYPE2':['LINEAR ',''], 'CRVAL1': [wave[0],''],  'CRPIX1':[1.,''] , 'CDELT1' : [1,''] ,\
              'WAT0_001':['system=equispec',''], 'WAT1_001':['wtype=linear label=Wavelength units=Angstroms',''], 'WAT2_001':['wtype=linear',''],\
              'APNUM1' : ['1 1 1 1','']}
        updateheader(_output1+'_1d.fits',0,dict)
        fileoutput.append(_output1+'_1d.fits')
        if show:
            raw_input('press enter to save and close figures')
        for i in range(3):
            fig = plt.figure(i+1)
            fig.suptitle('{}, ext={}'.format(obsid, jj))
            pp.savefig(fig)
            plt.close(fig)
    #Use scombine to merge all spectra?
    _input=','.join(fileoutput)
    _out=output+'_merge.1d.fits'
    os.system('rm '+output+'_merge.1d.fits')
    pp.close()
    from pyraf import iraf
    from iraf import imred, specred
    specred()
    iraf.scombine(_input,_out)

    if directory:
            os.chdir(current)

#    out_fits = fits.PrimaryHDU(header=_hd,data=_data)
#    out_fits.writeto(_output+'_2d.fits', clobber=True, output_verify='fix')
#    plt.ion()
#    plt.plot(w,net)
#    plt.clf()
#    plt.plot(w,flux)
#    reg=glob.glob(string.split(_infile,'sw')[0]+'/*reg')[0]
#    xpa('file '+_infile)
#    xpa('regions '+reg)
#    xpa('zoom to fit')
#    raw_input('ddd')



if __name__ == "__main__":
    parser = ArgumentParser(description=description)
    parser.add_argument('obsid', type=str, 
                        help='observation ID number')
    parser.add_argument("-o", "--output",dest="output",default='None',type=str,
                        help='-o  output [] \t [%(default)s]')
    parser.add_argument("-D", "--directory",dest="directory",default='./',type=str,
                        help='-D  directory [] \t [%(default)s]')

    parser.add_argument("-a", "--angle",dest="angle",default='',type=str,
                        help='-a  angle [] \t [%(default)s]')
    parser.add_argument("-r", "--ra",dest="ra",default='',type=str,
                        help='-r  ra [] \t [%(default)s]')
    parser.add_argument("-d", "--dec",dest="dec",default='',type=str,
                        help='-d  dec [] \t [%(default)s]')
    parser.add_argument("-s", "--slit",dest="slit",default='',type=str,
                        help='-s  slit [] \t [%(default)s]')
    parser.add_argument("-c", "--calfile", dest='calfile', default='', type=str,
                        help='-c calfile [] \t [%(default)s]')
    parser.add_argument("--show", dest='show', action='store_true',
                        help='--show \t [%(default)s], if set interactive plots are displayed')
    parser.add_argument('-v', '--verbose', dest='chatter', default=0, type=int,
                        help='-v verbose [] \t [%(default)s]')

    args = parser.parse_args()
    args.ra = float(args.ra)
    args.dec = float(args.dec)
    if (not args.ra) or (not args.dec): 
        parser.print_help()    

    if len(args.calfile)==0:
         args.calfile=None
    if len(args.angle)==0:
         args.angle = None
    else:
        args.angle = float(args.angle)
    main(args.obsid, output=args.output, directory=args.directory, angle=args.angle, ra=args.ra, dec=args.dec, 
         slit=args.slit, show=args.show, calfile=args.calfile, chatter=args.chatter)
    #_ra= args[1]
    #_dec= args[2]
    
