#!/usr/bin/env python

import os,string,re,sys
from astropy.io import fits as pyfits
#import pyfits
import pylab
import pylab as plt
from numpy import interp as ninterp
import numpy as np
import glob
from optparse import OptionParser
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
    from pyfits import open as opp
    print 'update header'
    try:
        imm=opp(image,mode='update')
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
            imm=opp(image,mode='update')
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

if __name__ == "__main__":
    parser = OptionParser(usage=usage,description=description, version="%prog agnkey")
    parser.add_option("-o", "--output",dest="output",default='None',type="str",
                      help='-o  output [] \t [%default]')
    parser.add_option("-D", "--directory",dest="directory",default='./',type="str",
                      help='-d  directory [] \t [%default]')

    parser.add_option("-a", "--angle",dest="angle",default='',type="str",
                      help='-a  angle [] \t [%default]')
    parser.add_option("-r", "--ra",dest="ra",default='',type="str",
                      help='-r  ra [] \t [%default]')
    parser.add_option("-d", "--dec",dest="dec",default='',type="str",
                      help='-d  dec [] \t [%default]')
    parser.add_option("-s", "--slit",dest="slit",default='',type="str",
                      help='-s  slit [] \t [%default]')
    
    option,args = parser.parse_args()
    _output=option.output
    _angle=option.angle
    _ra=option.ra
    _dec=option.dec
    _slit=option.slit
    _directory=option.directory
    _show=True

    if len(args)<1 : 
        sys.argv.append('--help')
    if not _ra:
        sys.argv.append('--help')
    if not _dec:
        sys.argv.append('--help')

    option,args = parser.parse_args()
    obsid = args[0]
    #_ra= args[1]
    #_dec= args[2]
    if not _output: 
        _output = re.sub('.img','',obsid)+'_out'


    current=os.getcwd()    
    
    if _directory:
        print 'change directory'
        os.chdir(_directory)

    try:
        _file=glob.glob('sw'+obsid+'*dt*')[0]
        f=pyfits.open(_file)
        ext=[]
        for i in range(2,len(f)): 
            ext.append(i)
    except:
        ext=[2]        
        
    if _show:
        pylab.ion()

    from uvotpy import uvotgetspec
    from uvotpy.uvotgetspec import getSpec
    fileoutput=[]
    os.system('rm '+_output+'*.pha')
    for jj in ext:
        _output1=_output+'_'+str(jj)

        if _slit:     
            uvotgetspec.trackwidth=float(_slit)

        uvotgetspec.getSpec(float(_ra),float(_dec),obsid, jj, fit_second=False, clobber=False,outfile=_output1,optimal_extraction=False,interactive=True,\
                background_lower=[20, 40], background_upper=[10, 30])
        #,spextwidth=13)
        _output1=re.sub('.pha','',glob.glob(_output1+'*pha')[0])
        f = pyfits.open(_output1+'.pha')
        d = f[2].data
        h = f[2].header
        w = d.field('LAMBDA') 
        net=d.field('NETRATE') 
        flux=d.field('FLUX')
        fluxerr=d.field('FLUXERR')
        wave=np.arange(w[0],w[-1],1)
        fluxnew=rebin_spec(w, flux, wave)
        os.system('rm '+_output1+'_1d.fits')
        out_fits = pyfits.PrimaryHDU(header=h,data=fluxnew)
        out_fits.writeto(_output1+'_1d.fits', clobber=True, output_verify='fix')
        dict={'CTYPE1':['LINEAR',''], 'CTYPE2':['LINEAR ',''], 'CRVAL1': [wave[0],''],  'CRPIX1':[1.,''] , 'CDELT1' : [1,''] ,\
              'WAT0_001':['system=equispec',''], 'WAT1_001':['wtype=linear label=Wavelength units=Angstroms',''], 'WAT2_001':['wtype=linear',''],\
              'APNUM1' : ['1 1 1 1','']}
        updateheader(_output1+'_1d.fits',0,dict)
        fileoutput.append(_output1+'_1d.fits')
        if _show:  raw_input('go on')

    _input=','.join(fileoutput)
    _out=_output+'_merge.1d.fits'
    os.system('rm '+_output+'_merge.1d.fits')

    from pyraf import iraf
    iraf.specred()
    iraf.scombine(_input,_out)

    if _directory:
            os.chdir(current)

#    out_fits = pyfits.PrimaryHDU(header=_hd,data=_data)
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
