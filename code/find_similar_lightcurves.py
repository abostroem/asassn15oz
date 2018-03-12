#Built in Python
import os
import sys
import glob

#Standard Packages
from astropy.io import ascii
from astropy import table
from astropy.time import Time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d

#Installed for this project
import extinction

#Mine
import visualization
get_ipython().magic('matplotlib inline')

import define_filters
sys.path.append('/Users/bostroem/Desktop/research/not_my_code/stefano')
from alice2 import leggifile as readfile

LC_DIR = '/Users/bostroem/Desktop/research/lightcurves/'
FIG_DIR = '../figures'
ALICE_TABLE_DIR = '../alice_tables'

def calc_extinction(A, obs_band, A_err=0., Rv=3.1, A_band='V', Rb=4.1):
    '''
    Use Cardelli et al (1989) extinction law to calculate extinction
    in magnitudes and error in extinction
    Inputs:
        A: extinction in magnitudes in filter=band - default is V band
        A_err: error in extinction
        Rv: default 3.1
        band: filter that extinction is given in. Right now only V and B are supported
    '''
    assert (A_band=='V') or (A_band=='B'), 'Error: only B and V band are supported'
    
    if A_band == 'B':
        Av = (Rv/Rb)*A
        Av_err = (Rv/Rb)*A_err
    else:
        Av = A
        Av_err = A_err
        
    bandpar = define_filters.define_filters()
    assert obs_band in bandpar.keys(), 'obs_band {} not in defined list: {}'.format(obs_band, bandpar.keys())
    cenwave = float(bandpar[obs_band][2])
    A_lam = extinction.ccm89(np.array([cenwave]), Av, Rv) 
    if Av != 0:
        A_lam_err = (A_lam/Av) * Av_err
    else:
        A_lam_err = np.array([0.0])
    return A_lam, A_lam_err
    
def calc_absmag(mag, dist_mod, ext_band, mag_err, dist_mod_err, 
                 obs_band, A_gal=0, A_gal_err=0, A_host=0, A_host_err=0, 
                 Rv=3.1, Rb=4.1):
    '''
    Define the absolute magnitude of a single number or an array
    applying extinction if needed
    
    Input:
        mag: apparent magnitudes - value or list
        dist_mod: distance modulus - value
        ext_band: filter that extinction is given in (only B and V supported)
        mag_err: err in apparent magnitudes
        dist_mod_err: error in distance modulus
        A_gal: galactic extinction in magnitudes
        A_gal_err: error in galactic extinction in magnitudes
        A_host: host galaxy extinction in magnitudes
        A_host_err: error in host galaxy extinction in magnitudes
        Rv: R for V band
        Rb: R for B band
    
    '''
    A_lam_gal, A_lam_gal_err = calc_extinction(A_gal, obs_band, A_err=A_gal_err, Rv=Rv, A_band=ext_band, Rb=Rb)
    A_lam_host, A_lam_host_err = calc_extinction(A_host, obs_band, A_err=A_host_err, Rv=Rv, A_band=ext_band, Rb=Rb)
    abs_mag = mag - dist_mod - A_lam_gal - A_lam_host
    abs_mag_err_stat = mag_err
    abs_mag_err_sys = np.sqrt(dist_mod_err**2 + A_lam_gal_err**2 + A_lam_host_err)
    
    return abs_mag, abs_mag_err_stat, abs_mag_err_sys

def build_sn_phot_table(LC_DIR, band):
    flist = glob.glob(os.path.join(LC_DIR, '*.dat'))
    tbdata = table.Table(names=['sn', 'type', 
                                'jd_expl', 
                                'jdmax', 'jdmax_err',
                                'magmax', 'magmax_err',
                                'abi', 'abi_err',
                                'abg', 'abg_err',
                                'absmagmax', 'absmagmax_err_stat', 'absmagmax_err_sys'], 
                         dtype=['S15', 'S4',
                                'f8', 
                                'f8', 'f8',
                                'f8', 'f8',
                                'f8', 'f8',
                                'f8', 'f8',
                                'f8', 'f8', 'f8'])
    for ifile in flist:
        if len(os.path.basename(ifile).split('_'))==1: #skip anything with comments in the title
            sndata = readfile(os.path.splitext(ifile)[0])
            if (band in sndata['mag'].keys()):
                if len(sndata['mag'][band]) > 0:
                    if sndata['magmax'][band] < 40:
                        absmagmax, absmagmax_err_stat, absmagmax_err_sys = calc_absmag(sndata['magmax'][band], 
                                                                               sndata['mu'], 
                                                                               'B', 
                                                                               sndata['magmax_err'][band],
                                                                               sndata['mu_err'],
                                                                               band, 
                                                                               A_gal=sndata['abg'],
                                                                               A_gal_err=sndata['abg_err'],
                                                                               A_host=sndata['abi'], 
                                                                               A_host_err=sndata['abi_err'])
                        if float(sndata['jdmax'][band])<100000: #Convert mjd to jd for consistency
                            sndata['jd_expl']=Time(float(sndata['jd_expl']),format='mjd').jd
                            sndata['jdmax'][band]=Time(float(sndata['jdmax'][band]),format='mjd').jd
                            sndata['jdmax_err'][band]=Time(float(sndata['jdmax_err'][band]),format='mjd').jd

                        irow = [sndata['sn'], sndata['sntype'], 
                                float(sndata['jd_expl']),
                                float(sndata['jdmax'][band]), float(sndata['jdmax_err'][band]), 
                                sndata['magmax'][band], sndata['magmax_err'][band],
                                sndata['abi'], sndata['abi_err'],
                                sndata['abg'], sndata['abg_err'],
                                absmagmax[0], absmagmax_err_stat, absmagmax_err_sys[0]]
                        tbdata.add_row(irow)
    if band.upper() == band:
        #Macs don't distinguish between case - so R was overwritting r
        tbdata.write(os.path.join(ALICE_TABLE_DIR, 'alice_sn_database_{}{}.csv'.format(band, band)), overwrite=True)
    else:
        tbdata.write(os.path.join(ALICE_TABLE_DIR, 'alice_sn_database_{}.csv'.format(band)), overwrite=True)

def calc_chisq(days1, days2, spec1, spec2, tlim=None):
    '''
    Calculate the reduced chisq of 2 spectra. These are unevenly sampled. 
    The second set of data is interpolated to the first set of data
    '''
    if tlim is not None:
        tlim_indx1 = days1<=tlim
        tlim_indx2 = days2<=tlim
        days1 = days1[tlim_indx1]
        spec1 = spec1[tlim_indx1]
        days2 = days2[tlim_indx2]
        spec2 = spec2[tlim_indx2]
    if len(days2) > 2 and len(days1) > 2:
        overlap_indx = (days1>=(days2[0])) & (days1 <= (days2[-1]))
        dof = len(days1[overlap_indx])
        if dof < 2:
            reduced_chisq=100
        else:
            interp_func = interp1d(days2, spec2)
            spec2_interp = interp_func(days1[overlap_indx])
            chisq = np.sum((spec1[overlap_indx] - spec2_interp)**2) #think about whether this should include errors, and how to interpolate?
            reduced_chisq = chisq/dof
    else:
        reduced_chisq=100
    return reduced_chisq

def normalize_data(sndata, band):
    '''
    Divide by the maximum magnitude from a light cruve
    Subtract off the date of the maximum from the days of an interpolated light curve
    '''
    galactic_ext, err = calc_extinction(sndata['abg'], band, A_band='B')
    host_ext, err = calc_extinction(sndata['abi'], band, A_band='B')
    
    norm_mag = (np.array(sndata['mag'][band])-galactic_ext-host_ext)/ \
                (np.min(np.float_(np.array(sndata['mag'][band])-galactic_ext-host_ext)))
    try:
        norm_day = np.array(sndata['jd'][band]) - float(sndata['jd_exp'])
    except KeyError:
        norm_day = np.array(sndata['jd'][band]) - float(sndata['jd_expl'])
    return norm_day, norm_mag    

def calc_overlap(day1, day2):
    '''
    Calculate whether 2 arrays overlap. If they do overlap, calculate the number of points
    day2 has contained within the limits of day1
    '''
    
    day1.sort()
    day2.sort()
    if (day1[0]>day2[-1]) or (day2[0]>day1[-1]):
        overlap = False
    else:
        overlap = True

    n_overlap = len(day2[(day2 > day1[0]) & (day2 < day1[-1])])

        
    return overlap, n_overlap

def make_chisq_table(band, comp_sn, tlim=None):
    '''
    Calculate the reduced chisq of a set of SN compared to a single SN
    Returned table is sorted with the smallest values first
    '''
    comp_sndata = readfile(os.path.join(LC_DIR, comp_sn))
    comp_day, comp_mag = normalize_data(comp_sndata, band)
    if band.upper() == band:
        tbdata = ascii.read(os.path.join(ALICE_TABLE_DIR, 'alice_sn_database_{}{}.csv'.format(band, band)))
    else:
        tbdata = ascii.read(os.path.join(ALICE_TABLE_DIR, 'alice_sn_database_{}.csv'.format(band)))
    rchisq_tbdata = table.Table(names=['sn', 'reduced_chisq','npts', 'overlap', 'n_overlap'], 
                                dtype=['S20', 'f8', 'i8', '?', 'i8'])
    for sn in tbdata['sn'][tbdata['type']=='II']:
        if sn != comp_sn:
            try:
                isndata = readfile(os.path.join(LC_DIR, sn))
                iday, imag = normalize_data(isndata, band)
                rchisq1 = calc_chisq(iday, comp_day, imag, comp_mag, tlim=tlim)
                rchisq2 = calc_chisq(comp_day, iday, comp_mag, imag, tlim=tlim)
                rchisq = (rchisq1+rchisq2)/2.
                overlap, n_overlap = calc_overlap(comp_day, iday)
                rchisq_tbdata.add_row([sn, rchisq, len(imag), overlap, n_overlap])
            except FileNotFoundError:
                print("can't open file for {}".format(sn))
    rchisq_tbdata.sort(['reduced_chisq'])
    return rchisq_tbdata

def plot_bestfit(comp_sn, rchisq_tbdata, band, ncomp=5, npts=5):
    comp_sndata = readfile(os.path.join(LC_DIR, comp_sn))
    comp_day, comp_mag = normalize_data(comp_sndata, band)
    fig = plt.figure(figsize=[8,8])
    number_found=0
    for i, irow in enumerate(rchisq_tbdata):
        if number_found == ncomp:
            break
        if (irow['n_overlap'] >= npts):
            number_found+=1
            print('Comparing SN {} to {}, reduced chisq = {:.2E}'.format(irow['sn'],
                                                                     comp_sn,
                                                                     irow['reduced_chisq']))
            isndata = readfile(os.path.join(LC_DIR, str(irow['sn'], 'utf-8')))
            iday, imag = normalize_data(isndata, band)
            ax = fig.add_subplot(ncomp, 1, number_found)
            ax.plot(comp_day, comp_mag, 'o', label=comp_sn)
            ax.plot(iday, imag, 'o', label=str(irow['sn'], 'utf-8'))
            #ax.set_title(str(irow['sn'], 'utf-8'))
            ax.legend(loc='best')
            ax.set_ylim(ax.get_ylim()[::-1])
            ax.set_xlim(comp_day[0]-10, comp_day[-1]+20)

def find_best_fit_sn(rchisq_tbdata, ncomp=5, npts=5):
    number_found=0
    sn_list = []
    for i, irow in enumerate(rchisq_tbdata):
        if number_found == ncomp:
            break
        if (irow['n_overlap'] >= npts):
            number_found+=1
            sn_list.append(str(irow['sn'], 'utf-8'))
    return sn_list

def build_bestfit_dict(comp_sn, npts=10, ncomp=5, tlim=None):
    comp_sndata = readfile(os.path.join(LC_DIR, comp_sn))
    best_fit_dict = {}
    for band in comp_sndata['mag'].keys():
        rchisq_tbdata = make_chisq_table(band, comp_sn, tlim=tlim)
        sn_list = find_best_fit_sn(rchisq_tbdata, npts=npts)
        for sn in sn_list:
            if sn in best_fit_dict.keys():
                best_fit_dict[sn]+=1
            else:
                best_fit_dict[sn]=1
    return best_fit_dict

def build_bestfit_table(best_fit_dict):
    bestfit_tbdata = table.Table(names=['sn', 'npts'], dtype=['S20', 'i8'])
    for sn in best_fit_dict.keys():
        bestfit_tbdata.add_row([sn, best_fit_dict[sn]])
    bestfit_tbdata.sort(['npts'])
    return bestfit_tbdata

def plot_bestfit_all_filters(bestfit_tbdata, comp_sn, ncomp=5):
    pp = PdfPages(os.path.join(FIG_DIR, 'bestfit_lc_{}.pdf'.format(comp_sn)))
    sndata = readfile(os.path.join(LC_DIR, comp_sn))
    nplots = len(sndata['mag'].keys())
    for irow in bestfit_tbdata[-ncomp:][::-1]:
        fig = plt.figure(figsize=[8,8])
        print('Reduced Chisq for {}'.format(irow['sn']))
        for i, band in enumerate(sndata['mag'].keys()):
            ax = fig.add_subplot(np.ceil(nplots/2), 2, i+1)
            comp_day, comp_mag = normalize_data(sndata, band)
            ax.plot(comp_day, comp_mag, 'o', label=comp_sn)
            isndata = readfile(os.path.join(LC_DIR, str(irow['sn'], 'utf-8')))
            if band in isndata['mag'].keys():
                iday, imag = normalize_data(isndata, band)
                ax.plot(iday, imag, 'o', label=str(irow['sn'], 'utf-8'))
                reduced_chisq1 = calc_chisq(comp_day, iday, comp_mag, imag)
                reduced_chisq2 = calc_chisq(iday, comp_day, imag, comp_mag)
                reduced_chisq = (reduced_chisq1 + reduced_chisq2)/2.
                print('\t{}:  {:.2E}'.format(band, reduced_chisq))
            ax.legend(loc='best')
            ax.set_ylim(ax.get_ylim()[::-1])
            ax.set_ylabel('{} Mag'.format(band))
            ax.set_xlim(comp_day[0]-10, comp_day[-1]+20)
        fig.savefig(pp, format='pdf')
        plt.close()
    pp.close()

def check_snname_matches_filename():
    flist = glob.glob(os.path.join(LC_DIR, '*.dat'))
    for ifile in flist:
        if len(os.path.basename(ifile).split('_'))==1:
            sndata = readfile(os.path.splitext(ifile)[0])
            snname = sndata['sn']
            if snname != os.path.splitext(os.path.basename(ifile))[0]:
                print(snname, ifile)

if __name__ == "__main__":
    bandpar = define_filters.define_filters()
    for band in bandpar.keys():
        build_sn_phot_table(LC_DIR, band)
        
    bestfit_dict = build_bestfit_dict('asassn15oz')
    bestfit_tbdata = build_bestfit_table(bestfit_dict)
    plot_bestfit_all_filters(bestfit_tbdata, 'asassn15oz', ncomp=5)
    print(bestfit_tbdata[-10:])
