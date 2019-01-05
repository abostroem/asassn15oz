import math
import glob
import os
import sys

import numpy as np

from astropy.io import ascii as asc
import astropy.units as u
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib as mpl

from utilities_az import visualization, supernova

#plt.style.use(['seaborn-paper', 'az-paper-onecol'])

DARK_DIR = '/Users/bostroem/dark'
FIG_DIR = '../figures/'
SAVE_DIR = '../snec_models'

snec_models = os.path.join(DARK_DIR,'SNEC/snec_models/')

def calc_abs_mag(app_mag, dist_mod, A_mw, A_host=0, dist_mod_err=0, app_mag_err=0, A_err_mw=0, A_err_host=0):
    '''
    Calculate the absolute magnitude in a given filter given the apparent magnitude, distance modulus, and extinction
    app_mag: arr or float
        array of apparent magnitudes
    dist_mod: float
        distance modulus
    A_mw: float
        Magnitudes of extinction in filter due to the Milky Way
    A_host: float
        Magnitudes of extinction in filter due to host; default=0
    dist_mod_err: float
        error in the distance modulus; default=0
    app_mag_err: arr or float
        array of error values for apparent magnitude (should be same size as app_mag); default=0
    A_err_mw: float
        error on A_mw (in magnitudes); default=0
    A_err_host: float
        err on A_host (in magnitudes); default=0
    '''
    abs_mag = app_mag - dist_mod - A_mw - A_host
    abs_mag_err = np.sqrt(app_mag_err**2 + dist_mod_err**2 + A_err_mw**2 + A_err_host**2)
    return abs_mag, abs_mag_err

class LightCurve2(object):
    '''
    Retrieve information from sndavis database
    
    Absolute magnitudes are corrected for extinction, apparent magnitude is not
    '''
    def __init__(self, name):
        self.name=name
        self.db, self.cursor = connect_to_sndavis.get_cursor()
        self.id = self.get_sn_id()
        if self.id is not None:
            self.cursor.execute('SELECT * FROM idsupernovae WHERE id = {}'.format(self.id))
            result = self.cursor.fetchone()
            self.ra = result['ra0']
            self.dec = result['dec0']
            self.type = result['sntype'] #TODO - update this to look at the type table and make a human readable type
            self.dist_mod = result['mu']
            self.dist_mod_err = result['muerr']
            self.ebv_mw = result['ebvg']
            self.ebv_mw_err = result['ebvgerr']
            self.ebv_host = result['ebvi'] 
            self.ebv_host_err = result['ebvierr']
            self.quality = result['quality']
            self.jdexpl = result['jdexpl']
            self.jdexpl_err = result['jdexplerr']
            self.jd = {}
            self.apparent_mag = {}
            self.apparent_mag_err = {}
            self.phase = {}
            self.abs_mag = {}
            self.abs_mag_err = {}
            self.slopes = {'s1':{}, 's2':{}, 's50':{}, 'tail':{}, 
                            's1_range':{}, 's2_range':{}, 's50_range':{}, 'tail_range':{},
                            's1_err':{}, 's2_err':{}, 's50_err':{}, 'tail_err':{}}
            self.A_host = {}
            self.A_err_host = {}
            self.A_mw = {}
            self.A_err_mw = {}
        else:
            print('SN {} not found in database'.format(self.name))
            #TODO add an absolute magnitude error
             
    def get_sn_id(self):
        self.cursor.execute('SELECT targetid FROM supernovanames WHERE name = {}'.format('"{}"'.format(self.name)))
        result = self.cursor.fetchone()
        if result is None:
            return None
        return result['targetid']
 
    def get_photometry(self, band='all'):
        '''
        Get apparent magnitude of a specific band
        '''
        if band == 'all':
            self.cursor.execute("SELECT DISTINCT BINARY(filter) FROM photometry WHERE targetid={}".format(self.id))
            results = self.cursor.fetchall()
            bands = [iband['BINARY(filter)'].decode('utf-8') for iband in results]

        elif len(band)==1: #Single filter
            bands = [band]
        else: #multiple filters
            bands = band
        for ifilter in bands:   
            A_host_band, A_err_host_band = spec.calc_extinction(self.ebv_host, ifilter)
            A_mw_band, A_err_mw_band = spec.calc_extinction(self.ebv_mw, ifilter)
            self.cursor.execute("SELECT jd, mag, magerr FROM photometry WHERE targetid={} AND filter=BINARY('{}')".format(self.id, ifilter))
            results = self.cursor.fetchall()
            jd = []
            mag = []
            mag_err = []
            for irow in results:
                jd.append(float(irow['jd']))
                mag.append(float(irow['mag']))
                mag_err.append(float(irow['magerr']))
            self.jd[ifilter] = np.array(jd)
            self.apparent_mag[ifilter] = np.array(mag)
            self.apparent_mag_err[ifilter] = np.array(mag_err)
            self.phase[ifilter] = self.jd[ifilter] - self.jdexpl
            self.A_host[ifilter] = A_host_band
            self.A_err_host[ifilter] = A_err_host_band
            self.A_mw[ifilter] = A_mw_band
            self.A_err_mw[ifilter] = A_err_mw_band
            
    def get_abs_mag(self, band='all'):
        '''
        Calcualte absolute magnitude
        '''
        print('Calculating Absolute Magntidue with Extinction')
        if band == 'all':
            bands = self.apparent_mag.keys()
        elif len(band)==1:
            bands = [band]
        else:
            bands = band
        if self.dist_mod_err is None:
            self.dist_mod_err = 0
        for iband in bands:
            if self.ebv_host is None:
                self.abs_mag[iband], self.abs_mag_err[iband] = calc_abs_mag(self.apparent_mag[iband], 
                                                   self.dist_mod, 
                                                   self.A_mw[iband], 
                                                   app_mag_err=self.apparent_mag_err[iband], 
                                                   dist_mod_err=self.dist_mod_err)
            else:
                self.abs_mag[iband], self.abs_mag_err[iband] = calc_abs_mag(self.apparent_mag[iband], 
                                                                            self.dist_mod, 
                                                                            self.A_mw[iband], 
                                                                            A_host=self.A_host[iband],
                                                                            app_mag_err=self.apparent_mag_err[iband], 
                                                                            dist_mod_err=self.dist_mod_err,
                                                                            A_err_mw=self.A_err_mw[iband],
                                                                            A_err_host=self.A_err_host[iband])
    
    def get_slope(self,slope_type, band='V'):
        '''
        '''
        self.cursor.execute('SELECT * FROM snslope WHERE  targetid={} AND \
                                                         slopetype="{}" AND \
                                                         filter=BINARY("{}")'.format(self.id, slope_type, band))
        result = self.cursor.fetchone()
        if result is not None:
            self.slopes['{}'.format(slope_type)][band]       = result['slope']
            self.slopes['{}_err'.format(slope_type)][band]   = result['slopeerr']
            self.slopes['{}_range'.format(slope_type)][band] = (result['tstart'], result['tstop'])
        else:
            self.slopes['{}'.format(slope_type)][band]       = None
            self.slopes['{}_err'.format(slope_type)][band]   = None
            self.slopes['{}_range'.format(slope_type)][band] = [None, None]





def prepare_model_data(model_dir):
    model_mag_tbdata = asc.read(os.path.join(model_dir,'magnitudes.dat'),
                                names=['time', 'temp', 'weird_mag', 
                                       'u', 'g', 'r', 'i', 'z', 'U', 
                                       'B', 'V', 'R', 'I'])
    #Observers call t_explosion the time when the explosion is first visible, therefore
    #the models should use this same conventions
    time_breakout = get_breakout_time(model_dir)
    if time_breakout is not None:
        model_mag_tbdata['time'] = ((model_mag_tbdata['time']-time_breakout)*u.second).to(u.day).value
    else:
        model_mag_tbdata=None
    return model_mag_tbdata
    
def get_breakout_time(model_dir):
    ofile = open(os.path.join(model_dir, 'info.dat'), 'r')
    all_lines = ofile.readlines()
    if len(all_lines)>6:
        time_breakout = float((all_lines[5].split('=')[1]).strip('seconds\n'))
    else: #SN never got to breakout
        time_breakout = None
    return time_breakout
                 
def calculate_bolometric_chisq():
    with open('missing_snec_models_bolo.txt', 'w') as missing_ofile:
        with open('chisq_bolo.csv', 'w') as chisq_ofile:
            for i, ini_mass in enumerate(ni_mass):
                for j, jenergy in enumerate(energies):
                    for m, mmass in enumerate(masses):
                        for k, kval in enumerate(Kvalues):
                            for r, rradius in enumerate(radii):
                                model_dir = os.path.join(snec_models, 
                                             'Ni_mass_{:1.4f}'.format(ini_mass),
                                             'Ni_mixing_{:1.1f}'.format(5.0),
                                             'M{:2.1f}'.format(mmass),
                                             'E_{:1.3f}'.format(jenergy),
                                             'K_{:2.1f}'.format(kval), 
                                             'R_{}'.format(rradius),
                                             'Data')
                                if os.path.exists(os.path.join(model_dir, 'lum_observed.dat')):
                                    model_tbdata = asc.read(os.path.join(model_dir, 'lum_observed.dat'), names=['phase', 'luminosity'])
                                    breakout_time = get_breakout_time(model_dir)
                                    if breakout_time is None:
                                        missing_ofile.write('breakout time is None for {}\n'.format(model_dir))
                                        chisq_bolo[i, j, m, :, k, r] = np.nan
                                        for toff in time_offsets:
                                            chisq_ofile.write('{},{},{},{},{},{},{}'.format(ini_mass, jenergy, mmass, kval, rradius, toff, np.nan))
                                    else:
                                        model_tbdata['phase_obs'] = ((model_tbdata['phase']-breakout_time)*u.s).to(u.day).value
                                        for t, toff in enumerate(time_offsets):
                                            model_lum_interp = np.interp(sn15oz_phot_tbdata['phase'], model_tbdata['phase_obs']+toff, model_tbdata['luminosity'])
                                            chisq_bolo[i, j, m, t, k, r] = np.sum((sn15oz_phot_tbdata['luminosity'] - model_lum_interp)**2/model_lum_interp)
                                            chisq_ofile.write('{},{},{},{},{},{},{}'.format(ini_mass, jenergy, mmass, kval, rradius, toff, chisq_bolo[i, j, m, t, k, r]))
                                else:
                                    missing_ofile.write('No file found for {}\n'.format(model_dir))
                                    chisq_bolo[i,j,m,:,k,r] = np.nan
                                    for toff in time_offsets:
                                            chisq_ofile.write('{},{},{},{},{},{},{}'.format(ini_mass, jenergy, mmass, kval, rradius, toff, np.nan))
    ofile_pickle = open(os.path.join(SAVE_DIR, 'chisq_bolo.pkl'), 'w')
    ofile_pickle.dump(chisq_bolo)
    ofile_pickel.close()
    min_indx_base_mod = np.where(chisq_bolo == np.nanmin(chisq_bolo))
    best_ni_mass_indx = min_indx_base_mod[0][0]
    best_energy_indx = min_indx_base_mod[1][0]
    best_mass_indx = min_indx_base_mod[2][0]
    best_time_indx = min_indx_base_mod[3][0]
    best_k_indx = min_indx_base_mod[4][0]
    best_r_indx = min_indx_base_mod[5][0]

    best_ni_mass = ni_mass[best_ni_mass_indx]
    best_ni_mix = 5.0
    best_energy = energies[best_energy_indx]
    best_mass = masses[best_mass_indx]
    best_time_offset = time_offsets[best_time_indx]
    best_Kvalue = Kvalues[best_k_indx]
    best_radius = radii[best_r_indx]

    ofile_write = open(os.path.join(SAVE_DIR, 'best_model_bolo.txt'), 'w')
    ofile_write.write(('best model parameters:\n mass: {}\n energy: {} \n K: {} \n R: {} \n t_offset: {} \n Ni mass: {}'.format(best_mass, best_energy, best_Kvalue, best_radius, best_time_offset, best_ni_mass)))
    
    R_phot = 967*u.Rsun
    M_csm_opti = 4*math.pi*(best_Kvalue*10**17*u.g/u.cm)*((best_radius*u.Rsun).to(u.cm) - R_phot.to(u.cm))
    ofile_write.write((M_csm_opti.to(u.Msun)))
    ofile_write.write(('M_CSM = {:2.2}'.format(M_csm_opti.to(u.Msun))))
    ofile_write.close()
    
    fig = plt.figure()
    fig.subplotpars.update(top=0.9)
    ax = fig.add_subplot(1,1,1)
    ax.set_title('log(ChiSq for Mass vs Energy)')
    im = ax.imshow(np.log(chisq_bolo[best_ni_mass_indx,:, :, best_time_indx, best_k_indx, best_r_indx ]), interpolation='nearest')
    ax.set_xlabel('Mass (Msun)')
    ax.set_ylabel('Energy (x10^51 ergs)')
    ax.set_xticks(np.arange(len(masses)))
    ax.set_xticklabels(list(masses))
    ax.set_yticks(np.arange(len(energies)))
    ax.set_yticklabels(list(energies))
    #ax.plot(5, 3, 'r*')
    #ax.plot(4, 2, 'w*')
    fig.colorbar(im)
    plt.savefig(os.path.join(FIG_DIR, 'chisq_Lbolo_MvsE.pdf'))
    
    fig = plt.figure()
    fig.subplotpars.update(top=0.9)
    ax = fig.add_subplot(1,1,1)
    ax.set_title('log(ChiSq for CSM parameters)')
    im = ax.imshow(np.log10(chisq_bolo[best_ni_mass_indx,best_energy_indx, best_mass_indx, best_time_indx, :, : ]), interpolation='nearest')
    ax.set_xlabel('Radius (100 Rsun)')
    ax.set_ylabel('K (x10^17 g/cm)')
    ax.set_xticks(np.arange(0,7))
    ax.set_xticklabels(np.int_(radii/100))
    ax.set_yticks(np.arange(0, 7))
    ax.set_yticklabels(Kvalues)
    #ax.plot(3, 0, 'r*')
    #ax.plot(1, 4, 'w*')
    fig.colorbar(im)
    plt.savefig(os.path.join(FIG_DIR, 'chisq_Lbolo_KvsR.pdf'))
    
    
    
def calc_color_chisq():
    sn15oz = supernova.LightCurve2('asassn-15oz')
    sn15oz.get_photometry()
    sn15oz.get_abs_mag()

    first_indx = {}
    second_indx = {}
    for ifilter in filters_1:
        first_indx[ifilter] = sn15oz.phase[ifilter]<end_uv
        if (first_indx[ifilter]==False).all():
            print('Warning: no data in first phase for filter {}'.format(ifilter))
    for ifilter in filters_2:
        second_indx[ifilter] = (sn15oz.phase[ifilter]>=end_uv) & (sn15oz.phase[ifilter]<=end_plateau)
        if (second_indx[ifilter]==False).all():
            print('Warning: no data in second phase for filter {}'.format(ifilter))    
    first_obs_array =  np.hstack([sn15oz.abs_mag[ifilter][first_indx[ifilter]] for ifilter in filters_1])
    second_obs_array = np.hstack([sn15oz.abs_mag[ifilter][second_indx[ifilter]] for ifilter in filters_2])
    norm1=None
    norm2=None
    norm=None
    norm = len(first_obs_array.flatten())+len(second_obs_array.flatten())
    if norm1:
        del norm1
        del norm2
    #norm1 = len(filters_1)
    #norm2 = len(filters2)
    #if norm:
        #del norm
    chisq_color = np.zeros((len(ni_mass), len(energies), len(masses), len(time_offsets), len(Kvalues), len(radii)))
    with open('missing_snec_models_color.txt', 'w') as missing_ofile:
        with open('chisq_color.csv', 'w') as chisq_ofile:
            for i, ini_mass in enumerate(ni_mass):
                for j, jenergy in enumerate(energies):
                    for m, mmass in enumerate(masses):
                        for k, kval in enumerate(Kvalues):
                            for r, rradius in enumerate(radii):
                                model_dir = os.path.join(snec_models, 
                                             'Ni_mass_{:1.4f}'.format(ini_mass),
                                             'Ni_mixing_{:1.1f}'.format(5.0),
                                             'M{:2.1f}'.format(mmass),
                                             'E_{:1.3f}'.format(jenergy),
                                             'K_{:2.1f}'.format(kval), 
                                             'R_{}'.format(rradius),
                                             'Data')
                                if os.path.exists(os.path.join(model_dir, 'magnitudes.dat')):
                                    model_tbdata = prepare_model_data(model_dir)
                                    if model_tbdata is None:
                                        missing_ofile.write('breakout time is None for {}\n'.format(model_dir))
                                        chisq_color[i,j,m,:,k,r] = np.nan
                                    else:
                                        first_lc = np.hstack([np.interp(sn15oz.phase[ifilter][first_indx[ifilter]], model_tbdata['time'], model_tbdata[ifilter]) for ifilter in filters_1])
                                        second_lc = np.hstack([np.interp(sn15oz.phase[ifilter][second_indx[ifilter]], model_tbdata['time'], model_tbdata[ifilter]) for ifilter in filters_2])
                                        for t, toff in enumerate(time_offsets):
                                            first_chisq = np.sum((first_obs_array - first_lc)**2/first_lc)
                                            second_chisq = np.sum((second_obs_array - second_lc)**2/second_lc)
                                            if norm1:
                                            
                                                ichi = first_chisq/norm1 + second_chisq/norm2
                                            elif norm:
                                                ichi = (first_chisq + second_chisq)/norm
                                            else:
                                                print('neither norm nor norm1 were defined')
                                                break
                                            chisq_color[i, j, m, t, k, r] = ichi
                                            chisq_ofile.write('{},{},{},{},{},{},{}'.format(ini_mass, jenergy, mmass, kval, rradius, toff, chisq_color[i, j, m, t, k, r]))
                                else:
                                    missing_ofile.write('No file found for {}\n'.format(model_dir))
                                    chisq_color[i,j,m,:,k,r] = np.nan
                                    for toff in time_offsets:
                                            chisq_ofile.write('{},{},{},{},{},{},{}'.format(ini_mass, jenergy, mmass, kval, rradius, toff, np.nan))
    ofile_pickle = open(os.path.join(PKL_DIR, 'chisq_color.pkl'))
    ofile_pickle.dump(chisq_color)
    ofile_pickel.close()
    min_indx_base_mod = np.where(chisq_color == np.nanmin(chisq_color))
    best_ni_mass_indx = min_indx_base_mod[0][0]
    best_energy_indx = min_indx_base_mod[1][0]
    best_mass_indx = min_indx_base_mod[2][0]
    best_time_indx = min_indx_base_mod[3][0]
    best_k_indx = min_indx_base_mod[4][0]
    best_r_indx = min_indx_base_mod[5][0]

    best_ni_mass = ni_mass[best_ni_mass_indx]
    best_ni_mix = 5.0
    best_energy = energies[best_energy_indx]
    best_mass = masses[best_mass_indx]
    best_time_offset = time_offsets[best_time_indx]
    best_Kvalue = Kvalues[best_k_indx]
    best_radius = radii[best_r_indx]
    
    ofile_write = open(os.path.join(SAVE_DIR, 'best_model_color.txt'), 'w')
    ofile_write.write(('best model parameters:\n mass: {}\n energy: {} \n K: {} \n R: {} \n t_offset: {} \n Ni mass: {}'.format(best_mass, best_energy, best_Kvalue, best_radius, best_time_offset, best_ni_mass)))
    ofile_write.write(R_phot = 967*u.Rsun)
    ofile_write.write(M_csm_opti = 4*math.pi*(best_Kvalue*10**17*u.g/u.cm)*((best_radius*u.Rsun).to(u.cm) - R_phot.to(u.cm)))
    ofile_write.close()

    fig = plt.figure()
    fig.subplotpars.update(top=0.9)
    ax = fig.add_subplot(1,1,1)
    ax.set_title('log(ChiSq for Mass vs Energy)')
    im = ax.imshow(chisq_color[best_ni_mass_indx,:, :, best_time_indx, best_k_indx, best_r_indx ], interpolation='nearest')
    ax.set_xlabel('Mass (Msun)')
    ax.set_ylabel('Energy (x10^51 ergs)')
    ax.set_xticks(np.arange(len(masses)))
    ax.set_xticklabels(list(masses))
    ax.set_yticks(np.arange(len(energies)))
    ax.set_yticklabels(list(energies))
    #ax.plot(5, 3, 'r*')
    #ax.plot(4, 2, 'w*')
    fig.colorbar(im)
    plt.savefig(os.path.join(FIG_DIR, 'chisq_color_MvE.pdf'))
    
    fig = plt.figure()
    fig.subplotpars.update(top=0.9)
    ax = fig.add_subplot(1,1,1)
    ax.set_title('log(ChiSq for CSM parameters)')
    im = ax.imshow(np.log10(chisq_color[best_ni_mass_indx,best_energy_indx, best_mass_indx, best_time_indx, :, : ]), interpolation='nearest')
    ax.set_xlabel('Radius (100 Rsun)')
    ax.set_ylabel('K (x10^17 g/cm)')
    ax.set_xticks(np.arange(0,7))
    ax.set_xticklabels(np.int_(radii/100))
    ax.set_yticks(np.arange(0, 7))
    ax.set_yticklabels(Kvalues)
    #ax.plot(3, 0, 'r*')
    #ax.plot(1, 4, 'w*')
    fig.colorbar(im)
    plt.savefig(os.path.join(FIG_DIR, 'chisq_color_KvsR.pdf'))

if __name__ == "__main__":
    filters_1 = ['U', 'B', 'V', 'R', 'I']
    filters_2 = ['V', 'R', 'I']
    end_uv = 30
    end_plateau = 150

    ni_mass =   np.array([0.083, 0.0965, 0.11])
    energies =  np.array([0.5, 0.8, 1.1, 1.4, 1.7, 2.0])
    masses =    np.array([11, 13, 15, 16, 17, 18, 21])
    ni_mixing = np.array([5.0])
    time_offsets = np.arange(-4, 4, 1)
    Kvalues =   np.array([0, 10, 20, 30, 35, 40, 50, 60])
    radii =     np.array([0, 1500, 1800, 2100, 2400, 2700, 3000, 3300])
    calc_color_chisq()
    calculate_bolometric_chisq()
    