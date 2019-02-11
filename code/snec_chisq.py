import math
import glob
import os
import sys
import datetime
import pickle

import numpy as np

from astropy.io import ascii as asc
import astropy.units as u
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib as mpl

from utilities_az import visualization, supernova

#plt.style.use(['seaborn-paper', 'az-paper-onecol'])

DARK_DIR = '/dark'
FIG_DIR = '../figures/'
SAVE_DIR = '../snec_models'

snec_models = os.path.join(DARK_DIR,'SNEC/snec_models/')

def prepare_model_data(model_dir, type='mag'):
    '''
    Read in SNEC models, correct for breakout time (model --> observer time) and convert
    from seconds to days
    
    Parameters:
    -----------
    model_dir: str
        name of model directory
    type: str
        type of model (mag [magnitudes] or lum[luminosity])
        
    Returns:
        astropy table with first column time (representing supernova phase in days) and either
        magnitudes in the other columns or luminosity
    '''
    if type=='mag':
        model_tbdata = asc.read(os.path.join(model_dir,'magnitudes.dat'),
                                names=['time', 'temp', 'weird_mag', 
                                       'u', 'g', 'r', 'i', 'z', 'U', 
                                       'B', 'V', 'R', 'I'])
    elif type=='lum':
        model_tbdata = asc.read(os.path.join(model_dir, 'lum_observed.dat'), names=['time', 'luminosity'])
    else:
        print('Type must be "lum" or "mag", {} not understood')
        sys.exit()
    #Observers call t_explosion the time when the explosion is first visible, therefore
    #the models should use this same conventions
    time_breakout = get_breakout_time(model_dir)
    if time_breakout is not None:
        model_tbdata['time'] = ((model_tbdata['time']-time_breakout)*u.second).to(u.day).value
    else:
        model_tbdata=None
    return model_tbdata
    
def get_breakout_time(model_dir):
    '''
    Parses SNEC output to get the breakout time in seconds. If model did not explode, then 
    None is returned
    
    Parameters:
    ------------
    model_dir: str
        name of model directory
        
    Returns:
    ---------
    the breakout time in seconds
    '''
    ofile = open(os.path.join(model_dir, 'info.dat'), 'r')
    all_lines = ofile.readlines()
    if len(all_lines)>6:
        time_breakout = float((all_lines[5].split('=')[1]).strip('seconds\n'))
    else: #SN never got to breakout
        time_breakout = None
    return time_breakout
    
def calc_chisq(data, model, error=0.1):
    '''
    Find the chi-square value given an array of data and model values
    
    Parameters:
    -------------
    data: array
        array of values representing the data
    model: array
        array of values representing the model
    error: flt or array
        if error is an array then it represents the error on each data point
        if error is a float then it will be used to calculate the error on the data (error *data)
    
    Returns:
    ---------
    the chi-square parameter (not reduced)
    '''
    if isinstance(error, float):
        error = error*data
    chisq = np.sum((data - model)**2/error**2)#uncertainty space
    if np.isfinite(chisq) == False:
        print('chisq is not finite')
        import pdb;pdb.set_trace()
    elif chisq < 0:
        print('chisq < 0')
        import pdb; pdb.set_trace()
    elif chisq == 0:
        print('chisq = 0')
        import pdb; pdb.set_trace()
    return chisq
def make_EM_plots(chisq, substr, best_mod):    
    fig = plt.figure()
    fig.subplotpars.update(top=0.9)
    ax = fig.add_subplot(1,1,1)
    ax.set_title('log(ChiSq for Mass vs Energy)')
    im = ax.imshow(np.log(chisq[best_mod.best_ni_mass_indx,:, :, best_mod.best_time_indx, best_mod.best_k_indx, best_mod.best_r_indx ]), interpolation='nearest')
    ax.set_xlabel('Mass (Msun)')
    ax.set_ylabel('Energy (x10^51 ergs)')
    ax.set_xticks(np.arange(len(masses)))
    ax.set_xticklabels(list(masses))
    ax.set_yticks(np.arange(len(energies)))
    ax.set_yticklabels(list(energies))
    fig.colorbar(im)
    plt.savefig(os.path.join(FIG_DIR, 'chisq_{}_MvsE.pdf'.format(substr)))
def make_KR_plots(chisq, substr, best_mod):    
    fig = plt.figure()
    fig.subplotpars.update(top=0.9)
    ax = fig.add_subplot(1,1,1)
    ax.set_title('log(ChiSq for CSM parameters)')
    im = ax.imshow(np.log10(chisq[best_mod.best_ni_mass_indx,best_mod.best_energy_indx, best_mod.best_mass_indx, best_mod.best_time_indx, :, : ]), interpolation='nearest')
    ax.set_xlabel('Radius (100 Rsun)')
    ax.set_ylabel('K (x10^17 g/cm)')
    ax.set_xticks(np.arange(len(radii)))
    ax.set_xticklabels(list(np.int_(radii/100)))
    ax.set_yticks(np.arange(len(Kvalues)))
    ax.set_yticklabels(list(Kvalues))
    #ax.plot(3, 0, 'r*')
    #ax.plot(1, 4, 'w*')
    fig.colorbar(im)
    plt.savefig(os.path.join(FIG_DIR, 'chisq_{}_KvsR.pdf'.format(substr)))

class best_model(object):
    def __init__(self,chisq):
        self.min_indx_base_mod = np.where(chisq == np.nanmin(chisq))
        self.best_ni_mass_indx = self.min_indx_base_mod[0][0]
        self.best_energy_indx =  self.min_indx_base_mod[1][0]
        self.best_mass_indx =    self.min_indx_base_mod[2][0]
        self.best_time_indx =    self.min_indx_base_mod[3][0]
        self.best_k_indx =       self.min_indx_base_mod[4][0]
        self.best_r_indx =       self.min_indx_base_mod[5][0]

        self.best_ni_mass = ni_mass[self.best_ni_mass_indx]
        self.best_ni_mix = 5.0
        self.best_energy = energies[self.best_energy_indx]
        self.best_mass = masses[self.best_mass_indx]
        self.best_time_offset = time_offsets[self.best_time_indx]
        self.best_Kvalue = Kvalues[self.best_k_indx]
        self.best_radius = radii[self.best_r_indx]

def make_identical_fake_data(model_dir, output_filename):
    model_tbdata = prepare_model_data(model_dir, type='lum')
    model_tbdata.write(output_filename, overwrite=True, format='ascii.no_header', delimiter=',')
    
def make_fake_data(model_dir, output_filename):
    model_tbdata = prepare_model_data(model_dir, type='lum')
    model_tbdata['luminosity'] = model_tbdata['luminosity']+np.random.randn(len(model_tbdata['luminosity']))
    model_tbdata.write(output_filename, overwrite=True, format='ascii.no_header', delimiter=',')
                 
def calculate_bolometric_chisq(obs_data_filename, end_plateau, log=True):
    #read in bolometric luminosity
    sn_tbdata = asc.read(obs_data_filename, names=['phase', 'luminosity'])
    if log==True:
        #Convert logL to luminosity
        sn_tbdata['luminosity'] = 10**sn_tbdata['logL']
    #Remove nebular phase from calculations
    sn_phot_tbdata = sn_tbdata[sn_tbdata['phase']<=end_plateau]
    #Create empty array for all chi-squares
    chisq_bolo = np.zeros((len(ni_mass), len(energies), len(masses), len(time_offsets), len(Kvalues), len(radii)))
    #file to record models for which a chi-square isn't calculated    
    with open(os.path.join(SAVE_DIR,'missing_snec_models_bolo.txt'), 'w') as missing_ofile:
        #file to record the parameters and chi-square for each loop 
        with open(os.path.join(SAVE_DIR,'chisq_bolo.csv'), 'w') as chisq_ofile:
            #Loop through parameters
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
                                #Check that a model was calculated for a set of parameters
                                if os.path.exists(os.path.join(model_dir, 'lum_observed.dat')):
                                    model_tbdata = prepare_model_data(model_dir, type='lum')
                                    #Check that model exploded
                                    if model_tbdata is None:
                                        missing_ofile.write('breakout time is None for {}\n'.format(model_dir))
                                        chisq_bolo[i, j, m, :, k, r] = np.nan
                                        #loop over time offsets representing uncertainty in explosion epoch
                                        for toff in time_offsets:
                                            chisq_ofile.write('{},{},{},{},{},{},{}\n'.format(ini_mass, jenergy, mmass, kval, rradius, toff, np.nan))
                                    else:
                                        #loop over time offsets representing uncertainty in explosion epoch
                                        for t, toff in enumerate(time_offsets):
                                            #Interpolate model to observed cadence
                                            model_lum_interp = np.interp(sn_phot_tbdata['phase']+toff, model_tbdata['time'], model_tbdata['luminosity'])                                           #Calculates chisquare, store in array, write to output files
                                            chisq_bolo[i, j, m, t, k, r] = calc_chisq(sn_phot_tbdata['luminosity'],model_lum_interp)
                                            chisq_ofile.write('{},{},{},{},{},{},{}\n'.format(ini_mass, jenergy, mmass, kval, rradius, toff, chisq_bolo[i, j, m, t, k, r]))
                                else:
                                    missing_ofile.write('No file found for {}\n'.format(model_dir))
                                    chisq_bolo[i,j,m,:,k,r] = np.nan
                                    #loop over time offsets representing uncertainty in explosion epoch
                                    for toff in time_offsets:
                                            chisq_ofile.write('{},{},{},{},{},{},{}\n'.format(ini_mass, jenergy, mmass, kval, rradius, toff, np.nan))
    ofile_pickle = open(os.path.join(SAVE_DIR, 'chisq_bolo.pkl'), 'wb')
    pickle.dump(chisq_bolo, ofile_pickle)
    ofile_pickle.close()
    min_indx_base_mod = np.where(chisq_bolo == np.nanmin(chisq_bolo))
    best_bolo_model = best_model(chisq_bolo)

    ofile_write = open(os.path.join(SAVE_DIR, 'best_model_bolo.txt'), 'w')
    ofile_write.write(('best model parameters:\n mass: {}\n energy: \
                        {} \n K: {} \n R: {} \n t_offset: {} \n Ni mass: {}\n'.format(best_bolo_model.best_mass, 
                                                                                      best_bolo_model.best_energy, 
                                                                                      best_bolo_model.best_Kvalue, 
                                                                                      best_bolo_model.best_radius, 
                                                                                      best_bolo_model.best_time_offset, 
                                                                                      best_bolo_model.best_ni_mass)))
    
    R_phot = 967*u.Rsun
    M_csm_opti = 4*math.pi*(best_bolo_model.best_Kvalue*10**17*u.g/u.cm)*((best_bolo_model.best_radius*u.Rsun).to(u.cm) - R_phot.to(u.cm))
    ofile_write.write(('\nM_CSM = {:2.2}'.format(M_csm_opti.to(u.Msun))))
    ofile_write.close()
    
    make_EM_plots(chisq_bolo, 'Lbolo', best_bolo_model)
    make_KR_plots(chisq_bolo, 'Lbolo', best_bolo_model)
    
    
    
    
def calc_color_chisq():
    #Get absolute magnitudes from SNDAVIS database
    sn15oz = supernova.LightCurve2('asassn-15oz')
    sn15oz.get_photometry()
    sn15oz.get_abs_mag()
    #Break into two time segments (to apply different filters to each)
    first_indx = {}
    second_indx = {}
    #Create dictionary of arrays of where photometry is in the first time frame
    for ifilter in filters_1:
        first_indx[ifilter] = sn15oz.phase[ifilter]<end_uv
        #check that time segment has some data
        if (first_indx[ifilter]==False).all():
            print('Warning: no data in first phase for filter {}'.format(ifilter))
    #create a dictionary for array for where photometry is in the second time frame and not nebular
    for ifilter in filters_2:
        second_indx[ifilter] = (sn15oz.phase[ifilter]>=end_uv) & (sn15oz.phase[ifilter]<=end_plateau)
        #Check that time segment has some data
        if (second_indx[ifilter]==False).all():
            print('Warning: no data in second phase for filter {}'.format(ifilter))    
    #build array of observed magnitude values for first time segment for all filters in first time segment
    first_obs_array =  np.hstack([sn15oz.abs_mag[ifilter][first_indx[ifilter]] for ifilter in filters_1])
    #build array of observed magnitude values for second time segment for all filters in second time segment
    second_obs_array = np.hstack([sn15oz.abs_mag[ifilter][second_indx[ifilter]] for ifilter in filters_2])
    #Define normalization arrays (for reduced chisq)
    first_err_array = np.hstack([sn15oz.abs_mag_err[ifilter][first_indx[ifilter]] for ifilter in filters_1])
    second_err_array = np.hstack([sn15oz.abs_mag_err[ifilter][second_indx[ifilter]] for ifilter in filters_2])
    norm1=None
    norm2=None
    norm=None
    norm = len(first_obs_array.flatten())+len(second_obs_array.flatten()) #Total number of points in chisquare
    if norm1:
        del norm1
        del norm2
    #norm1 = len(filters_1)
    #norm2 = len(filters2)
    #if norm:
        #del norm
    #Create chisq array to populate
    chisq_color = np.zeros((len(ni_mass), len(energies), len(masses), len(time_offsets), len(Kvalues), len(radii)))
    ##Create file for missing models
    with open(os.path.join(SAVE_DIR,'missing_snec_models_color.txt'), 'w') as missing_ofile:
        #Create a file for models and chisq values
        with open(os.path.join(SAVE_DIR,'chisq_color.csv'), 'w') as chisq_ofile:
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
                                #Check that model exists
                                if os.path.exists(os.path.join(model_dir, 'magnitudes.dat')):
                                    model_tbdata = prepare_model_data(model_dir)
                                    #Check that model finished running
                                    if model_tbdata is None:
                                        missing_ofile.write('breakout time is None for {}\n'.format(model_dir))
                                        chisq_color[i,j,m,:,k,r] = np.nan
                                    #Calculate chi-square
                                    else:
                                        for t, toff in enumerate(time_offsets):
                                            #Interpolate model values to observed values
                                            first_lc = np.hstack(np.array([np.interp(sn15oz.phase[ifilter][first_indx[ifilter]]+toff, model_tbdata['time'], model_tbdata[ifilter]) for ifilter in filters_1]))
                                            second_lc = np.hstack(np.array([np.interp(sn15oz.phase[ifilter][second_indx[ifilter]]+toff, model_tbdata['time'], model_tbdata[ifilter]) for ifilter in filters_2]))
                                            #Calcuate chisquare for the first and second part of the light curve
                                            first_chisq = calc_chisq(first_obs_array, first_lc, first_err_array)
                                            second_chisq = calc_chisq(second_obs_array, second_lc, second_err_array)
                                            #Calculate a reduced chisquare
                                            if norm1:
                                                ichi = np.abs(first_chisq/norm1 + second_chisq/norm2)
                                            elif norm:
                                                ichi = np.abs((first_chisq + second_chisq)/norm) #magnitudes are negative, so this is negative w/o abs()
                                            else:
                                                print('neither norm nor norm1 were defined')
                                                break
                                            chisq_color[i, j, m, t, k, r] = ichi
                                            chisq_ofile.write('{},{},{},{},{},{},{}\n'.format(ini_mass, jenergy, mmass, kval, rradius, toff, chisq_color[i, j, m, t, k, r]))
                                else:
                                    missing_ofile.write('No file found for {}\n'.format(model_dir))
                                    chisq_color[i,j,m,:,k,r] = np.nan
                                    for toff in time_offsets:
                                            chisq_ofile.write('{},{},{},{},{},{},{}\n'.format(ini_mass, jenergy, mmass, kval, rradius, toff, np.nan))
    #Save chisquare
    ofile_pickle = open(os.path.join(SAVE_DIR, 'chisq_color.pkl'), 'wb')
    pickle.dump(chisq_color, ofile_pickle)
    ofile_pickle.close()

    best_color_model = best_model(chisq_color) 
    
    ofile_write = open(os.path.join(SAVE_DIR, 'best_model_color.txt'), 'w')
    ofile_write.write(('best model parameters:\n mass: {}\n energy: \
                        {} \n K: {} \n R: {} \n t_offset: {} \n Ni mass: {}'.format(best_color_model.best_mass, 
                                                                                    best_color_model.best_energy, 
                                                                                    best_color_model.best_Kvalue, 
                                                                                    best_color_model.best_radius, 
                                                                                    best_color_model.best_time_offset, 
                                                                                    best_color_model.best_ni_mass)))
    R_phot = 967*u.Rsun
    M_csm_opti = 4*math.pi*(best_color_model.best_Kvalue*10**17*u.g/u.cm)*((best_color_model.best_radius*u.Rsun).to(u.cm) - R_phot.to(u.cm))
    ofile_write.write(('\nM_CSM = {:2.2}'.format(M_csm_opti.to(u.Msun))))
    ofile_write.close()
    make_EM_plots(chisq_color, 'color', best_color_model)
    make_KR_plots(chisq_color, 'color', best_color_model)

if __name__ == "__main__":
    #filters_1 = ['U', 'B', 'V', 'R', 'I']
    #filters_2 = ['V', 'R', 'I']
    filters_1 = [ 'U', 'V', 'g', 'r', 'i']
    filters_2 = ['g', 'r', 'i']
    end_uv = 20
    end_plateau = 150
    #Testing section
    #ni_mass =   np.array([0.083])
    #energies =  np.array([0.5, 0.8])
    #masses =    np.array([11, 13])
    #ni_mixing = np.array([5.0])
    #time_offsets = np.arange(-4, 4, 1)
    #Kvalues =   np.array([0, 10, 20])
    #radii =     np.array([0, 1500, 1800])
    #make_fake_data(os.path.join(snec_models, 
    #                                         'Ni_mass_{:1.4f}'.format(0.083),
    #                                         'Ni_mixing_{:1.1f}'.format(5.0),
    #                                         'M{:2.1f}'.format(11),
    #                                         'E_{:1.3f}'.format(0.5),
    #                                         'K_{:2.1f}'.format(10), 
    #                                         'R_{}'.format(1500),
    #                                         'Data'),
    #                        'identical_fake_lum.csv')
    #obs_data_filename = 'identical_fake_lum.csv'
    ni_mass =   np.array([0.083, 0.0965, 0.11])
    energies =  np.array([0.5, 0.8, 1.1, 1.4, 1.7, 2.0])
    masses =    np.array([11, 13, 15, 16, 17, 18, 21])
    ni_mixing = np.array([5.0])
    time_offsets = np.arange(-4, 4, 1)
    Kvalues =   np.array([0, 10, 20, 30, 35, 40, 50, 60])
    radii =     np.array([0, 1500, 1800, 2100, 2400, 2700, 3000, 3300])
    start = datetime.datetime.now()
    calc_color_chisq()
    end = datetime.datetime.now()
    print('color calc took {}'.format(end-start))
    start = datetime.datetime.now()
    obs_data_filename = '../data/bolometric.txt'
    calculate_bolometric_chisq(obs_data_filename, end_plateau, log=False)
    end = datetime.datetime.now()

    print('Bolo luminosity calc took {}'.format(((end-start)*u.second).to(u.minute)))
    
