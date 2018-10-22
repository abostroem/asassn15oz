 
# coding: utf-8

# Creates: 
#     * example_fits.pdf

 
# In[1]:


import os
import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import yaml
import line_analysis_BSNIP
import fit_spectral_lines


 
 
# In[2]:


plt.style.use(['seaborn-paper','az-paper-twocol'])


 
 
# In[3]:


FIG_DIR = '.'
OUTPUT_DIR = '../../data/line_info'
input_file_Sc = '../../data/line_vel_analysis_input.yml'
input_file_Ca = os.path.join(OUTPUT_DIR, 'CaII_input.yaml')


 
 
# In[4]:


filename_Ca = 'asassn-15oz_20150924_redblu_123847.580.fits'
dir_Ca = '../../data/spectra/lco'

filename_Sc = 'asassn15oz_20151006_redblu_101906.800.fits'
dir_Sc = dir_Ca


 
 
# In[5]:


spectrum_Ca = line_analysis_BSNIP.read_iraf_spectrum(os.path.join(dir_Ca, filename_Ca))
spectrum_Ca.__setattr__('filename', filename_Ca)
interactive_fig, ax1, ax2 = fit_spectral_lines.plot_spectrum(spectrum_Ca)  
ax1, x_cont, y_cont, flux_norm = fit_spectral_lines.continuum_normalization(ax1, spectrum_Ca, 
                                                                            input_filename=input_file_Ca, 
                                                                            interactive=False, absorption=True)
continuum_l = fit_spectral_lines.build_continuum_object(spectrum_Ca, x_cont[0], y_cont[0])
continuum_r = fit_spectral_lines.build_continuum_object(spectrum_Ca, x_cont[1], y_cont[1])
fit_wave_Ca = np.arange(x_cont.min(), x_cont.max()+1, 0.01)
#-----------------------
#Fit the feature
#-----------------------
(x1, y1), (x2, y2) = fit_spectral_lines.get_id_fit_region_from_file(spectrum_Ca.filename, input_file_Ca)
#Select the normalized fit spectrum
line_indx_Ca = (spectrum_Ca.wave<max(x1, x2)) & (spectrum_Ca.wave > min(x1, x2))
line_wave_Ca = spectrum_Ca.wave[line_indx_Ca]
line_flux_Ca = flux_norm[line_indx_Ca]
#Select the line centers
center_list = fit_spectral_lines.get_line_centers_from_file(spectrum_Ca.filename, input_file_Ca)
#Select the fitting function
fit_type = fit_spectral_lines.get_fit_type_from_file(spectrum_Ca.filename, input_file_Ca)
fit_Ca, lines, args, ax1, ax2 = fit_spectral_lines.fit_feature(line_wave_Ca, line_flux_Ca, fit_wave_Ca, fit_type, center_list, ax1, ax2, 
                                                               continuum_l, continuum_r,
                                                            offsets=[44, 164], fixed_offset=True, similar_widths=True)
fit_edge_l = fit_spectral_lines.build_continuum_object(spectrum_Ca, x1, y1)
fit_edge_r = fit_spectral_lines.build_continuum_object(spectrum_Ca, x2, y2)
plt.close()


 
 
# In[6]:


iline = 'ScII5662'
ofile = open(input_file_Sc, 'r')
feature_dict = yaml.load(ofile)[iline]


#Read in spectrum
spectrum_Sc = line_analysis_BSNIP.read_iraf_spectrum(os.path.join(dir_Sc, filename_Sc))
#Remove CR and other large deviations
smooth_flux = line_analysis_BSNIP.smooth_signal(spectrum_Sc.flux, 
                            feature_dict['smooth_param']['width'], 
                            feature_dict['smooth_param']['deg'])
edge_results = line_analysis_BSNIP.find_edges(spectrum_Sc, feature_dict, smooth_flux, os.path.join(dir_Sc, filename_Sc))
if edge_results is not None:
    err_binsize, blue_edge_indx, red_edge_indx, wcenter_indx, continuum_l_Sc, continuum_r_Sc = edge_results
    #Calculate the most common velocity
    wcenter = spectrum_Sc.wave[wcenter_indx]
    vel_fit = line_analysis_BSNIP.find_velocity(spectrum_Sc.wave, smooth_flux, spectrum_Sc.error, wcenter, continuum_l_Sc, continuum_r_Sc, err_binsize)
    #Find the error in the pseudo equivalent width
    line_indx_Sc = np.arange(len(spectrum_Sc.wave))[(spectrum_Sc.wave>=continuum_l_Sc.wave)&(spectrum_Sc.wave<=continuum_r_Sc.wave)]
    min_indx_Sc = int(np.floor(line_indx_Sc[0]-err_binsize/2))
    max_indx_Sc = int(np.ceil(line_indx_Sc[-1]+err_binsize/2+1))
    continuum_extended = line_analysis_BSNIP.calc_continuum(spectrum_Sc.wave[min_indx_Sc:max_indx_Sc], continuum_l_Sc, continuum_r_Sc)
    flux_var = line_analysis_BSNIP.calc_flux_variance(spectrum_Sc.flux[min_indx_Sc:max_indx_Sc]-continuum_extended,
                                vel_fit(spectrum_Sc.wave[min_indx_Sc:max_indx_Sc]), err_binsize) #These don't include the errors from the continuum subtraction yet; ok for EW calc
    if len(flux_var) > len(spectrum_Sc.flux[line_indx_Sc]):
        flux_var = flux_var[1:-1]
    continuum_Sc = line_analysis_BSNIP.calc_continuum(spectrum_Sc.wave[line_indx_Sc], continuum_l_Sc, continuum_r_Sc)
    continuum_var = line_analysis_BSNIP.calc_continuum_variance(spectrum_Sc.wave[line_indx_Sc], continuum_l_Sc, continuum_r_Sc)
    delta_wave = np.median(spectrum_Sc.wave[1:]-spectrum_Sc.wave[:-1])
    #Find the velocity error
    vel_err = line_analysis_BSNIP.calc_velocity_error(spectrum_Sc.wave[line_indx_Sc], spectrum_Sc.flux[line_indx_Sc], vel_fit, continuum=continuum_Sc)
    vel_min_Sc = spectrum_Sc.wave[line_indx_Sc][np.argmin(vel_fit(spectrum_Sc.wave[line_indx_Sc]))]


 
 
# In[7]:


fig = plt.figure()

left_Ca, width_Ca = 0.1, 0.4
bottom_Ca, height_Ca = 0.17, 0.8

left_Sc, width_Sc = 0.59, 0.4
bottom_Sc, height_Sc = 0.17, 0.8

ax_Ca = plt.axes([left_Ca, bottom_Ca, width_Ca, height_Ca])
ax_Sc = plt.axes([left_Sc, bottom_Sc, width_Sc, height_Sc])

min_indx = int(np.where(spectrum_Ca.wave == spectrum_Ca.wave[line_indx_Ca][0])[0])
max_indx = int(np.where(spectrum_Ca.wave == spectrum_Ca.wave[line_indx_Ca][-1])[0])
wave_Ca = spectrum_Ca.wave[min_indx-100:max_indx+100]

l0 = ax_Ca.axvspan(fit_edge_l.wave, fit_edge_r.wave, 
              ymin=ax_Ca.get_ylim()[0]*10**15, ymax = ax_Ca.get_ylim()[1]*10**15, 
              alpha=0.1, color='k', label='Fit region')
continuum_Ca = np.polyval(np.polyfit([continuum_l.wave, continuum_r.wave], [continuum_l.flux, continuum_r.flux], 1), wave_Ca)
l1=ax_Ca.plot(wave_Ca, spectrum_Ca.flux[min_indx-100:max_indx+100]*10**15, label='Spectrum')

for ifit in fit_Ca[1:]:
    l2=ax_Ca.plot(wave_Ca, (fit_Ca[0](wave_Ca) - ifit(wave_Ca))*continuum_Ca*10**15, ':', label='Individual fit')
l3=ax_Ca.plot(wave_Ca, fit_Ca(wave_Ca)*continuum_Ca*10**15, label='Combined fit')

l4 = ax_Ca.plot([continuum_l.wave, continuum_r.wave], 
               [continuum_l.flux*10**15, continuum_r.flux*10**15],
               marker='o', ls='--', label='continuum')

ax_Ca.set_xlabel(r'Wavelength ($\rm \AA$)')
ax_Ca.set_ylabel(r'Flux ($\rm x10^{-15} erg/cm{^2}/s/\AA$)')

min_indx_Sc = int(np.where(spectrum_Sc.wave == spectrum_Sc.wave[line_indx_Sc][0])[0])
max_indx_Sc = int(np.where(spectrum_Sc.wave == spectrum_Sc.wave[line_indx_Sc][-1])[0])
wave_Sc = spectrum_Sc.wave[min_indx_Sc-15:max_indx_Sc+15]

ax_Sc.axvspan(continuum_l_Sc.wave, continuum_r_Sc.wave, 
              ymin=ax_Sc.get_ylim()[0]*10**15*10**15, ymax = ax_Sc.get_ylim()[1]*10**15*10**15, 
              alpha=0.1, color='k', label='Fit region')
ax_Sc.plot(wave_Sc, spectrum_Sc.flux[min_indx_Sc-15:max_indx_Sc+15]*10**15, label='Spectrum', color=l1[0].get_color())
ax_Sc.plot(spectrum_Sc.wave[line_indx_Sc], (vel_fit(spectrum_Sc.wave[line_indx_Sc])+continuum_Sc)*10**15, label='Fit', color=l3[0].get_color())
ax_Sc.plot([continuum_l_Sc.wave, continuum_r_Sc.wave], np.array([continuum_l_Sc.flux, continuum_r_Sc.flux])*10**15, marker='o', ls='--', label='continuum', color=l4[0].get_color())
ax_Sc.set_xlabel(r'Wavelength ($\rm\AA$)')

ax_Ca.legend([l1[0], l2[0], l3[0], l4[0],l0],[l1[0].get_label(), l2[0].get_label(), l3[0].get_label(), l4[0].get_label(),l0.get_label()],  
             loc='best', fontsize='small', framealpha=1.0,
             bbox_to_anchor=(0.6, -0.04, 0.2, .43))
ax_Sc.legend(loc='best', fontsize='small', framealpha=1.0,
             bbox_to_anchor=(0.56, -0.05, 0.2, .43))
plt.savefig('example_fits.pdf')


 