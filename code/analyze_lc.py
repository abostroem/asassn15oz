
import sys
sys.path.append('/Users/bostroem/Desktop/research/not_my_code/SNEC-1.01')
import numpy as np

import supernova as sn
import chisq_analysis

energies = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.5]
masses = [12, 15, 17, 19]
ni_mixing = [5.0]
time_offsets = np.arange(-4, 4, 1)
Kvalues = [1.0E17, 10.0E17, 20.0E17, 30.0E17, 40.0E17, 50.0E17, 60.0E17]
Kvalues_str = ['1.0E17', '10.0E17', '20.0E17', '30.0E17', '40.0E17', '50.0E17', '60.0E17']
radii = [700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]

snec_models = '/Users/bostroem/Desktop/research/snec_models_pregrid'
snname = 'asassn15oz'
S2_start = 50
S2_end = 88  #Vary this parameter
#time_offsets = np.array([0])


sn15oz = sn.LightCurve2(snname)
sn15oz.get_photometry()
sn15oz.get_abs_mag()
    
snec_15oz = chisq_analysis.SnecAnalysis(snname, snec_models, S2_start, S2_end, 
                 ni_mixing, masses, energies, time_offsets, 
                 Kvalues=Kvalues, Kvalues_str = Kvalues_str, radii=radii, fig_dir='../figures')
    
snec_15oz.calc_chisq_base_model(sn15oz)
snec_15oz.plot_chisq_base_mod()
snec_15oz.calc_chisq_wind(sn15oz)
snec_15oz.plot_chisq_csm_mod()
snec_15oz.plot_lightcurve(sn15oz)