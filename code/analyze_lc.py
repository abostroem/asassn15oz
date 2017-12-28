
import sys
sys.path.append('/Users/bostroem/Desktop/research/not_my_code/SNEC-1.01')
import numpy as np

import supernova as sn
import chisq_analysis

ni_mass = [0.083, 0.0965, 0.11]
energies = [0.5, 0.8, 1.1, 14, 1.7, 2.0]
masses = [11, 13, 15, 16, 17, 18, 21]
ni_mixing = [5.0]
time_offsets = np.arange(-4, 4, 1)
Kvalues = [10, 20, 30, 35, 40, 50, 60]
radii = [1500, 1800, 2100, 2400, 2700, 3000, 3300]

snec_models = '/Users/bostroem/Desktop/research/not_my_code/SNEC-1.01/snec_models/'
snname = 'asassn15oz'
S2_start = 50
S2_end = 88  #Vary this parameter
#time_offsets = np.array([0])


#sn15oz = sn.LightCurve2(snname)
#sn15oz.get_photometry()
#sn15oz.get_abs_mag()
   
 
snec_15oz = chisq_analysis.SnecAnalysis(snname, snec_models, S2_start, S2_end, 
                 ni_mass, ni_mixing, masses, energies, time_offsets, 
                 Kvalues, radii, fig_dir='../figures')
    
snec_15oz.read_chisq()
snec_15oz.plot_2D('Progenitor Mass', 'Explosion Energy')

input('Enter to close')