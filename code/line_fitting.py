import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import fits
from astropy.time import Time
from astropy.table import Table, Column
import yaml

from spectroscopy import calc_wavelength
from fit_spectral_lines2 import fit_feature
import line_analysis_BSNIP


SPEC_DIR_LCO = '../data/spectra/lco'
SPEC_DIR_EFOSC = '../data/spectra/EFOSC'
FIG_DIR = '../figures'
OUTPUT_DIR = '../data/line_info'

append = False

line_list = {'FeII_multi':{'num_components':2, 
                           'file_indx': np.int_(np.arange(3,14))},
             'HA-cachito':{'num_components':2,
                            'file_indx': np.int_(np.arange(0,14))},
             'HB':{'num_components':1,
                            'file_indx': np.int_(np.arange(0,14))},
             'NaI':{'num_components':1,
                            'file_indx': np.int_(np.arange(9,14))}, 
             'OI':{'num_components': 1,
                            'file_indx': np.int_(np.arange(5,14))}
             } 

spectrum_filelist = ['asassn15oz_20150904_redblu_122216.314.fits', #0
                     'asassn-15oz_20150906_redblu_105042.698a.fits', #1 The other spectrum from this night is with a worse sensitivity function --> lower S/N
                     'asassn15oz_20150907_redblu_123835.277.fits', #2
                     'asassn15oz_20150911_redblu_105336.349.fits', #3
                     'asassn15oz_20150916_redblu_120911.274.fits', #4
                     'asassn-15oz_20150920_redblu_135034.512.fits', #5
                     'asassn-15oz_20150924_redblu_123847.580.fits', #6
                     'asassn-15oz_20150930_redblu_122858.217.fits', #7
                     'tASASSN-15oz_20151003_Gr13_Free_slit1.0_57720_1_e.fits', #8
                     'asassn15oz_20151006_redblu_101906.800.fits', #9
                     'asassn15oz_20151014_redblu_112918.305.fits', #10
                     'asassn-15oz_20151025_redblu_102221.833.fits', #11
                     'tASAS-SN_15oz_20151107_Gr13_Free_slit1.5_57723_1_e.fits', #12
                     'tASAS-SN_15oz_20151118_Gr13_Free_slit1.0_57723_1_e.fits'] #13

for iline in line_list.keys():
    with PdfPages(os.path.join(FIG_DIR, '{}_analysis.pdf'.format(iline))) as pdf:
        print('***************** Now Fitting {} *****************'.format(iline))
        line_table = Table(names=['line', 'date'], dtype=('S15', 'S25'))
        for component in range(line_list[iline]['num_components']):
            line_table.add_column(Column(name='vel{}'.format(component), dtype='f8'))
            line_table.add_column(Column(name='vel_err_left_{}'.format(component), dtype='f8'))
            line_table.add_column(Column(name='vel_err_right_{}'.format(component), dtype='f8'))
            line_table.add_column(Column(name='vel_pew_{}'.format(component), dtype='f8'))
            line_table.add_column(Column(name='vel_pew_err{}'.format(component), dtype='f8'))
        for spectrum_filename in np.array(spectrum_filelist)[line_list[iline]['file_indx']]:
            print(spectrum_filename)
            if os.path.exists(os.path.join(SPEC_DIR_LCO, spectrum_filename)):
                filename = os.path.join(SPEC_DIR_LCO, spectrum_filename)
                binsize=20
            else:
                filename = os.path.join(SPEC_DIR_EFOSC,spectrum_filename)
                binsize=5
            spectrum = line_analysis_BSNIP.read_iraf_spectrum(filename)
            min_list, pew_list, fig = fit_feature(spectrum, iline, binsize, absorption=True, similar_widths=True, input_filename=os.path.join(OUTPUT_DIR, '{}_input.yaml'.format(iline)))
            if filename.endswith('.fits'):
                date = fits.getval(filename, 'date-obs', 0)
            else:
                date = input('enter date for file: {}'.format(ifile))

            irow = [iline, date]
            for imin_component, ipew_component in zip(min_list, pew_list):
                [irow.append(icomponent) for icomponent in [imin_component[0], imin_component[1], imin_component[2], ipew_component[0], np.sqrt(ipew_component[1])]]
            line_table.add_row(irow)
            pdf.savefig(fig)
            plt.close()
    

    if os.path.exists(os.path.join(OUTPUT_DIR,'{}.tab'.format(iline))):
        if append is True:
            old_file = Table.read(os.path.join(OUTPUT_DIR,'{}.tab'.format(iline)), format='ascii.tab')
            new_table = vstack([old_file, line_table])
            new_table.write(os.path.join(OUTPUT_DIR,'{}.tab'.format(iline)), overwrite=True, format='ascii.tab')
        else:
            line_table.write(os.path.join(OUTPUT_DIR, '{}.tab'.format(iline)), overwrite=True, format='ascii.tab')
    else:
        line_table.write(os.path.join(OUTPUT_DIR, '{}.tab'.format(iline)), format='ascii.tab')   
    