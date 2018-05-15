import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import numpy as np

from astropy.table import Table, vstack
from astropy.io import fits
from astropy.io import ascii as asc

import yaml

import os
import line_analysis_BSNIP



FIG_DIR = '../figures'
OUTPUT_DIR = '../data/line_info'
input_file = '../data/line_vel_analysis_input.yml'
#plt.ioff()
append = False
ofile = open(input_file)
feature_dict = yaml.load(ofile)
for iline in feature_dict.keys():
    line_table = Table(names=['line', 'date', 'velocity', 'vel_err_left', 'vel_err_right', 'pEW', 'pEW_err'], dtype=('S15', 'S15', 'f8', 'f8', 'f8', 'f8', 'f8'))
    with PdfPages(os.path.join(FIG_DIR, '{}_analysis.pdf'.format(iline))) as pdf:
        for ifile in feature_dict[iline]['files']:
            for idir_key in feature_dict[iline]['directories'].keys():
                filename = os.path.join(feature_dict[iline]['directories'][idir_key], ifile)
                if os.path.exists(filename):
                    found = True
                    break
            if not found:
                print('Unable to locate file {}'.format(ifile))
            else:
                if filename.endswith('.fits'):
                    date = fits.getval(filename, 'date-obs', 0)
                else:
                    date = input('enter date for file: {}'.format(ifile))
                line_return = line_analysis_BSNIP.characterize_line(feature_dict[iline], filename, visualization_level=2)
                if line_return is not None:
                    pew, pew_var, vel_min, vel_err, fig = line_return
                    line_table.add_row([iline, date, vel_min, vel_err[0], vel_err[1], pew, np.sqrt(pew_var)])
                    pdf.savefig(fig)
                    plt.close(fig)
    if os.path.exists(os.path.join(OUTPUT_DIR,'{}.tab'.format(iline))):
        if append is True:
            old_file = Table.read(os.path.join(OUTPUT_DIR,'{}.tab'.format(iline)), format='ascii.tab')
            new_table = vstack([old_file, line_table])
            new_table.write(os.path.join(OUTPUT_DIR,'{}.tab'.format(iline)), overwrite=True, format='ascii.tab')
        else:
            line_table.write(os.path.join(OUTPUT_DIR, '{}.tab'.format(iline)), overwrite=True, format='ascii.tab')
    else:
        line_table.write(os.path.join(OUTPUT_DIR, '{}.tab'.format(iline)), format='ascii.tab')
plt.ion()