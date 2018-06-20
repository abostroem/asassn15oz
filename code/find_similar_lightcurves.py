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
import matplotlib
matplotlib.style.use('seaborn-colorblind')

from scipy.interpolate import interp1d

#Installed for this project
import extinction

#Mine
import visualization
#get_ipython().magic('matplotlib inline')
import connect_to_sndavis
import define_filters
import supernova

FIG_DIR = '../figures'

def build_sn_list():
    db, cursor = connect_to_sndavis.get_cursor()
    # 6 = II, 17=IIP-like, 18=IIL-like
    query_str = '{} {} {} {} {} {} {} {} {} {}'.format(
                'SELECT DISTINCT idsupernovae.`id`, sntype,name, slope, slopetype', 
                'FROM idsupernovae', 
                'JOIN supernovanames', 
                'ON idsupernovae.`id` = supernovanames.`targetid`',
                'JOIN snslope',
                'ON idsupernovae.`id` = snslope.`targetid` ',
                'WHERE (sntype = 6 ',
                'OR sntype = 17',
                'OR sntype = 18)',
                "AND slopetype = 's50';")
    query = cursor.execute(query_str)
    results = cursor.fetchall()
    id = []
    name = []
    slope = []
    for idict in results:
        id.append(idict['id'])
        name.append(idict['name'])
        slope.append(idict['slope'])
    tbdata = table.Table([id, name, slope], names = ['id', 'name', 's50'])
    return tbdata

def find_closest_slope(tbdata):
    slope_diff = np.abs((tbdata['s50'] - tbdata['s50'][tbdata['name']=='ASASSN-15oz']))
    indx = np.argsort(slope_diff)
    return indx

def compare_sn(snname1, snname2, rank, band='all', sn2_phase_offset = 0):
    sn1 = supernova.LightCurve2(snname1)
    sn1.get_photometry(band=band)
    sn2 = supernova.LightCurve2(snname2)
    sn2.get_photometry(band=band)
    common_bands = set(sn1.apparent_mag.keys())&(sn2.apparent_mag.keys())
    fig = plt.figure(figsize=[8.5, 11])
    for plot_num, iband in enumerate(common_bands):
        if plot_num == 5:
            plt.savefig(os.path.join(FIG_DIR, 'similar_lc_{}_2.pdf'.format(sn2.name)))
            plt.close()
            fig = plt.figure(figsize=[8.5, 11])
        ax = fig.add_subplot(3, 2, plot_num%6+1)
        ax.plot(sn1.phase[iband], sn1.apparent_mag[iband]/sn1.apparent_mag[iband].min(), 'o', label = sn1.name, markersize=2)
        ax.plot(sn2.phase[iband]+sn2_phase_offset, sn2.apparent_mag[iband]/sn2.apparent_mag[iband].min(), 's', label = sn2.name, markersize=2)
        ax.set_title('{}, {} band, rank={}'.format(sn2.name, iband, rank), fontsize='small')
        ax.set_ylim(ymax=np.max((sn1.apparent_mag[iband])[sn1.phase[iband] < 100])/sn1.apparent_mag[iband].min()+0.05)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlim(0, 100)
        ax.legend(loc='best', fontsize='xx-small')
    fig.tight_layout()
    
    plt.savefig(os.path.join(FIG_DIR, 'similar_lc_{}.pdf'.format(sn2.name)))
    plt.close()

if __name__ == "__main__":
    num_sn = 10
    tbdata = build_sn_list()
    best_match_indx = find_closest_slope(tbdata)
    print(tbdata[best_match_indx[1:num_sn+1]])
    for rank, sn_indx in enumerate(best_match_indx[1:num_sn+1]): #start at 1 b/c the SN is always the best match to itself
        compare_sn('ASASSN-15oz', tbdata['name'][sn_indx], rank+1)
    compare_sn('ASASSN-15oz', '2016zb', 1, sn2_phase_offset = 8)



