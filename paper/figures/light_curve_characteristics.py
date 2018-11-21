 
# coding: utf-8

# Creates:
#     * slopes.tex
#     * slope_compare.pdf
#     * abs_mag_vs_slope_sample.pdf

 
# In[1]:


import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
from astropy.time import Time
from astropy.table import Table
from utilities_az import connect_to_sndavis, supernova


 
 
# In[2]:


plt.style.use(['seaborn-paper','az-paper-onecol'])


 
 
# In[3]:


db, cursor = connect_to_sndavis.get_cursor()


 
# # Light Curve Characteristics

 
# In[4]:


peak_mag_query_str = 'SELECT mag, magerr, jd FROM snpeakmagnitude WHERE targetid=322 AND filter="V"'
cursor.execute(peak_mag_query_str)
results = cursor.fetchone()
sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz.get_photometry(band='V')
sn15oz.apparent_mag['V'] = np.append(sn15oz.apparent_mag['V'],results['mag'])
sn15oz.apparent_mag_err['V'] = np.append(sn15oz.apparent_mag_err['V'],results['magerr'])
sn15oz.get_abs_mag(band='V')


 
 
# In[5]:


print('Peak Abs Mag = {} +/- {} on {}'.format(sn15oz.abs_mag['V'][-1], sn15oz.abs_mag_err['V'][-1], results['jd']-sn15oz.jdexpl))
print('first obs: {}'.format(sn15oz.jd['V'][0]-sn15oz.jdexpl))


 
 
# In[6]:


sn15oz.abs_mag['V'][:-1][(sn15oz.jd['V']-sn15oz.jdexpl)>50][0]


 
# # Create the table for paper 

 
# In[7]:


sn15oz = supernova.LightCurve2('asassn-15oz')
sn15oz.get_photometry(band='V')


 
# ## s1 slope

 
# In[8]:


query_str_s1 = 'SELECT DISTINCT slope, slopeerr, tstart, tstop, filter, yvalue, yvalueerr FROM snslope WHERE targetid = 322 AND slopetype = "s1"'
cursor.execute(query_str_s1)
results_s1 = cursor.fetchall()[0]


 
# ## s2 Slope

 
# In[9]:


query_str_s2 = 'SELECT DISTINCT slope, slopeerr, tstart, tstop, filter, yvalue, yvalueerr  FROM snslope WHERE targetid = 322 AND slopetype = "s2"'
cursor.execute(query_str_s2)
results_s2 = cursor.fetchall()[0]


 
# ## s50 Slope

 
# In[10]:


query_str_s50 = 'SELECT DISTINCT slope, slopeerr, tstart, tstop, filter, yvalue, yvalueerr  FROM snslope WHERE targetid = 322 AND slopetype = "s50"'
cursor.execute(query_str_s50)
results_s50 = cursor.fetchall()[0]


 
# ## Tail Slope

 
# In[11]:


query_str_stail = 'SELECT DISTINCT slope, slopeerr, tstart, tstop, filter, yvalue, yvalueerr  FROM snslope WHERE targetid = 322 AND slopetype = "tail"'
cursor.execute(query_str_stail)
results_stail = cursor.fetchall()[0]


 
# ## Latex output

 
# In[12]:


slope_type = ['$s_1$', '$s_2$', '$s_{50V}$', '$s_{tail}$']
slope = [results_s1['slope'], results_s2['slope'], results_s50['slope'], results_stail['slope']]
slope_err = [results_s1['slopeerr'], results_s2['slopeerr'], results_s50['slopeerr'], results_stail['slopeerr']]
jd_start = np.array([results_s1['tstart'], results_s2['tstart'], results_s50['tstart'], results_stail['tstart']])
jd_end = np.array([results_s1['tstop'], results_s2['tstop'], results_s50['tstop'], results_stail['tstop']])
phase_start = jd_start -sn15oz.jdexpl
phase_end = jd_end - sn15oz.jdexpl

tbdata = Table([slope_type, np.array(slope)*50, np.array(slope_err)*50, phase_start, phase_end], 
               names = ['Slope Type', r'Slope (mag $\rm{(50 days)^{-1}}$)', r'Slope Error (mag $\rm{(50 days)^{-1}}$)', 'Start Phase (day)', 'End Phase (day)'])

tbdata[r'Slope (mag $\rm{(50 days)^{-1}}$)'].format = '1.2f'
tbdata[r'Slope Error (mag $\rm{(50 days)^{-1}}$)'].format='1.3f'
tbdata['Start Phase (day)'].format = '2.1f'
tbdata['End Phase (day)'].format = '2.1f'


 
 
# In[13]:


tbdata.write('../slopes.tex', format='latex',
            latexdict={'preamble':r'\centering',
                       'caption':r'The best-fit slope to the V-band light curve of ASASSN-15oz measured between the start and end phase listed in the table. All slopes are measured in units of magnitudes per 50 days.\label{tab:slope}',
                       'data_start':r'\hline',
                       'data_end': r'\hline',
                       'header_start':r'\hline'},
            overwrite=True)

ofile = open('../slopes.tex', 'r')
all_lines = ofile.readlines()
ofile.close()
ofile = open('../slopes.tex', 'w')
all_lines[0] = all_lines[0].replace('table', 'table*')
all_lines[-1] = all_lines[-1].replace('table', 'table*')
for iline in all_lines:
    ofile.write(iline)
ofile.close()


 
# # Plot fits on LC

 
# In[14]:


def get_model(result, sn):
    time = np.arange(result['tstart'], result['tstop']+1)
    phase = time - sn.jdexpl
    app_mag = result['slope']*time + result['yvalue']
    return phase, app_mag

phase_s50, app_mag_s50 = get_model(results_s50, sn15oz)
phase_s1, app_mag_s1 = get_model(results_s1, sn15oz)
phase_s2, app_mag_s2 = get_model(results_s2, sn15oz)
phase_stail, app_mag_stail = get_model(results_stail, sn15oz)



fig = plt.figure()
fig.subplotpars.update(left=0.175)
ax = fig.add_subplot(1,1,1)
ax.errorbar(sn15oz.phase['V'], sn15oz.apparent_mag['V'], sn15oz.apparent_mag_err['V'], fmt='o', capsize=1.0, 
            zorder=0, elinewidth=0.5, markeredgecolor='LightGray', markerfacecolor='LightGray', ecolor='LightGray',
           label='$V$-band light curve')
ax.plot(phase_stail, app_mag_stail, lw=2.0, label=r'$\rm s_{tail}$')
ax.plot(phase_s1, app_mag_s1,   lw=2.0, ls='-', label=r'$\rm s_{1}$', zorder=3)
ax.plot(phase_s50, app_mag_s50, lw=3.0, label=r'$\rm s_{50}$')
ax.plot(phase_s2, app_mag_s2,   lw=2.0, ls='-', label=r'$\rm s_{2}$')

ax.legend()
ax.invert_yaxis()
ax.set_ylabel('Apparent $V$-band magnitude')
ax.set_xlabel('Phase (day)')
plt.savefig('slopes.pdf')


 
# # Place in Context; Make figure for paper of S_tail vs S_50v

 
# In[15]:


query_str = "SELECT supernovanames.name, snslope.targetid,slope, slopeerr, tstart, tstop, filter, slopetype   FROM snslope JOIN idsupernovae ON idsupernovae.id=snslope.targetid   JOIN supernovanames ON idsupernovae.id=supernovanames.targetid   WHERE snslope.filter='V' AND slopetype='s50' AND (sntype IN (6, 17, 18))"
cursor.execute(query_str)
results = cursor.fetchall()
s50vall = []
s50vallerr = []

for idict in results:
    s50vall.append(idict['slope'])
    s50vallerr.append(idict['slopeerr'])
s50vall = np.array(s50vall)
s50vallerr = np.array(s50vallerr)


 
 
# In[16]:


query_str = "SELECT supernovanames.name, snslope.targetid,slope, slopeerr, tstart, tstop, filter, slopetype   FROM snslope JOIN idsupernovae ON idsupernovae.id=snslope.targetid   JOIN supernovanames ON idsupernovae.id=supernovanames.targetid   WHERE snslope.filter='V' AND slopetype='tail' AND (sntype IN (6, 17, 18))"
cursor.execute(query_str)
results = cursor.fetchall()
stailall = []
stailallerr = []

for idict in results:
    stailall.append(idict['slope'])
    stailallerr.append(idict['slopeerr'])
    
stailall = np.array(stailall)
stailallerr = np.array(stailallerr)


 
# ## Get all objects with both a photometric and a tail slope

 
# In[17]:


query_str = "SELECT supernovanames.name, snslope.targetid,slope, slopeerr, tstart, tstop, filter, slopetype   FROM snslope JOIN idsupernovae ON idsupernovae.id=snslope.targetid   JOIN supernovanames ON idsupernovae.id=supernovanames.targetid   WHERE idsupernovae.id IN (SELECT idsupernovae.id      FROM snslope JOIN idsupernovae ON idsupernovae.id=snslope.targetid      WHERE (sntype IN (6, 17, 18)) AND      (snslope.filter='V') AND      ((snslope.slopetype='s50') OR (snslope.slopetype='tail'))      GROUP BY snslope.targetid      HAVING (COUNT(snslope.slopetype) > 1))      AND (snslope.filter='V') AND      ((snslope.slopetype='s50') OR (snslope.slopetype='tail'))      ORDER BY idsupernovae.id, snslope.slopetype"
cursor.execute(query_str)
results = cursor.fetchall()


 
 
# In[18]:


s50v = []
s50verr = []
stail = []
stailerr = []

for sndict50, sndict_tail in zip(results[::2], results[1::2]):
    assert sndict50['slopetype'] == 's50'
    assert sndict_tail['slopetype'] =='tail'
    assert sndict50['name'] == sndict_tail['name']
    s50v.append(sndict50['slope'])
    s50verr.append(sndict50['slopeerr'])
    stail.append(sndict_tail['slope'])
    stailerr.append(sndict_tail['slopeerr'])
s50v = np.array(s50v)
s50verr = np.array(s50verr)
stail = np.array(stail)
stailerr = np.array(stailerr)


 
# ## Plot

 
# In[19]:


from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()  

histogram_bins = np.arange(0, 2.5, 0.1)

left, width = 0.22, 0.55
bottom, height = 0.2, 0.6
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
fig = plt.figure()
fig.subplotpars.update(left=0.5)
axCenter = fig.add_axes(rect_scatter)
axHistx = fig.add_axes(rect_histx)
axHisty = fig.add_axes(rect_histy)

axCenter.axvline(0.5, color='grey', ls=':')
axCenter.text(0.6, 1.05, 'IIL')
axCenter.text(0.2, 1.05, 'IIP')
axCenter.axhline(0.98/2, color='grey', ls='--')
axCenter.text(1.1, 0.26, 'Complete \nTrapping', fontsize=8)

l1 = axCenter.errorbar( s50v*50, stail*50, xerr = s50verr*50, yerr= stailerr*50, fmt='o', label='SN Sample')
l2 = axCenter.errorbar(results_s50['slope']*50, results_stail['slope']*50,   fmt='s',
                  xerr=results_s50['slopeerr']*50, yerr=results_stail['slopeerr']*50,
                  label='ASASSN-15oz', markeredgecolor='k')

axHistx.hist(s50v*50, bins=histogram_bins, edgecolor='k', orientation='vertical')
n50, bins50, patches50 = axHistx.hist(s50vall*50, bins=histogram_bins, edgecolor='k', orientation='vertical', histtype='step', hatch='/')

axHisty.hist(stail*50, bins=histogram_bins, edgecolor='k', orientation='horizontal', label='Both slopes')
ntail, binstail, patchestail = axHisty.hist(stailall*50, bins=histogram_bins, edgecolor='k', orientation='horizontal', label='At least one slope', histtype='step', hatch='/')

# no labels
axCenter.set_ylim(binstail[:-1][ntail>0][0],binstail[1:][ntail>0][-1])
axCenter.set_xlim(bins50[:-1][n50>0][0],bins50[1:][n50>0][-1])


axHisty.yaxis.set_major_formatter(nullfmt)
axHisty.set_yticks(axCenter.get_yticks())
axHisty.set_ylim(axCenter.get_ylim())
axHistx.xaxis.set_major_formatter(nullfmt)
axHistx.set_xticks(axCenter.get_xticks())
axHistx.set_xlim(axCenter.get_xlim())

axHistx.legend(loc='best')

axCenter.set_ylabel(r"$\rm s_{tail}$ (mag/50day)")
axCenter.set_xlabel(r"$\rm s_{50V}$ (mag/50day)")
axCenter.legend([l1[0], l2[0]], ['SN Sample', 'ASASSN-15oz'], loc=(0.05, 0.8))#bbox_to_anchor=(0.25,0.75, 0.2, 0.2))
plt.savefig('slope_compare.pdf')
print(len(s50v))


 
# # Make LC comparison picture

 
# In[23]:


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot([], [], label=None)
sn_list = ['asassn-15oz','1999em', '2004et', '2013by', '2013ej', '2005cs']
labels = ['ASASSN-15oz', 'SN1999em', 'SN2004et', 'SN2013by', 'SN2013ej', 'SN2005cs']
for snname, snlabel in zip(sn_list, labels):
    sn = supernova.LightCurve2(snname)
    sn.get_photometry()
    sn.get_abs_mag()
    indx = np.argsort(sn.phase['V'])
    ax.plot(sn.phase['V'][indx], sn.abs_mag['V'][indx], '.', ls='-', label=snlabel)
    if snname == 'asassn-15oz':
        ax.set_ylim(np.min(sn.abs_mag['V']), np.max(sn.abs_mag['V'])+1)
        ax.set_xlim(0, 250)
    
ax.invert_yaxis()
ax.legend()
ax.set_xlabel('Phase (day)')
ax.set_ylabel(r'Absolute $V$-band magnitude')
plt.savefig('lc_comparison.pdf')


 
# # Make Figure for Paper of S_50v vs Abs Mag

 
# In[25]:


query_str = "SELECT DISTINCT supernovanames.name, snslope.targetid,slope, slopeerr, slopetype, mag, magerr, mu, muerr, ebvg, ebvgerr, ebvi, ebvierr FROM snslope JOIN idsupernovae ON idsupernovae.id=snslope.targetid JOIN supernovanames ON idsupernovae.id=supernovanames.targetid JOIN snpeakmagnitude ON idsupernovae.id=snpeakmagnitude.targetid WHERE (sntype IN (6, 17, 18)) AND (slopetype='s50') AND (snslope.filter='V') AND (snpeakmagnitude.filter='V') ORDER BY snslope.targetid"
cursor.execute(query_str)
results = cursor.fetchall()


 
 
# In[26]:


plt.style.use(['seaborn-paper','az-paper-onecol'])
name = []
slope = []
slopeerr = []
absmag = []
absmagerr = []
sn15oz_dict = {}
special_sn_dict = {}
for idict in results:
    if (idict['mu'] is not None) and (idict['muerr'] is not None):
        name.append(idict['name'])
        slope.append(idict['slope'])
        slopeerr.append(idict['slopeerr'])
        if idict['ebvi'] is None:
            idict['ebvi'] = 0
        iabs_mag, iabs_mag_err = supernova.calc_abs_mag(idict['mag'], idict['mu'], idict['ebvg'], 'V', 
                                         host_ebv=idict['ebvi'], dist_mod_err=idict['muerr'], app_mag_err=idict['magerr'])
        absmag.append(iabs_mag)
        absmagerr.append(iabs_mag_err)
        if idict['targetid'] == 322: #ASASSN-15oz
            sn15oz_dict = idict
            sn15oz_dict['absmag'] = iabs_mag
            sn15oz_dict['absmagerr'] = iabs_mag_err
        if idict['name'] in sn_list[1:]:
            special_sn_dict[idict['name']] = {}
            special_sn_dict[idict['name']]['absmag'] = iabs_mag
            special_sn_dict[idict['name']]['absmagerr'] = iabs_mag_err
            special_sn_dict[idict['name']]['slope'] = idict['slope']
            special_sn_dict[idict['name']]['slopeerr'] = idict['slopeerr']
name = np.array(name)
slope = np.array(slope)
slopeerr = np.array(slopeerr)
absmag = np.array(absmag)
absmagerr = np.array(absmagerr)
print(len(slope))
    
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.axvline(0.5, color='gray', ls=':')
ax.text(0.1,-21, 'IIP')
ax.text(0.6,-21, 'IIL')
ax.errorbar(slope*50, absmag, xerr=slopeerr*50, yerr=absmagerr, fmt='o', label='SN Sample')
for snname, snlabel in zip(sn_list, labels):
    if snname == 'asassn-15oz':
        ax.errorbar(sn15oz_dict['slope']*50, sn15oz_dict['absmag'], 
            xerr=sn15oz_dict['slopeerr']*50, yerr=sn15oz_dict['absmagerr'],
           fmt='s', label=snlabel, markeredgecolor='k')
    else:
        ax.errorbar(special_sn_dict[snname]['slope']*50, special_sn_dict[snname]['absmag'], 
            xerr=special_sn_dict[snname]['slopeerr']*50, yerr=special_sn_dict[snname]['absmagerr'],
           fmt='o', label=snlabel, markeredgecolor='w')


ax.invert_yaxis()
ax.set_xlabel(r'$\rm s_{50V}$ (mag/50day)')
ax.set_ylabel('Absolute $V$-band magnitude')
ax.legend(loc='lower right')

plt.savefig('abs_mag_vs_slope_sample.pdf')


 
# ## Combine figures 

 
# In[29]:


plt.style.use(['seaborn-paper','az-paper-twocol'])
name = []
slope = []
slopeerr = []
absmag = []
absmagerr = []
sn15oz_dict = {}
special_sn_dict = {}
for idict in results:
    if (idict['mu'] is not None) and (idict['muerr'] is not None):
        name.append(idict['name'])
        slope.append(idict['slope'])
        slopeerr.append(idict['slopeerr'])
        if idict['ebvi'] is None:
            idict['ebvi'] = 0
        iabs_mag, iabs_mag_err = supernova.calc_abs_mag(idict['mag'], idict['mu'], idict['ebvg'], 'V', 
                                         host_ebv=idict['ebvi'], dist_mod_err=idict['muerr'], app_mag_err=idict['magerr'])
        absmag.append(iabs_mag)
        absmagerr.append(iabs_mag_err)
        if idict['targetid'] == 322: #ASASSN-15oz
            sn15oz_dict = idict
            sn15oz_dict['absmag'] = iabs_mag
            sn15oz_dict['absmagerr'] = iabs_mag_err
        if idict['name'] in sn_list[1:]:
            special_sn_dict[idict['name']] = {}
            special_sn_dict[idict['name']]['absmag'] = iabs_mag
            special_sn_dict[idict['name']]['absmagerr'] = iabs_mag_err
            special_sn_dict[idict['name']]['slope'] = idict['slope']
            special_sn_dict[idict['name']]['slopeerr'] = idict['slopeerr']
name = np.array(name)
slope = np.array(slope)
slopeerr = np.array(slopeerr)
absmag = np.array(absmag)
absmagerr = np.array(absmagerr)
fig = plt.figure()
fig.subplotpars.update(wspace=0, left=0.1, right=0.975)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax2.axvline(0.5, color='gray', ls=':')
ax2.text(0.1,-21, 'IIP')
ax2.text(0.6,-21, 'IIL')
ax2.errorbar(slope*50, absmag, xerr=slopeerr*50, yerr=absmagerr, fmt='o', label='SN Sample')
for snname, snlabel in zip(sn_list, labels):
    if snname == 'asassn-15oz':
        ax2.errorbar(sn15oz_dict['slope']*50, sn15oz_dict['absmag'], 
            xerr=sn15oz_dict['slopeerr']*50, yerr=sn15oz_dict['absmagerr'],
           fmt='s', label=snlabel, markeredgecolor='k')
    else:
        ax2.errorbar(special_sn_dict[snname]['slope']*50, special_sn_dict[snname]['absmag'], 
            xerr=special_sn_dict[snname]['slopeerr']*50, yerr=special_sn_dict[snname]['absmagerr'],
           fmt='o', label=snlabel, markeredgecolor='Gray')

ax2.set_yticklabels([])
ax2.set_ylim(-22.5, -9.5)
ax2.invert_yaxis()
ax2.set_xlabel(r'$\rm s_{50V}$ (mag/50day)')
ax2.legend(loc='lower right')

ax1.plot([], [], label=None)
sn_list = ['asassn-15oz','1999em', '2004et', '2013by', '2013ej', '2005cs']
labels = ['ASASSN-15oz', 'SN1999em', 'SN2004et', 'SN2013by', 'SN2013ej', 'SN2005cs']
for snname, snlabel in zip(sn_list, labels):
    sn = supernova.LightCurve2(snname)
    sn.get_photometry()
    sn.get_abs_mag()
    indx = np.argsort(sn.phase['V'])
    ax1.plot(sn.phase['V'][indx], sn.abs_mag['V'][indx], '.', ls='-', label=snlabel)
    if snname == 'asassn-15oz':
        #ax1.set_ylim(np.min(sn.abs_mag['V']), np.max(sn.abs_mag['V']))
        ax1.set_xlim(0, 250)
    
ax1.set_ylim(-22.5, -9.5)
ax1.invert_yaxis()
ax1.legend()
ax1.set_xlabel('Phase (day)')
ax1.set_ylabel(r'Absolute $V$-band magnitude')
plt.savefig('lc_comp_abs_slope.pdf')


 