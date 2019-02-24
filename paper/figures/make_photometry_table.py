 
# coding: utf-8

# Creates a table of all imaging observations in the database for paper:
# * lc_obs.tex

 
# In[1]:


import numpy as np
from astropy import table
from astropy.table import Table
from astropy.time import Time

from utilities_az import supernova, connect_to_sndavis


 
 
# In[2]:


db, cursor = connect_to_sndavis.get_cursor()


 
 
# In[3]:


sn15oz = supernova.LightCurve2('asassn-15oz')


 
 
# In[4]:


query_str = '''
SELECT DISTINCT mag, magerr, BINARY(filter), jd, source 
FROM photometry 
WHERE targetid = 322
ORDER BY jd'''


 
 
# In[5]:


cursor.execute(query_str)
results = cursor.fetchall()


 
 
# In[6]:


loc_dict = {
1: {'short':'OGG 2m', 'long': 'Haleakala Observatory - 2m'}, #ogg2m001-fs02   2m0
2: {'short':'COJ 2m', 'long': 'Siding Springs Observatory - 2m'},  #coj2m002-fs03
3: {'short':'COJ 1m', 'long': 'Siding Springs Observatory - 1m'},  #coj1m003-kb71
4: {'short':'LSC 1m', 'long': 'CTIO - Region IV'}, #lsc1m004-kb77 1m0
5: {'short':'LSC 1m', 'long': 'CTIO - Region IV'}, #lsc1m005-kb78 1m0
8: {'short':'ELP 1m', 'long': 'McDonald Observatory - 1m'},#elp1m008-kb74 1m0
9: {'short':'LSC 1m', 'long': 'CTIO - Region IV'}, #lsc1m009-fl03 1m0
10: {'short': 'CPT 1m', 'long': 'SAAO - Sutherland Facilities - 1m'}, #cpt1m010
11: {'short': 'COJ 1m', 'long': 'Siding Springs Observatory - 1m'}, #coj1m011-kb05 1m0
12: {'short': 'CPT 1m', 'long': 'SAAO - Sutherland Facilities - 1m'}, #cpt1m012-kb75 1m0
13: {'short': 'CPT 1m', 'long': 'SAAO - Sutherland Facilities - 1m '},#cpt1m013-kb76 1m0
88: {'short': 'Swift', 'long': 'Swift'}, #Swift
-88: {'short': 'Swift', 'long': 'Swift'}} #Swift; non-detections in the DB have negative sources


 
 
# In[7]:


band = []
jd = []
date = []
mag = []
mag_err = []
source = []
phase = []
for iresult in results:
    #rename the Swift filters to their proper names
    ifilter = iresult['BINARY(filter)']
    if ifilter in [b'us', b'bs', b'vs']:
        band.append(ifilter.decode('utf-8')[0])
    elif ifilter in [b'uw1', b'uw2', b'um2']:
        band.append('uv'+ifilter.decode('utf-8')[1:])
    else:
        band.append(ifilter)
    jd.append(iresult['jd'])
    mag.append(iresult['mag'])
    mag_err.append(iresult['magerr'])
    if ifilter in [b'J', b'H', b'K']:
        source.append('NTT')
    else:
        source.append(loc_dict[iresult['source']]['short'])
    date.append(Time(iresult['jd'], format='jd', out_subfmt='date').iso)
    phase.append((Time(iresult['jd'], format='jd') - Time(sn15oz.jdexpl, format='jd')).value)
tbdata = Table([date, jd, phase, mag, mag_err, band, source], 
               names=['Date-Obs','JD', 'Phase (Day)',
                      'Apparent Magnitude', 
                      'Apparent Magnitude Error', 
                      'Filter', 
                      'Source'])
tbdata.sort(keys=['JD', 'Filter'])


 
 
# In[8]:


#tbdata.write('../lc_obs.tex', format='aastex', 
#             formats={'JD':'%8.2f', 
#                      'Phase (Day)':'%4.2f',
#                      'Apparent Magnitude':'%2.2f',
#                      'Apparent Magnitude Error': '%1.2f'}, overwrite=True,
#            latexdict={'preamble':r'\centering',
#                       'caption':r'Imaging Observations of ASASSN-15oz.\label{tab:LcObs}',
#                       'data_start':r'\hline'})


 
 
# In[9]:


tbdata_short = tbdata[0:5].copy()
tbdata_short.write('../lc_obs_short.tex', format='latex', 
             formats={'JD':'%8.2f', 
                      'Phase (Day)':'%4.2f',
                      'Apparent Magnitude':'%2.2f',
                      'Apparent Magnitude Error': '%1.2f'}, overwrite=True,
            latexdict={'preamble':r'\centering',
                       'caption':r'Sample of Imaging Observations of ASASSN-15oz. Full table available on-line.\label{tab:LcObs}',
                       'data_start':r'\hline',
                       'data_end':r'\hline',
                       'header_start':r'\hline',
                       'tabletype': 'table*'})


 
 
# In[10]:


#tbdata.write('../lc_obs.dat',  overwrite=True, format='ascii')
#ofile = open('../lc_obs.dat', 'r')
#all_lines = ofile.readlines()
#ofile.close()
#header = '''#Photometric observations of ASASSN-15oz.
##Columns:
##Date-Obs: (str) Human readable date of observation
##JD: (float) Julian Date of observation
##Phase: (float) Days since explosion, where explosion is defined as {}
##Apparent Magnitude: (float)
##Apparent Magntidue Error: (float)
##Filter: (str) Filter used for observation
##Source: (str) Observatory used to take the data. OGG, COJ, LSC,  ELP, and CPT are Las Cumbres Observatory Telescopes.\n
#'''.format(sn15oz.jdexpl)
#ofile = open('../asassn15oz_lc_obs.dat', 'w')
#ofile.write(header)
#for iline in all_lines[1:]:
#    ofile.write(iline)
#ofile.close()


 
 
# In[11]:


tbdata.write('../lc_obs.csv',  overwrite=True)
ofile = open('../lc_obs.csv', 'r')
all_lines = ofile.readlines()
ofile.close()
header = '''#Photometric observations of ASASSN-15oz.
#Columns:
#Date-Obs: (str) Human readable date of observation
#JD: (float) Julian Date of observation
#Phase: (float) Days since explosion, where explosion is defined as {}
#Apparent Magnitude: (float)
#Apparent Magntidue Error: (float)
#Filter: (str) Filter used for observation
#Source: (str) Observatory used to take the data. OGG, COJ, LSC,  ELP, and CPT are Las Cumbres Observatory Telescopes.\n
'''.format(sn15oz.jdexpl)
ofile = open('../asassn15oz_lc_obs.csv', 'w')
ofile.write(header)
for iline in all_lines[1:]:
    ofile.write(iline)
ofile.close()


 