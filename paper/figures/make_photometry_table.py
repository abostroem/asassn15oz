 
# coding: utf-8

# Creates a table of all imaging observations in the database for paper:
# * lc_obs.tex

 
# In[10]:


import numpy as np
from astropy import table
from astropy.table import Table
from astropy.time import Time

from utilities_az import supernova, connect_to_sndavis


 
 
# In[11]:


db, cursor = connect_to_sndavis.get_cursor()


 
 
# In[12]:


sn15oz = supernova.LightCurve2('asassn-15oz')


 
 
# In[13]:


query_str = '''
SELECT DISTINCT mag, magerr, BINARY(filter), jd, source 
FROM photometry 
WHERE targetid = 322
ORDER BY jd'''


 
 
# In[14]:


cursor.execute(query_str)
results = cursor.fetchall()


 
 
# In[15]:


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
13: {'short': 'CPT 1m', 'long': 'SAAO - Sutherland Facilities - 1m '},
88: {'short': 'Swift', 'long': 'Swift'}} #cpt1m013-kb76 1m0


 
 
# In[18]:


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
               names=['Date-Obs','JD', 'Phase',
                      'Apparent Magnitude', 
                      'Apparent Magnitude Error', 
                      'Filter', 
                      'Source'])
tbdata.sort(keys=['JD', 'Filter'])


 
 
# In[19]:


tbdata.write('../lc_obs.tex', format='aastex', 
             formats={'JD':'%8.2f', 
                      'Phase':'%4.2f',
                      'Apparent Magnitude':'%2.2f',
                      'Apparent Magnitude Error': '%1.2f'}, overwrite=True,
            latexdict={'preamble':r'\centering',
                       'caption':r'Imaging Observations of ASASSN-15oz.\label{tab:LcObs}',
                       'data_start':r'\hline'})


 