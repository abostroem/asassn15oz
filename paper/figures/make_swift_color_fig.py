 
# coding: utf-8

# Creates:
#     * uv_slope_comp.pdf

 
# In[1]:


import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
#get_ipython().run_line_magic('matplotlib', 'inline')

from utilities_az import supernova,visualization,define_filters,connect_to_sndavis


 
 
# In[2]:


plt.style.use(['seaborn-paper', 'az-paper-onecol'])


 
 
# In[3]:


sn15oz = supernova.LightCurve2('ASASSN-15oz')
sn15oz.get_photometry()


 
# # Look at lots of supernovae

 
# In[4]:


db, cursor = connect_to_sndavis.get_cursor()


 
# ## Select all SN with uw2 and vs data

 
# In[5]:


sql_query = '''SELECT DISTINCT name, slope
FROM photometry
JOIN supernovanames ON photometry.targetid=supernovanames.targetid
JOIN snslope ON photometry.targetid=snslope.targetid
WHERE photometry.filter='uw2' AND slopetype='s50';'''
cursor.execute(sql_query)
results = cursor.fetchall()
snname_uw2=[]
snslope_uw2=[]
for i in results:
    snname_uw2.append(i['name'])
    snslope_uw2.append(i['slope'])
snslope_uw2 = np.array(snslope_uw2)


 
 
# In[6]:


sql_query = '''SELECT DISTINCT name
FROM photometry
JOIN supernovanames ON photometry.targetid=supernovanames.targetid
WHERE filter='vs';'''
cursor.execute(sql_query)
results = cursor.fetchall()
snname_vs=[]
for i in results:
    snname_vs.append(i['name'])


 
 
# In[7]:


snlist = list(set(snname_uw2) & set(snname_vs))
max_slope = np.max(snslope_uw2)*100+0.001


 
 
# In[8]:


cm_rainbow = visualization.make_rainbow_cm()


 
 
# In[14]:


fig = plt.figure()
fig.subplotpars.update(left=0.18, right=0.95, top=0.85)
gs1 = mpl.gridspec.GridSpec(9, 1)
ax1 = plt.subplot(gs1[1:10,0])
ax2 = plt.subplot(gs1[0,0])
snlist.sort()
for sn in snlist:
    if sn == 'ASASSN-15oz' or sn == '2009kr':
        continue
    slope_color=cm_rainbow(float(snslope_uw2[np.array(snname_uw2)==sn]*100/max_slope))
    snobj = supernova.LightCurve2(sn)
    snobj.get_photometry()
    if sn.startswith('ASA'):
        label = 'ASASSN-{}'.format(sn.split('-')[1])
    else:
        label=sn
    vs_interp = np.interp(snobj.jd['uw2'], snobj.jd['vs'], snobj.apparent_mag['vs'])
    jderr = np.ones(len(snobj.jd['uw2']))*snobj.jdexpl_err
    #ax1.errorbar(snobj.jd['uw2']-snobj.jdexpl, snobj.apparent_mag['uw2']-vs_interp, xerr=jderr, fmt='o', c=slope_color, label=label, ecolor='k', elinewidth=0.5)
    ax1.errorbar(snobj.jd['uw2']-snobj.jdexpl, snobj.apparent_mag['uw2']-vs_interp, xerr=jderr, fmt='.', c=slope_color, label=label, elinewidth=0.5)
vs_interp = np.interp(sn15oz.jd['uw2'], sn15oz.jd['vs'], sn15oz.apparent_mag['vs'])
slope_color=cm_rainbow(float(snslope_uw2[np.array(snname_uw2)=='ASASSN-15oz']*100/max_slope))
jderr_15oz = np.ones(len(sn15oz.jd['uw2']))*sn15oz.jdexpl_err
ax1.errorbar(sn15oz.jd['uw2']-sn15oz.jdexpl, sn15oz.apparent_mag['uw2']-vs_interp, xerr=jderr_15oz, fmt='s', c=slope_color, elinewidth=0.5, ecolor='k',label='ASASSN-15oz', mec='k', zorder=15)
ax1.set_xlim(0,42)
ax1.set_xlabel('Phase (day)')
ax1.set_ylabel('UVW2-V (mag)')
ax1.legend(ncol=2, bbox_to_anchor=[0.31, -0.1, 0.2, 0.5], framealpha=1)
#ax1.legend( ncol=1, framealpha=1.0, loc='lower right')
ax1.set_ylim(-2.5, 5.75)
norm = mpl.colors.Normalize(vmin=0, vmax=1.0/100*max_slope*50)
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm_rainbow,
                                norm=norm, orientation='horizontal')
ax2.tick_params(bottom=False,top=True, which='major',
               labelbottom=False, labeltop=True, pad=0.5)
ax1.text(7, 8.0,'Slope (mag/50day)')
plt.savefig('uv_slope_comp.pdf')


 