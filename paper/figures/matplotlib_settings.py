import matplotlib as mpl
from cycler import cycler
import numpy as np
from matplotlib import pyplot as plt

cols = [(0,0,0)]
for x in np.linspace(0,1, 254):
    rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
    gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
    bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
    cols.append((rcol, gcol, bcol))
cols.append((1,1,1))
cm_rainbow = mpl.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)


plt.style.use('seaborn-paper')
mpl.rcParams['axes.prop_cycle'      ]= cycler('color', ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 
        '#CC6677', '#882255', '#AA4499', '#661100',  '#AA4466','#4477AA'])
mpl.rcParams['legend.fontsize'      ]= 7.0
mpl.rcParams['figure.figsize'       ]= [3.5, 3.5]
mpl.rcParams['figure.subplot.bottom']= 0.175
#mpl.rcParams['figure.subplot.hspace']= 0.2
mpl.rcParams['figure.subplot.left'  ]= 0.225
mpl.rcParams['figure.subplot.right' ]= 0.95
mpl.rcParams['figure.subplot.top'   ]= 0.95
#mpl.rcParams['figure.subplot.wspace']= 0.245
mpl.rcParams['lines.markersize'     ]= 8
mpl.rcParams['savefig.pad_inches'   ]= 0
mpl.rcParams['image.cmap'           ]= 'viridis'
mpl.rcParams['lines.linewidth'      ]= 1.5
mpl.rcParams['xtick.minor.visible'  ]= True
mpl.rcParams['ytick.minor.visible'  ]= True
#mpl.rcParams['font.size'            ]= 20.0
mpl.rcParams['axes.linewidth'       ]= 2.0 
mpl.rcParams['axes.labelsize'       ]= 12.0
mpl.rcParams['xtick.labelsize'      ]= 12.0
mpl.rcParams['ytick.labelsize'      ]= 12.0
mpl.rcParams['axes.axisbelow'       ]= False 
mpl.rcParams['xtick.top'            ]= True
mpl.rcParams['xtick.major.size'     ]= 8.0
mpl.rcParams['xtick.minor.size'     ]= 4.0
mpl.rcParams['xtick.major.width'    ]= 1.2
mpl.rcParams['xtick.minor.width'    ]= 0.6
mpl.rcParams['xtick.major.pad'      ]= 8.0
mpl.rcParams['xtick.minor.pad'      ]= 8.0
mpl.rcParams['xtick.direction'      ]= 'in'

mpl.rcParams['ytick.right'          ]= True 
mpl.rcParams['ytick.major.size'     ]= 8.0  
mpl.rcParams['ytick.minor.size'     ]= 4.0  
mpl.rcParams['ytick.major.width'    ]= 1.2  
mpl.rcParams['ytick.minor.width'    ]= 0.6  
mpl.rcParams['ytick.major.pad'      ]= 5.0  
mpl.rcParams['ytick.minor.pad'      ]= 5.0  
mpl.rcParams['ytick.direction'      ]= 'in'   

mpl.rcParams['legend.framealpha'    ]= 0.7  
#mpl.rcParams['legend.edgecolor'     ]= 0.6 
mpl.rcParams['legend.borderpad'     ]= 0.25 
mpl.rcParams['legend.labelspacing'  ]= 0.25 
mpl.rcParams['legend.handlelength'  ]= 1.0  
mpl.rcParams['legend.handletextpad' ]= 0.3  
mpl.rcParams['legend.borderaxespad' ]= 1.1  
#mpl.rcParams['figure.subplot.left'  ]= 0.185  
#mpl.rcParams['figure.subplot.right' ]= 0.955  
#mpl.rcParams['figure.subplot.bottom']= 0.17   
#mpl.rcParams['figure.subplot.top'   ]= 0.965  
#mpl.rcParams['figure.subplot.wspace']= 0.03   
#mpl.rcParams['figure.subplot.hspace']= 0.03     

