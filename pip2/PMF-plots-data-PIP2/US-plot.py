import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import matplotlib as mpl

import matplotlib
import math
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from scipy.stats import kde

matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams['axes.linewidth'] = 3.5
matplotlib.rcParams['xtick.major.width'] = 2.5
matplotlib.rcParams['ytick.major.width'] = 2.5
matplotlib.rcParams['xtick.minor.width'] = 2.0
import matplotlib.pyplot as plt
font = {'family' : 'serif',
        'serif' : ['Liberation Serif'],
        'size'   : 30}
matplotlib.rc('font', **font)
import numpy as np

plt.figure(figsize=(10,7))

rescolors = ['orange', 'blue']
req_systems = ["00","00M"]

counter = 0
for i in req_systems:
    us_data = np.loadtxt(i+"-US-100BS")
    us_data[:,1] = us_data[:,1]/ 4.184 # kcal correction on profile
    us_data[:,2] = us_data[:,2]/ 4.184 # kcal correction on errors
    plt.plot(us_data[:,0],us_data[:,1],color=rescolors[counter], linewidth =7, label ='PI2'+i)
    plt.fill_between(us_data[:,0],us_data[:,1]-us_data[:,2],us_data[:,1]+us_data[:,2],color =rescolors[counter],alpha = 0.5)
    counter = counter + 1

plt.xticks([3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5])
plt.yticks([0,-4,-8,-12,-16,-20])

plt.xlabel('CV (nm) ')
plt.ylabel('PMF (kcal/mol)')
plt.legend(loc= "lower right", prop={'size': 25})
plt.savefig("pmf-plot-00.png",dpi  = 300)
