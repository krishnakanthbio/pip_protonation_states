#	evaluating del G as a function of pH for PIP2 mutant system
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 5 July 2023

#  Load the required packages

import numpy as np
import matplotlib.pylab as plt
#from pymaxent import *
import matplotlib.colors as mcolors
plt.rcParams['font.family'] = 'comfortaa'

pip2_labels = ["4'5'","4'5","45'","45"]
pka_mean = np.array([6.25, 6.44, 7.60, 7.41])
pka_sd = np.array([0.04, 0.05,0.04, 0.04])
my_colors = list(mcolors.TABLEAU_COLORS.values())

def my_ph(x,pka_array):
	pka1, pka2, pka4,pka5 = pka_array
	conc_4 = 10**(x-pka1)
	conc_5 = 10**(x-pka2)
	conc_45 = conc_4*10**(x-pka4) + conc_5*10**(x-pka5)
	f_sum = 1+conc_4+conc_5+conc_45
	f = [conc_45/f_sum,conc_4/f_sum,conc_5/f_sum]
	f = f+[1-sum(f)]
	return(np.array(f))



req_systems = ["00","01","10","11"] 

sim_del_gs_wt = []
sim_del_gs_wt_sd = []
for i in req_systems:
	us_data = np.loadtxt("PMF-plots-data-PIP2/"+i+"-US-100BS")#,skiprows = 17)
	us_data[:,1] = us_data[:,1]*1e3 # kJ to J
	my_y = us_data [:,1] - np.mean(us_data[-10:-1,1])
	sim_del_gs_wt.append(-1.0*np.min(my_y))
	sim_del_gs_wt_sd.append(us_data[np.where(my_y==np.min(my_y))[0][0],2]*1e3)
	
sim_del_gs_mutant = []
sim_del_gs_mutant_sd = []
for i in req_systems:
	us_data = np.loadtxt("PMF-plots-data-PIP2/"+i+"M-US-100BS")#,skiprows = 17)
	us_data[:,1] = us_data[:,1]*1e3 # kJ to J
	my_y = us_data [:,1] - np.mean(us_data[-10:-1,1])
	sim_del_gs_mutant.append(-1.0*np.min(my_y))
	sim_del_gs_mutant_sd.append(us_data[np.where(my_y==np.min(my_y))[0][0],2]*1e3)	

sim_del_gs_wt = np.round(np.array(sim_del_gs_wt),0)
sim_del_gs_wt_sd = np.round(np.array(sim_del_gs_wt_sd),0)
sim_del_gs_mutant = np.round(np.array(sim_del_gs_mutant),0)
sim_del_gs_mutant_sd = np.round(np.array(sim_del_gs_mutant_sd),0)


ph_range = np.linspace(4,10, 3000)
f = np.zeros([3000,4])
f_u = np.zeros([3000,4])	
f_l = np.zeros([3000,4])
ph_del_g_wt = []
ph_del_g_wt_sd = []
ph_del_g_mutant = []
ph_del_g_mutant_sd = []

for i in range(3000):
	f[i,:] = np.round(my_ph(ph_range[i], pka_mean),3)
	f_u[i,:] = my_ph(ph_range[i], pka_mean+pka_sd)
	f_l[i,:] = my_ph(ph_range[i], pka_mean-pka_sd)
	f_sd = np.max(np.row_stack([f[i,:],f_l[i,:],f_u[i,:]]),axis=0)-f[i,:]
	f_sd = np.round(f_sd,3)
	ph_del_g_wt.append(np.sum(np.round(f[i,:]*sim_del_gs_wt,0)))
	ph_del_g_wt_sd.append(np.sum(np.round(sim_del_gs_wt_sd*np.round(f[i,:],3) + sim_del_gs_wt*np.round(f_sd,3),0)))	
	ph_del_g_mutant.append(np.sum(np.round(f[i,:]*sim_del_gs_mutant,0)))
	ph_del_g_mutant_sd.append(np.sum(np.round(sim_del_gs_mutant_sd*np.round(f[i,:],3) + sim_del_gs_mutant*np.round(f_sd,3),0)))		


for j in range(4):
	plt.plot(ph_range,f[:,j], label = pip2_labels[j], color = my_colors[j])
	plt.fill_between(ph_range,f_l[:,j],f_u[:,j],  color = my_colors[j], alpha=0.5)

ax = plt.gca()
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
plt.xlabel("pH")
plt.ylabel("Ionization State Fraction")
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.tight_layout()
plt.savefig("fractions.pdf")
plt.close()



plt.plot(ph_range, np.array(ph_del_g_wt)/4184, label ="PIP$_2$ - Wild Type",color ="C2")
plt.fill_between(ph_range,(np.array(ph_del_g_wt)+np.array(ph_del_g_wt_sd))/4184,(np.array(ph_del_g_wt)-np.array(ph_del_g_wt_sd))/4184,  color = "C2", alpha=0.2)
plt.plot(ph_range, np.array(ph_del_g_mutant)/4184, label ="PIP$_2$ - Mutant",color ="C3")
plt.fill_between(ph_range,(np.array(ph_del_g_mutant)+np.array(ph_del_g_mutant_sd))/4184,(np.array(ph_del_g_mutant)-np.array(ph_del_g_mutant_sd))/4184,  color = "C3", alpha=0.2)


ax = plt.gca()
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
plt.xlabel("pH")
plt.ylabel("$\Delta G^{sim}_{binding}$"+" "+"$(kcal\:mol^{-1})$")
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

plt.tight_layout()

plt.savefig("pip2_mutant_wild_type_delg_ph.pdf")
