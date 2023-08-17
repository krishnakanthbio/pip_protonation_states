#	evaluating del G as a function of pH for PIP3 mutant system
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 11 July 2023

#  Load the required packages

import numpy as np
import matplotlib.pylab as plt
import matplotlib.colors as mcolors
plt.rcParams['font.family'] = 'comfortaa'

pip3_labels = ["3'4'5'","3'4'5", "3'45'","34'5'",  "3'45", "34'5", "345'", "345"]
pka_mean = np.array([6.48, 6.10, 6.52, 7.80, 7.11, 8.19, 8.08, 7.08, 7.67, 9.14, 9.83, 9.24])
pka_sd = np.array([0.08, 0.07,0.07, 0.09, 0.05, 0.1, 0.06, 0.1, 0.1, 0.1, 0.1, 0.2])

my_colors = list(mcolors.TABLEAU_COLORS.values())

def my_ph(x,pka_array):
	pka1, pka2, pka3, pka7, pka8, pka9, pka10, pka11, pka12, pka16, pka17, pka18 = pka_array 
	conc_3 = 10**(x-pka1)
	conc_4 = 10**(x-pka2)
	conc_5 = 10**(x-pka3)
	conc_34 = conc_3*10**(x-pka7) + conc_4*10**(x-pka9)
	conc_35 = conc_3*10**(x-pka8) + conc_5*10**(x-pka11)
	conc_45 = conc_4*10**(x-pka10) + conc_5*10**(x-pka12)
	conc_345 = conc_34*10**(x-pka16) + conc_35*10**(x-pka17) + conc_45*10**(x-pka18)
	f_sum = 1+conc_3+conc_4+conc_5+conc_34+conc_35+conc_45+conc_345
	f = [conc_345/f_sum, conc_34/f_sum,conc_35/f_sum,conc_45/f_sum,conc_3/f_sum, conc_4/f_sum,conc_5/f_sum ]
	f = f+[1-sum(f)]
	return (np.array(f))



req_systems = ["000","001","010","100","011","101","110","111"] 

sim_del_gs_wt = []
sim_del_gs_wt_sd = []
for i in req_systems:
	us_data = np.loadtxt("PMF-plots-data-PIP3/"+i+".txt",skiprows = 17)
	us_data[:,1] = us_data[:,1]*1e3 # kJ to J
	my_y = us_data [:,1] - np.mean(us_data[-10:-1,1])
	sim_del_gs_wt.append(-1.0*np.min(my_y))
	sim_del_gs_wt_sd.append(us_data[np.where(my_y==np.min(my_y))[0][0],2]*1e3)
	
sim_del_gs_mutant = []
sim_del_gs_mutant_sd = []
for i in req_systems:
	us_data = np.loadtxt("PMF-plots-data-PIP3/"+i+"M.txt")#,skiprows = 17)
	us_data[:,1] = us_data[:,1]*1e3 # kJ to J
	my_y = us_data [:,1] - np.mean(us_data[-10:-1,1])
	sim_del_gs_mutant.append(-1.0*np.min(my_y))
	sim_del_gs_mutant_sd.append(us_data[np.where(my_y==np.min(my_y))[0][0],2]*1e3)	

sim_del_gs_wt = np.round(np.array(sim_del_gs_wt),0)
sim_del_gs_wt_sd = np.round(np.array(sim_del_gs_wt_sd),0)
sim_del_gs_mutant = np.round(np.array(sim_del_gs_mutant),0)
sim_del_gs_mutant_sd = np.round(np.array(sim_del_gs_mutant_sd),0)


ph_range = np.linspace(4,10, 3000)
f = np.zeros([3000,8])
f_u = np.zeros([3000,8])	
f_l = np.zeros([3000,8])
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


for j in range(8):
	plt.plot(ph_range,f[:,j], label = pip3_labels[j], color = my_colors[j])
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



plt.plot(ph_range, np.array(ph_del_g_wt)/4184, label ="PIP$_3$ - Wild Type",color ="C2")
plt.fill_between(ph_range,(np.array(ph_del_g_wt)+np.array(ph_del_g_wt_sd))/4184,(np.array(ph_del_g_wt)-np.array(ph_del_g_wt_sd))/4184,  color = "C2", alpha=0.2)
plt.plot(ph_range, np.array(ph_del_g_mutant)/4184, label ="PIP$_3$ - Mutant",color ="C3")
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

plt.savefig("pip3_mutant_wild_type_delg_ph.pdf")
