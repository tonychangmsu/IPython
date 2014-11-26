#Title: PFT_ACI.py
#Original Author: Benjamin Poulter 
#Translation Author: Tony Chang
#Date: 11/20/2014
#Abstract: Solves the A-Ci curve for 3 different limiting factors, (1) RuBP limiting zone, 
#(2) RuBP regen limited, (3) TPU limited. Requires the user to determine the points where
# the plant switches to different limiting agents.
#
#Dependencies: python 3.x, numpy, pandas, scipy

#March 2011
#Follows Sharkey et al.
#Slightly different results from Excel Solver - RMSE is lower with solver than with R

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
  
def Ac_func(x, Vcmax, Rd):
	Ci_pa_rubp = x[0]
	Roa = x[1]
	Kc = x[2]
	O = x[3]
	Ko = x[4]
	return(Vcmax*((Ci_pa_rubp-Roa)/(Ci_pa_rubp+Kc*(1+O/Ko)))-Rd)

def Aj_func(x, J):
	Ci_pa = x[0]
	Roa = x[1]
	Rd = x[2]
	return((J*((Ci_pa-Roa)/(4*Ci_pa + 8 *Roa))-Rd))
	
def At_func(x, TPU):
	Rd = x[0]
	return(3*TPU-Rd)

def A_curves(A, Ci_ppm, Ci_thresh, T= 28, Patm = 101, O = 21):
#function for building all the parameters for the A-Ci curves under the different ranges
#input requires A and Ci values, as well as distinction of where the initial Ci_threshhold value is for Pa
#optional inputs for Patm (kPa), T (C), and O (kPa)
#output the parameters for each of the Michelis-Menton curves

	#Michelis menton parameters
	#Ci_thresh = 30		#30 kPa play with this term to get the best fit....
	#T = 28			#C
	#Could have as variables
	#Patm = 101			#kPa
	#O = 21 			#kPa

	#Michelis menton parameters
	R = 0.008314		#
	Kc_c = 35.9774		#
	Ko_c = 12.3772		#
	Kroa_c = 11.187		#
	Vcmax_c = 26.355	#
	Rd_c = 18.715		#

	Kc_delH = 80.99		#
	Ko_delH = 23.72		#
	Kroa_delH = 24.46	#
	Vcmax_delH = 65.33	#
	Rd_delH = 46.39		#

	#
	Ci_pa = Ci_ppm*Patm*0.001
	Ci_pa_rubp = Ci_pa[Ci_pa < Ci_thresh]
	A_rubp = A[Ci_pa < Ci_thresh]

	#
	Kc = np.exp(Kc_c-Kc_delH/(R*(T+273.15)))
	Ko = np.exp(Ko_c-Ko_delH/(R*(T+273.15)))
	Roa = np.exp(Kroa_c-Kroa_delH/(R*(T+273.15)))
	Vcmax_adj = np.exp(Vcmax_c-Vcmax_delH/(R*(T+273.15)))
	Rd_adj = np.exp(Rd_c-Rd_delH/(R*(T+273.15)))

	#
	xdata = np.array([Ci_pa_rubp, Roa, Kc, Ko, O])
	ydata = A_rubp
	pdata = (0,0) #initialize Vcmax = 0 and Rd = 0

	Vcmax_mod = curve_fit(Ac_func, xdata, ydata, p0=pdata)
	Vcmax = Vcmax_mod[0][0]
	Rd = Vcmax_mod[0][1]
	Vcmax25 = Vcmax / Vcmax_adj
	Rd25 = Rd / Rd_adj

	#new xdata for the Aj fit for points beyond pressure threshold
	Ci_pa_ps = Ci_pa[Ci_pa > Ci_thresh]
	A_ps = A[Ci_pa > Ci_thresh]
	xdata_J = np.array([Ci_pa_ps, Roa, Rd])
	ydata_J = A_ps
	J_mod = curve_fit(Aj_func, xdata_J, ydata_J, p0 = 0)
	J = J_mod[0][0]

	#now build a function to determine the TPU limited portion and optimize it.
	TPU_thresh = 60
	Ci_pa_ts = Ci_pa[Ci_pa > TPU_thresh]
	A_ts = A[Ci_pa > TPU_thresh]
	xdata_TPU = np.array([Rd])
	ydata_ts = A_ts
	TPU_mod = curve_fit(At_func, xdata_TPU, ydata_ts, p0 = 0)
	TPU = TPU_mod[0][0]

	xy = pd.DataFrame(np.array([Ci_pa, O, Vcmax, Vcmax25, Rd, Rd25, Kc, Ko, Roa, J, TPU]))
	out = xy.T
	out.columns = ['Ci_pa', 'O', 'VcmaxT', 'Vcmax25', 'Rd', 'Rd25', 'Kc', 'Ko', 'Roa', 'J', 'TPU']
	return(out)

#=============================MAIN==============================
#test dataset
#A =  np.array([-3.27,3.49,10.8,17.8,23.5,27.6,30.5,31.7,33,33.3,33.3])
#Ci_ppm = np.array([20.7,77.5,135,197,266,344,428,517,606,698,791])

#try for a different dataset?
leafdata = pd.read_csv('D:\\CHANG\\PhD_Material\\MSU_Course_Material\\Fall_2014\\BIOE_532\\Project\\CHANG_PILLET_Datasheet_of_Leaf_Area.csv')
picoLA = leafdata['Leaf Area (m^2)'][leafdata.Sample =='PICOC4']
pialLA = leafdata['Leaf Area (m^2)'][leafdata.Sample =='PIALC1']
la = [picoLA.iloc[0], pialLA.iloc[0]]
data = pd.read_csv('D:\\CHANG\\PhD_Material\\MSU_Course_Material\\Fall_2014\\BIOE_532\\Project\\CHANG_PILLET_Experiment_data_summary_11212014.csv')

pico = data[data.Sample=='PICOC4']
pial = data[data.Sample=='PIALC1']

species =[pico, pial]
treatments = ['7C', '17C', '27C', 'Amb']
colors = ['green', 'purple', 'red', 'blue']
species_names = ['PICO', 'PIAL']

plt.rcParams['figure.figsize'] = 20,16
for i in range(len(species)):
	for t in range(len(treatments)-1): #do not include ambient case
		ax = plt.subplot2grid((len(species),(len(treatments))-1),(i,t))
		A = species[i].Photo[species[i].Treatment==treatments[t]] #assimilation value
		Ci_ppm = species[i].CO2S[species[i].Treatment==treatments[t]] #C_i values

		ind = np.array(np.argsort(Ci_ppm))
		A = A.iloc[ind]
		Ci_ppm = Ci_ppm.iloc[ind]

		Ci_thresh = 100	
		out = A_curves(A, Ci_ppm, Ci_thresh)
		Vcmax = out.VcmaxT[0]
		Rd = out.Rd[0]
		Ci_pa = out.Ci_pa[0]
		O =out.O[0]
		Roa = out.Roa[0]
		Kc = out.Kc[0]
		Ko = out.Ko[0]
		J = out.J[0]
		TPU = out.TPU[0]

		Ac =  (Vcmax*((Ci_pa-Roa)/(Ci_pa+Kc*(1+O/Ko)))-Rd) * la[i] 
		Aj = (J*((Ci_pa-Roa)/(4*Ci_pa + 8 *Roa))-Rd) *la[i]
		At = (np.ones(len(Ci_pa))*3* TPU - Rd) * la[i]
		A = A * la[i]
		plt.scatter(Ci_pa, A, label ='A obs')
		plt.plot(Ci_pa, Ac, color = 'red', label = 'RuBisCo')
		#plt.plot(Ci_pa, Aj, color = 'green', label = 'RuBP Regen') 
		#plt.plot(Ci_pa, At, color = 'purple', label = 'TPU')
		plt.grid()
		plt.title('%s %s'%(species_names[i], treatments[t]))
		plt.legend(loc = 'lower right')
		plt.xticks(np.arange(0,180,20),np.arange(0,1800,200))
		plt.xlabel("$C_i$ ($CO_2$ $Ppm$)")
		plt.ylabel("$A$ ($\mu mol$ $s^{-1}$)")

#things to do (1. remove outlier data which is faulty, 2.diagnose where the pressure threshold zone is for each species, 3. normalize all values by their leaf area)
#compare the absolute values...
'''
#plt.savefig('D:\\chang\\phd_material\\msu_course_material\\fall_2014\\bioe_532\\project\\aci_plot_species.png', bbox_inches='tight')
#solve for the full curve?
plt.rcParams['figure.figsize'] = 20,14
treatments = ['Amb', '7C', '17C', '27C']
colors= ['blue', 'red']
m = ['x','o']
Vcmax_array = [[],[]]
J_array = [[],[]]
Rd_array= [[],[]]
for t in range(len(treatments)): #do not include ambient case
	#ax = plt.subplot2grid((len(species)-1,(len(treatments))-1),(0,t))
	ax1 = plt.subplot(2,2,t+1)
	ax2= ax1.twinx()
	for i in range(len(species)):
		A = species[i].Photo[species[i].Treatment==treatments[t]] #assimilation value
		Ci_ppm = species[i].CO2S[species[i].Treatment==treatments[t]] #C_i values
		cond = species[i].Cond[species[i].Treatment==treatments[t]] * la[i]
		ind = np.array(np.argsort(Ci_ppm))
		A = A.iloc[ind]
		Ci_ppm = Ci_ppm.iloc[ind]

		Ci_thresh = 120	
		out = A_curves(A, Ci_ppm, Ci_thresh)
		Vcmax = out.VcmaxT[0]
		Rd = out.Rd[0]
		Ci_pa = out.Ci_pa[0]
		O =out.O[0]
		Roa = out.Roa[0]
		Kc = out.Kc[0]
		Ko = out.Ko[0]
		J = out.J[0]
		TPU = out.TPU[0]
		Vcmax_array[i].append(Vcmax)
		J_array[i].append(J)
		Rd_array[i].append(Rd)

		Ac =  (Vcmax*((Ci_pa-Roa)/(Ci_pa+Kc*(1+O/Ko)))-Rd) * la[i] 
		Aj = (J*((Ci_pa-Roa)/(4*Ci_pa + 8 *Roa))-Rd) *la[i]
		At = (np.ones(len(Ci_pa))*3* TPU - Rd) * la[i]
		A = A * la[i]
		ax1.scatter(Ci_pa, A, marker = m[i],color = colors[i]) #label ='%s  obs'%(species_names[i])
		ax1.plot(Ci_pa, Ac, color = colors[i], label = '$A_n$ $%s$'%(species_names[i]))
		
		ax2.plot(Ci_pa, cond, color = colors[i], ls = '--', label = '$g_c$ $%s$'%(species_names[i]))
		#plt.plot(Ci_pa, Aj, color = 'green', label = 'RuBP Regen') 
		#plt.plot(Ci_pa, At, color = 'purple', label = 'TPU')
		
		plt.title('Treatment: %s'%(treatments[t]), fontsize = 18)
		ax1.legend(loc = 'upper left')
		ax2.legend(loc = 'upper right')
		plt.xticks(np.arange(0,180,20),np.arange(0,1800,200))
		ax1.set_xlabel("$C_i$ ($CO_2$ $Ppm$)")
	ax1.set_ylabel("$A$ ($\mu mol$ $s^{-1}$)", fontsize = 14)
	ax2.set_ylabel('$g_c$ ($mmol$ $s^{-1}$)', fontsize = 14)
	plt.grid()
	plt.tight_layout()
plt.savefig('D:\\chang\\phd_material\\msu_course_material\\fall_2014\\bioe_532\\project\\curve_gc_plot_species.png', bbox_inches='tight')

#bar plots of Vcmax
Vcmax_array = np.array(Vcmax_array)
J_array = np.array(J_array)
Rd_array = np.array(Rd_array)
 
ind = np.arange(0,len(treatments)/2, 0.5)
width = 0.2
plt.rcParams['figure.figsize'] = 8,5
fig, ax = plt.subplots()
rects1 = ax.bar(ind, Vcmax_array[0]*la[0], width, color='b', alpha = 0.7)
rects2 = ax.bar(ind+width, Vcmax_array[1]*la[1], width, color='r', alpha = 0.7)
ax.set_ylabel('$V_{cmax}$ ($\mu mol$ $s^{-1}$)', fontsize = 18)
ax.set_title('$V_{cmax}$ by species and treatment', fontsize=15)
ax.set_xticks(ind+width)
ax.set_xticklabels( ('$Amb$', '$7^oC$', '$17^oC$', '$27^oC$'), fontsize =14)
ax.legend( (rects1[0], rects2[0]), ('PICO', 'PIAL') , loc= 'upper left')
plt.grid()
#plt.savefig('D:\\chang\\phd_material\\msu_course_material\\fall_2014\\bioe_532\\project\\vcmax_species.png', bbox_inches='tight')

#bar plots of Rd

ind = np.arange(0,len(treatments)/2, 0.5)
width = 0.2
plt.rcParams['figure.figsize'] = 8,5
fig, ax = plt.subplots()
rects1 = ax.bar(ind, Rd_array[0]*la[0], width, color='b', alpha = 0.7)
rects2 = ax.bar(ind+width, Rd_array[1]*la[1], width, color='r', alpha = 0.7)
ax.set_ylabel('$R_d$ ($\mu mol$ $s^{-1}$)', fontsize = 18)
ax.set_title('$R_d$ by species and treatment', fontsize=15)
ax.set_xticks(ind+width)
ax.set_xticklabels( ('$Amb$', '$7^oC$', '$17^oC$', '$27^oC$'), fontsize =14)
ax.legend( (rects1[0], rects2[0]), ('PICO', 'PIAL') , loc= 'upper left')
plt.grid()
#plt.savefig('D:\\chang\\phd_material\\msu_course_material\\fall_2014\\bioe_532\\project\\Rd_species.png', bbox_inches='tight')

#Finally, plotting the absolute An for the different treatments by 400 ppm and 1200pm, doesn't really make that much sense. We can see the values from the scatter plot...
'''