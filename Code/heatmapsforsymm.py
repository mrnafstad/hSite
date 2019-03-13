"""
cd c:\Users\Halvor\Documents\Master\Kode\
python heatmapsforsymm.py

"""
#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys
from matplotlib import ticker, cm

def placementmass( massarr, farr ):

	maxima = []
	di = []
	rvir = []
	for f_R0 in farr:
		d_vir = []
		d_ta = []
		d_c = []
		rvirrta = []
		avirata = []
		variable = []
		avir = []
		acoll = []
		delta_c = []
		for M in massarr:
			filename = "Numbers\Chameleon\M%.15ff_R0%.15f.txt" % (np.log10(M), np.log10(abs(1-f_R0)))
			file = open(filename, "r")

			i = 0

			for j in file.readlines():
				try:
					value = float(j)
				except ValueError:
					break
				if i == 0:
					variable.append(value)
					i += 1
				elif i == 1:
					d_vir.append(value)
					i += 1
				elif i == 2:
					d_ta.append(value)
					i += 1
				elif i == 3:
					d_c.append(value)
					i += 1

				elif i == 4:
					rvirrta.append(value)
					i += 1
				elif i == 5:
					avirata.append(value)
					i += 1
				elif i == 6:
					avir.append(value)
					i += 1

				elif ( i == 7 ):
					acoll.append(value)
					i += 1

				elif i == 8:
					delta_c.append(value)
					i = 0
					

			file.close()

		maxima.append(np.argmax(rvirrta))
		di.append(delta_c)
		rvir.append(np.argmin(d_vir))

	actmax = np.argmax(rvir)
	minf = farr[actmax]
	minmass = massarr[maxima[actmax]]
	delta_is = di[actmax]

	return minmass, minf, delta_is, actmax, maxima



masses = [1e14] #np.logspace(10, 18, 30)
beta_list = [0.5] #np.flip(np.linspace(0, 1.5, 30), 0)
L_list = np.logspace(-1, 1, 30)
z_ssb_list = np.flip(np.linspace(0, 3, 30), 0)

X, Y = np.meshgrid(L_list, z_ssb_list)

d_vir = np.zeros(((len(L_list)), len(L_list)))
d_ta = np.zeros(((len(L_list)), len(L_list)))
d_c = np.zeros(((len(L_list)), len(L_list)))
rvirrta = np.zeros(((len(L_list)), len(L_list)))
avirata = np.zeros(((len(L_list)), len(L_list)))
#variable = np.zeros((range(len(masses))), range(len(f_R0_list)))
avir = np.zeros(((len(L_list)), len(L_list)))
acoll = np.zeros(((len(L_list)), len(L_list)))
delta_c = np.zeros(((len(L_list)), len(L_list)))

"""
Mlocalmin, f_R0, delta_is, minima, massloc = placementmass(masses, f_R0_list)

print massloc[minima], minima
print masses[massloc[minima]], f_R0_list[minima]
print Mlocalmin, f_R0
minmin = minima 
minplus = minima 
Mcomp1 = masses[massloc[minima] - 3]
f_R0_comp1 = f_R0_list[minplus]
Mcomp2 = masses[massloc[minima] + 3]
f_R0_comp2 = f_R0_list[minmin]
"""

k = 0
for M in masses:
	for beta in beta_list:	
		for L in L_list:
	
			t = 0
			for z_ssb in z_ssb_list:		
		
				filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(M), z_ssb, L, beta)
			
				file = open(filename, "r")

				i = 0

				for j in file.readlines():
					try:
						value = float(j)
					except ValueError:
						break
					if i == 0:
						#variable[k,t] = value
						i += 1
					elif i == 1:
						d_vir[k,t] = 1 - value/295.618133959
						i += 1
					elif i == 2:
						d_ta[k,t] = 1 - value/7.098509211
						i += 1
					elif i == 3:
						d_c[k,t] = 1 - value/381.001857758968
						i += 1
					elif i == 4:
						rvirrta[k,t] = 1 - value/0.481722136884160
						i += 1
					elif i == 5:
						avirata[k,t] = 1 - value/1.6697492018
						i += 1
					elif i == 6:
						avir[k,t] = 1 - value/0.9188909889262
						i += 1
					elif ( i == 7 ):
						acoll[k,t] = 1 - value/0.9999907897
						i += 1
					elif i == 8:
						delta_c[k,t] = 1 - value/0.0005904
						i = 0				

				file.close()
				t += 1
			k += 1


lev = np.linspace(np.min(delta_c), np.max(delta_c), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(delta_c), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
mpl.xscale("log")
#mpl.yscale("log")
#mpl.xlabel(r"$\mathrm{h}^{-1}$$\mathrm{M}_{\odot}$ ")
mpl.xlabel(r"$\mathrm{L}$ $\mathrm{h}^{-1}\mathrm{Mpc}$", rotation = 0, labelpad = 20)
#mpl.ylabel(r"$\beta$", rotation = 0, labelpad = 10)
mpl.ylabel(r"$z_{\mathrm{ssb}}$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\delta_{\mathrm{i, symm}}}{\delta_{\mathrm{i}, \Lambda\mathrm{CDM}}}$")
#mpl.savefig("Figures\Heatmaps\Symm_z_delta_i.pdf", bbox_inches = "tight")
mpl.savefig("Figures\Heatmaps\Symm_L_z_delta_i.pdf", bbox_inches = "tight")
mpl.clf()


lev = np.linspace(np.min(d_ta), np.max(d_ta), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(d_ta), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
mpl.xscale("log")
#mpl.xlabel(r"$\mathrm{h}^{-1}$$\mathrm{M}_{\odot}$ ")
mpl.xlabel(r"$\mathrm{L}$ $\mathrm{h}^{-1}\mathrm{Mpc}$", rotation = 0, labelpad = 20)
#mpl.ylabel(r"$\beta$", rotation = 0, labelpad = 10)
#mpl.yscale("log")
mpl.ylabel(r"$z_{\mathrm{ssb}}$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\Delta_{\mathrm{ta, symm}}}{\Delta_{\mathrm{ta}, \Lambda\mathrm{CDM}}}$")
#mpl.savefig("Figures\Heatmaps\Symm_z_Delta_ta.pdf", bbox_inches = "tight")
mpl.savefig("Figures\Heatmaps\Symm_L_z_Delta_ta.pdf", bbox_inches = "tight")
mpl.clf()



lev = np.linspace(np.min(d_vir), np.max(d_vir), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(d_vir), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
mpl.xscale("log")
#mpl.ylabel(r"$\beta$", rotation = 0, labelpad = 10)
mpl.xlabel(r"$\mathrm{L}$ $\mathrm{h}^{-1}\mathrm{Mpc}$", rotation = 0, labelpad = 20)
#mpl.xlabel(r"$\mathrm{h}^{-1}$$\mathrm{M}_{\odot}$ ")
#mpl.yscale("log")
mpl.ylabel(r"$z_{\mathrm{ssb}}$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\Delta_{\mathrm{vir, symm}}}{\Delta_{\mathrm{vir}, \Lambda\mathrm{CDM}}}$")
#mpl.savefig("Figures\Heatmaps\Symm_z_Delta_vir.pdf", bbox_inches = "tight")
mpl.savefig("Figures\Heatmaps\Symm_L_z_Delta_vir.pdf", bbox_inches = "tight")
mpl.clf()



lev = np.linspace(np.min(d_c), np.max(d_c), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(d_c), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
mpl.xscale("log")
#mpl.yscale("log")
#mpl.xlabel(r"$\mathrm{h}^{-1}$$\mathrm{M}_{\odot}$ ")
mpl.xlabel(r"$\mathrm{L}$ $\mathrm{h}^{-1}\mathrm{Mpc}$", rotation = 0, labelpad = 20)
#mpl.ylabel(r"$\beta$", rotation = 0, labelpad = 10)
mpl.ylabel(r"$z_{\mathrm{ssb}}$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\Delta_{\mathrm{c, symm}}}{\Delta_{\mathrm{c}, \Lambda\mathrm{CDM}}}$")
#mpl.savefig("Figures\Heatmaps\Symm_z_Delta_c.pdf", bbox_inches = "tight")
mpl.savefig("Figures\Heatmaps\Symm_L_z_Delta_c.pdf", bbox_inches = "tight")
mpl.clf()
#

lev = np.linspace(np.min(rvirrta), np.max(rvirrta), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(rvirrta), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
#mpl.axhline(y = abs(1 - f_R0_list[minima]), color = "r", linestyle = "--")
#mpl.axvline(x = Mlocalmin, color = "r", linestyle = "--")
mpl.xscale("log")
#mpl.yscale("log")
#mpl.xlabel(r"$\mathrm{h}^{-1} $$\mathrm{M}_{\odot}$ ")
mpl.xlabel(r"$\mathrm{L}$ $\mathrm{h}^{-1}\mathrm{Mpc}$", rotation = 0, labelpad = 20)
#mpl.ylabel(r"$\beta$", rotation = 0, labelpad = 10)
mpl.ylabel(r"$z_{\mathrm{ssb}}$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\tilde{R}_{\mathrm{vir,symm}}/\tilde{R}_{\mathrm{ta,symm}}}{\tilde{R}_{\mathrm{vir},\Lambda\mathrm{CDM}}/\tilde{R}_{\mathrm{ta},\Lambda \mathrm{CDM}}}$")
#mpl.savefig("Figures\Heatmaps\Symm_z_rvirrta.pdf", bbox_inches = "tight")
mpl.savefig("Figures\Heatmaps\Symm_L_z_rvirrta.pdf", bbox_inches = "tight")
mpl.clf()


"""
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, avirata)
cbar = fig.colorbar(cs)
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel(r"$M_{\odot} h^{-1}$")
mpl.suptitle(r"$\frac{a_{vir}}{a_{ta}}_{\Lambda \mathrm{CDM}} - \frac{a_{vir}}{a_{ta}}_{cham}$")
mpl.savefig("Figures\Heatmaps\Cham_avirata.pdf", bbox_inches = "tight")
mpl.clf()



fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, avir)
cbar = fig.colorbar(cs)
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel(r"$M_{\odot} h^{-1}$")
mpl.ylabel(r"$|f_{\mathcal{R}_0} - 1|$", rotation = 0, labelpad = 20)
mpl.title(r"$a_{vir,\Lambda \mathrm{CDM}} - a_{vir, cham}$")
mpl.savefig("Figures\Heatmaps\Cham_avir.pdf", bbox_inches = "tight")
mpl.clf()

"""

fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(acoll))
cbar = fig.colorbar(cs)
mpl.xscale("log")
#mpl.yscale("log")
#mpl.xlabel(r"$\mathrm{h}^{-1}$$\mathrm{M}_{\odot}$ ")
mpl.xlabel(r"$\mathrm{L}$ $\mathrm{h}^{-1}\mathrm{Mpc}$", rotation = 0, labelpad = 20)
#mpl.ylabel(r"$\beta$", rotation = 0, labelpad = 10)
mpl.ylabel(r"$z_{\mathrm{ssb}}$", rotation = 0, labelpad = 20)
mpl.title(r"$1 - a_{\mathrm{c, symm}}$")
#mpl.savefig("Figures\Heatmaps\Symm_z_acoll.pdf", bbox_inches = "tight")
mpl.savefig("Figures\Heatmaps\Symm_L_z_acoll.pdf", bbox_inches = "tight")
mpl.clf()
