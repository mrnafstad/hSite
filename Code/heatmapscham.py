"""
cd c:\Users\Halvor\Documents\Master\Kode\
python heatmapproduction.py

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
		rvir.append(rvirrta[np.argmax(rvirrta)])

	actmax = np.argmax(rvir)
	minf = farr[maxima[actmax]]
	minmass = massarr[actmax]
	delta_is = di[actmax]

	return minmass, minf, delta_is, actmax, maxima



masses = np.logspace(10, 18, 30)
f_R0_list = np.flip(np.logspace(-12, -3, 30) + 1, 0)

X, Y = np.meshgrid(masses, abs(f_R0_list - 1))

d_vir = np.zeros(((len(masses)), len(f_R0_list)))
d_ta = np.zeros(((len(masses)), len(f_R0_list)))
d_c = np.zeros(((len(masses)), len(f_R0_list)))
rvirrta = np.zeros(((len(masses)), len(f_R0_list)))
avirata = np.zeros(((len(masses)), len(f_R0_list)))
#variable = np.zeros((range(len(masses))), range(len(f_R0_list)))
avir = np.zeros(((len(masses)), len(f_R0_list)))
acoll = np.zeros(((len(masses)), len(f_R0_list)))
delta_c = np.zeros(((len(masses)), len(f_R0_list)))

k = 0

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

for M in masses:
	t = 0
	for f_R0 in f_R0_list:
		
		filename = "Numbers\Chameleon\M%.15ff_R0%.15f.txt" % (np.log10(M), np.log10(abs(1-f_R0)))
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
mpl.yscale("log")
mpl.xlabel(r"$ \mathrm{h}^{-1}\mathrm{M}_{\odot}$")
mpl.ylabel(r"$|f_{\mathcal{R}_0} - 1|$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\delta_{\mathrm{i, cham}}}{\delta_{\mathrm{i}, \Lambda\mathrm{CDM}}}$")
mpl.savefig("Figures\Heatmaps\Cham_delta_i.pdf", bbox_inches = "tight")
mpl.clf()


lev = np.linspace(np.min(d_ta), np.max(d_ta), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(d_ta), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel(r"$\mathrm{h}^{-1}\mathrm{M}_{\odot} $")
mpl.ylabel(r"$|f_{\mathcal{R}_0} - 1|$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\Delta_{\mathrm{ta, cham}}}{\Delta_{\mathrm{ta}, \Lambda\mathrm{CDM}}}$")
mpl.savefig("Figures\Heatmaps\Cham_Delta_ta.pdf", bbox_inches = "tight")
mpl.clf()



lev = np.linspace(np.min(d_vir), np.max(d_vir), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(d_vir), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel(r"$\mathrm{h}^{-1}\mathrm{M}_{\odot} $")
mpl.ylabel(r"$|f_{\mathcal{R}_0} - 1|$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\Delta_{\mathrm{vir, cham}}}{\Delta_{\mathrm{vir}, \Lambda\mathrm{CDM}}}$")
mpl.savefig("Figures\Heatmaps\Cham_Delta_vir.pdf", bbox_inches = "tight")
mpl.clf()



lev = np.linspace(np.min(d_c), np.max(d_c), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(d_c), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[massloc[minima]]], [abs(f_R0_list[minima] - 1)],'ro', mfc='none', label = "Global extrema")
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel(r"$\mathrm{h}^{-1}\mathrm{M}_{\odot} $")
mpl.ylabel(r"$|f_{\mathcal{R}_0} - 1|$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\Delta_{\mathrm{c, cham}}}{\Delta_{\mathrm{c}, \Lambda\mathrm{CDM}}}$")
mpl.savefig("Figures\Heatmaps\Cham_Delta_c.pdf", bbox_inches = "tight")
mpl.clf()
#

lev = np.linspace(np.min(rvirrta), np.max(rvirrta), 20)
fig, ax = mpl.subplots()
cs = ax.contourf(X, Y, np.transpose(rvirrta), lev)
cbar = fig.colorbar(cs)
#mpl.plot([masses[minima]], [abs(f_R0_list[massloc[minima]] - 1)],'ro', mfc='none', label = "Global extrema")
#mpl.axhline(y = abs(1 - f_R0_list[minima]), color = "r", linestyle = "--")
#mpl.axvline(x = Mlocalmin, color = "r", linestyle = "--")
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel(r"$\mathrm{h}^{-1}\mathrm{M}_{\odot} $")
mpl.ylabel(r"$|f_{\mathcal{R}_0} - 1|$", rotation = 0, labelpad = 20)
mpl.suptitle(r"$1 - \frac{\tilde{R}_{\mathrm{vir,cham}}/\tilde{R}_{\mathrm{ta,cham}}}{\tilde{R}_{\mathrm{vir},\Lambda\mathrm{CDM}}/\tilde{R}_{\mathrm{ta},\Lambda \mathrm{CDM}}}$")
mpl.savefig("Figures\Heatmaps\Cham_rvirrta.pdf", bbox_inches = "tight")
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
mpl.yscale("log")
mpl.xlabel(r"$\mathrm{h}^{-1}\mathrm{M}_{\odot} $")
mpl.ylabel(r"$|f_{\mathcal{R}_0} - 1|$", rotation = 0, labelpad = 20)
mpl.title(r"$1 - a_{\mathrm{c, cham}}$")
mpl.savefig("Figures\Heatmaps\Cham_acoll.pdf", bbox_inches = "tight")
mpl.clf()
