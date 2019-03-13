#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys
from matplotlib import ticker, cm

def plac(z_ssb, L, beta):
	d_vir = []
	d_ta = []
	d_c = []
	rvirrta = []
	avirata = []
	variable = []
	avir = []
	acoll = []
	delta_c = []
				
	for M in masses:	
		filename = "Numbers\Symmetron\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(M), z_ssb, L, beta)
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

	return d_vir, d_ta, d_c, rvirrta, avirata, avir, acoll, delta_c

masses = np.logspace(10, 18, 30)
beta_list = [0.5, 1]
L_list = [0.5, 1.0, 1.5]
z_ssb_list = [0.5, 1.0, 2, 3]
linestyles = ["-.", "--", ":"]
colors1 = ["grey", "black"]	
colors = ["yellow", "forestgreen", "limegreen", "seagreen"]
tit = ["$\Delta_{vir}", "$\Delta_{ta}", "$\Delta_{c}", r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}", r"$\frac{a_{vir}}{a_{ta}}", "$a_{vir}", "$a_c", "$\delta_i"]

y = [295.62, 7.099, 381.00, 0.4817, 1.66881, 0.916595229171737, 1.0, 0.0005904]

t = 1
for L in L_list:
	k = 0
	for z_ssb in z_ssb_list:
		small_beta_nested = plac(z_ssb, L, beta_list[0])
		big_beta_nested = plac(z_ssb, L, beta_list[1])

		for i in range(8):
			mpl.figure(i+1)
			if k == 0:
				mpl.plot(masses, small_beta_nested[i], linestyle = "-.", color = colors1[0], linewidth = 0.75, label = r"$\beta = 0.5$")
				mpl.plot(masses, big_beta_nested[i],  linestyle = "-", color = colors1[1], linewidth = 0.75, label = r"$\beta = 1.0$")
			else:
				mpl.plot(masses, small_beta_nested[i], linestyle = "-.", color = colors1[0], linewidth = 0.75)
				mpl.plot(masses, big_beta_nested[i],  linestyle = "-", color = colors1[1], linewidth = 0.75)
			mpl.fill_between(masses, big_beta_nested[i], small_beta_nested[i], facecolor = colors[k], alpha = 0.5, label = r"$z_{ssb} = %.1f$" % z_ssb, rasterized = True)

		k += 1
	t += 1
	for i in range(8):
		t = 1
		titlestring = tit[i] + ", \quad L = %.1f \mathrm{Mpc} h^{-1}$" % L

		mpl.figure(i+1)
		mpl.axhline(y[i], color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.suptitle(titlestring)
		mpl.legend()
		mpl.xscale("log")
		mpl.xlabel("$M_{\odot}h^{-1}$")
		mpl.ylabel(tit[i]+"$", rotation = 0, labelpad = 10)
		mpl.savefig("Figures\Symmetron\Shadingvar%.0fL%.0f.pdf" % (i, L), bbox_inches = "tight")
		mpl.clf()
	mpl.close("all")
