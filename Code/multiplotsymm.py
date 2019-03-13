"""
cd c:\Users\Halvor\Documents\Master\Kode\
python multiplotssymm.py

"""
#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys

def placementmass( massarr, z_ssb_arr, L_arr, beta_arr ):
	label1 = r"$\mathrm{h}^{-1} \mathrm{M}_{\odot}$"

	if len(z_ssb_arr) > 1:
		L = L_arr[0]
		beta =beta_arr[0]
		for z_ssb in z_ssb_arr:	
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
				filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(M), z_ssb, L, beta)
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


			seccondvar = r"$z_{\mathrm{ssb}} = %.1f$ " %  z_ssb

			mpl.figure(1)
			mpl.plot(variable, d_vir, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{vir}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")
			
			mpl.figure(2)
			mpl.plot(variable, d_ta, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{ta}$", fontsize = 15, rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(3)
			mpl.plot(variable, rvirrta, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(4)
			mpl.plot(variable, avirata, linewidth = 0.75, label = seccondvar)
			
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\frac{a_{vir}}{a_{ta}}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(5)
			mpl.plot(variable, avir, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			
			mpl.ylabel(r"$a_{vir}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(6)
			mpl.plot(variable, d_c, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{c}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")

			mpl.figure(7)
			mpl.plot(variable, abs(np.array(acoll)-1), linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$|a_c - 1|$", rotation = 0, labelpad = 20)
			mpl.xscale("log")


			mpl.figure(8)
			mpl.plot(variable, delta_c, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\delta_i$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")
			mpl.yscale("log")

		mpl.figure(1)
		mpl.axhline(y = 295.618, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{h}^{-1} \mathrm{Mpc}$, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Zdvir.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(2)
		mpl.axhline(y = 7.085, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{h}^{-1} \mathrm{Mpc}$, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Zdta.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(3)
		mpl.axhline(y = 0.481722, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{h}^{-1} \mathrm{Mpc}$, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Zrvirrta.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(4)
		mpl.axhline(y = 1.66974, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{h}^{-1} \mathrm{Mpc}$, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Zavirata.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(5)
		mpl.axhline(y = 0.91889, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{Mpc}$, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Zavir.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(6)
		mpl.axhline(y = 381.00185776, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{h}^{-1}\mathrm{Mpc} $, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Zdelta_c.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(7)
		#mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{h}^{-1}\mathrm{Mpc} $, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Za_c.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(8)
		mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$L = %.1f \mathrm{h}^{-1}\mathrm{Mpc} $, $\beta = %.1f$" % (L, beta))
		mpl.savefig("Figures\Symmetron\Zdelta_i.pdf", bbox_inches = "tight")
		mpl.clf()



	elif len(L_arr) > 1:
		z_ssb = z_ssb_arr[0]
		beta = beta_arr[0]
		for L in L_arr:	
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
				filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(M), z_ssb, L, beta)
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


			seccondvar = r"$L = %.1f \mathrm{Mpc}\mathrm{h}^{-1}$ " %  L

			mpl.figure(1)
			mpl.plot(variable, d_vir, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{vir}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")
			
			mpl.figure(2)
			mpl.plot(variable, d_ta, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{ta}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")

			mpl.figure(3)
			mpl.plot(variable, rvirrta, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(4)
			mpl.plot(variable, avirata, linewidth = 0.75, label = seccondvar)
			
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\frac{a_{vir}}{a_{ta}}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(5)
			mpl.plot(variable, avir, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			
			mpl.ylabel(r"$a_{vir}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(6)
			mpl.plot(variable, d_c, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{c}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")

			mpl.figure(7)
			mpl.plot(variable, abs(np.array(acoll)-1), linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$|a_c - 1|$", rotation = 0, labelpad = 20)
			mpl.xscale("log")


			mpl.figure(8)
			mpl.plot(variable, delta_c, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\delta_i$", rotation = 0, labelpad = 20)
			mpl.xscale("log")
			mpl.yscale("log")

		mpl.figure(1)
		mpl.axhline(y = 295.6181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $\beta = %.1f$" % (z_ssb, beta))
		mpl.savefig("Figures\Symmetron\Ldvir.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(2)
		mpl.axhline(y = 7.0985, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend(loc = 2)
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $\beta = %.1f $" % (z_ssb, beta))
		mpl.savefig("Figures\Symmetron\Ldta.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(3)
		mpl.axhline(y = 0.481722, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.savefig("Figures\Symmetron\Lrvirrta.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(4)
		mpl.axhline(y = 1.6697, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $\beta = %.1f $" % (z_ssb, beta))
		mpl.savefig("Figures\Symmetron\Lavirata.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(5)
		mpl.axhline(y = 0.91889, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $\beta = %.1f $" % (z_ssb, beta))
		mpl.savefig("Figures\Symmetron\Lavir.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(6)
		mpl.axhline(y = 381.00186, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $\beta = %.1f$" % (z_ssb, beta))
		mpl.savefig("Figures\Symmetron\Ldelta_c.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(7)
		#mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $\beta = %.1f $" % (z_ssb, beta))
		mpl.savefig("Figures\Symmetron\La_c.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(8)
		mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $\beta = %.1f $" % (z_ssb, beta))
		mpl.savefig("Figures\Symmetron\Ldelta_i.pdf", bbox_inches = "tight")
		mpl.clf()




	elif len(beta_arr) > 1:
		z_ssb = z_ssb_arr[0]
		L = L_arr[0]
		for beta in beta_arr:	
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
				filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(M), z_ssb, L, beta)
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


			seccondvar = r"$\beta= %.1f$ " %  beta

			mpl.figure(1)
			mpl.plot(variable, d_vir, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{vir}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")
			
			mpl.figure(2)
			mpl.plot(variable, d_ta, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{ta}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")

			mpl.figure(3)
			mpl.plot(variable, rvirrta, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(4)
			mpl.plot(variable, avirata, linewidth = 0.75, label = seccondvar)
			
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\frac{a_{vir}}{a_{ta}}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(5)
			mpl.plot(variable, avir, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			
			mpl.ylabel(r"$a_{vir}$", rotation = 0, labelpad = 20)
			mpl.xscale("log")

			mpl.figure(6)
			mpl.plot(variable, d_c, linewidth = 0.75, label = seccondvar)
			mpl.legend(fontsize = 15)
			mpl.xlabel(label1, fontsize = 15)
			mpl.ylabel(r"$\Delta_{c}$", rotation = 0, labelpad = 20, fontsize = 15)
			mpl.xscale("log")

			mpl.figure(7)
			mpl.plot(variable, abs(np.array(acoll)-1), linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$|a_c - 1|$", rotation = 0, labelpad = 20)
			mpl.xscale("log")


			mpl.figure(8)
			mpl.plot(variable, delta_c, linewidth = 0.75, label = seccondvar)
			mpl.legend()
			mpl.xlabel(label1)
			mpl.ylabel(r"$\delta_i$", rotation = 0, labelpad = 20)
			mpl.xscale("log")
			mpl.yscale("log")

		mpl.figure(1)
		mpl.axhline(y = 295.6181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f \mathrm{h}^{-1} \mathrm{Mpc}$" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betadvir.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(2)
		mpl.axhline(y = 7.0985, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f  \mathrm{h}^{-1}\mathrm{Mpc}$" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betadta.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(3)
		mpl.axhline(y = 0.481722, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f \mathrm{h}^{-1}\mathrm{Mpc} $" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betarvirrta.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(4)
		mpl.axhline(y = 1.6697, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f  \mathrm{h}^{-1}\mathrm{Mpc}$" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betaavirata.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(5)
		mpl.axhline(y = 0.918889, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f\mathrm{h}^{-1} \mathrm{Mpc} $" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betaavir.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(6)
		mpl.axhline(y = 381.0018577, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f \mathrm{h}^{-1}\mathrm{Mpc} $" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betadelta_c.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(7)
		#mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f \mathrm{h}^{-1}\mathrm{Mpc} $" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betaa_c.pdf", bbox_inches = "tight")
		mpl.clf()

		mpl.figure(8)
		mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
		mpl.legend()
		mpl.title(r"$z_{\mathrm{ssb}} = %.1f$, $L = %.1f \mathrm{h}^{-1} \mathrm{Mpc} $" % (z_ssb, L))
		mpl.savefig("Figures\Symmetron\Betadelta_i.pdf", bbox_inches = "tight")
		mpl.clf()
	return 

masses = np.logspace(10, 18, 30)
beta_list = np.linspace(0, 1.5, 30)
L_list = np.logspace(-1, 1, 30)
z_ssb_list = np.linspace(0, 3, 30)

placementmass(masses, [z_ssb_list[5], z_ssb_list[10], z_ssb_list[20], z_ssb_list[-1]], [1], [0.5])
placementmass(masses, [1], [1], [beta_list[5], beta_list[10], beta_list[20]])
placementmass(masses, [1], [L_list[0], L_list[15], L_list[-1]], [0.5])