"""
cd c:\Users\Halvor\Documents\Master\Kode\
python multiplots.py

"""
#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys

def placementmass( label1, massarr, farr ):

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


		seccondvar = r"$|1 - f_{\mathcal{R}_0}| = 10^{%.0f}$ " %  np.log10(abs(1 - f_R0))

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
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")

		mpl.figure(4)
		mpl.plot(variable, avirata, linewidth = 0.75, label = seccondvar)
		
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$\frac{a_{vir}}{a_{ta}}$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")

		mpl.figure(5)
		mpl.plot(variable, avir, linewidth = 0.75, label = seccondvar)
		mpl.legend(fontsize = 15)
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
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$|a_c - 1|$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")


		mpl.figure(8)
		mpl.plot(variable, delta_c, linewidth = 0.75, label = seccondvar)
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$\delta_i$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")
		mpl.yscale("log")

	mpl.figure(1)
	mpl.axhline(y = 294.605, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmaindvir.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(2)
	mpl.axhline(y = 7.09, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmaindta.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(3)
	mpl.axhline(y = 0.48181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmainrvirrta.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(4)
	mpl.axhline(y = 1.66881, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmainavirata.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(5)
	mpl.axhline(y = 0.916595229171737, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmainavir.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(6)
	mpl.axhline(y = 382.562, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmaindelta_c.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(7)
	#mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmaina_c.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(8)
	mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Massmaindelta_i.pdf", bbox_inches = "tight")
	mpl.clf()

	return 

def placementF( label1, massarr, farr ):

	for M in massarr:	
		d_vir = []
		d_ta = []
		d_c = []
		rvirrta = []
		avirata = []
		variable = []
		avir = []
		acoll = []
		delta_c = []
		for f_R0 in farr:
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


		seccondvar = r"$M = 10^{%.0f} M_{\odot} h^{-1}$" % np.log10(M)

		mpl.figure(1)
		mpl.plot(abs(farr - 1), d_vir, linewidth = 0.75, label = seccondvar)
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$\Delta_{vir}$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")
		
		mpl.figure(2)
		mpl.plot(abs(farr - 1), d_ta, linewidth = 0.75, label = seccondvar)
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$\Delta_{ta}$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")

		mpl.figure(3)
		mpl.plot(abs(farr - 1), rvirrta, linewidth = 0.75, label = seccondvar)
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")

		mpl.figure(4)
		mpl.plot(abs(farr - 1), avirata, linewidth = 0.75, label = seccondvar)
		
		mpl.legend()
		mpl.xlabel(label1)
		mpl.ylabel(r"$\frac{a_{vir}}{a_{ta}}$", rotation = 0, labelpad = 20)
		mpl.xscale("log")

		mpl.figure(5)
		mpl.plot(abs(farr - 1), avir, linewidth = 0.75, label = seccondvar)
		mpl.legend()
		mpl.xlabel(label1)
		
		mpl.ylabel(r"$a_{vir}$", rotation = 0, labelpad = 20)
		mpl.xscale("log")

		mpl.figure(6)
		mpl.plot(abs(farr - 1), d_c, linewidth = 0.75, label = seccondvar)
		mpl.legend(fontsize = 15)
		mpl.xlabel(label1, fontsize = 15)
		mpl.ylabel(r"$\Delta_{c}$", rotation = 0, labelpad = 20, fontsize = 15)
		mpl.xscale("log")

		mpl.figure(7)
		mpl.plot(abs(farr - 1), abs(np.array(acoll)-1), linewidth = 0.75, label = seccondvar)
		mpl.legend()
		mpl.xlabel(label1)
		mpl.ylabel(r"$|a_c - 1|$", rotation = 0, labelpad = 20)
		mpl.xscale("log")


		mpl.figure(8)
		mpl.plot(abs(farr - 1), delta_c, linewidth = 0.75, label = seccondvar)
		mpl.legend()
		mpl.xlabel(label1)
		mpl.ylabel(r"$\delta_i$", rotation = 0, labelpad = 20)
		mpl.xscale("log")
		mpl.yscale("log")

	mpl.figure(1)
	mpl.axhline(y = 294.605, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmaindvir.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(2)
	mpl.axhline(y = 7.09, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmaindta.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(3)
	mpl.axhline(y = 0.48181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmainrvirrta.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(4)
	mpl.axhline(y = 1.66881, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmainavirata.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(5)
	mpl.axhline(y = 0.916595229171737, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmainavir.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(6)
	mpl.axhline(y = 382.562, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmaindelta_c.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(7)
	#mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmaina_c.pdf", bbox_inches = "tight")
	mpl.clf()

	mpl.figure(8)
	mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.legend()
	mpl.savefig("Figures\Gback\Fmaindelta_i.pdf", bbox_inches = "tight")
	mpl.clf()

	return 


masses = np.logspace(10, 18, 30)
f_R0_list = np.logspace(-12, -3, 30) + 1

F_for_masses = [f_R0_list[12], f_R0_list[19], f_R0_list[25], f_R0_list[29]]
M_for_curvatures = [masses[0], masses[9], masses[19], masses[29]]


placementmass(r"$\mathrm{h}^{-1} \mathrm{M}_{\odot}$", masses, F_for_masses)
placementF(r"$|f_{\mathcal{R}_0} - 1|$", M_for_curvatures, f_R0_list)

