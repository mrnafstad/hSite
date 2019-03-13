"""
cd c:\Users\Halvor\Documents\Master\Kode\
python inspectsymm.py

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

#in natural (Planck) units:
H_0 = 1.491e-32 #9.369e-32										#eV 
G = 1.0												#1

Msun = 1.989e30*3.755e-21 #5.977e-22 							#eV-1
M_pl = 4.341e-9*3.755e-21 #5.977e-22  
#M_pl = 1/np.sqrt(8*np.pi)
GtimeM_sun = G*Msun
rho_c0 = 3.*H_0**2*M_pl**2# /(8.*np.pi*G)

def gammaBack(Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, perturbation):
	return 1 

def d_rho(Omega_m0, Lambda, delta_i, R):
	return Omega_m0*(1 + delta_i)/R**3 								#This!
	

def LCDMacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, delta_rho, n ):
	return (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/E(Omega_m0, Lambda, a, 0, 0, 0, 0, 0, 0, 0 )

def ELCDMnorad(Omega_m0, Lambda, a, gamma, beta, s, t, u, v, w):
	return Omega_m0/a**3 + Lambda

def chameleonacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, delta_rho, n ):
	return  ( Lambda*r - gamma( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, True )*Omega_m0*(1 + delta_i)/(2*r**2) + \
		3./2.*drdy*Omega_m0/a**3)/E(Omega_m0, Lambda, a, gamma, 0, 0, 0, 0, 0, 0 )

def Echamnorad( Omega_m0, Lambda, a, gamma, beta, n, M, f_R0, delta_rho, r ):
	return gamma( Omega_m0, Lambda, 0, beta, n, M, f_R0, delta_rho, a, a, False )*Omega_m0/a**3 + Lambda

def Phi_N( r, M, Omega_m0, delta_i ):
	Phi_N = ((G*M*H_0)**2*Omega_m0*(1 + delta_i)/2)**(1./3.)/r
	return Phi_N

def gamma_symm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, delta_rho, r, a, perturbation ):
	#In comparison to chameleon, g -> beta, z_ssb -> f_R0, L -> n notationally

	deltaR_overR = DRORsymm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, r, a, delta_rho, perturbation )

	
	Rs = 1 - deltaR_overR
	return 1 + 2*beta**2*( 1 - Rs**3)
def phi_symm(rho, rho_ssb, L, M, beta):
	return 2*L**2/M_pl*beta*rho_ssb*np.sqrt(1 - rho/rho_ssb)

def DRORsymm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, r, a, delta_rho, perturbation ):

	if beta**2 <= 0:
		return 0

	else:
		
		mu = 1/(np.sqrt(2)*L)
		a_ssb = 1./float(1. + z_ssb)
		
		M_phisquare = 2*Omega_m0*rho_c0*L**2/a_ssb**3
		l = M_pl**2 * a_ssb**3/(4*Omega_m0*rho_c0*beta**2*L**2)
		
		delta = delta_rho(Omega_m0, Lambda, delta_i, r)
		rho_ssb = rho_c0*Omega_m0/a_ssb**3
		rho_b = rho_c0*delta_rho(Omega_m0, Lambda, 0, a)
		rho_p = rho_c0*delta

		Phi = Phi_N(r, M, Omega_m0, delta_i)
		if (rho_p < rho_ssb and perturbation):
			phi_c = phi_symm(rho_p, rho_ssb, L, M, beta)
		else:
			phi_c = 0

		if rho_b < rho_ssb:
			phi_inf = phi_symm(rho_b, rho_ssb, L, M, beta)
		else:
			phi_inf = 0

			#phi_c = phi_symm(rho_p, rho_ssb, L, M, beta)
		d = abs(phi_c - phi_inf)/(6*beta*Phi*M_pl)

		if d < 0:
			d = 0
		elif d > 1:
			d = 1


		return d

def placementmass( massarr, beta_arr, L, z_ssb ):

	maxima = []
	di = []
	rvir = []
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

		maxima.append(np.argmin(d_vir))
		di.append(delta_c)
		rvir.append(rvirrta[np.argmin(d_vir)])

	actmax = np.argmin(rvir)
	#minf = beta_arr[maxima[actmax]]
	minmass = massarr[actmax]
	delta_is = di[actmax]

	return maxima, delta_is

def d_i_arr( mass, beta_arr, L_arr, z_ssb_arr ):
	d_vir1 = []
	d_ta1 = []
	d_c1 = []
	rvirrta1 = []
	avirata1 = []
	#variable1 = []
	avir1 = []
	acoll1 = []
	delta_c1 = []
	for L in L_arr:
		filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(mass), 1, L, 0.5)
		file = open(filename, "r")
		i = 0
		for j in file.readlines():
			try:
				value = float(j)
			except ValueError:
				break
			if i == 0:
				#variable.append(value)
				i += 1
			elif i == 1:
				d_vir1.append(value)
				i += 1
			elif i == 2:
				d_ta1.append(value)
				i += 1
			elif i == 3:
				d_c1.append(value)
				i += 1
			elif i == 4:
				rvirrta1.append(value)
				i += 1
			elif i == 5:
				avirata1.append(value)
				i += 1
			elif i == 6:
				avir1.append(value)
				i += 1

			elif ( i == 7 ):
				acoll1.append(value)
				i += 1

			elif i == 8:
				delta_c1.append(value)
				i = 0
					

		file.close()

	d_vir2 = []
	d_ta2 = []
	d_c2 = []
	rvirrta2 = []
	avirata2 = []
	#variable2 = []
	avir2 = []
	acoll2 = []
	delta_c2 = []
	for z_ssb in z_ssb_arr:
		filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(mass), z_ssb, 1, 0.5)
		file = open(filename, "r")
		i = 0
		for j in file.readlines():
			try:
				value = float(j)
			except ValueError:
				break
			if i == 0:
				#variable.append(value)
				i += 1
			elif i == 1:
				d_vir2.append(value)
				i += 1
			elif i == 2:
				d_ta2.append(value)
				i += 1
			elif i == 3:
				d_c2.append(value)
				i += 1
			elif i == 4:
				rvirrta2.append(value)
				i += 1
			elif i == 5:
				avirata2.append(value)
				i += 1
			elif i == 6:
				avir2.append(value)
				i += 1

			elif ( i == 7 ):
				acoll2.append(value)
				i += 1

			elif i == 8:
				delta_c2.append(value)
				i = 0
					

		file.close()
	d_vir3 = []
	d_ta3 = []
	d_c3 = []
	rvirrta3 = []
	avirata3 = []
	#variable2 = []
	avir3 = []
	acoll3 = []
	delta_c3 = []
	for beta in beta_arr:
		filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(mass), 1, 1, beta)
		file = open(filename, "r")
		i = 0
		for j in file.readlines():
			try:
				value = float(j)
			except ValueError:
				break
			if i == 0:
				#variable.append(value)
				i += 1
			elif i == 1:
				d_vir3.append(value)
				i += 1
			elif i == 2:
				d_ta3.append(value)
				i += 1
			elif i == 3:
				d_c3.append(value)
				i += 1
			elif i == 4:
				rvirrta3.append(value)
				i += 1
			elif i == 5:
				avirata3.append(value)
				i += 1
			elif i == 6:
				avir3.append(value)
				i += 1

			elif ( i == 7 ):
				acoll3.append(value)
				i += 1

			elif i == 8:
				delta_c3.append(value)
				i = 0
					

		file.close()
	return delta_c1, delta_c2, delta_c3

N = 5000000
y0 = np.log(1e-4)
y = np.linspace(y0, -1e-15, N)
a = np.exp(y)
Omega_m0 = 0.25
Lambda = 0.75

def figs(M, z_ssb, beta, L, labelstring, lcol, delta_i):

	rho_ssb = rho_c0*Omega_m0*(1 + z_ssb)**3
	delta_rho_ssb = Omega_m0*(1 + z_ssb)**3
	r = open("Numbers\Symmetron\Sym2\RadiusM%.15fz_ssb%.15fL%.2fbeta%.2f.txt" % (np.log10(M), z_ssb, L, beta), "r")
	rad = np.array(r.readlines(), dtype = np.float32)
	r.close()
	ta = np.argmax(rad)

	mpl.figure(1)
	mpl.plot(a, rad, linestyle = "-", color = lcol, linewidth = 0.75, label = labelstring)
	mpl.plot([a[ta]], [rad[ta]], "o", color = lcol, markersize = 5, label = r"$a_{\mathrm{ta}} = %.4f$" % a[ta])
	L = L*3.09e22/1.97e-7/0.7
	M = M*Msun/0.7

	geff = np.zeros(len(a))
	gback = np.zeros(len(a))

	phi = np.zeros(len(a))
	phi_b = np.zeros(len(a))

	Phi = Phi_N( rad, M*Msun/0.7, Omega_m0, delta_i )

	Delta = d_rho(Omega_m0, Lambda, delta_i, rad)
	Delta_b = d_rho(Omega_m0, Lambda, 0, a)
	
	rho_b = rho_c0*d_rho(Omega_m0, Lambda, 0, a)
	rho = rho_c0*d_rho(Omega_m0, Lambda, delta_i, rad)

	for i in range(len(a)):
		geff[i] = gamma_symm(Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, d_rho, rad[i], a[i], True) - 1
		
		if Delta[i]/delta_rho_ssb < 1:
			phi[i] = phi_symm(rho[i], rho_ssb, L, M, beta)
			

		if Delta_b[i]/delta_rho_ssb < 1:
			phi_b[i] = phi_symm(rho_b[i], rho_ssb, L, M, beta)
			
	gmax_symm = 2*beta**2
	phissymm = [max(phi), max(phi_b)]

	geff_max = max(geff)
	geff_max_ind = np.argmax(geff)
	
	mpl.figure(2)
	mpl.plot(a, geff, linestyle = "-", color = lcol, linewidth = 0.75, label = labelstring)
	mpl.plot([a[geff_max_ind]], [geff_max], "o", color = lcol, markersize = 5, label = r"$a_{\mathrm{max}} = %.3f$" % a[geff_max_ind])

	mpl.figure(3)
	mpl.plot(a, phi/max(phissymm), linestyle = "-", color = lcol, linewidth = 0.75, label = r"$\phi_{\mathrm{c}}$: " + labelstring)
	mpl.plot(a, phi_b/max(phissymm), linestyle = "--", color = lcol, linewidth = 0.75, label = r"$\phi_{\mathrm{\infty}}$: " + labelstring)

	return

def savefigs( titlestring, filename ):
	mpl.figure(1)
	mpl.legend()
	mpl.title(titlestring)
	mpl.xlabel("a")
	mpl.ylabel(r"$\tilde{R}$", rotation = 0, labelpad = 10)
	mpl.savefig("Figures\Symmetron\Radforcomp_"+filename + ".pdf",addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
	mpl.clf()

	mpl.figure(2)
	mpl.legend()
	mpl.title(titlestring)
	mpl.xlabel("a")
	mpl.ylabel(r"$\gamma_{eff}$", rotation = 0, labelpad = 10)
	mpl.savefig("Figures\Symmetron\Gammaforcomp_"+filename + ".pdf",addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
	mpl.clf()

	mpl.figure(3)
	mpl.legend()
	mpl.title(titlestring)
	mpl.xlabel("a")
	mpl.ylabel(r"$\frac{\phi}{\phi_{\mathrm{max}}}$", rotation = 0, labelpad = 10)
	mpl.savefig("Figures\Symmetron\Fieldforcomp_"+filename + ".pdf",addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
	mpl.clf()
	return


masses = np.logspace(10, 18, 30)
beta_list = np.linspace(0, 1.5, 30)
L_list = np.logspace(-1, 1, 30)
z_ssb_list = np.linspace(0, 3, 30)

beta_ind = np.abs(beta_list - 0.5).argmin()
beta_close = beta_list[beta_ind]
beta_min = beta_list[beta_ind - 5]
beta_plus = beta_list[beta_ind + 5]

L_ind = np.abs(L_list - 1).argmin()
L_close = L_list[L_ind]
L_var_min = L_list[L_ind - 5]
L_var_plus = L_list[L_ind + 5]

z_ssb_ind = np.abs((z_ssb_list - 1)).argmin()
z_ssb_close = z_ssb_list[z_ssb_ind]
z_ssb_min = z_ssb_list[z_ssb_ind - 5]
z_ssb_plus   = z_ssb_list[z_ssb_ind + 10]

maxim, deltas = placementmass(masses, [beta_close], 1, 1)

maxmass = masses[maxim]

delta_i_L, delta_i_z, delta_i_beta = d_i_arr(maxmass, beta_list, L_list, z_ssb_list)

figs(maxmass, 1, 0.5, L_close, r"$L = %.2f \mathrm{Mpc} \mathrm{h}^{-1}$" % L_close, "red", delta_i_L[L_ind])
figs(maxmass, 1, 0.5, L_var_min, r"$L = %.2f \mathrm{Mpc} \mathrm{h}^{-1}$" % L_var_min, "blue", delta_i_L[L_ind - 5])
figs(maxmass, 1, 0.5, L_var_plus, r"$L = %.2f \mathrm{Mpc} \mathrm{h}^{-1}$" % L_var_plus, "green", delta_i_L[L_ind + 5])

savefigs(r"$M = 10^{%.2f}  \mathrm{M}_{\odot} \mathrm{h}^{-1}$, $\beta = 0.5$, $z_{\mathrm{ssb}} = 1$" % np.log10(maxmass), "Lvar")

figs(maxmass, z_ssb_close, 0.5, 1, r"$z_{\mathrm{ssb}} = %.2f$" % z_ssb_close, "red", delta_i_z[z_ssb_ind])
figs(maxmass, z_ssb_min, 0.5, 1, r"$z_{\mathrm{ssb}} = %.2f$" % z_ssb_min, "blue", delta_i_z[z_ssb_ind - 5])
figs(maxmass, z_ssb_plus, 0.5, 1, r"$z_{\mathrm{ssb}} = %.2f$" % z_ssb_plus, "green", delta_i_z[z_ssb_ind + 10])

savefigs(r"$M = 10^{%.2f} \mathrm{M}_{\odot} \mathrm{h}^{-1}$, $\beta = 0.5$, $L = 1 \mathrm{Mpc} \mathrm{h}^{-1}$" % np.log10(maxmass), "Zvar")

figs(maxmass, 1, beta_close, 1, r"$\beta = %.2f$" % beta_close, "red", delta_i_beta[beta_ind])
figs(maxmass, 1, beta_min, 1, r"$\beta = %.2f$" % beta_min, "blue", delta_i_beta[beta_ind - 5])
figs(maxmass, 1, beta_plus, 1, r"$\beta = %.2f$" % beta_plus, "green", delta_i_beta[beta_ind + 5])

savefigs(r"$M = 10^{%.2f} \mathrm{M}_{\odot} \mathrm{h}^{-1}$, $z_{\mathrm{ssb}} = 1$, $L = 1 \mathrm{Mpc} \mathrm{h}^{-1}$" % np.log10(maxmass), "betavar")
