"""
cd c:\Users\Halvor\Documents\Master\Kode\
python inspectcham.py

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

def phi_cham( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho):
	if delta_i == 0:
		delta = delta_rho(Omega_m0, Lambda, 0, a)
	else:
		delta = delta_rho(Omega_m0, Lambda, delta_i, r)

	phi = M_pl*abs(1 - f_R0)/(2*beta)*( (Omega_m0 + 4*Lambda)/(delta + 4*Lambda) )**(n + 1)

	#phi = M_pl*abs(1-f_R0)/(2*beta)*( (n + 1)**2 *(Omega_m0 + 4*Lambda)/(delta + 4*Lambda) )**(1/(n**2 + n + 1))

	return phi

def DROR( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho, perturbation ):

	if beta**2 <= 0:
		return 0

	else:
		if perturbation:			
			Phi = Phi_N(r, M, Omega_m0, delta_i)
			phi_c = phi_cham(Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho)
			phi_inf = phi_cham(Omega_m0, Lambda, 0, beta, n, M, f_R0, r, a, delta_rho)
			d = abs( phi_inf - phi_c )/(6*beta*M_pl*Phi)


		else:

			phi_c = 0
			Phi = Phi_N(a, M, Omega_m0, 0)
			phi_inf = phi_cham(Omega_m0, Lambda, 0, beta, n, M, f_R0, a, a, delta_rho)
			d = abs( phi_inf - phi_c )/(6*beta*M_pl*Phi)
			
			#d = 1
		
		try:
			for i in range(len(a)):
				if d[i] < 0:
					d[i] = 0
				elif d[i] > 1:
					d[i] = 1
		except:
			if d < 0:
				d = 0
			elif d > 1:
				d = 1

		
		return d



def gammaHuSawicki1( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, perturbation ):
	
	deltRoverR = DROR( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho, perturbation )

	Rs = 1 - deltRoverR
	return 1 + 2*beta**2*( 1 - Rs**3)

def placementmass( massarr, farr ):

	maxima = []
	di = []
	rvir = []
	d_i_s = []
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
		d_i_s.append(delta_c)

		maxima.append(np.argmax(rvirrta))
		di.append(delta_c)
		rvir.append(rvirrta[np.argmax(rvirrta)])

	actmax = np.argmax(rvir)
	minf = farr[maxima[actmax]]
	minmass = massarr[actmax]
	delta_is = di[actmax]

	return d_i_s


masses = np.logspace(10, 18, 30)
f_R0_list = np.flip(np.logspace(-12, -3, 30) + 1, 0)

f_R0 = f_R0_list[25]


delta_is = placementmass(masses, f_R0_list)
"""
minmin = minima 
minplus = minima 
Mcomp1 = masses[minima]
f_R0_comp1 = f_R0_list[massloc[minplus] - 3]
Mcomp2 = masses[minima]
f_R0_comp2 = f_R0_list[massloc[minmin] + 3]
"""
#print Mlocalmin, f_R0

f_R0_loc = np.abs(f_R0_list - ( 1+ 1.62e-6)).argmin()
f_R0 = f_R0_list[f_R0_loc]
f_R0_comp1 = f_R0_list[f_R0_loc - 3]
f_R0_comp2 = f_R0_list[f_R0_loc + 3]

massloc = np.abs(masses - 5.754e12).argmin()
Mlocalmin = masses[massloc]
Mcomp1 = masses[massloc]
Mcomp2 = masses[massloc]

delta_i1 = delta_is[f_R0_loc][massloc]
delta_i2 = delta_is[f_R0_loc - 3][massloc]
delta_i3 = delta_is[f_R0_loc + 3][massloc]


rlocmin = open("Numbers\Chameleon\RadiusM%.15ff_R0%.15f.txt" % (np.log10(Mlocalmin), np.log10(abs(1-f_R0))), "r")
radlocmin = np.array(rlocmin.readlines(), dtype = np.float32)
rlocmin.close()
rcom1 = open("Numbers\Chameleon\RadiusM%.15ff_R0%.15f.txt" % (np.log10(Mcomp1), np.log10(abs(1-f_R0_comp1))), "r")
radcom1 = np.array(rcom1.readlines(), dtype = np.float32)
rcom1.close()
rcom2 = open("Numbers\Chameleon\RadiusM%.15ff_R0%.15f.txt" % (np.log10(Mcomp2), np.log10(abs(1-f_R0_comp2))), "r")
radcom2 = np.array(rcom2.readlines(), dtype = np.float32)
rcom2.close()

N = 5000000
y0 = np.log(1e-4)
y = np.linspace(y0, -1e-15, N)
a = np.exp(y)
betacham = 1./np.sqrt(6)
ncham = 1.0
Omega_m0 = 0.25
Lambda = 0.75

"""
titlestring = r"$|f_{\mathcal{R}_0} - 1| = 10^{%.2f}$" % np.log10(abs(f_R0-1))
labelstring1 = r"$\mathrm{M} = 10^{%.2f}\mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mlocalmin)
labelstring2 = r"$\mathrm{M} = 10^{%.2f}\mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mcomp1)
labelstring3 = r"$\mathrm{M} = 10^{%.2f}\mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mcomp2)

"""
titlestring = r"$\mathrm{M} = 10^{%.2f}\mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mlocalmin)
labelstring1 = r"$|f_{\mathcal{R}_0} - 1| = 10^{%.2f}$" % np.log10(abs(f_R0-1))
labelstring2 = r"$|f_{\mathcal{R}_0} - 1| = 10^{%.2f}$" % np.log10(abs(f_R0_comp1-1))
labelstring3 = r"$|f_{\mathcal{R}_0} - 1| = 10^{%.2f}$" % np.log10(abs(f_R0_comp2-1))


mpl.figure(1)
mpl.plot(a, radlocmin, "b", linewidth = 0.75, label = labelstring1)
mpl.plot(a, radcom1, "r", linewidth = 0.75, label =  labelstring2)
mpl.plot(a, radcom2, "g", linewidth = 0.75, label =  labelstring3)

mpl.plot(a[np.argmax(radlocmin)], radlocmin[np.argmax(radlocmin)], "bo", markersize = 5, label = r"$a_{\mathrm{ta}} = %.4f$" % a[np.argmax(radlocmin)])
mpl.plot(a[np.argmax(radcom1)], radcom1[np.argmax(radcom1)], "ro", markersize = 5, label = r"$a_{\mathrm{ta}} = %.4f$" % a[np.argmax(radcom1)])
mpl.plot(a[np.argmax(radcom2)], radcom2[np.argmax(radcom2)], "go", markersize = 5, label = r"$a_{\mathrm{ta}} = %.4f$" % a[np.argmax(radcom2)])

mpl.legend()
mpl.xlabel("a")
mpl.ylabel(r"$\tilde{R}$", rotation = 0, labelpad = 10)
mpl.title(titlestring)
print "--"
mpl.savefig("Figures\Chameleon\Radforcomp_fvar2.pdf", addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
mpl.clf()

geff_locmin = np.zeros(len(a))
geff_comp1 = np.zeros(len(a))
geff_comp2 = np.zeros(len(a))

gback = np.zeros(len(a))

phi_locmin = np.zeros(len(a))
phi_comp1 = np.zeros(len(a))
phi_comp2 = np.zeros(len(a))

phi_blocmin = np.zeros(len(a))
phi_bcomp1 = np.zeros(len(a))
phi_bcomp2 = np.zeros(len(a))

Phi_newtcham_locmin = Phi_N( radlocmin, Mlocalmin*Msun/0.7, Omega_m0, delta_i1)
Phi_newtcham_comp1 = Phi_N( radcom1, Mcomp1*Msun/0.7, Omega_m0, delta_i2 )
Phi_newtcham_comp2 = Phi_N( radcom2, Mcomp2*Msun/0.7, Omega_m0, delta_i3 )

dlocmin = d_rho(Omega_m0, Lambda, delta_i1, radlocmin)
dcom1 = d_rho(Omega_m0, Lambda, delta_i2, radcom1)
dcom2 = d_rho(Omega_m0, Lambda, delta_i3, radcom2)

curve1 = np.zeros(len(a))
curve2 = np.zeros(len(a))
curve3 = np.zeros(len(a))

first1 = first2 = first3 = True


for i in range(len(a)):
	geff_locmin[i] = gammaHuSawicki1(Omega_m0, Lambda, delta_i1, betacham, ncham, Mlocalmin*Msun/0.7, f_R0, d_rho, radlocmin[i], a[i], True) - 1
	geff_comp1[i] = gammaHuSawicki1(Omega_m0, Lambda, delta_i2, betacham, ncham, Mcomp1*Msun/0.7, f_R0_comp1, d_rho, radcom1[i], a[i], True) - 1
	geff_comp2[i] = gammaHuSawicki1(Omega_m0, Lambda, delta_i3, betacham, ncham, Mcomp2*Msun/0.7, f_R0_comp2, d_rho, radcom2[i], a[i], True) - 1

	if ( geff_locmin[i] > 1e-4 and first1 ):
		gammanon1 = i
		first1 = False
	if ( geff_comp1[i] > 1e-4 and first2 ):
		first2 = False
		gammanon2 = i
	if ( geff_comp2[i] > 1e-4 and first3 ):
		first3 = False
		gammanon3 = i
	
	gback[i] = gammaHuSawicki1( Omega_m0, Lambda, 0, 0, 0, 0, 0, 0, a, a, False ) - 1

	phi_locmin[i] = phi_cham(Omega_m0, Lambda, delta_i1, betacham, ncham, Mlocalmin*Msun/0.7, f_R0, radlocmin[i], a[i], d_rho)
	phi_comp1[i] = phi_cham(Omega_m0, Lambda, delta_i2, betacham, ncham, Mcomp1*Msun/0.7, f_R0_comp1, radcom1[i], a[i], d_rho)
	phi_comp2[i] = phi_cham(Omega_m0, Lambda, delta_i3, betacham, ncham, Mcomp2*Msun/0.7, f_R0_comp2, radcom2[i], a[i], d_rho)

	phi_blocmin[i] = phi_cham(Omega_m0, Lambda, 0, betacham, ncham, Mlocalmin*Msun/0.7, f_R0, a[i], a[i], d_rho)
	phi_bcomp1[i] = phi_cham(Omega_m0, Lambda, 0, betacham, ncham, Mcomp1*Msun/0.7, f_R0_comp1, a[i], a[i], d_rho)
	phi_bcomp2[i] = phi_cham(Omega_m0, Lambda, 0, betacham, ncham, Mcomp2*Msun/0.7, f_R0_comp2, a[i], a[i], d_rho)

	curve1[i] = np.sqrt(M_pl*(f_R0 - 1)/(2*betacham*phi_locmin[i]))
	curve2[i] = np.sqrt(M_pl*(f_R0 - 1)/(2*betacham*phi_comp1[i]))
	curve3[i] = np.sqrt(M_pl*(f_R0 - 1)/(2*betacham*phi_comp2[i]))

mpl.plot(a, geff_locmin, "b", linewidth = 0.75, label = labelstring1)
mpl.plot(a, geff_comp1, "r", linewidth = 0.75, label = labelstring2)
mpl.plot(a, geff_comp2, "g", linewidth = 0.75, label = labelstring3)

mpl.plot(a[np.argmax(geff_locmin)], geff_locmin[np.argmax(geff_locmin)], "bo", markersize = 5, label = r"$a_{\mathrm{max}} = %.4f$" % a[np.argmax(geff_locmin)])
mpl.plot(a[np.argmax(geff_comp1)], geff_comp1[np.argmax(geff_comp1)], "ro", markersize = 5, label = r"$a_{\mathrm{max}} = %.4f$" % a[np.argmax(geff_comp1)])
mpl.plot(a[np.argmax(geff_comp2)], geff_comp2[np.argmax(geff_comp2)], "go", markersize = 5, label = r"$a_{\mathrm{max}} = %.4f$" % a[np.argmax(geff_comp2)])

mpl.plot(a[gammanon1], geff_locmin[gammanon1], "bx", markersize = 5, label = r"$a_{\mathrm{screening}} = %.4f$" % a[gammanon1])
mpl.plot(a[gammanon2], geff_comp1[gammanon2], "rx", markersize = 5, label = r"$a_{\mathrm{screening}} = %.4f$" % a[gammanon2])
mpl.plot(a[gammanon3], geff_comp2[gammanon3], "gx", markersize = 5, label = r"$a_{\mathrm{screening}} = %.4f$" % a[gammanon3])


#mpl.plot(a, gback, "c", linewidth = 0.75, label =  r"$background$")
mpl.legend()
mpl.xlabel("a")
mpl.title(titlestring)
mpl.ylabel(r"$\gamma_{eff}$", rotation = 0, labelpad = 10)
mpl.savefig("Figures\Chameleon\Gammaforcomp_fvar3.pdf",addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
mpl.clf()

"""
mpl.plot(a, dlocmin, "b--", linewidth = 0.75, label = labelstring1)
mpl.plot(a, dcom1, "r-", linewidth = 0.75, label = labelstring2)
mpl.plot(a, dcom2, "g:", linewidth = 0.75, label = labelstring3)
mpl.legend()
mpl.xlabel("a")
mpl.title(titlestring)
mpl.yscale("log")
mpl.ylabel(r"$\Delta_{\rho}$", rotation = 0, labelpad = 10)
mpl.savefig("Figures\Chameleon\Deltacomp_fvar2.pdf", bbox_inches = "tight")
mpl.clf()

mpl.plot(a, curve1, "b--", linewidth = 0.75, label = labelstring1)
mpl.plot(a, curve2, "r-", linewidth = 0.75, label = labelstring2)
mpl.plot(a, curve3, "g:", linewidth = 0.75, label = labelstring3)
mpl.legend()
mpl.xlabel("a")
mpl.title(titlestring)
mpl.yscale("log")
mpl.ylabel(r"$\frac{\mathcal{R}}{\mathcal{R}_0}$", rotation = 0, labelpad = 10)
mpl.savefig("Figures\Chameleon\Curvaturescalar_fvar2.pdf",addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
mpl.clf()

mpl.plot(a, phi_locmin,"b--", linewidth = 0.75, label = r"$\phi_c:$" + labelstring1)
mpl.plot(a, phi_blocmin,"b", linewidth = 0.75, label = r"$\phi_{\infty}:$" + labelstring1)
mpl.plot(a, phi_comp1,"r--", linewidth = 0.75, label = r"$\phi_c:$" + labelstring2)
mpl.plot(a, phi_bcomp1,"r", linewidth = 0.75, label = r"$\phi_{\infty}:$" + labelstring2)
mpl.plot(a, phi_comp2,"g--", linewidth = 0.75, label = r"$\phi_c:$" + labelstring3)
mpl.plot(a, phi_bcomp2,"g", linewidth = 0.75, label = r"$\phi_{\infty}:$" + labelstring3)

mpl.legend()
mpl.xlabel("a")
mpl.yscale("log")
mpl.ylabel(r"$\phi$", rotation = 0, labelpad = 10)
mpl.title(titlestring)
mpl.savefig("Figures\Chameleon\Phisforcomp_fvar2.pdf", addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
mpl.clf()
"""
"""
mpl.plot(a, Phi_newtcham_locmin, "b--", linewidth = 0.75, label = r"$M_{min} = 10^{%.2f} \mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mlocalmin))
mpl.plot(a, Phi_newtcham_comp1, "r-", linewidth = 0.75, label = r"$M = 10^{%.2f} \mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mcomp1))
mpl.plot(a, Phi_newtcham_comp2, "g:", linewidth = 0.75, label = r"$M = 10^{%.2f} \mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mcomp2))
mpl.legend()
mpl.xlabel("a")
mpl.yscale("log")
mpl.title(r"$\mathrm{M} = 10^{%.2f}\mathrm{M}_{\odot}\mathrm{h}^{-1}$" % np.log10(Mlocalmin))
mpl.ylabel(r"$\Phi_{N}$", rotation = 0, labelpad = 10)
mpl.savefig("Figures\Chameleon\Newtpotforcomp2.pdf", bbox_inches = "tight")
mpl.clf()
"""