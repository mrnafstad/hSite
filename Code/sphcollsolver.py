"""
cd c:\Users\Halvor\Documents\Master\Kode\
python betternewvir.py

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
import time

from numba import cuda


"""
H_0 = 3.44e-35				#eV
G = 6.71e-57				#eV^-2
M_pl = 2.44e27				#eV
rho_c0 = 3*H_0**2/(8*np.pi*G)
Msun = 2e30*5.61e26
GtimeM_sun = G*Msun
"""
#in natural (Planck) units:
H_0 = 1.491e-32 #9.369e-32										#eV 
G = 1.0												#1

Msun = 1.989e30*3.755e-21 #5.977e-22 							#eV-1
M_pl = 4.341e-9*3.755e-21 #5.977e-22  
#M_pl = 1/np.sqrt(8*np.pi)
GtimeM_sun = G*Msun
rho_c0 = 3.*H_0**2*M_pl**2# /(8.*np.pi*G)
"""

hbar = 6.58e-16
c = 3e8
G = 6.71e-39										#GeV^-2, revisit the value of G in NU
Msun = 1.99e30*5.61e26								#GeV, should be 1.12e66
M_pl = 2.44e18										#GeV
H_0 = 70/(3.09e19*5.068e-39)/1.519e24							#GeV
rho_c0 = 3.*H_0**2*M_pl**2
GtimeM_sun = G*Msun
"""
print H_0, G, rho_c0, Msun, M_pl   

def vircheck(R, y, Omega_m0, Lambda, delta_i, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file):

	a = np.exp(y)

	rad = R[:, 0]
	drdy = R[:, 1]
	ddrddy = np.zeros(len(a))
	T = np.zeros(len(a))
	W = np.zeros(len(a))
	
	for i in range(len(a)):
		ddrddy[i] = acceleration( Omega_m0, Lambda, delta_i, a[i], rad[i], drdy[i], E, gamma, beta, M, f_R0, delta_rho, n )
		T[i] = kinetic(drdy[i])
		W[i] = potential(rad[i], drdy[i], ddrddy[i], a[i], E, Omega_m0, Lambda, gamma, beta, delta_rho, M, delta_i, n, f_R0)
	s = j = np.argmax(rad)
	W_ta = potential(rad[s], drdy[s], ddrddy[s], a[s], E, Omega_m0, Lambda, gamma, beta, delta_rho, M, delta_i, n, f_R0 )
	

	

	vir = False
	
	while ( j < len(rad) - 1 and rad[j] >= 0):

		if abs( 2*T[j] + W[j] ) <= 1e-4:

			i = k = j
			vir = True

		elif ( 2*T[j] < abs(W[j]) ):

			i = k = j
			vir = True

		j += 1
	
	if not vir:
		i = j
		print "It did not work"

	t = 0
	while rad[t] > 0:
		t += 1

	if vir:
		while i < len(rad):
			rad[i] = rad[k]
			i += 1


		odensity = (Omega_m0*(1+delta_i)/rad[k]**3)/(Omega_m0/a[k]**3)
		odensitymax = (Omega_m0*(1+delta_i)/rad[s]**3)/(Omega_m0/a[s]**3)
		odensitycoll = (Omega_m0*(1+delta_i)/rad[k]**3)/(Omega_m0/a[t]**3)

		file.write("{:1.15f} \n".format(odensity))
		file.write("{:1.15f} \n".format(odensitymax))
		file.write("{:1.15f} \n".format(odensitycoll))

		rviroverrta = rad[k]/rad[s]
		aviroverata = a[k]/a[s]

		file.write("{:1.15f} \n".format(rviroverrta))
		file.write("{:1.15f} \n".format(aviroverata))
		file.write("{:1.15f} \n".format(a[k]))
		file.write("{:1.15f} \n".format(a[t]))
	else:
		file.write("These parameters do not lead to collapse. \n")
	return rad, T, W, s, k


def kinetic( drdy ):
	return 3./10.*drdy**2

def potential( radius, dotR, ddotR, a, E, Omega_m0, Lambda, gamma, beta, delta_rho, M, delta_i, n, f_R0 ):
	if gamma == 0:
		W = 3./5.*radius*(-3*Omega_m0/(2*a**3*E(Omega_m0, Lambda, a, gamma, beta, n, M, f_R0, delta_rho, radius))*dotR + ddotR)

	else:
		W = 3./5.*radius*(-3*gammaBack( Omega_m0, Lambda, 0, beta, n, M, f_R0, delta_rho, a, a, True )*Omega_m0/(2*a**3*E(Omega_m0, Lambda, a, gamma, beta, n, M, f_R0, delta_rho, radius))*dotR + ddotR)
	return W


def r( x, y, acceleration, Omega_m0, Lambda, delta_i, E, gamma, beta, M, f_R0, delta_rho, n ):
	#input function to be used by odeint. Generic form for a seccond order ode
	rad = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = acceleration( Omega_m0, Lambda, delta_i, a, rad, drdy, E, gamma, beta, M, f_R0, delta_rho, n )	

	return rr



def findcoll(tolerance, acceleration, model, E, gamma, beta, M, f_R0, delta_rho, n, y0, plot, file):

	#Set initial conditions and time array
	N = 5000000

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	a = np.exp(y)

	tol = np.log(1/(1+tolerance))

	#Set density parameters according to argument "model"
	if model == "EdS":
		Omega_m0 = 1.0
		Lambda = 0.0
	else:
		Omega_m0 = 0.25
		Lambda = 0.75

	collmax = False 
	collmin = False 
	collmid = False

	colltime_max = 10

	delta_max = 0.001
	delta_min = 0.00000000001
	j = 0
	c = 0

	while ( abs(colltime_max) >= abs(tol) ):

		#set bisection point
		delta_mid = (delta_max + delta_min)/2.0

		# solve for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_max, E, gamma, beta, M, f_R0, delta_rho, n) )


		#solve for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_min, E, gamma, beta, M, f_R0, delta_rho, n) )


		#solve for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_mid, E, gamma, beta, M, f_R0, delta_rho, n) )

		for i in range(len(radiusmax[:,0])):
			if radiusmax[i,0] <= 0:
				colltime_max = y[i]
				collmax = True
				#assert np.isnan(radiusmax[i,0]), "radius is %.3e, delta_max = %.3e" % (radiusmax[i,0], delta_max)
				#print "max"
				break
			else:
				collmax = False

		for i in range(len(radiusmid[:,0])):
			if radiusmid[i,0] <= 0:
				colltime_mid = y[i]
				collmid = True
				#print "mid"
				break
			else:
				collmid = False

		for i in range(len(radiusmin[:,0])):
			if radiusmin[i,0] <= 0:
				colltime_min = y[i]
				collmin = True
				#print "min"
				break
			else:
				collmin = False


		if ( collmax and collmid ):
			#check wether deltamax and deltamid gives collapse
			#sets new delta_max to delta_mid for next iteration
			delta_max = delta_mid


		elif ( not collmid ):
			#checks wether deltamax gives collapse, but deltamid does not
			#sets new deltamin to deltamid for next iteration
			delta_min = delta_mid

		elif ( not collmax):
			print "Increase delta_max"
			break

		if ( collmin ):
		
			print "Well damn.."
			#delta_max = 1e-15
			delta_next = delta_min*1e-7
			nextmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_next, E, gamma, beta, M, f_R0, delta_rho, n) )
			for i in range(len(radiusmin[:,0])):
				if radiusmin[i,0] <= 0:
					nextmintime = y[i]
					nextcoll = True
					#print "min"
					break
				else:
					nextcoll = False

			if abs(nextmintime) > abs(colltime_min):
				delta_min *= 2

			elif ( nextcoll == False or abs(nextmintime) < colltime_min): 
				delta_min = delta_next
				c += 1
				if c > 10:
					print colltime_min
					break
			
			

		#set x_coll, which is returned from findcoll()
		x_coll = colltime_max
		"""
		collmax = False 
		collmin = False 
		collmid = False
		"""
		

		j += 1

		if j > 200:
			print "This may be infinite"
			break


	ct = np.exp(-x_coll) -1
	
	#file.write("{:1.7e} \n".format(ct))
	print delta_max

	R = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_max, E, gamma, beta, M, f_R0, delta_rho, n ) )
	k = np.argmin(R[:,0])
	#print k, R[k,0]


	#mpl.plot(y, R[:,0], linewidth = 0.75, label = model)
	

	rvir, T, W, ta_loc, vir_loc = vircheck( R, y, Omega_m0, Lambda, delta_max, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file )
	
	if plot:
		if model == "Chameleon":
			mpl.plot(a, rvir, "b-", linewidth = 0.75, label = "Cham, $M = 10^{%.0f} \mathrm{M}_{\odot} \mathrm{h}^{-1}$, \n $|f_{\mathcal{R}_0} - 1 |= 10^{%.0f}$, \n $\delta_\mathrm{i}$ = %1.5e" % (np.log10(M/Msun*0.7), np.log10(abs(f_R0 - 1)), delta_max))
			mpl.plot([a[ta_loc]], [rvir[ta_loc]], "bo", markersize = 5, label = r"$\mathrm{Cham}_{\mathrm{ta}}$")
			mpl.plot([a[vir_loc]], [rvir[vir_loc]], "b.", markersize = 5, label = r"$\mathrm{Cham}_{\mathrm{vir}}$")
		if model == "Symmetron":
			mpl.plot(a, rvir, "r-", linewidth = 0.75, label = "Symm, $M = 10^{%.0f} \mathrm{M}_{\odot} \mathrm{h}^{-1}$, \n $z_{\mathrm{ssb}} = %.1f, L = %.0f \mathrm{Mpc} \mathrm{h}^{-1}$, \n $\delta_\mathrm{i}$ = %1.5e" % (np.log10(M/Msun*0.7), f_R0, n/3.09e22*1.97e-7*0.7, delta_max))
			mpl.plot([a[ta_loc]], [rvir[ta_loc]], "ro", markersize = 5, label = r"$\mathrm{Symm}_{\mathrm{ta}}$")
			mpl.plot([a[vir_loc]], [rvir[vir_loc]], "r.", markersize = 5, label = r"$\mathrm{Symm}_{\mathrm{vir}}$")
		elif model == "LCDM":
			mpl.plot(a, rvir, "k-", linewidth = 0.75, label = r"$\Lambda\mathrm{CDM}, \delta_\mathrm{i}$ = %1.5e" % delta_max)
			mpl.plot([a[ta_loc]], [rvir[ta_loc]], "ko", markersize = 5, label = r"$\Lambda\mathrm{CDM}_{\mathrm{ta}}$")
			mpl.plot([a[vir_loc]], [rvir[vir_loc]], "k.", markersize = 5, label = r"$\Lambda\mathrm{CDM}_{\mathrm{vir}}$")
		elif model == "EdS":
			mpl.plot(a, rvir, "g-", linewidth = 0.75, label = r"EdS, $\delta_\mathrm{i}$ = %1.5e" % delta_max)
			mpl.plot([a[ta_loc]], [rvir[ta_loc]], "go", markersize = 5, label = r"$\mathrm{EdS}_{\mathrm{ta}}$")
			mpl.plot([a[vir_loc]], [rvir[vir_loc]], "g.", markersize = 5, label = r"$\mathrm{EdS}_{\mathrm{vir}}$")
	file.write("{:1.9f} \n".format(delta_max))
	if plot:
		return T, W, y, delta_max, ta_loc, vir_loc
	else:
		return rvir, delta_max

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
	#assert r < 0, "Unphysical radius in Phi_N for %.2e" % delta_i
	#Phi_N = ( (GtimeM_sun*M*H_0)**(2./3.)*(Omega_m0*(1 + delta_i)/2.)**(1./3.)/r )
	#Phi_N = (M/(np.sqrt(3)*np.pi))**(2./3.)*(2*rho_c0*Omega_m0*(1 + delta_i))**(1./3.)/(M_pl**2*r)
	#Phi_N = ( 2*M*np.sqrt(Omega_m0*(1 + delta_i))/H_0**2 )**(2./3.)*rho_c0/(6*M_pl**2*r)
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

		


first = True
values = open("Numbers\Vals.txt", "w")
def gammaHuSawicki1( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, perturbation ):
	
	deltRoverR = DROR( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho, perturbation )

	Rs = 1 - deltRoverR
	return 1 + 2*beta**2*( 1 - Rs**3)


def gamma_symm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, delta_rho, r, a, perturbation ):
	#In comparison to chameleon, g -> beta, z_ssb -> f_R0, L -> n notationally

	deltaR_overR = DRORsymm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, r, a, delta_rho, perturbation )

	
	Rs = 1 - deltaR_overR
	return 1 + 2*beta**2*( 1 - Rs**3)

def controll( model1, model2, E1, E2, Gamma, beta, M, f_R0, n, delta_rho, delta_i ):

	N = 5000000
	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)
	a = np.exp(y)

	f1 = odeint( r, [r0, drdx0], y, args = ( model1, Omega_m0, Lambda, delta_i, E1, Gamma, beta, M, f_R0, delta_rho, n ) )
	f2 = odeint( r, [r0, drdx0], y, args = ( model2, Omega_m0, Lambda, delta_i, E2, Gamma, beta, M, f_R0, delta_rho, n ) )

	diff = f1[:,0] - f2[:,0]

	mpl.plot(a, diff, "-.", linewidth = 0.75, label = r"$R_{\Lambda CDM} - R_{Cham}$ with $\beta$ = %.2f" % beta)
	mpl.xlabel("a")
	mpl.ylabel("Diff")
	mpl.xscale("log")
	return

def placement( filename, label1, seccondvar, rootfinding ):

	file = open(filename, "r")
	d_vir = []
	d_ta = []
	d_c = []
	rvirrta = []
	avirata = []
	variable = []
	avir = []
	acoll = []
	if rootfinding:
		delta_c = []
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
		elif ( i == 7 and not rootfinding ):
			acoll.append(value)
			i = 0

		elif ( i == 7 and rootfinding):
			acoll.append(value)
			i += 1

		elif i == 8:
			delta_c.append(value)
			i = 0
		

	file.close()
	
	mpl.figure(1)
	mpl.plot(variable, d_vir, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\Delta_{vir}$", rotation = 0)
	mpl.xscale("log")
	
	mpl.figure(2)
	mpl.plot(variable, d_ta, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\Delta_{ta}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(3)
	mpl.plot(variable, rvirrta, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(4)
	mpl.plot(variable, avirata, linewidth = 0.75, label = seccondvar)
	
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\frac{a_{vir}}{a_{ta}}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(5)
	mpl.plot(variable, avir, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	
	mpl.ylabel(r"$a_{vir}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(6)
	mpl.plot(variable, d_c, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\Delta_{c}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(7)
	mpl.plot(variable, acoll, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$a_c$", rotation = 0)
	mpl.xscale("log")

	if rootfinding:
		mpl.figure(8)
		mpl.plot(variable, delta_c, linewidth = 0.75, label = seccondvar)
		mpl.legend()
		mpl.xlabel(label1)
		mpl.ylabel(r"$\delta_i$", rotation = 0)
		mpl.xscale("log")
		mpl.yscale("log")

	return 


def makefigs(Name):
	mpl.figure(1)
	mpl.xlabel("a")
	mpl.ylabel(r"$\tilde{R}$", rotation = 0, labelpad = 10)
	mpl.legend()
	mpl.savefig("Figures/Evolution" + Name, addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")

	mpl.figure(2)
	mpl.yscale("log")
	mpl.xlabel("a")
	mpl.ylabel(r"$\bar{E}$", rotation = 0, labelpad = 10)
	mpl.legend()
	mpl.savefig("Figures\Energies" + Name, addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
	mpl.clf()


	mpl.figure(3)
	#mpl.yscale("log")
	mpl.xlabel("a")
	mpl.ylabel(r"$\frac{\bar{K}}{|\bar{W}|}$", rotation = 0, labelpad = 10)
	mpl.legend()
	mpl.savefig("Figures\Relative_Energies" + Name, addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
	mpl.clf()


tolerance = 0.0001
y0 = np.log(1e-4)
a_i = np.exp(y0)
N = 5000000
y = np.linspace(y0, -1e-15, N)
r0 = np.exp(y0)
drdx0 = np.exp(y0)
a = np.exp(y)



#findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, y0)
#findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, y0)


M = 1e14*Msun/0.7


Omega_m0 = 0.25
Lambda = 0.75

print "Working on LCDM"
fileLCDM = open("Numbers\LCDMviracc.txt", "w")
mpl.figure(1)
T1, W1, y, d_LCDM, LCDM_ta, LCDM_vir = findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, 0, 0, 0, 0, 0, 0, y0, True, fileLCDM)
fileLCDM.close()

mpl.figure(2)
mpl.plot(a, T1, "k-", linewidth = 0.75, label = r"$\bar{K}_{\Lambda \mathrm{CDM}}$")
mpl.plot(a, abs(W1), "k-.", linewidth = 0.75, label = r"$|\bar{W}_{\Lambda \mathrm{CDM}}|$")
mpl.axvline(x = a[LCDM_ta], ymin = 0.2, ymax = 0.7, color = "k", linestyle = "--", linewidth = 0.75,label = r"$\Lambda\mathrm{CDM}_{\mathrm{ta}}$")
mpl.axvline(x = a[LCDM_vir], ymin = 0.2, ymax = 0.7, color = "k", linestyle = ":", linewidth = 0.75, label = r"$\Lambda\mathrm{CDM}_{\mathrm{vir}}$")

mpl.figure(3)
relLCDM = T1/abs(W1)
mpl.plot(a, relLCDM, "k-", linewidth = 0.75, label = r"$\frac{\bar{K}_{\Lambda \mathrm{CDM}}}{|\bar{W}_{\Lambda \mathrm{CDM}}|}$")
mpl.axvline(x = a[LCDM_ta], ymin = relLCDM[LCDM_ta] - 0.1, ymax = relLCDM[LCDM_ta] + 0.1, color = "k", linestyle = "--", linewidth = 0.75,label = r"$\Lambda\mathrm{CDM}_{\mathrm{ta}}$")
mpl.axvline(x = a[LCDM_vir], ymin = relLCDM[LCDM_vir] - 0.1, ymax = relLCDM[LCDM_vir] + 0.1, color = "k", linestyle = ":", linewidth = 0.75, label = r"$\Lambda\mathrm{CDM}_{\mathrm{vir}}$")




print "Working on EdS"
fileEdS = open("Numbers\EdSvirracc.txt", "w")
mpl.figure(1)
T2, W2, y, d_EdS, EdS_ta, EdS_vir = findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, 0, 0, 0, 0, 0, 0, y0, True, fileEdS)
fileEdS.close()

mpl.figure(2)
mpl.plot(a, T2, "g-", linewidth = 0.75, label = r"$\bar{K}_{\mathrm{EdS}}$")
mpl.plot(a, abs(W2), "g-.", linewidth = 0.75, label = r"$|\bar{W}_{\mathrm{EdS}}|$")
mpl.axvline(x = a[EdS_ta], ymin = 0.2, ymax = 0.7, color = "c", linestyle = "--", linewidth = 0.75, label = r"$\mathrm{EdS}_{\mathrm{ta}}$")
mpl.axvline(x = a[EdS_vir], ymin = 0.2, ymax = 0.7, color = "c", linestyle = ":", linewidth = 0.75, label = r"$\mathrm{EdS}_{\mathrm{vir}}$")

mpl.figure(3)
relEdS = T2/abs(W2)
mpl.plot(a, relEdS, "g-", linewidth = 0.75, label = r"$\frac{\bar{K}_{\mathrm{EdS}}}{|\bar{W}_{\mathrm{EdS}}|}$")
mpl.axvline(x = a[EdS_ta], ymin = relEdS[EdS_ta] - 0.1, ymax = relEdS[EdS_ta] + 0.1, color = "g", linestyle = "--", linewidth = 0.75, label = r"$\mathrm{EdS}_{\mathrm{ta}}$")
mpl.axvline(x = a[EdS_vir], ymin = relEdS[EdS_vir] - 0.1, ymax = relEdS[EdS_vir] + 0.1, color = "g", linestyle = ":", linewidth = 0.75, label = r"$\mathrm{EdS}_{\mathrm{vir}}$")

makefigs("_EdS_LCDM.pdf")




beta = 1/np.sqrt(6)
f_R0 = 1 + 1e-3
n = 1.0
filecham = open("Numbers\Chameleon10.txt", "w")
print "Working on chameleon"
mpl.figure(1)
Tcham, Wcham, y, d_cham, cham_ta, cham_vir = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M, f_R0, d_rho, 1.0, y0, True, filecham)
filecham.close()

mpl.figure(2)
mpl.plot(a, Tcham, "b-", linewidth = 0.75, label = r"$\bar{K}_{\mathrm{Cham}}$")
mpl.plot(a, abs(Wcham), "b-.", linewidth = 0.75, label = r"$|\bar{W}_{\mathrm{Cham}}|$")
mpl.axvline(x = a[cham_ta], ymin = 0.2, ymax = 0.7, color = "b", linestyle = "--", linewidth = 0.75,label = r"$\mathrm{Cham}_{\mathrm{ta}}$")
mpl.axvline(x = a[cham_vir], ymin = 0.2, ymax = 0.7, color = "b", linestyle = ":", linewidth = 0.75, label = r"$\mathrm{Cham}_{\mathrm{vir}}$")
mpl.plot(a, T1, "k-", linewidth = 0.75, label = r"$\bar{K}_{\Lambda \mathrm{CDM}}$")
mpl.plot(a, abs(W1), "k-.", linewidth = 0.75, label = r"$|\bar{W}_{\Lambda \mathrm{CDM}}|$")
mpl.axvline(x = a[LCDM_ta], ymin = 0.2, ymax = 0.7, color = "k", linestyle = "--", linewidth = 0.75,label = r"$\Lambda\mathrm{CDM}_{\mathrm{ta}}$")
mpl.axvline(x = a[LCDM_vir], ymin = 0.2, ymax = 0.7, color = "k", linestyle = ":", linewidth = 0.75, label = r"$\Lambda\mathrm{CDM}_{\mathrm{vir}}$")
mpl.plot(a, T2, "g-", linewidth = 0.75, label = r"$\bar{K}_{\mathrm{EdS}}$")
mpl.plot(a, abs(W2), "g-.", linewidth = 0.75, label = r"$|\bar{W}_{\mathrm{EdS}}|$")

mpl.figure(3)
relCham = Tcham/abs(Wcham)
mpl.plot(a, relCham, "b-", linewidth = 0.75, label = r"$\frac{\bar{K}_{\mathrm{Cham}}}{|\bar{W}_{\mathrm{cham}}|}$")
mpl.plot(a, relEdS, "g-", linewidth = 0.75, label = r"$\frac{\bar{K}_{\mathrm{EdS}}}{|\bar{W}_{\mathrm{EdS}}|}$")
mpl.plot(a, relLCDM, "k-", linewidth = 0.75, label = r"$\frac{\bar{K}_{\Lambda \mathrm{CDM}}}{|\bar{W}_{\Lambda \mathrm{CDM}}|}$")

mpl.axvline(x = a[LCDM_ta], ymin = relLCDM[LCDM_ta] - 0.1, ymax = relLCDM[LCDM_ta] + 0.1, color = "k", linestyle = "--", linewidth = 0.75,label = r"$\Lambda\mathrm{CDM}_{\mathrm{ta}}$")
mpl.axvline(x = a[LCDM_vir], ymin = relLCDM[LCDM_vir] - 0.1, ymax = relLCDM[LCDM_vir] + 0.1, color = "k", linestyle = ":", linewidth = 0.75, label = r"$\Lambda\mathrm{CDM}_{\mathrm{vir}}$")
mpl.axvline(x = a[cham_ta], ymin = relCham[cham_ta] - 0.1, ymax = relCham[cham_ta] + 0.1, color = "b", linestyle = "--", linewidth = 0.75,label = r"$\mathrm{Cham}_{\mathrm{ta}}$")
mpl.axvline(x = a[cham_vir], ymin = relCham[cham_vir] - 0.1, ymax = relCham[cham_vir] + 0.1, color = "b", linestyle = ":", linewidth = 0.75, label = r"$\mathrm{Cham}_{\mathrm{vir}}$")

makefigs("_EdS_LCDM_Cham.pdf")




print "On Symmetron"
beta = .5
z_ssb = 1.0
L = 1.

filesymm = open("Numbers\Symmetronacc.txt", "w")
mpl.figure(1)
Tsymm, Wsymm, y, d_symm, sym_ta, sym_vir = findcoll(tolerance, chameleonacc, "Symmetron", Echamnorad, gamma_symm, beta, M, z_ssb, d_rho, L*3.09e22/1.97e-7/0.7, y0, True, filesymm)
filesymm.close()

mpl.figure(2)
mpl.plot(a, Tsymm, "r-", linewidth = 0.75, label = r"$\bar{K}_{\mathrm{Symm}}$")
mpl.plot(a, abs(Wsymm), "r-.", linewidth = 0.75, label = r"$|\bar{W}_{\mathrm{Symm}}|$")
mpl.axvline(x = a[sym_ta], ymin = 0.2, ymax = 0.7, color = "r", linestyle = "--", linewidth = 0.75,label = r"$\mathrm{Symm}_{\mathrm{ta}}$")
mpl.axvline(x = a[sym_vir], ymin = 0.2, ymax = 0.7, color = "r", linestyle = ":", linewidth = 0.75, label = r"$\mathrm{Symm}_{\mathrm{vir}}$")
mpl.plot(a, T1, "k-", linewidth = 0.75, label = r"$\bar{K}_{\Lambda \mathrm{CDM}}$")
mpl.plot(a, abs(W1), "k-.", linewidth = 0.75, label = r"$|\bar{W}_{\Lambda \mathrm{CDM}}|$")
mpl.axvline(x = a[LCDM_ta], ymin = 0.2, ymax = 0.7, color = "k", linestyle = "--", linewidth = 0.75,label = r"$\Lambda\mathrm{CDM}_{\mathrm{ta}}$")
mpl.axvline(x = a[LCDM_vir], ymin = 0.2, ymax = 0.7, color = "k", linestyle = ":", linewidth = 0.75, label = r"$\Lambda\mathrm{CDM}_{\mathrm{vir}}$")
mpl.plot(a, T2, "g-", linewidth = 0.75, label = r"$\bar{K}_{\mathrm{EdS}}$")
mpl.plot(a, abs(W2), "g-.", linewidth = 0.75, label = r"$|\bar{W}_{\mathrm{EdS}}|$")


mpl.figure(3)
relsymm = Tsymm/abs(Wsymm)
mpl.plot(a, relsymm, "r-.", linewidth = 0.75, label = r"$\frac{\bar{K}_{\mathrm{Cham}}}{|\bar{W}_{\mathrm{cham}}|}$")
mpl.plot(a, relEdS, "g:", linewidth = 0.75, label = r"$\frac{\bar{K}_{\mathrm{EdS}}}{|\bar{W}_{\mathrm{EdS}}|}$")
mpl.plot(a, relLCDM, "k-", linewidth = 0.75, label = r"$\frac{\bar{K}_{\Lambda \mathrm{CDM}}}{|\bar{W}_{\Lambda \mathrm{CDM}}|}$")

mpl.axvline(x = a[LCDM_ta], ymin = relLCDM[LCDM_ta] - 0.1, ymax = relLCDM[LCDM_ta] + 0.1, color = "k", linestyle = "--", linewidth = 0.75,label = r"$\Lambda\mathrm{CDM}_{\mathrm{ta}}$")
mpl.axvline(x = a[LCDM_vir], ymin = relLCDM[LCDM_vir] - 0.1, ymax = relLCDM[LCDM_vir] + 0.1, color = "k", linestyle = ":", linewidth = 0.75, label = r"$\Lambda\mathrm{CDM}_{\mathrm{vir}}$")
mpl.axvline(x = a[sym_ta], ymin = relsymm[sym_ta] - 0.1, ymax = relsymm[sym_ta] + 0.1, color = "r", linestyle = "--", linewidth = 0.75,label = r"$\mathrm{Symm}_{\mathrm{ta}}$")
mpl.axvline(x = a[sym_vir], ymin = relsymm[sym_vir] - 0.1, ymax = relsymm[sym_vir] + 0.1, color = "r", linestyle = ":", linewidth = 0.75, label = r"$\mathrm{Symm}_{\mathrm{vir}}$")

makefigs("_EdS_LCDM_Cham_Symm.pdf")


"""

masses = np.logspace(10, 18, 30)
f_R0_list = np.logspace(-12, -3, 30) + 1
First = True
diff = 0
for M in masses:
	for f_R0 in f_R0_list:
		
		filename = "Numbers\Chameleon\M%.15ff_R0%.15f.txt" % (np.log10(M), np.log10(abs(1-f_R0)))
	
		try:
			file = open("Numbers\Chameleon\RadiusM%.15ff_R0%.15f.txt" % (np.log10(M), np.log10(abs(1-f_R0))), "r")
			#f = open(filename, "r")
			print "file exists"

		except IOError:
			t1 = time.time()
			Mloop = open(filename, "w")

			t21 = time.time()
			print "1 - f_R0 = 10^ %.15f" %  np.log10(abs(1 - f_R0))
			Mloop.write("{:1.1f} \n".format(M))
			R, d_cham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M*Msun/0.7, f_R0, d_rho, 1.0, y0, False, Mloop)
			Mloop.close()
			file = open("Numbers\Chameleon\RadiusM%.15ff_R0%.15f.txt" % (np.log10(M), np.log10(abs(1-f_R0))), "w")
			for i in range(len(R)):
				file.write("{:1.15f} \n".format(R[i]))
			file.close()
			t2 = time.time()
			diff += abs(t2 - t1)
			if First:
				ETT = diff*len(masses)*len(f_R0_list)
				First = False
			ETRM = (ETT - diff)/60.
			ETRH = ETRM/60.
			print "Estimated time remaining: %.5e minutes or %.5f hours" % (ETRM, ETRH)




First = True
diff = 0
masses = [1e14] #np.logspace(10, 18, 30)
beta_list = [0.5]# np.linspace(0, 1.5, 30)
L_list = np.logspace(-1, 1, 30)
z_ssb_list = np.linspace(0, 3, 30)
t = 1
for M in masses:
	for z_ssb in z_ssb_list:
		for beta in beta_list:
			for L in L_list:
		
				filename = "Numbers\Symmetron\Sym2\M%.15fz_ssb%.3fL%.3fbeta%.3f.txt" % (np.log10(M), z_ssb, L, beta)
			
				try:
					f = open(filename, "r")
					print "file exists"
					t += 1
					

				except IOError:
					t1 = time.time()
					Mloop = open(filename, "w")

					t21 = time.time()
					print "{:.15f}   {:1.2f}   {:1.2f}   {:1.2f}".format(np.log10(M), z_ssb, L, beta)
					Mloop.write("{:1.1f} \n".format(M))
					R, d_cham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gamma_symm, beta, M*Msun/0.7, z_ssb, d_rho, L*3.09e22/1.97e-7/0.7, y0, False, Mloop)
					Mloop.close()
					file = open("Numbers\Symmetron\Sym2\RadiusM%.15fz_ssb%.15fL%.2fbeta%.2f.txt" % (np.log10(M), z_ssb, L, beta), "w")
					for i in range(len(R)):
						file.write("{:1.15f} \n".format(R[i]))
					file.close()
					t2 = time.time()
					diff += abs(t2 - t1)
					k = diff/t
					if First:
						ETT = diff*len(masses)*len(z_ssb_list)*len(L_list)*len(beta_list)
						First = False
					ETA = k*len(masses)*len(z_ssb_list)*len(L_list)*len(beta_list)
					ETRM = (ETT - diff)/60.
					ETRH = ETRM/60.
					print "Estimated time remaining: %.5e minutes or %.5f hours" % (ETRM, ETRH)
					print "ETA: %.5f" % ((ETA-diff)/60/60)
					t += 1
"""