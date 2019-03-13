"""
cd c:\Users\Halvor\Documents\Master\Kode\
python combinedfields.py

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
	return rad, T, W


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
	

	rvir, T, W = vircheck( R, y, Omega_m0, Lambda, delta_max, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file )
	
	if plot:
		if model == "Chameleon":
			mpl.plot(a, rvir, "-.", linewidth = 0.75, label = "Chameleon, $M = 10^{%.0f} M_{\odot} h^{-1}$, \n$f_{R0} - 1 = 10^{%.0f}$, $\delta_i$ = %1.5e" % (np.log10(M/Msun*0.7), np.log10(abs(f_R0 - 1)), delta_max))
			#mpl.plot(y, R[:,0], "-.", linewidth = 0.75, label = "Chameleon no vir")
			#Rmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_min, E, gamma, beta, M, f_R0, delta_rho, n ) )
			#rvirmin, Tmin, Wmin = vircheck( Rmin, y, Omega_m0, Lambda, delta_min, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file )
			#mpl.plot(y, rvirmin, "-.", linewidth = 0.75, label = r"Chameleon, M = %1.0e, $f_{R0} - 1$ = %1.0e, $\delta_i$ = %1.5e" % (M, f_R0 - 1, delta_min))
		if model == "Symmetron":
			mpl.plot(a, rvir, "-", linewidth = 0.75, label = "Symmetron, $M = 10^{%.0f} M_{\odot} h^{-1}$, \n$z_{ssb} = %.1f, L = %.0f Mpc h^{-1}, \delta_i$ = %1.5e" % (np.log10(M/Msun*0.7), f_R0, n/3.09e22*1.97e-7*0.7, delta_max))
		elif model == "LCDM":
			mpl.plot(a, rvir, "--", linewidth = 0.75, label = r"$\Lambda$CDM, $\delta_i$ = %1.5e" % delta_max)
		elif model == "EdS":
			mpl.plot(a, rvir, ":", linewidth = 0.75, label = r"EdS$ \delta_i$ = %1.5e" % delta_max)
	file.write("{:1.7f} \n".format(delta_max))
	if plot:
		return T, W, y, delta_max
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

		


def gammaHuSawicki1( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, perturbation ):
	
	deltRoverR = DROR( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho, perturbation )

	Rs = 1 - deltRoverR
	return 1 + 2*beta**2*( 1 - Rs**3)


def gamma_symm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, delta_rho, r, a, perturbation ):
	#In comparison to chameleon, g -> beta, z_ssb -> f_R0, L -> n notationally

	deltaR_overR = DRORsymm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, r, a, delta_rho, perturbation )

	
	Rs = 1 - deltaR_overR
	return 1 + 2*beta**2*( 1 - Rs**3)

tolerance = 0.0001
y0 = np.log(1e-4)

a_i = np.exp(y0)
N = 5000000
y = np.linspace(y0, -1e-15, N)
a = np.exp(y)
r0 = np.exp(y0)
drdx0 = np.exp(y0)
Omega_m0 = 0.25
Lambda = 0.75

betacham = 1./np.sqrt(6)
ncham = 1.0
f_R0 = 1e-3 + 1

betasymm = 0.5
z_ssb = 1.0
L = 3.09e22/1.97e-7/0.7
rho_ssb = rho_c0*Omega_m0*(1 + z_ssb)**3
delta_rho_ssb = Omega_m0*(1 + z_ssb)**3

M = 1e14*Msun/0.7

try:
	f1 = open("Numbers\Chameleon\Forfigs.txt", "r")
	file1 = open("Numbers\Chameleon\Comb.txt", "r")
	res1 = file1.readlines()
	r = f1.readlines()
	rcham = np.array(r, dtype=np.float32)
	d_cham = float(res1[-1])
	f1.close()
	file1.close()

	file2 = open("Numbers\Symmetron\Comb.txt", "r")
	res2 = file2.readlines()
	d_symm = float(res2[-1])
	f2 = open("Numbers\Symmetron\Forfigs.txt", "r")
	r = f2.readlines()
	rsymm =	np.array(r, dtype=np.float32)
	f2.close()
	file2.close()

except IOError:
	file1 = open("Numbers\Chameleon\Comb.txt", "w")
	rcham, d_cham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, betacham, M, f_R0, d_rho, ncham, y0, False, file1)
	file1.close()

	file2 = open("Numbers\Symmetron\Comb.txt", "w")
	rsymm, d_symm = findcoll(tolerance, chameleonacc, "Symmetron", Echamnorad, gamma_symm, betasymm, M, z_ssb, d_rho, L, y0, False, file2)
	file2.close()

	f1 = open("Numbers\Chameleon\Forfigs.txt", "w")
	f2 = open("Numbers\Symmetron\Forfigs.txt", "w")
	for i in range(len(rsymm)):
		f1.write("{:1.15f} \n".format(rcham[i]))
		f2.write("{:1.15f} \n".format(rsymm[i]))
	f1.close()
	f2.close()



rhocham = rho_c0*d_rho(Omega_m0, Lambda, d_cham, rcham)
rhosymm = rho_c0*d_rho(Omega_m0, Lambda, d_symm, rsymm)
delta_rhocham = d_rho(Omega_m0, Lambda, d_cham, rcham)
delta_rhosymm = d_rho(Omega_m0, Lambda, d_symm, rsymm)

delta_rho_b = d_rho(Omega_m0, Lambda, 0, a)

rho_b = rho_c0*d_rho(Omega_m0, Lambda, 0, a)

delta_chamdiff = np.zeros(len(a))
delta_symmdiff = np.zeros(len(a))

d_cont_cham = np.zeros(len(a))
d_cont_symm = np.zeros(len(a))

geff_cham = np.zeros(len(a))
geff_symm = np.zeros(len(a))

phicham = np.zeros(len(a))
phisymm = np.zeros(len(a))

phi_bcham = np.zeros(len(a))
phi_bsymm = np.zeros(len(a))

Phi_newtcham = Phi_N( rcham, M, Omega_m0, d_cham )
Phi_newtsymm = Phi_N( rsymm, M, Omega_m0, d_symm )
first1 = True
first2 = True

for i in range(len(a)):
	geff_cham[i] = gammaHuSawicki1(Omega_m0, Lambda, d_cham, betacham, ncham, M, f_R0, d_rho, rcham[i], a[i], True) - 1
	geff_symm[i] = gamma_symm(Omega_m0, Lambda, d_symm, betasymm, L, M, z_ssb, d_rho, rsymm[i], a[i], True) - 1
	
	phicham[i] = phi_cham(Omega_m0, Lambda, d_cham, betacham, ncham, M, f_R0, rcham[i], a[i], d_rho)
	if delta_rhosymm[i]/delta_rho_ssb < 1:
		phisymm[i] = phi_symm(rhosymm[i], rho_ssb, L, M, betasymm)
		if first1:
			print 1, a[i], rhosymm[i], rho_ssb
			first1 = False


	phi_bcham[i] = phi_cham(Omega_m0, Lambda, 0, betacham, ncham, M, f_R0, a[i], a[i], d_rho)
	if delta_rho_b[i]/delta_rho_ssb < 1:
		phi_bsymm[i] = phi_symm(rho_b[i], rho_ssb, L, M, betasymm)
		if first2:
			print 2, a[i], rho_b[i], rho_ssb
			first2 = False


	delta_chamdiff[i] =  delta_rhocham[i] - delta_rho_b[i] 
	delta_symmdiff[i] =  delta_rhocham[i] - delta_rho_b[i]

	d_cont_symm[i] = delta_rhosymm[i]/delta_rho_b[i]
	d_cont_cham[i] = delta_rhocham[i]/delta_rho_b[i]

gmax_cham = 2*betacham**2
gmax_symm = 2*betasymm**2

phischam = [max(phicham), max(phi_bcham)]
phissymm = [max(phisymm), max(phi_bsymm)]


mpl.plot(a, geff_cham, linewidth = 0.75)
#mpl.plot(a, geff_symm, linewidth = 0.75, label = "Symmetron")
#mpl.legend()
mpl.xlabel("a", fontsize = 15)
mpl.ylabel(r"$\gamma_{\mathrm{eff}}$", rotation = 0, labelpad = 15, fontsize = 15)
mpl.savefig("Figures\Chameleon\Geffspec.pdf", bbox_inches = "tight")
mpl.clf()


mpl.plot(a, geff_symm, linewidth = 0.75)
#mpl.legend()
mpl.xlabel("a", fontsize = 15)
mpl.ylabel(r"$\gamma_{\mathrm{eff}}$", rotation = 0, labelpad = 15, fontsize = 15)
mpl.savefig("Figures\Symmetron\Geffspec.pdf", bbox_inches = "tight")
mpl.clf()
"""

mpl.plot(d_cont_cham, geff_cham/gmax_cham, linewidth = 0.75, label = "Chameleon")

mpl.xlabel(r"$\delta$")
mpl.xscale("log")
mpl.ylabel(r"$\frac{\gamma_{eff}}{\gamma_{max}}$", rotation = 0, labelpad = 15)
mpl.savefig("Figures\Chameleon\Geffrelcont.pdf", bbox_inches = "tight")
mpl.clf()

mpl.plot(d_cont_symm, geff_symm/gmax_symm, linewidth = 0.75, label = "Symmetron")
mpl.xscale("log")
mpl.legend()
mpl.xlabel(r"$\delta$")
mpl.ylabel(r"$\frac{\gamma_{eff}}{\gamma_{max}}$", rotation = 0, labelpad = 15)
mpl.savefig("Figures\Geffrelcont.pdf", bbox_inches = "tight")
mpl.clf()
"""
mpl.plot(a, phicham/max(phischam), linewidth = 0.75, label = r"$\phi_{\mathrm{c}}$")
mpl.plot(a, phi_bcham/max(phischam), linewidth = 0.75, label = r"$\phi_{\infty}$")
mpl.legend(fontsize = 15)
mpl.xlabel("a", fontsize = 15)
mpl.ylabel(r"$\frac{\phi }{ \phi_{\mathrm{max}}}$", rotation = 0, labelpad = 15, fontsize =  15)
#mpl.yscale("log")
mpl.savefig("Figures\Chameleon\scalarfield.pdf", bbox_inches = "tight")
mpl.clf()

mpl.plot(a, phisymm/max(phissymm), linewidth = 0.75, label = r"$\phi_{\mathrm{c}}$")
mpl.plot(a, phi_bsymm/max(phissymm), linewidth = 0.75, label = r"$\phi_{\infty}$")
mpl.legend(fontsize = 15)
mpl.xlabel("a", fontsize = 15)
mpl.ylabel(r"$\frac{\phi }{ \phi_{\mathrm{max}}}$", rotation = 0, labelpad = 15, fontsize = 15)
#mpl.yscale("log")
mpl.savefig("Figures\Symmetron\scalarfield.pdf", bbox_inches = "tight")
mpl.clf()

"""

mpl.plot(Phi_newtcham, geff_cham/gmax_cham, linewidth = 0.75, label = "Chameleon")
mpl.xlabel(r"$\Phi_N$")
mpl.xscale("log")
mpl.yscale("log")
mpl.ylabel(r"$\frac{\gamma_{eff}}{\gamma_{max}}$", rotation = 0, labelpad = 15)
mpl.savefig("Figures\Chameleon\Phi_N_to_gamma.pdf", bbox_inches = "tight")
mpl.clf()

mpl.plot(Phi_newtsymm, geff_symm/gmax_symm, linewidth = 0.75, label = "Chameleon")
mpl.xlabel(r"$\Phi_N$")
mpl.xscale("log")
mpl.yscale("log")
mpl.ylabel(r"$\frac{\gamma_{eff}}{\gamma_{max}}$", rotation = 0, labelpad = 15)
mpl.savefig("Figures\Symmetron\Phi_N_to_gamma.pdf", bbox_inches = "tight")
mpl.clf()


mpl.plot(a, Phi_newtcham, linewidth = 0.75, label = "Chameleon")
mpl.plot(a, Phi_newtsymm, linewidth = 0.75, label = "Symmetron")
mpl.legend()
mpl.xlabel("a")
mpl.ylabel(r"$\Phi_N$")
mpl.yscale("log")
mpl.savefig("Figures\Phi_N_comb.pdf", bbox_inches = "tight")
mpl.clf()

mpl.plot(a, Phi_newtcham, linewidth = 0.75, label = "Symmetron")
mpl.xlabel("a")
mpl.ylabel(r"$\Phi_N$")
mpl.yscale("log")
mpl.savefig("Figures\Chameleon\Phi_N_spec.pdf", bbox_inches = "tight")
mpl.clf()

mpl.plot(a, Phi_newtsymm, linewidth = 0.75, label = "Symmetron")
mpl.xlabel("a")
mpl.ylabel(r"$\Phi_N$")
mpl.yscale("log")
mpl.savefig("Figures\Symmetron\Phi_N_spec.pdf", bbox_inches = "tight")
mpl.clf()
"""