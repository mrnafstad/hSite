from matplotlib.pylab import *
from numpy import *

def bolus( bolustime, units ):
	time  = bolustime*4
	#print time
	for i in range(16):
		b = time + i
		bolus_in_body[b] += Humalog( units, i )
		#print in_body[b]

def basal( basaltime, units ):
	time  = basaltime
	print time
	for i in range(24):

		basal_in_body[i] += Insulatard( units, i )
		#in_body[b] += Insulatard( units, i )

def total( basal_in_body, bolus_in_body ):
	for i in range(len(bolus_in_body)):
		if i/4 == int:
			in_body = basal_in_body[i]
		
		

def Humalog( units, time ):
	humalog_profile = array([ 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.425, 0.35, 0.275, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.0])/48.8
	a = units*humalog_profile[time]
	#print time/4.
	return a

def Insulatard( units, time):
	insulatard_profile = array([1.0, 0.99, 0.95, 0.85, 0.8, 0.7, 0.6, 0.5, 0.45, 0.4, 0.3, 0.25,\
	 0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025, 0.0125, 0.01, 0.005, 0.0 ])/87.175
	a = units*insulatard_profile[time]
	return a



basal_in_body = zeros(36)
bolus_in_body = zeros(144)
in_body = zeros(144)		#for 24 hours
real_time = linspace(0, 36 , 144)

basaltime = linspace(0, 36, 36)

basal( 8, 10)
bolus( 8, 2 )
bolus( 13, 1)
bolus( 17, 2)
bolus(20, 2)
bolus(23, 1)
print sum(basal_in_body)
print sum(bolus_in_body)
#print in_body
#print sum(in_body)
plot( real_time, bolus_in_body, "r-o")
hold("on")
plot( basaltime, basal_in_body, "b-o")
legend(["Bolus Insulin", "Basal Insulin"])
title("Insulin left from injenction")
show()
hold("off")