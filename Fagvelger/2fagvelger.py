print "----------------"*5
print """
fagvelger 2.0; you can now take notes, and the program is in general faster and more intuitive, programingwise. There will probably be an update wich allows you
to create a plot to show how you have worked over time. fagvelger now uses dictionaries and functions! The only file you need to create is the timeliste.txt file,
everything else is taken care of by the program. This "issue" will be solved in due time.
"""
"""
Need to update to be more copatible with pomodoros, i.e. 25 minutes bulks of time. Don't need the time option, always one pomodoro.
More important to input worktime!
"""

print "----------------"*5
import numpy as np
import sys
import matplotlib.pylab as mpl

try:
	numbers = np.array(np.loadtxt("timeliste.txt", skiprows=0, unpack=True))
except IOError:
	print "Create a file called timeliste.txt on the form %.1e %.1e %.1e" % (0,0,0)
	sys.exit(1)

#Create lists and dictionaries

CME = {"Type": "FYS5190", "Total": numbers[0], "Latest": 0, "Notes": "i"}
KF = {"Type": "Special syllabus", "Total": numbers[1], "Latest": 0, "Notes": "i"}
EF = {"Type": "Master thesis", "Total": numbers[2], "Latest": 0, "Notes": "i"}
C = [CME, KF, EF]

#Define various functions

def time():
	#generates random time
	k = np.random.randint(1,4)
	t = 30*k
	return t

def cchooser(t):
	#Generates random course
	k = np.random.randint(0,2)
	c = C[k]
	print "----------------"*5
	c["Notes"] = raw_input("Notes for %s ? \n" % c["Type"])
	print "----------------"*5
	c["Latest"] = t
	c["Total"] += c["Latest"]
	return c

def cupdater(c, t):
	#Updates specified course
	print "----------------"*5
	c["Notes"] = raw_input("Notes for %s ? \n" % c["Type"])
	print "----------------"*5
	c["Latest"] = t
	c["Total"] += c["Latest"]
	return c

def saver(c):
	#Saves info
	numbers = [CME["Total"], KF["Total"], EF["Total"]]
	np.savetxt("timeliste.txt", numbers)
	note = open("Notes\%snotes.txt" %c["Type"], "a")
	note.write("%s \n" % c["Notes"])
	note.close()
	latest = open("Steps\%ssteps.txt" %c["Type"], "a")
	latest.write("%.2e" % c["Latest"])
	latest.close()
	print "----------------"*5
	print "Work with %s for %.0f minutes" % (c["Type"], t)
	print "----------------"*5

def weeksave():
	#Saves weekly info
	week = raw_input("Week #:")
	weekly = open("Weeks\hours.txt", "r")
	lines = weekly.readlines()

	klas = []	
	print type(klas)
	for line in lines:
		hrs = line.split(" ")
		for i in range(len(hrs)):
			hrs[i] = float(hrs[i])
			print hrs
		klas.append(hrs)
	k = [np.linalg.norm(klas[:0]), np.linalg.norm(klas[:1]), np.linalg.norm(klas[:2])]

	weekly = open("Weeks\hours.txt", "a")
	numb = [CME["Total"] - k[0], KF["Total"] - k[1], EF["Total"] - k[2]]
	print numb
	weekly.write("%.2f %.2f %.2f \n" % (numb[0], numb[1], numb[2]))
	weekly.close()

def plotter():
	#Creates a bar plot
	n = 3
	ind = np.arange(n)
	fig = mpl.figure()
	ax = fig.add_subplot(111)

	rects1 = ax.bar(ind, [CME["Total"], KF["Total"], EF["Total"]] , 0.35, color = "black", yerr=0)

	ax.set_xlim(-0.35, len(ind))
	ax.set_ylim(0,int(max(numbers))+2)
	ax.set_ylabel("Hours")
	ax.set_title("Hours of work")
	xTickMarks = ["%s" % CME["Type"], "%s" % KF["Type"], "%s" % EF["Type"]]
	ax.set_xticks(ind+0.35)
	xtickNames = ax.set_xticklabels(xTickMarks)
	mpl.setp(xtickNames, rotation=45, fontsize=10)
	mpl.show()

#Do various test to determine time spent and course worked on

if len(sys.argv) == 1:
	#Generates all at random
	t = time()
	c = cchooser(t)
	saver(c)

elif len(sys.argv) == 2:
	#Generates random time for specific course
	c = C[int(sys.argv[1])]
	t = time()
	c = cupdater(c, t)
	saver(c)

elif len(sys.argv) == 3 and float(sys.argv[2]) > 0:
	#Generates with specific time
	t = float(sys.argv[2])
	if int(sys.argv[1]) >= 3:
		#Random course
		c = cchooser(t)
	else:
		#Specific course
		c = C[int(sys.argv[1])]
		c = cupdater(c, t)
	saver(c)

elif len(sys.argv) == 4:
	#Generates with specific time
	t = float(sys.argv[2])
	if int(sys.argv[1]) >= 3:
		#Random course
		c = cchooser(t)
	else:
		#Specific course
		c = C[int(sys.argv[1])]
		c = cupdater(c, t)
	saver(c)
	plotter()

#elif int(sys.argv[3]) == 1:
#	weeksave()