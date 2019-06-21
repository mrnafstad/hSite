from random import randint
import numpy as np
import sys


try:
	numbers = np.array(np.loadtxt("timeliste.txt", skiprows=0, unpack=True))
	full = np.array(np.loadtxt("full.txt", skiprows=0, unpack=True))
except IOError:
	print "Create a file named timeliste.txt looking like this:\n %.0f %2.0f %2.0f \n and try again" % (0,0,0)
	sys.exit(1)

fulllist = [np.array[full[0]], np.array[full[1]], np.array[full[2]]]
number = [numbers[0], numbers[1], numbers[2]]
print "2014.01.11", "---"*10, "New features:"
print "----------------"*10
print """You can add or remove hours in the same way as you could decide the alloted time yourself. Just do this: 
\n Terminal>python fagvelger.py 1 2 \n \n where 1 is the number of half hours you want to add / remove, 
and 2 corresponds to adding to Ast1100, writing 3 instead will add to Mat1120 and 4 will add to Fys1120.
\n \n If however 2 == 1, you will get a bar plot of your hours,
so you dont have to check the file to know!
"""
print "----------------"*10

ri = randint(0,3)
t = randint(1,4)
time = 30*t

if len(sys.argv) >= 2 and float(sys.argv[1]) >= 1 or len(sys.argv) >= 2 and float(sys.argv[1]) <= 1:
	time = float(sys.argv[1])*30

if len(sys.argv) >= 3 and float(sys.argv[2]) > 1:
		numbers[float(sys.argv[2])-2] = numbers[float(sys.argv[2])-2] + float(sys.argv[1])*30/60
		print "The workload has been added / removed"

if len(sys.argv) >= 4 and float(sys.argv[3]) == 1:
	import matplotlib.pylab as mpl
	fulllist = np.append(fulllist, number, 1)
	np.savetxt("full.txt", fulllist)
	print fulllist
	k = range(len(fulllist[0]))
	t = np.linspace(0, k, k)
	mpl.plot(fulllist[0], t)
	mpl.hold("on")
	mpl.plot(fulllist[1], t)
	mpl.plot(fulllist[2], t)
	mpl.show()


else:
	if ri == 0:
		print "\n Do Ast1100! a.k.a ast, a.k.a awesome.."
		numbers[0] = float(numbers[0]) + time/60.
		k = True
	elif ri == 1:
		print "\n Do Mat1120! a.k.a linalg, a.k.a awesome when you get it.."
		numbers[1] = float(numbers[1]) + time/60.
		k = True
	elif ri == 2:
		print "\n Do Fys1120! a.k.a read book instead of going to class, a.k.a actually cool just silly course.."
		numbers[2] = float(numbers[2]) + time/60.
		k = True
	elif ri == 3:
		print "\n Dont work! Work out, do some gaming, sleep, or whatever.."
		k = False

	if len(sys.argv) == 2 and float(sys.argv[1]) >= 1:
		print "\n Work for %.1f minutes" % (float(sys.argv[1])*30)
	elif len(sys.argv) == 2 and float(sys.argv[1]) < 1 or ri == 3:
		print "\n Do this for %.0f minutes" % time
	elif k:
		print "\n Work for %.0f minutes" % time


np.savetxt("timeliste.txt",numbers)

if len(sys.argv) == 3 and float(sys.argv[2]) == 1:
	import matplotlib.pyplot as mpl
	n = 3
	ind = np.arange(n)
	fig = mpl.figure()
	ax = fig.add_subplot(111)

	rects1 = ax.bar(ind, numbers, 0.35, color = "black", yerr=0)

	ax.set_xlim(-0.35, len(ind))
	ax.set_ylim(0,int(max(numbers))+2)
	ax.set_ylabel("Hours")
	ax.set_title("Hours of work")
	xTickMarks = ["Ast 1100", "Mat 1120", "Fys 1120"]
	ax.set_xticks(ind+0.35)
	xtickNames = ax.set_xticklabels(xTickMarks)
	mpl.setp(xtickNames, rotation=45, fontsize=10)
	mpl.show()	