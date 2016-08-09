#!/usr/bin/env python
# Created by Daniele Silvestro on 27/08/2010. => dsilvestro@senckenberg.de 
import sys
import os
import os.path
import platform
from math import *
from numpy import *
import random
import time
import csv
try: 
	from scipy.special import gamma
	from scipy.special import gammainc
	from scipy import integrate
except: pass # if run==0: print "\n Warning: scipy library not found, the program will run with some limitations!"
try: from multiprocessing import Pool
except: 
	print "\n Error: multiprocessing library not found!"
	quit()
from likelihood import *
from prior import *
from files import *
from help import *
from mcmc_tdi	import MCMC_clade as MCMC
from mcmc_tdi import temp_gen

###################################
if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate - Clade\x07")
version= "BayesRate 1.5 build 20121117\n		    Clade-specific analysis"
print "\n		%s\n\tMarginal Likelihood via Thermodynamic Integration\n" % version

while True: 
	try:
		#input_tree=raw_input("> input tree-file (NEXUS): ")
		trees_ingroup, input_file, name_file, path_dir=read_tree_nexus()
		if 2>1: # input_tree in ('h','help', ''): help("nexus")
			part_file=raw_input("> load partition file: ")
			if part_file in ('h','help', ''): help("part_fileTDI")
			else: 
				runs=int(cmd("> number scaling classes", 6, "runs"))
				try: BETA=int(cmd("> temperatures distribution:\n\t0) uniform\n\t1) beta (alpha=0.3)", "", "runs"))
				except: 
					print "Assuming beta..."
					BETA=1
				break
	except(KeyboardInterrupt): quit()
	except: pass

def startMCMC(i):
	t1 = time.clock()
	MCMC([trees_ingroup, input_file, name_file, path_dir, part_file, BETA], i, runs, 0)
	if i==0: print "\n MCMC execution time: {0:.2f} min".format((time.clock()-t1)/60.)
	if runs>1 and i==0: print " finishing remote chains..." 
bflogs="%sBFlogs" % (path_dir)
try: os.mkdir(bflogs) 
except(OSError): pass

pool = Pool(runs)
pool.map(startMCMC, range(0,runs)) 

part_file=part_file.replace("\ ",'sp')
part_file=part_file.replace(' ','')
part_file=part_file.replace('sp',' ')
f = open(part_file, "U")
L=f.readlines()
for l, line in enumerate(L):
	line= line.split()
	if "output_file" in line: 
		output_file=line[1]
		break

out_list=list()
for i in range(0, runs): 
	out_file = "%sBFlogs/%s_sum_%s.txt" % (path_dir,output_file, i)
	out_list.append(out_file)
print " summarizing %s runs..." % (runs)

temperatures= temp_gen(BETA, runs)
marginal_likelihood(runs, out_list, path_dir, name_file, output_file, temperatures)

print "\n"
quit()

