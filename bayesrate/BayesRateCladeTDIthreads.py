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
import thread
from likelihood import *
from prior import *
from files import *
from help import *
from mcmc_tdi import MCMC_clade as MCMC
###################################
if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate - Clade\x07")
version= "BayesRate 1.5 build 20121117\n		    Clade-specific analysis"
if __name__=="__main__":
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

r=list()

def startMCMC(part_file, i, thread):
	#t1 = time.clock()
	arg=list()
	r.append(MCMC([trees_ingroup, input_file, name_file, path_dir, part_file, BETA], i, runs, 0))
	#if i==0: print "\n MCMC execution time: {0:.2f} min".format((time.clock()-t1)/60.)
	if runs>1 and i==0: print " finishing remote chains..." 
	if len(r)==runs:
		part_file, input_file_, name_file_, path_di_r=fix_path(part_file)
		f = open(part_file, 'r')
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

		import mcmc_tdi
		temperatures=mcmc_tdi.temp_gen(BETA, runs)
		marginal_likelihood(runs, out_list, path_dir, name_file, output_file, temperatures)

		print "\n Use 'Ctrl-C' to quit.\n"


def startMCMCseq(part_file):
	sequential=1
	for i in range(0, runs):
		MCMC([trees_ingroup, input_file, name_file, path_dir, part_file, BETA], i, runs, sequential)
	
	part_file, input_file_, name_file_, path_di_r=fix_path(part_file)
	f = open(part_file, 'r')
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

	import mcmc_tdi
	temperatures=mcmc_tdi.temp_gen(BETA, runs)
	marginal_likelihood(runs, out_list, path_dir, name_file, output_file, temperatures)


bflogs="%sBFlogs" % (path_dir)
try: os.mkdir(bflogs) 
except(OSError): pass

if __name__=="__main__":
	done=0
	for n in range(0,runs):
		thread.start_new_thread(startMCMC,(part_file, n,"thread"))
	try:
		while True:pass
	except: quit()
else: startMCMCseq(part_file)






