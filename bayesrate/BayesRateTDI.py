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
try: from multiprocessing import Pool
except: 
	print "\n Error: multiprocessing library not found!"
	quit()
from likelihood import *
from prior import *
from files import *
from mcmc_tdi import *
###################################
if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate-TDI\x07")
version= "BayesRate TDI 1.5 build 20121117"
print "\n		%s\n\tMarginal Likelihood via Thermodynamic Integration\n" % version

####################################	


K, input_file, name_file, path_dir= load_treefile()
if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate - TDI - %s\x07" % input_file)

print "\t%s trees found" % (len(K))
# DEFAULTS
global bd, exp_M, rj, categ
runs,iteration,rep,burnin_trees,mod_times,mod_rates,sample_freq,rand_tree,print_freq,missing_taxa,sample_fraction,fixed_times,fix_k,BF,bd,rj\
= 6,100000,100,0,.5,1,100,0,10000,0,1,list(),0,0,0,0
prior_shift, datasets=list(), 0
lam_rA, lam_mA = 0, 0
print "> please select the model of diversification: \n\t0) pure-birth\n\t1) birth death\
\n\t2) birth-death bothvar\n\t3) n-rate pure-birth\
\n\t4) models with taxon sampling"
while True:
	model=raw_input("> ") 
	try: 
		model=int(model)
		break
	except:
		if model=='h' or model=='help': help("model")
if model<0 or model>4: model=0
if model==0: bd, exp_M, rj, categ = 0, 0, 0, 1 # PB1
if model==1: bd, exp_M, rj, categ = 1, 0, 0, 1 # BD1
if model==2: 
	bd, exp_M, rj, categ = 1, 1, 0, 1 # BOTHVAR
	fix_k=cmd("> model of diversification: \n\t0) variable speciation and extinction (default)\n\t1) constant speciation\n\t2) constant extinction\n\t3) no extinction", "", "bothvar")

if model==3:								   # PBn
	bd, exp_M, rj, categ = 0, 0, 0, 1
	categ=cmd("> number of rates", 2, "n_rates")
	if categ<2: 
		print "\tassuming 2 rates..."
		categ=2
	categ=int(categ)
if model==4: 
	bd=cmd("> model of diversification (1): \n\t0) pure-birth\n\t1) birth death", "", "pb_bd")
	if bd != 0 and bd != 1: bd=1
	exp_M, rj, categ = 0, 0, 1
	sample_fraction=cmd("> sampling fraction", 1, "sample_fraction")

rep=int(cmd("> number of trees", rep, "rep"))
if rep> len(K): 
	print "\trunning on %s trees..." % (len(K))
	rep=len(K)
iteration=int(cmd("> number per-tree MCMC iterations", iteration, "iteration"))
runs=int(cmd("> number scaling classes (chains)", runs, "runs"))
try: BETA=int(cmd("> temperatures distribution:\n\t0) uniform\n\t1) beta (alpha=0.3)", "", "runs"))
except: 
	print "Assuming beta..."
	BETA=1

# ADVANCED
burnin=10000
if bd==1: rj, mod_times, mod_rates=0, 0, .66
while True:
	print """\n> Select option or type 'run' to start the analysis:
	1) MCMC settings
	2) Set priors (rates)
	3) Set constraints (shift times)
"""
	opt = raw_input("> ")
	try: 
		if opt in ('h', 'help'): help("adv_settings")
		opt=int(opt)
		if opt<1 or opt>3: pass
	except(NameError, ValueError): pass

	if opt==1:
		# " MCMC SETTINGS:"
		burnin_trees=cmd("> tree burnin fraction", (burnin_trees), "tree_burnin")
		starting_tree= 0
		q="> tree selection (%s):\n\t0) from %s to %s\n\t1) regularly sampled\n\t2) random" % (rand_tree,starting_tree+1,starting_tree+rep)
		rand_tree= cmd(q, "", "tree_selection")
		if rand_tree=="": rand_tree=0
		burnin=int(cmd("> MCMC burnin generations", burnin, "MCMC_burnin"))
		print_freq=int(cmd("> log-to-screen frequency", print_freq, "print_freq"))
		sample_freq=int(cmd("> sampling frequency", sample_freq, "sample_freq"))
		mod_rates=cmd("> update rates frequency", mod_rates, "mod_rates")
		if model==3:
			mod_times=cmd("> update shift-times frequency", mod_times, "mod_times")
			if mod_times==0:
				for i in range(1,categ): fixed_times.append(input("> please enter shift time (%s/%s): " % (i, categ-1)))

	# "\n PRIOR SETTINGS:"
	if opt==2:
		while True: 
			try: 
				pr=raw_input("> prior on net diversification (0): \n\t0) uniform (default)\n\t1) exponential\n> ")
				if pr in ('h', 'help'): help("prior_r")
				pr=int(pr)
				break
			except(ValueError): pass
		if pr==0: lam_rA = 0
		else: lam_rA=cmd("> lambda parameter", 1, "lam_r")
		if model==1 or model==2 or model==4:
			while True: 
				try: 
					pr=raw_input("> prior on extinction fraction (0): \n\t0) uniform (default)\n\t1) beta\n> ")
					if pr in ('h', 'help'): help("prior_m")
					pr=int(pr)
					break
				except(ValueError): pass
			if pr==0: lam_mA =0	
			elif pr==1: 
				a=(cmd("> alpha parameter", 1, "lam_m_a"))
				b=(cmd("> beta parameter", 1, "lam_m_b"))
				lam_mA = [a,b]
				print lam_mA, a, b
				
	if opt==3:
		if mod_times>0 and model==3:
			print "> constraints on shift times (uniform prior):"
			for i in range(1,categ): 
				q1= "\tlower bound of rate-shift no. %s of %s" % (i, categ-1)
				t1=cmd(q1, 0, "shift_1")
				prior_shift.append(t1)
				q2= "\tupper bound  of rate-shift no. %s of %s" % (i, categ-1)
				prior_shift.append(cmd(q2, t1, "shift_2"))
		else: print "\toption not available!"

	# RUN
	if opt in ('run', 'r', 'mcmc'): break

if exp_M==1: bd=1
prior_shift=list(sort(prior_shift))

def startMCMC(i):
	global bd, exp_M, rj, categ
	t1 = time.clock()
	MCMC(K,i, runs, iteration, burnin, categ, rep, mod_times, mod_rates, rj, sample_freq, rand_tree, bd, exp_M, model,\
	lam_rA, lam_mA, path_dir, name_file, fixed_times, prior_shift, sample_fraction, print_freq, fix_k,0,BETA)
	if i==0: print "\n MCMC execution time: {0:.2f} min".format((time.clock()-t1)/60.)
	if runs>1 and i==0: print " finishing remote chains..." 
bflogs="%sBFlogs" % (path_dir)
try: os.mkdir(bflogs) 
except(OSError): pass
pool = Pool(runs)
pool.map(startMCMC, range(0,runs)) 
out_list=list()
for i in range(0, runs): 
	if bd==1: MM="BD"
	if bd==0: MM="PB%s" % (categ)
	out_file = "%sBFlogs/%s_sum_%s_%s.txt" % (path_dir,name_file, MM, i)
	#out_file="%s_%s" % (name_file, i)
	out_list.append(out_file)
#print " results saved in files: %s\n" % (out_list)
print " summarizing %s runs (%s trees)..." % (runs, rep)
import mcmc_tdi
temperatures=mcmc_tdi.temp_gen(BETA, runs)
marginal_likelihood(runs, out_list, path_dir, name_file, MM, temperatures)

print "\n"
quit()