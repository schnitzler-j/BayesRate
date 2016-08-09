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
from likelihood import *
from prior import *
from help import *
from files import *
from mcmc_tdi import *

###################################
if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate\x07")
version= "BayesRate 1.4 build 20120212"
#print "\n		%s\n" % version

def MCMC():

	K, input_file, name_file, path_dir= load_treefile()
	if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate - %s\x07" % input_file)

	print "\t%s trees found" % (len(K))
	# DEFAULTS
	global bd, exp_M, categ
	runs,iteration,rep,burnin_trees,mod_times,mod_rates,sample_freq,rand_tree,print_freq,sample_fraction,fixed_times,fix_k,BF,bd\
	= 1,100000,100,0,.5,.5,100,0,10000,1,list(),0,0,0
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
	if model<0 or model>5: model=0
	if model==0: bd, exp_M, categ = 0, 0, 1 # PB1
	if model==1: bd, exp_M, categ = 1, 0, 1 # BD1
	if model==2: 
		bd, exp_M, categ = 1, 1, 1 # BOTHVAR
		fix_k=cmd("> model of diversification: \n\t0) variable speciation and extinction (default)\n\t1) constant speciation\n\t2) constant extinction\n\t3) no extinction", "", "bothvar")

	if model==3:								   # PBn
		bd, exp_M, categ = 0, 0, 1
		categ=cmd("> number of rates", 2, "n_rates")
		if categ<2: 
			print "\tassuming 2 rates..."
			categ=2
		categ=int(categ)
	if model==4: 
		bd=cmd("> model of diversification (1): \n\t0) pure-birth\n\t1) birth death", "", "pb_bd")
		if bd != 0 and bd != 1: bd=1
		exp_M, categ = 0, 1
		sample_fraction=cmd("> sampling fraction", 1, "sample_fraction")

	if model==5: 
		#bd=cmd("> model of diversification (1): \n\t0) pure-birth\n\t1) birth death", "", "pb_bd")
		#if bd != 0 and bd != 1: bd=1
		exp_M, categ = 0, 1
		bd=1
		sample_fraction=cmd("> sampling fraction", 1, "sample_fraction")


	rep=int(cmd("> number of trees", rep, "rep"))
	if rep> len(K): 
		print "\trunning on %s trees..." % (len(K))
		rep=len(K)
	iteration=int(cmd("> number per-tree MCMC iterations", iteration, "iterations"))
	# ADVANCED
	burnin=10000
	if bd==1: mod_times, mod_rates= 0, .66
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
			b=1-burnin_trees
			starting_tree= int(len(K)-(len(K)*b)) # first 25% of tree are excluded
			q="> tree selection (%s):\n\t0) from %s to %s\n\t1) regularly sampled\n\t2) random" % (rand_tree,starting_tree+1,starting_tree+rep)
			rand_tree= cmd(q, "", "tree_selection")
			if rand_tree=="": rand_tree=0
			#except: pass
			burnin=int(cmd("> MCMC burnin generations", burnin, "MCMC_burnin"))
			print_freq=int(cmd("> log-to-screen frequency", print_freq, "print_freq"))
			sample_freq=int(cmd("> sampling frequency", sample_freq, "sample_freq"))
			mod_rates=cmd("> update rates frequency", mod_rates, "mod_rates")
			if model==3:
				mod_times=cmd("> update shift-times frequency", mod_times, "mod_times")
				if mod_times==0:
					for i in range(1,categ): 
						q="> please enter shift time (%s/%s)" % (i, categ-1)
						fixed_times.append(cmd(q, ""))

		# "\n PRIOR SETTINGS:"
		if opt==2:
			while True: 
				try: 
					pr=raw_input("> prior on net diversification: \n\t0) uniform (default)\n\t1) exponential\n> ")
					if pr in ('h', 'help'): help("prior_r")
					pr=int(pr)
					break
				except(ValueError): pass
			if pr==0: lam_rA = 0
			else: lam_rA=cmd("> lambda parameter", 1, "lam_r")
			if model==1 or model==2 or model==4:
				while True: 
					try: 
						pr=raw_input("> prior on extinction fraction: \n\t0) uniform (default)\n\t1) beta\n> ")
						if pr in ('h', 'help'): help("prior_m")
						pr=int(pr)
						break
					except(ValueError): pass
				if pr==0: lam_mA =0	
				elif pr==1: 
					a=(cmd("> alpha parameter", 1, "lam_m_a"))
					b=(cmd("> beta parameter", 1, "lam_m_b"))
					lam_mA = [a,b]
					#print lam_mA, a, b
					
		if opt==3:
			if mod_times>0 and model==3:
				print "> constraints on shift times (uniform prior):"
				for i in range(1,categ): 
					q1= "\tlower bound of rate-shift no. %s of %s" % (i, categ-1)
					t1=cmd(q1, 0, "shift_1")
					prior_shift.append(t1)
					q2= "\tupper bound of rate-shift no. %s of %s" % (i, categ-1)
					prior_shift.append(cmd(q2, t1, "shift_2"))
			else: print "\toption not available!"

		# RUN
		if opt in ('run', 'r', 'mcmc'): break




	if exp_M==1: bd=1
	prior_shift=list(sort(prior_shift))

	#else: startMCMC()
	run=0
	if bd==1: MM="BD"
	if bd==0: MM="PB%s" % (categ)
	out = "%s%s_sum_%s.txt" % (path_dir, name_file, MM)
	outfile = open(out , "wb") 
	out_log = "%s%s_%s.log" % (path_dir, name_file, MM)
	logfile = open(out_log , "wb")
	if bd==0 and categ>1:
		head ="iteration\tposterior\tprior\tlikelihood\t"
		for i in range(categ): head += "sp.rate_%s\t" % (i+1)
		head+="\tstd_r\tstd_t\talpha(r)\tE(ext.frac)\tmodel\tk.par\tz.par\ttree"
	elif exp_M==0: head="iteration\tposterior\tprior\tlikelihood\tspeciation.rate\textinction.rate\tnet.diversification\textinction.fraction\tstd_r\tstd_t\talpha(r)\tE(ext.frac)\tmodel\tk.par\tz.par\ttree"
	else: head="iteration\tposterior\tprior\tlikelihood\tspeciation.rate\textinction.rate\tnet.diversification\textinction.fraction\tstd_r\tstd_t\talpha(r)\tE(ext.frac)\tlam0\tk.par\tz.par\ttree"
	wlog=csv.writer(logfile, delimiter='	')
	head=head.split('\t')
	wlog.writerow(head)
	logfile.flush()

	acc=PostA=0.
	#list_times, list_rates=list(),list()

	try: rtree_list, mean_root, ex_trees, roots = list_trees(K,rep, rand_tree, burnin_trees)
	except: rtree_list, mean_root, ex_trees, roots = list_trees(K,rep, 0, 0)
	timesA, ratesA, mratesA, std_rA, std_tA, kA, zA=initial_params(mean_root, categ, fixed_times, prior_shift)
	std_mA=std_m=std_rA
	model_analysis=model
	if bd==0 and categ==1: model="Yule process - constant rate"
	if bd==0 and categ>1: model="Yule process - %s rates (%s)" % (categ, prior_shift)
	if bd==1: model="Birth-Death process - constant rates"
	if exp_M==1: model="Birth-Death process - BOTHVAR"
	if model_analysis==4: model= "%s - taxon sampling(%s)" % (model, sample_fraction)
	if BF==1: doing_bf="(performing BF test)"
	else: doing_bf=""
	out_1= "\n file: %s\n %s trees (%s excluded, %s sampled)\n %s nodes, mean root age: %s"\
	% (name_file, len(K), ex_trees, rep, len(K[0]), mean_root)
	out_2= "\n %s iterations (%s burnin), sampling-freq: %s, remote chains: %s %s" % (iteration, burnin, sample_freq, runs-1, doing_bf)
	out_3= "\n times-update: %s, rates-update: %s" % (mod_times, mod_rates)
	out_4= "\n starting time frames: %s l-rates: %s m-rates: %s" % (timesA, ratesA, mratesA)
	if exp_M==1: out_4= "\n l-rates: %s m-rates: %s k: %s z: %s" % (ratesA, mratesA, kA, zA)
	out_5= "\n sampled trees: %s" % (rtree_list)
	out_6= "\n model: %s\n" % (model)			
	out_mcmc = out_1+out_2+out_3+out_4+out_5+out_6
	output(outfile, out_mcmc, run, 1)
	tot_gen=0
	
	# log-files marginal rates, times
	if model_analysis==3:
		out_log = "%s%s_rates_%s.log" % (path_dir, name_file, MM)
		lograte = open(out_log , "wb")
		head="it\t"
		for i in range(int(min(roots))): head=head+ "birth.rate_%s\t" % i
		head=head+ "\n"
		lograte.writelines(head)
		out_log = "%s%s_t_shift_%s.log" % (path_dir, name_file, MM)
		logshift = open(out_log , "wb")
		head="it\t"
		for i in range(categ-1): head=head+ "shift_%s\t" % (i)
		head=head+"\n"
		logshift.writelines(head)
	IT=iteration
	for k in range(0, rep):
		rtree, root=rtree_list[k], (roots[k])
		it, num_tree = 0, k
		timesA, ratesA, mratesA, std_rA, std_tA, kA, zA=initial_params(root, categ, fixed_times, prior_shift)
		if run==0: print "\n sampling tree no. %s [%s] root age: %s" % (k+1, rtree, root)
		for it in range(IT): # MCMC generations
			if it > burnin: tot_gen=tot_gen+1
			rates, std_r=update_rate(ratesA, std_rA, mod_rates)
			rates.reverse()
			if bd==1: 
				mrates, std_m=update_mrate(mratesA, rates, std_mA, mod_rates)
				mrates.reverse()
			else: mrates=0		
			lam_r, lam_m=lam_rA, lam_mA

			if it==0: rates, mrates=ratesA, mratesA
			
			temp=list(K[rtree])
			temp.sort()
			temp.reverse()
			if len(timesA)>2:
				up=temp[1]
				lo=temp[len(temp)-1]
			else: up,lo = temp[0], 0
			#if it==0 and run==0:
			#	print K[rtree]
			#	print temp
			#	print up, lo

			if random.random() < mod_times and it>0 and len(timesA)>2: times=update_time(timesA, prior_shift, up, lo) # UPDATE TIMES
			else: times, std_t=timesA, std_tA

			times.sort() 
			times.reverse()

			prior_t, prior_l, prior_m=update_prior(root, up, lo, times, rates, mrates, rep, lam_r, lam_m, bd, prior_shift)

			if exp_M==1: 
				k, z, rates, mrates = update_BOTHVAR(ratesA, mratesA, kA, zA, mod_rates)
				if fix_k==1: k=0.00001
				if fix_k==2: z=10000
				if fix_k==3: z,mrates[0]=10000,0
			else: k, z, = 0,0

			if 2>1: #try:
				Prior=prior_l+prior_t+prior_m
				if Prior> -inf:
					if exp_M==0: lnL=update_likelihood(K, times, rates, mrates, rtree, bd, model_analysis, sample_fraction)
					else: lnL, lam0A=lnL_BOTHVAR(K, rates, k, mrates, z, rtree)
				else: lnL=-inf
			#except: Prior, lnL = (prior_l+prior_t+prior_m), -inf

			if lnL >=1 or lnL<1: pass
			else: lnL = -inf
			Post = lnL + Prior 
			if it ==0: lnLA, PriorA, PostA = lnL, Prior, Post
			if Post-PostA >= log(random.random()): # acceptance
				lnLA, PriorA, PostA = lnL, Prior, Post
				timesA, ratesA, mratesA = times, rates, mrates
				std_rA, std_tA, std_mA = std_r, std_t, std_m
				lam_rA, lam_mA = lam_r, lam_m
				kA, zA = k, z
				acc=acc+1.
			else: pass

			#if it > burnin and it % sample_freq==0: # sample all parameters
			#	list_times.append(timesA)
			#	list_rates.append(ratesA)
			if fabs(lnLA)<100 and it % print_freq==0 and run==0: print it, "\t{0:.2f}  ".format(float(lnLA)), "t:", timesA#lnL, Prior
			if fabs(lnLA)>=100 and it % print_freq==0 and run==0: print it, "\t{0:.2f} ".format(float(lnLA)), "t:", timesA #, lnL, Prior
			if it % print_freq==0 and run==0 and bd==0: print "\t\tl:", ratesA #, std_rA, std_tA, "P", lam_rA #, "\n\t", rates #marginal_rates
			if it % print_freq==0 and run==0 and bd==1: print "\t\tr:", ratesA, "\ta:", mratesA #, std_m, lam_mA
			if it % print_freq==0 and run==0 and exp_M==1: print "\t\tk:", kA, "z:", zA, "lam0:", lam0A
			# LOG CHAIN
			if it > burnin:
				if lam_m==0: prM=0    # UNIFORM
				else: 
					a=float(lam_m[0])
					b=float(lam_m[1])
					prM=a/(a+b)
				
				
				mean_Brate=0.
				mean_Drate=0.
				for i in range(0, len(ratesA)): mean_Brate=mean_Brate+ ratesA[i]*(timesA[i]-timesA[i+1])
				if bd ==1:
					for i in range(0, len(mratesA)): mean_Drate=mean_Drate+ mratesA[i]*(timesA[i]-timesA[i+1])
				if tot_gen % sample_freq==0 and exp_M==0 and model_analysis!=3: log_chain(logfile, wlog, [tot_gen, PostA, PriorA, lnLA, mean_Brate/root/(1-mean_Drate/root), \
				-mean_Brate/root*mean_Drate/root/(mean_Drate/root-1), mean_Brate/root, mean_Drate/root, std_rA, std_tA, lam_rA, prM, len(ratesA), kA, zA, num_tree]) # IF YOU SAMPLE r AND a
				#if tot_gen % sample_freq==0 and exp_M==0: log_chain(logfile, tot_gen, PostA, PriorA, lnLA, mean_Brate/root, mean_Drate/root, \
				#(mean_Brate/root-mean_Drate/root), ((mean_Drate/root)/(mean_Brate/root)), std_rA, std_tA, lam_rA, lam_mA, len(ratesA), kA, zA, num_tree)
	
				elif tot_gen % sample_freq==0 and exp_M==1: log_chain(logfile, wlog, [tot_gen, PostA, PriorA, lnLA, mean_Brate/root, mean_Drate/root, \
				(mean_Brate/root-mean_Drate/root), ((mean_Drate/root)/(mean_Brate/root)), std_rA, std_tA, lam_rA, prM, lam0A, kA, zA, num_tree])
				
				elif model_analysis==3 and tot_gen % sample_freq==0:
					all_args=[tot_gen, PostA, PriorA, lnLA]+ratesA+[std_rA, std_tA, lam_rA, prM, len(ratesA), kA, zA, num_tree]
					log_chain(logfile, wlog, all_args)
					list_t, list_r=timesA, ratesA
					marginal_rates(lograte, logshift, tot_gen, list_t, list_r, roots, run)

	out_mcmc ="\n Acceptance probability: {0:.4f}".format(acc/(IT*rep))
	output(outfile, out_mcmc, run, 1)
	print " The MCMC sample was saved to file: '%s_%s.log'" % (name_file, MM)
	if model_analysis==3:
		print " Marginal rates were saved to file: '%s_rates_%s.log'" % (name_file, MM)
		print " Posterior rate-shift times were saved to file: '%s_t_shift_%s.log'" % (name_file, MM)
	print "\n"
