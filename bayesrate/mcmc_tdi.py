#!/usr/bin/env python
# Created by Daniele Silvestro on 27/08/2010. => dsilvestro@senckenberg.de 
import sys
import os
import os.path
import platform
from numpy import *
import random
import time
import csv
try: from multiprocessing import Pool
except: print "\n Error: multiprocessing library not found!"
from likelihood import *
from prior import *
from files import *
from help import *
###################################
def temp_gen(BETA, runs):
	if BETA==0:
		temperatures=list() 
		C=runs-1
		d=tuple(range(C+1))
		for i in d: temperatures.append(i*(1./C))
		temperatures=list(temperatures)
	else:
		# TI # Xie et al. 2011; Baele et al. 2012
		K=runs-1.
		k=array(range(int(K+1))) #K+1 categories
		beta=k/K
		alpha=.3 # categories are beta distributed
		temperatures=list(beta**(1./alpha))
		temperatures.reverse()
	#print temperatures
	return temperatures

def MCMC_clade(arg, run, runs,sequential):

	[trees_ingroup, input_file, name_file, path_dir, part_file, BETA]=arg
		
	n_datasets, dataset, size_clades, name_clades, age_clades, partition, models, rhos, trees, iteration, sampling_freq, \
	burnin, output_file, special_rule, mod_rates, print_freq, std_updates, lam_r, lam_m \
	=parse_file(part_file, trees_ingroup, input_file, name_file, path_dir, run)
	
	temperatures=temp_gen(BETA, runs)
	temp=temperatures[run]

	iteration=1+iteration
	if mod_rates==0: mod_rates=2./(max(partition)+1)

	out = "%sBFlogs/%s_sum_%s.txt" % (path_dir, output_file, run)
	outfile = open(out , "wb") 
	
	out_log = "%sBFlogs/%s%s.log" % (path_dir, output_file, run)
	logfile = open(out_log , "wb") 
	head="it\tposterior\tprior\tlikelihood\t"
	for i in range(len(dataset)): head += "net_div_%s\text_frac_%s\t" % (name_clades[i], name_clades[i])
	for i in range(len(dataset)): head += "speciation_%s\textinction_%s\t" % (name_clades[i], name_clades[i])
	head += "tree"
	wlog=csv.writer(logfile, delimiter='	')
	head=head.split('\t')
	wlog.writerow(head)
	logfile.flush()
	
	out0= " partition sequence: %s" % partition
	out_info(outfile, out0, run)

	for i in range(len(dataset)):
		out1=  "\n clade %s: %s" % (i, name_clades[i])
		out2=  " %s taxa included, taxon sampling: %.3f" % (size_clades[i], rhos[i])
		mod=['Pure-Birth', 'Birth-Death']
		out3=  " mean age: %.3f, diversification model: %s" % (age_clades[i], mod[int(models[i])])
		out4=  " output file: %s%s" % (output_file, run)
		out_info(outfile, out1+out2+out3, run)				
	acc=0
	it=0
	out0= "\n %s iterations (%s burnin), sampling-freq: %s, temperature: %s" % (iteration, burnin, sampling_freq, temp)
	out1= "\n rates update frequency: %s\n" % (mod_rates)
	out_info(outfile, out0+out1, run)		
	
	for t in range(trees): # tree
		lambdasA, musA=list(), list()
		for i in range(0,n_datasets): # start from random rates at every tree
			j=random.random()         # to avoid nan, -inf with temperature=inf
			lambdasA.append(j)
			if models[i]==1: musA.append(random.uniform(0,j))
			else: musA.append(0)
	
		if run==0 and sequential==0: print "\n sampling tree no. %s" % (t+1)
		if sequential==1: print "\n sampling tree no. %s temperature: %s" % (t+1,temperatures[run])
		list_coldlnL = list()
		for n in range(iteration): # generations
			lambdas, std_r=update_rate(lambdasA, float(std_updates), float(mod_rates))
			mus, std_m = update_mrate(musA, lambdas, float(std_updates), float(mod_rates))
			j=0
			for i in range(len(dataset)): # partition constraints
				if i>0:
					if partition[i]!=partition[i-1]: j=i
					else: j=i-1
				lambdas[i]=lambdas[j]
				if models[i]==1: mus[i]=mus[j]
				else: mus[i]=0
			
			if special_rule==1: lambdas[1]=(lambdas[0]/(1.-mus[0]))*(1.-mus[1])                           # constrains speciation rates [BD, BD]
			if special_rule==2: lambdas[1], mus[0]=(lambdas[0]/(1.-mus[0])), 0                            # constrains speciation rates [BD, PB]
			if special_rule==3: lambdas[0], mus[1]=(lambdas[1]/(1.-mus[1])), 0                            # constrains speciation rates [BD, PB]
			if special_rule==4: mus[1]=mus[0]*lambdas[0]/(lambdas[1]-mus[0]*lambdas[1]+mus[0]*lambdas[0]) # constrains extinction rates [BD, BD]			

			lik,prior=0,0
			for d in range(n_datasets):
				K=dataset[d]
				x=array(K[t])
				l=[lambdas[d]]
				m=[mus[d]]
				prior += sum(update_prior(1, 1, 1, list(), l, m, 1, lam_r, lam_m, models[d], list()))
				if rhos[d]<1 and prior>-inf: lik=lik+lnL_BD_Yang(x, lambdas[d], mus[d], models[d], rhos[d])
				elif rhos[d]==1 and prior>-inf: lik=lik+lnL_BirthDeath(x, lambdas[d], mus[d], models[d])
				
			if lik >=1 or lik<1: Post=lik*(temperatures[run])+prior
			else: Post=-inf
			if n==0: PostA, likA, priorA=Post, lik, prior
			if Post-PostA >= log(random.random()): # acceptance
				lambdasA, musA, likA= lambdas, mus, lik
				PostA=Post
				priorA=prior
				acc=acc+1.
			if n % print_freq ==0: 
				if run==0 or sequential==1:
					print n, "\t {0:.2f}".format(likA)   #" %s \tlnL %s" % (n,PostA)
					print "\t r",lambdasA
					print "\t a", musA
			if n % sampling_freq ==0 and n> burnin: 
				out="%s\t%s\t%s\t%s\t" % (it, PostA, priorA, likA)
				for i in range(len(dataset)): out=out+ "%s\t%s\t" % (lambdasA[i], musA[i])
				for i in range(len(dataset)): out=out+ "%s\t%s\t" % (lambdasA[i]/(1-musA[i]), -musA[i]*lambdasA[i]/(musA[i]-1))
				list_coldlnL.append(likA)
				out=out+"%s" % (t)
				out=out.split('\t')
				log_chain(logfile,wlog,out)
			if n > burnin: it=it+1
		out_mcmc ="\n no. tree: %s\n mean(lnL): %s \n std(lnL): %s\n acceptance: %s" % (t+1, mean(list_coldlnL), std(list_coldlnL), acc/n)#, mean(lnLS))
		out_info(outfile, out_mcmc, 1)
		if run==0: "\n acceptance: {0:.4f}\n".format(acc/n)
	return i
	

def update_likelihood(K, times, rates, mrates, rtree, bd, model, sample_fraction):
	lnL=0
	times=list(sort(times))
	times.reverse()
	rates.reverse()
	if bd!=0: mrates.reverse()
	for i in range(0, len(times)-1):
		st1=times[i]
		st2=times[i+1]
		lrate_i=rates[i]
		if bd!=0: mrate_i=mrates[i]
		else: mrate_i=0
		if bd==0 and model !=4: up, lo, nv, t =filter_tree(K, rtree, st1, st2)
		else: 	
			x=K[rtree]
		#	effnodes= range(0, len(K[0]))
		#	x=x[effnodes]
		if bd==0 and model !=4: lnL= lnL + lnL_Kendall(up, lo, nv, lrate_i)
		elif model==4: lnL= lnL_BD_Yang(x, lrate_i, mrate_i, bd, sample_fraction)
		else: lnL= lnL_BirthDeath(x, lrate_i, mrate_i, bd)
	return lnL 

def log_shift(root, tot_gen, t, run, logshift):
	s="%s\t" % (tot_gen)
	for k in range(1, len(t)-1): s=s+"%s\t" % (t[k])
	output(logshift, s, run, 0)

def log_rates(l_marginal, tot_gen, root, run, lograte, roots):
	s="%s\t" % tot_gen
	for i in range(int(min(roots))): s=s+"%s\t" % l_marginal[i]
	output(lograte, s, run, 0)

def marginal_rates(lograte, logshift, tot_gen, t, r, roots, run):
	root=int(mean(roots))
	cat=1
	l_marginal=list()
	for i in range(0,root, cat): # cycles from 0 to the root
		for j in range(0, len(t)-1): # cycles over time frames
			if i < t[j] and i+cat > t[j+1]: l_marginal.append(r[j]) # diversification rates
	log_rates(l_marginal, tot_gen, root, run, lograte, roots)
	log_shift(roots, tot_gen, t, run, logshift)

def MCMC(K,run, runs, IT, burnin, categ, rep, mod_times, mod_rates, rj, sample_freq, rand_tree, \
	bd, exp_M, model_analysis, lam_rA, lam_mA, path_dir, name_file, fixed_times, prior_shift, sample_fraction, \
	print_freq, fix_k, sequential,BETA):
	temperatures=list() #=0.0, 0.33333333333333331, 0.66666666666666663, 1.0 # 1., 1.5, 3., 100000.
	C=runs-1
	d=tuple(range(C+1))
	for i in d: temperatures.append(i*(1./C))
	#temperatures=list(temperatures)
	temperatures=temp_gen(BETA, runs)

	if bd==1: MM="BD"
	if bd==0: MM="PB%s" % (categ)
	out = "%sBFlogs/%s_sum_%s_%s.txt" % (path_dir, name_file, MM, run)
	outfile = open(out , "wb") 
	out_log = "%sBFlogs/%s_%s_%s.log" % (path_dir, name_file, MM, run)
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
	lnLS= list()  # list_std_r, list_std_t,

	try: rtree_list, mean_root, ex_trees, roots = list_trees(K,rep, rand_tree, burnin_trees)
	except: rtree_list, mean_root, ex_trees, roots = list_trees(K,rep, 0, 0)
	timesA, ratesA, mratesA, std_rA, std_tA, kA, zA=initial_params(mean_root, categ, fixed_times, prior_shift)
	std_mA=std_m=std_rA
	if bd==0 and rj==0 and categ==1: model="Yule process - constant rate"
	if bd==0 and rj==0 and categ>1: model="Yule process - %s rates (%s)" % (categ, prior_shift)
	if bd==0 and rj>0: model="Yule process - variable rate"	
	if bd==1: model="Birth-Death process - constant rates"
	if exp_M==1: model="Birth-Death process - BOTHVAR"
	if model_analysis==4: model= "%s - taxon sampling(%s)" % (model, sample_fraction)
	out_1= "\n file: %s\n %s trees (%s excluded, %s sampled)\n %s nodes, mean root age: %s"\
	% (name_file, len(K), ex_trees, rep, len(K[0]), mean_root)
	out_2= "\n %s iterations (%s burnin), sampling-freq: %s, remote chains: %s" % (IT, burnin, sample_freq, runs-1)
	out_3= "\n times-update: %s, rates-update: %s" % (mod_times, mod_rates)
	out_4= "\n starting time frames: %s l-rates: %s m-rates: %s" % (timesA, ratesA, mratesA)
	if exp_M==1: out_4= "\n l-rates: %s m-rates: %s k: %s z: %s" % (ratesA, mratesA, kA, zA)
	out_6= "\n model: %s - temperatures: %s [%s]\n" % (model, temperatures, temperatures[run])			
	out_mcmc = out_1+out_2+out_3+out_4+out_6
	output(outfile, out_mcmc, run, 1)
	tot_gen=0
	
	# log-files marginal rates, times
	if model_analysis==3:
		out_log = "%sBFlogs/%s_rates_%s_%s.log" % (path_dir, name_file, MM, run)
		lograte = open(out_log , "wb")
		head="it\t"
		for i in range(int(min(roots))): head=head+ "birth.rate_%s\t" % i
		head=head+ "\n"
		lograte.writelines(head)
		out_log = "%sBFlogs/%s_t_shift_%s_%s.log" % (path_dir, name_file, MM, run)
		logshift = open(out_log , "wb")
		head="it\t"
		for i in range(categ-1): head=head+ "shift_%s\t" % (i)
		head=head+"\n"
		logshift.writelines(head)

	for k in range(0, rep):
		list_coldlnL, acc=list(), 0
		rtree, root=rtree_list[k], (roots[k])
		it, num_tree = 0, k
		timesA, ratesA, mratesA, std_rA, std_tA, kA, zA=initial_params(root, categ, fixed_times, prior_shift)
		if run==0 and sequential==0: print "\n sampling tree no. %s [%s] root age: %s" % (k+1, rtree, root)
		if sequential==1: print "\n sampling tree no. %s [%s] root age: %s temperature: %s" \
		% (k+1, rtree, root, temperatures[run])
		for it in range(IT+1): # MCMC generations
			missing=list()
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
			#print up, lo

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
				else: 
					lnL=-inf
			if lnL >=1 or lnL<1: Post = lnL*(temperatures[run])+Prior #*(temperatures[run])+Prior
			else:
				lnL, Post = -inf, -inf
			if it ==0: lnLA, PriorA, PostA = lnL, Prior, Post
			if Post-PostA >= log(random.random()): # acceptance
				lnLA, PriorA, PostA = lnL, Prior, Post
				timesA, ratesA, mratesA = times, rates, mrates
				std_rA, std_tA, std_mA = std_r, std_t, std_m
				lam_rA, lam_mA = lam_r, lam_m
				kA, zA = k, z
				acc=acc+1.
				missingA = missing
			else: pass

			if it > burnin and it % sample_freq==0: # sample all parameters
				lnLS.append(PostA)#/temperatures[run]) #lnLA*(temperatures[run])) # heated
				list_coldlnL.append(lnLA)			  # cold
			if it % print_freq==0 and run==0: print it, "\t{0:.4f}  ".format(float(lnLA)), timesA
			if it % print_freq==0 and sequential==1 and run>0: print it,"\t{0:.4f}  ".format(float(lnLA)), timesA
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

		out_mcmc ="\n no. tree: %s\n mean(lnL): %s \n std(lnL): %s\n acceptance: %s" % (num_tree, mean(list_coldlnL), std(list_coldlnL), acc/(IT+1.))#, mean(lnLS))
		output(outfile, out_mcmc, run, 0)
	return i
	######