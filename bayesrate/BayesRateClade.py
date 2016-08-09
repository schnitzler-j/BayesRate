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
except: pass # print "\n Warning: scipy library not found, the program will run with some limitations!"
from likelihood import *
from prior import *
from files import *
from help import *
###################################
#if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate - Clade\x07")
version= "BayesRate 1.4 build 20120212\n		   Clade-specific analysis"
#print "\n		%s\n" % version

def MCMC():
	#trees_ingroup, input_file, name_file, path_dir, taxa_names= load_file(0)
	tree_list, input_file, name_file, path_dir=read_tree_nexus()
	dataset, partition, rhos, models=list(), list(), list(), list()
	size_clades, age_clades, name_clades=list(),list(),list()
	std_updates=.05
	sampling_freq, print_freq, burnin, special_rule, mod_rates=100, 10000, 10000, 0, 0
	lam_r, lam_m=0,0
	if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate - %s\x07" % input_file)
	while True:
		try:
			n_datasets, dataset, size_clades, name_clades, age_clades, partition, models, rhos, trees, iteration, sampling_freq, \
			burnin, output_file, special_rule, mod_rates, print_freq, std_updates, lam_r, lam_m \
			=parse_file(0, tree_list, input_file, name_file, path_dir, 0)
			break
		except(ValueError,IOError): print "\tFile format not recognized."

	iteration=1+iteration
	if mod_rates==0: mod_rates=2./(max(partition)+1)
	
	out = "%s%s_sum.txt" % (path_dir, output_file)
	outfile = open(out , "wb") 
	
	out_log = "%s%s.log" % (path_dir, output_file)
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
	print models
	out_info(outfile, out0, 0)

	for i in range(len(dataset)):
		out1= "\n clade %s: %s" % (i, name_clades[i])
		out2= "\n %s taxa included, taxon sampling: %.3f" % (size_clades[i], rhos[i])
		mod=['Pure-Birth', 'Birth-Death']
		out3= "\n mean age: %.3f, diversification model: %s" % (age_clades[i], mod[int(models[i])])
		out_info(outfile, out1+out2+out3, 0)		
		lambdasA, musA=list(), list()
	
	out4= "\n\n output file: %s" % output_file
	out_info(outfile, out4, 0)
			
	for i in range(0,n_datasets): 
		j=random.random()
		lambdasA.append(j)
		if models[i]==1: musA.append(random.uniform(0,j))
		else: musA.append(0)
	it=0
	IT=0
	out0= "\n %s trees, %s iterations (%s burnin), sampling-freq: %s" % (trees, iteration, burnin, sampling_freq)
	out1= "\n rates update frequency: %s\n" % (mod_rates)
	out_info(outfile, out0+out1, 0)
	
	for t in range(trees): # trees
		acc=0
		print " tree number %s" % (t+1)
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
				
			if lik> -inf: Post=lik+prior
			else: Post=-inf
			
			if n==0: PostA, likA, priorA=Post, lik, prior
			if Post-PostA >= log(random.random()): # acceptance
				lambdasA, musA= lambdas, mus
				PostA, likA, priorA=Post, lik, prior
				acc=acc+1.
			if n % print_freq ==0: 
				print n, "\t {0:.2f}".format(likA)   #" %s \tlnL %s" % (n,PostA)
				print "\t r",lambdasA
				print "\t a", musA
			if n % sampling_freq ==0 and n> burnin: 
				IT += sampling_freq
				out="%s\t%s\t%s\t%s\t" % (IT, PostA, priorA, likA)
				list_coldlnL.append(likA)
				for i in range(len(dataset)): out=out+ "%s\t%s\t" % (lambdasA[i], musA[i])
				for i in range(len(dataset)): out=out+ "%s\t%s\t" % (lambdasA[i]/(1-musA[i]), -musA[i]*lambdasA[i]/(musA[i]-1))
				out=out+"%s" % (t)
				out=out.split('\t')
				wlog.writerow(out)
				logfile.flush()

			it=it+1

		out_mcmc ="\n no. tree: %s\n mean(lnL): %s \n std(lnL): %s\n acceptance: %s" % (t+1, mean(list_coldlnL), std(list_coldlnL), acc/n)#, mean(lnLS))
		out_info(outfile, out_mcmc, 1)
		print "\n acceptance: {0:.4f}\n".format(acc/n)




