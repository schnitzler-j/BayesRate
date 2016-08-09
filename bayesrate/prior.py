#!/usr/bin/env python
# Created by Daniele Silvestro on 27/08/2010. => dsilvestro@senckenberg.de 
import sys
import os
import os.path
import platform
from numpy import *
import random
from help import *
# priors, initial and sample values

# sample positive numbers in range [m, M]
# from normal distr. centered in i
def wiener_process(m, s, lowest, greatest):  # NORMAL VARIATION
	k= random.normalvariate(m, s)
	if k<lowest: k= random.uniform(lowest, m)
	if k>greatest: k= random.uniform(m, greatest)
	return k

def update_parameter(i, std, m, M, freq): 
	if random.random()<=freq:
		i=fabs(random.normalvariate(i,std))
		if i>m: i=fabs(i-m)+m
		if i>M: i=fabs(M-(i-M))
	return i
	
def update_rate(rates, std_r, mod_rates):
	new_rates=list()
	std_r=update_parameter(float(std_r), .5, .01, 2.5, mod_rates)
	for i in rates:
		i=update_parameter(i, std_r, .001, 50, mod_rates)
		new_rates.append(i)
	return new_rates, std_r

def update_mrate(mrates, lrates, std_m, mod_rates):
	new_mrates=list() # print lrates, mrates
	std_m=update_parameter(std_m, .5, .01, .3, mod_rates) 
	for i in mrates:
		i=update_parameter(i, std_m, 0, 1, mod_rates)
		new_mrates.append(i)
	return new_mrates, std_m

def update_BOTHVAR(lrates, mrates, kA, zA, mod_rates):
	lratesN, mratesN= list(), list()
	lratesN.append(update_parameter(lrates[0], .25, .001, 5, mod_rates))
	mratesN.append(update_parameter(mrates[0], .25, 0, lratesN[0], mod_rates))
	kN=update_parameter(kA, .1, .001, 5, mod_rates)
	zN=update_parameter(zA, random.uniform(.1, 10), .001, 50, mod_rates)
	return kN, zN, lratesN, mratesN

def update_time(times, prior_shift, UP, LO):
	times=list(sort(times))
	times.reverse()
	new_times=list()
	new_times=times
	ind=int(random.uniform(1.5, len(times)-1.5))
	up=times[ind-1]
	lo=times[ind+1]
	if ind-1==0: up=UP
	if ind+1==len(times): lo=LO
	
	if len(prior_shift)>1: 	
		j=ind-1
		new_times[ind]=random.uniform(max(new_times[ind]-.5,lo+.5, prior_shift[j*2-1]), min(new_times[ind]+.5,up-.5, prior_shift[j*2-2]))
	else: new_times[ind]=random.uniform(max(new_times[ind]-.5,lo+.5, .5), min(new_times[ind]+.5,up-.5, max(times)-.5))
			#		if times[i]>=prior_shift[j*2-1] and times[i]<=prior_shift[j*2-2]: pass #prior_t.append(log(1./(prior_shift[j*2-1]-prior_shift[j*2-2])))
		#		else: prior_t.append(-100)

	times=list(sort(new_times))
	return times

def update_prior(root, up, lo, times, lrates, mrates, rep, lam_r, lam_m, bd, prior_shift):
	prior_l, prior_m, prior_t=list(), list(), list()	
	#print lam_m
	# UNIFORM PRIOR
	if len(times)>2:
		for i in range(1,len(times)-1):
			#lo, up = min(times[1:len(times)-1]), max(times[1:len(times)-1])
			if times[i] <= lo or times[i]>=up: prior_t.append(-inf)
			#elif fabs(times[i]-times[i+1])<1 and times[i+1]>0 :prior_t.append(-inf)
			else: prior_t.append(log(1./root)) # UNIFORM
			
			if len(prior_shift)>1:
				j=i-1
				if times[i]<=prior_shift[j*2-1] and times[i]>=prior_shift[j*2-2]: pass #prior_t.append(log(1./(prior_shift[j*2-1]-prior_shift[j*2-2])))
				else: prior_t.append(-inf)
				#print times[i], prior_shift[j*2-1], prior_shift[j*2-2]
			
	else: prior_t.append(0)
	
	for i in lrates:
		if lam_r==0:
			if i <=0 or i>5: prior_l.append(-inf)
			else: prior_l.append(log(1./5))     # UNIFORM
		else: 
			lam_r=1./lam_r
			prior_l.append(log(lam_r)-lam_r*i)  # EXPONENTIAL (lambda=0 if uniform)

	if bd==1:
		for i in range(0, len(mrates)): 
			if mrates[i]>= 1 or mrates[i]<0: prior_m.append(-inf) #lrates[i]
			try:
				if len(lam_m)>2: lam=lam_m[i]
				else: lam=lam_m
				a=float(lam[0]) 
				b=float(lam[1])
				if b>0:                           # UN-NORMALIZED BETA 
					prior_m.append((a-1)*log(mrates[i])+(b-1)*log(mrates[i]))
					#log(gamma(a+b)/(gamma(a)*gamma(b))*mrates[i]**(a-1)*(1-mrates[i])**(b-1)))
				else:
					l=-log(1-.95)/a
					#print l, log(l)-l*i, a, i
					prior_m.append(log(l)-l*mrates[i])    # EXPONENTIAL (semi-PB)
			except: 
				#print lam_m
				lam=lam_m
				prior_m.append(log(1.))  # UNIFORM
			
			#elif lam==0 and len(lam)==1: prior_m.append(log(1.))  # UNIFORM
			#else: 
	return sum(prior_t), sum(prior_l), sum(prior_m)

def filter_tree(K, rtree, st1, st2):
	nvp=0
	x=K[rtree]
	effnodes= range(0, len(K[0]))
	x=x[effnodes]
	x=list(x)
	x=list(sort(x))
	x.reverse()

	def f(x): return x >= st1
	def l(x): return x >= st2
	def m(x): return x < st1 and x >= st2
	nv = list(filter(m, x))   # vector of node heights within the time-frame
	lo = len(filter(f, x)) +1 # number of lineages at the start (n. nodes + 1)
	up = len(filter(l, x)) +1 # number of lineages at the end

	nv=list(nv)
	#res = list()
	if st1 <= x[0]: nv.insert(0,st1)
	nv=subtract(nv, st2) # node heights scaled by the 'end' age
	#return nv # 
	return up, lo, nv, st1-st2

def initial_params(root, categ, fixed_times, prior_shift):
	times=list()
	rates=list()
	mrates=list()
	times.append(root)
	for i in range(categ):
		if i==0: times.append(0)
		elif len(prior_shift)>1: times.append(random.uniform(prior_shift[i*2-1], prior_shift[i*2-2]))
		else: times.append(random.uniform(1,root*.33))
		#min(root*.33,max(1., random.lognormvariate(2,2)))) # random.uniform(lo, up))

	for i in range(categ):
		#if i==0: times.append(0)
		#elif len(prior_shift)>1: times.append(random.uniform(prior_shift[i-1], prior_shift[i]))
		#else: times.append(min(root*.33,max(1., random.lognormvariate(2,2)))) # random.uniform(lo, up))
		lrate=random.random() #max(0.0001, fabs(random.gammavariate(.5, 1)))
		mrate=max(0, random.uniform(0, lrate-.1))
		rates.append(lrate)
		mrates.append(mrate)
	if sum(fixed_times) != 0: 
		lim=root,  0. 
		times = fixed_times + list(lim)
	times=list(sort(times))
	times.reverse()
	std_t=std_r=.1
	kA, zA = random.uniform(.0001, .2), random.uniform(1, 10)
	return times, rates, mrates, std_r, std_t, kA, zA

def list_trees(K, rep, rand_tree, burnin_trees):
	rtree_list=list()
	roots=list()
	b=1-burnin_trees
	rtree= int(len(K)-(len(K)*b)) # first 25% of tree are excluded
	for j in range(rep):
		if rand_tree==0: rtree_list.append(j+rtree)
		if rand_tree==1: rtree_list.append(j*((len(K)-rtree)//rep)+rtree)
		if rand_tree==2: rtree_list.append(int(random.uniform(rtree, len(K))))		
		x=K[j]
	#	print len(x), len(K[0])
	#	effnodes= range(0, len(K[0]))
	#	x=x[effnodes]
		roots.append(max(x))
	root=int(mean(roots)) # change it to integer for custom max age
	#root=int(min(roots))
	return rtree_list, root, int(len(K)-(len(K)*b)), roots

def cmd(arg1, default, H):
	try: 
		while True:
			#value=raw_input("%s (%s, 'h'): " % (arg1, default))
			if default=="": value=raw_input("%s\n> " % (arg1))
			else: value=raw_input("%s (%s): " % (arg1, default))
			if value=="help" or value=="h": help(H)
			else: 
				return float(value)
				break
	except(ValueError): 
		if default != "": print "\tAssuming %s" % (default)
		return default
