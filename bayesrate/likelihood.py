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
try: from scipy import integrate
except: pass
# BIRTH-DEATH AND PURE-BIRTH LIKELIHOOD FUNCTIONS

def lnL_Kendall(up, lo, nv, lrate):
	nv=array(nv)
	nvp=sum(nv[1:len(nv)])
	s1=sum(log(range(lo, up)))	
	s2=(up-lo)*log(lrate)		
	s3=-(lo*nv[0]+nvp)*lrate	
	return s1+s2+s3

def lnL_BirthDeath(nv, r, a, bd):
	nv=list(nv)
	nv.sort()
	nv.reverse()
	nv.insert(0, 0)
	s=len(nv)	
	t=array(nv)	
	if bd==0: r,a=r,0
	return (sum(log(range(1,s)))+((s-2)*log(r))+(r*sum(t[2:s+1]))+(s*log(1-a))-2*sum(log(exp(r*t[1:s+1])-a)))

def lnL_BD_Yang(x, r, a, bd, rho):
	l=r/(1-a)
	m= -a*r/(a-1)
	t=array(x)
	t1=max(t)
	s=len(x)+1
	if bd==0: m=0
	def p0(t): return 1 - (rho*(l-m) / (rho*l + (l*(1-rho) -m) * exp(-(l-m)*t)))
	def p1(t): return  rho*(l-m)**2*exp(-(l-m)*t)/(rho*l + (l*(1-rho) -m)*exp(-(l-m)*t))**2
	return log(s-2) + log(p1(max(t))/(1-p0(max(t)))) + sum(log(l*p1(t[0:len(t)])))

def lnL_BD_Yang_scale_total_length(x, r, a, bd, rho): # rescale the total length of the tree to its 
	t=array(x)
	s=len(x)+1
	scal=.5*s/sum(t)
	t=t*scal
	scal_rate=max(t)/max(x)
	l=(r/(1-a))/scal_rate
	m=(-a*r/(a-1))/scal_rate
	s=len(x)+1
	if bd==0: m=0
	def p0(t): return 1 - (rho*(l-m) / (rho*l + (l*(1-rho) -m) * exp(-(l-m)*t)))
	def p1(t): return  rho*(l-m)**2*exp(-(l-m)*t)/(rho*l + (l*(1-rho) -m)*exp(-(l-m)*t))**2
	return log(s-2) + log(p1(max(t))/(1-p0(max(t)))) + sum(log(l*p1(t[0:len(t)])))

def lnL_BD_lineage_fast(x, ancstate, rhos, lambdas, mus):
	r= array(lambdas)
	a= array(mus)
	L= r/(1-a)
	M= -a*r/(a-1)
	l=L[ancstate]
	m=M[ancstate]
	def p0(t,l,m,rho): return 1 - (rho*(l-m) / (rho*l + (l*(1-rho) -m) * exp(-(l-m)*t)))
	def p1(t,l,m,rho): return rho*(l-m)**2*exp(-(l-m)*t)/(rho*l + (l*(1-rho) -m)*exp(-(l-m)*t))**2
	clade_lik = sum(log(l[0:len(x)-1]*p1(x[0:len(x)-1],l[0:len(x)-1],m[0:len(x)-1],rhos[0:len(x)-1])))
	return (log(len(x)-1) + log(p1(max(x), mean(l), mean(m), mean(rhos))/(1-p0(max(x), mean(l), mean(m), mean(rhos)))) + clade_lik)

def lnL_BD_lineage_fast_spec(x, ancstate, rhos, lambdas, mus): # samples speciation rate and extinction fraction
	r= array(lambdas)
	a= array(mus)
	L= r
	M= r*a
	l=L[ancstate]
	m=M[ancstate]
	def p0(t,l,m,rho): return 1 - (rho*(l-m) / (rho*l + (l*(1-rho) -m) * exp(-(l-m)*t)))
	def p1(t,l,m,rho): return rho*(l-m)**2*exp(-(l-m)*t)/(rho*l + (l*(1-rho) -m)*exp(-(l-m)*t))**2
	clade_lik = sum(log(l[0:len(x)-1]*p1(x[0:len(x)-1],l[0:len(x)-1],m[0:len(x)-1],rhos[0:len(x)-1])))
	return (log(len(x)-1) + log(p1(max(x), mean(l), mean(m), mean(rhos))/(1-p0(max(x), mean(l), mean(m), mean(rhos)))) + clade_lik)

def lnL_BD_lineage_pred_div(x, ancstate, rhos, lambdas, mus):
	ancstate=array(ancstate)
	rhos=array(rhos)
	# values in state 0 are the training nodes
	rhos=rhos[ancstate==0]
	x=x[ancstate==0]
	ancstate=ancstate[ancstate==0]
	r= array(lambdas)
	a= array(mus)
	L= r/(1-a)
	M= -a*r/(a-1)
	l=L[ancstate]
	m=M[ancstate]
	def p0(t,l,m,rho): return 1 - (rho*(l-m) / (rho*l + (l*(1-rho) -m) * exp(-(l-m)*t)))
	def p1(t,l,m,rho): return rho*(l-m)**2*exp(-(l-m)*t)/(rho*l + (l*(1-rho) -m)*exp(-(l-m)*t))**2
	clade_lik = sum(log(l[0:len(x)-1]*p1(x[0:len(x)-1],l[0:len(x)-1],m[0:len(x)-1],rhos[0:len(x)-1])))
	return (log(len(x)-1) + log(p1(max(x), mean(l), mean(m), mean(rhos))/(1-p0(max(x), mean(l), mean(m), mean(rhos)))) + clade_lik)

def lnL_BOTHVAR(K, lrates, k, mrates, z, rtree):
	mu0 = mrates[0] #, mrates[0]/lrates[0]
	x=K[rtree]
	effnodes= range(0, len(K[0]))
	x=x[effnodes]
	x=list(sort(x))
	x.reverse()
	Tmax=max(x)
	lam0 = lrates[0] # + mu0*(1 - exp(-z * Tmax))) * exp(Tmax * k)
	x=array(x)
	x=list(Tmax-x)
	x.insert(0, x[0])
	t0=x[0]
	x=array(x)
	def lambdaFx(lam0, k, time): return lam0 * exp(-k * time)

	def rhoFxn(lam0, k, mu0, z, timelow, timehigh):
		return ((mu0 * timehigh) + (mu0/z) * exp(-z * timehigh) + (lam0/k) * \
		exp(-k * timehigh) - (mu0 * timelow) - (mu0/z) * \
		exp(-z * timelow) - (lam0/k) * exp(-k * timelow))

	def pTt(lam0, k, mu0, z, t0, Tmax): 
		def FXrho(x): return ((mu0 * x) + (mu0/z) * exp(-z * x) + (lam0/k) * exp(-k * x))
		def fx1(x): return ((1 - exp(-z * x)) * exp(FXrho(x) - FXrho(t0)))
		try:
			temp_fx= integrate.quad(fx1, t0, Tmax) #, limit=100)
			return (1/(1 + mu0 * temp_fx[0]))
		except: return -inf

	L3=0
	L1 = sum(log(range(1,len(x)))) + sum(log(lambdaFx(lam0, k, x[2:len(x)])))
	ptT = pTt(lam0, k, mu0, z, t0, Tmax)
	rho = rhoFxn(lam0, k, mu0, z, t0, Tmax)
	L2 = (2 * log(ptT) + 2 * rho)
	for i in range(2, len(x)): L3 = (L3 + 2 * log(pTt(lam0, k, mu0, z, x[i], Tmax)) + rhoFxn(lam0, k, mu0, z, x[i], Tmax))
	return L1+L2+L3, lam0
	
def marginal_likelihood(runs, names, path_dir, name_file, output_file, x):
	import numpy as np
	summ = "%s%s_%s_marginal.txt" % (path_dir, name_file, output_file)
	outfile = open(summ , "wb")
	f,sd,ac=list(),list(),list()
	for run in range(0, runs):
		acc, lnL, stds=list(), list(), list()
		input_file = names[run]
		infile = open(input_file, "U")
		L=infile.readlines()
		if run==0:
			for l, line in enumerate(L):
				if "no. tree" in line: break
				else: outfile.writelines(line)
		for l, line in enumerate(L):
			if "mean(lnL)" in line:
				lineL = line.split(':')
				lnL.append(float(lineL[1]))
			if "std(lnL)" in line:
				lineL = line.split(':')
				stds.append(float(lineL[1]))
			if "acceptance" in line: 
				lineL = line.split(':')
				acc.append(float(lineL[1]))
		f.append(mean(lnL))
		sd.append(mean(stds))
		ac.append(mean(acc))
		
	s= " mean log likelihoods: %s\n temperatures: %s\n std(lnL): %s\n acceptance: %s" % (f,x,sd,ac)
	print s
	outfile.writelines(s)

	
	mL=0
	dT=np.diff(x)
	if abs(dT[0]-dT[len(dT)-1]) < 0.001: # Bezier approximation
		for b in range(2, len(x)): mL += .5 * (x[b]-x[b-1]) * (f[b-1]+f[b])
		c0 = 1./5*f[0]+4./5*f[1]
		c1 = (x[1]*f[2] - x[2]*f[1])/(x[1]-x[2])
		mL += 1./20*((x[1]-x[0])*(f[0]+3*c0+6*c1+10*f[1]))

	else:   # integration by trapezoidal rule
		for i in range(len(x)-1): 
			M=((f[i]+f[i+1])/2.)*(x[i]-x[i+1]) 
			mL+=M
		
		
	sML= "\n\n Log Marginal Likelihood: %s" % (mL)
	print sML
	print """ 
 The marginal likelihood can be used to compare different analyses and to perform model selection
 and hypothesis testing via Bayes factors. 

 The log files of each individual chain used for thermodynamic integration were saved in a 'BFlogs'
 folder created in the directory of the input file. The marginal likelihood was saved to file: 
 '%s_%s_marginal.txt' \n""" % (name_file, output_file)
	outfile.writelines(sML)	


