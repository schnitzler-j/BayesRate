import random as rand
import os
import csv
from numpy import *
from help import *
from files import fix_path

def runBMA():
	print """
	Bayesian Model Averaging (BMA) generates a joint posterior distribution by resampling
	MCMC samples obtained under different models based on their relative probability. 

	Load the log files (from 'cold' MCMC chains) and enter the respective 
	marginal likelihood obtained through thermodynamic integration (use log
	units e.g. -124.65). Type 'BMA' to stop adding log-files and calculate 
	model relative probabilities."""

	infiles=list()
	marginal=list()

	while True:
		infile=raw_input("\n> input log file: ")
		try:
			if infile == "help" or infile=="h": raise(ValueError)
			if infile=="BMA" or infile=="bma": break
			else:
				infile, input_file, name_file, path_dir=fix_path(infile)
				infiles.append(infile)
				while True:
					try:
						marginal.append(input("> marginal likelihood: "))
						break
					except(KeyboardInterrupt): raise(SyntaxError)
					except: help("BMA_ml")
		except(ValueError): help("log_BMA")
		except: print " "


	print len(infiles), "models"
	marginal=array(marginal)
	rel_prob=exp(marginal-max(marginal))
	rel_prob=rel_prob/sum(rel_prob)
	for i in range(len(rel_prob)): print "relative probability model %s: %s" % (i, rel_prob[i])

	row_count = sum(1 for row in csv.reader( open(infiles[0]) ))

	while True:
		max_it=row_count-1
		rand_seq=list()
		for p in rel_prob: 
			rand_seq.append(int(round((max_it/max(rel_prob))*p)))
		print rand_seq, "\n"
		out_log = "bayesrate/logBMA.log" # % (path_dir)
		logfile = open(out_log , "wb") 
		allK=list()
		it=0
		threshold=input("\n> min retained probability: ")
		for i in range(len(infiles)):
			if rel_prob[i]>threshold:
				t = csv.reader(open(infiles[i], 'rb'), delimiter='\t', quotechar='|')
				#print max_it, rand_seq[i]
				inds=rand.sample(xrange(1,max_it+1), rand_seq[i])
				#else: inds=range(1,max_it)
				#print min(inds), max(inds)
				K=list()
				for row in t: K.append(row)
				try: K=array(K)
				except(ValueError): print inds
				head=K[0]
				K=K[inds]
				print "Model %s: %s samples" % (i,len(K))
				itt=0
				for j in K: 
					itt+=10
					j[0]=it+itt
					allK.append(j)
				it += 100
		
		
		w=csv.writer(logfile, delimiter='	', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		w.writerow(head)
		for j in allK: w.writerow(j)
			#logfile.writelines(head)
		logfile.close()
		stem_out=raw_input("\n> output file name: ")

		file1=file(out_log, 'U')
		out="%s%s_BMA.log" % (path_dir, stem_out)
		newfile = open(out, "wb") 
		L=file1.readlines()
		h=L[0]
		h=h.replace(' ', '\t')
		newfile.writelines(h)

		for i in range(1,len(L)):
			z=L[i].split()
			z[0]=int(i*100)
			for j in range(0,len(z)): 
				a="%s\t" % (z[j])
				newfile.writelines(a)
			newfile.writelines("\n")
		newfile.close()
		print "\n%s MCMC samples included, %s models with %s cumulative probability." \
		% (len(L), len(rel_prob[rel_prob>threshold]), round(sum(rel_prob[rel_prob>threshold]), 3))
		print "The combined posterior was saved as:\n%s%s_BMA.log" % (path_dir, stem_out)
		break
