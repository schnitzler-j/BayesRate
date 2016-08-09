#!/usr/bin/env python 
# Created by Daniele Silvestro on 13/06/2011. => daniele.silvestro@unil.ch 
import sys
import os
import os.path
import platform
self_path = os.getcwd()
from help import *
from prior import cmd as cmd_1

if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate\x07")
print """
                           BayesRate 1.65
                    build 20160326 Cross Platform

            Bayesian estimation of diversification rates

                Daniele Silvestro and Jan Schnitzler
				
           Biodiversity and Climate Research Centre (BiK-F)
                   Senckenberg Research Institute
                     Frankfurt am Main, Germany
                     daniele.silvestro@unil.ch
                    j.schnitzler@uni-leipzig.de

 Type 'help' for information on the commands that are available.
\n"""
#################################### LOAD LIBRARIES #####################################
def exit_error(msg): 
	print msg
	raw_input("> Press Enter to quit...") 
	quit()
print "> Importing required libraries:"
try: 
	import numpy 
	from numpy import *
	print "\tnumpy................loaded"
except: exit_error("\tnumpy...........Error: numpy library not found.\n\tYou can download numpy at: http://sourceforge.net/projects/numpy/files/ \n")

try:
	import dendropy
	print "\tdendropy.............loaded"	
except:
	exit_error("\tdendropy.............Error: dendropy library not found.\n\tYou can download dendropy at: \
http://packages.python.org/DendroPy/downloading.html#source-download-and-installation \n")

try: 
	from scipy.special import gamma
	from scipy.special import gammainc
	from scipy import integrate
	print "\tscipy................loaded\n"
except: print "\tscipy................Warning: scipy library not found.\n\tThe program will run with some limitations!\n"

from files import fix_path, read_tree_nexus,parse_bisse_file

#################################### MAIN FUNCTIONS #####################################
def menu_opt(model_ind):
	list_cmds=['mcmc_rtt', 'mcmc_clade', 'mcmc_SM','musse', 'tdi_rtt', 'tdi_clade', 'tdi_SM',\
	'sub_trees', 'remove_outgroup', 'BMA', 'plot_rates', 'plot_rtt', 'plot_rate_shift']
	model=list_cmds[model_ind-1]
	#print model
	
	if model=='mcmc_rtt': 
		import BayesRate1
		BayesRate1.MCMC()
	if model=='mcmc_clade': 
		analysis=int(cmd_1("> select analysis: \n\t0) clade specific rates\n\t1) key innovation test (Moore and Donoghue 2009)", "", "clade_key_inn"))
		if analysis==0:
			import BayesRateClade
			BayesRateClade.MCMC()
		if analysis==1: 
			import BayesRateCladeKey
			BayesRateCladeKey.MCMC()

	# MuSSE MENU
	if 'musse' in model:
		while True:
			infile1= raw_input("> tree-file (NEXUS): ")
			infile1, input_file1, name_file1, path_dir1=fix_path(infile1)
			infile2= raw_input("> trait file: ")
			infile2, input_file2, name_file2, path_dir2=fix_path(infile2)
			i= int(cmd_1("> select model: \n\t0) MuSSE\n\t1) GeoSSE\n\t2) ClaSSE", "", ""))
			SSE_models=["musse","geosse","classe"]
			SSE_model = SSE_models[i]
			tax_sampling= raw_input("""> Sampling fractions are provided for each character state as space-separated numbers.
E.g. for a MuSSE analysis with 3 character states > 0.5 0.25 0.7\n> taxon sampling: """)
			if tax_sampling in ('h','help'): help('tax_sampling_musse')
			else:
				burnin= cmd_1("> burnin", 100, "")
				ntrees= cmd_1("> number of trees: ", 100, "no_tree_subsample")
				iterations= cmd_1("> number of iterations: ", 1000000, "")
				sampling_frac= cmd_1("> sampling frequency: ", 100, "")
				print_frac= cmd_1("> log-to-screen frequency: ", 100, "")
				win_size= cmd_1("> window size: ", 0.05, "")
				prior= cmd_1("> exponential prior: ", 5, "")
				break
		args= "--b %s --t %s --i %s --s %s --p %s --d %s --r %s " % (int(burnin),int(ntrees),int(iterations),int(sampling_frac),int(print_frac),win_size,prior)
		cmd= """cd "%s/bayesrate/r_functions" &&Rscript BDiversitree1.r %s %s %s %s --rho "%s" """ % (self_path, infile1,infile2,SSE_model, args, tax_sampling)
		print cmd
		os.system(cmd)
	
	
	# BiSSE-BMA MENU
	if '_SM' in model:
		if 'tdi' in model: TDI=1
		else: TDI=0
		print "\nThe BiSSE-BMA method requires the R library Diversitree (Fitzjohn 2012) and its dependencies.\n"
		while True:
			infile1= raw_input("> tree-file (NEXUS): ")
			if infile1 in ('h','help'): help('nexus')
			else:
				infile1, input_file1, name_file1, path_dir1=fix_path(infile1)
				infile2= raw_input("> trait table: ")
				if infile2 in ('h','help'): help('mcmc_bisse_tbl')
				else:
					infile2, input_file2, name_file2, path_dir2=fix_path(infile2)
					S=parse_bisse_file(infile2)
					break
		cmd= """cd "%s/bayesrate/r_functions" &&Rscript BiSSEBMA.r "%s" "%s" "%s" %s %s""" % (self_path, infile1,infile2,path_dir1, TDI, S)
		#print cmd
		os.system(cmd)

	if 'tdi_rtt' in model or 'tdi_clade' in model:
		while True:	
			try: 
				th=int(cmd_1("> run options: \n\t0) sequential\n\t1) parallel (multi-thread)\n\t2) parallel (multiple processors)", "", "threads"))
				break
			except(ValueError): pass

		if th==0:
			thr="sequential"
			if model =='tdi_rtt': import BayesRateTDIthreads
			if model =='tdi_clade': import BayesRateCladeTDIthreads				
		else:
			if th==1: thr="threads"
			if th==2: thr=""
		
			if platform.system() != "Windows" and platform.system() != "Microsoft": dotslash,winD="python",""
			else: dotslash,winD="","/D"

			if model =='tdi_rtt' : cmd= """cd %s "%s/bayesrate" && %s "BayesRateTDI%s.py" """ % (winD,self_path, dotslash,thr)
			if model =='tdi_clade' : cmd= """cd %s "%s/bayesrate" && %s "BayesRateCladeTDI%s.py" """ % (winD, self_path, dotslash,thr)
			#if model =='tdi_SM': cmd="""cd %s/bayesrate &&BayesRateSMTDI%s.py""" % (self_path, thr)

			try:
				cmdfile = "%s/bayesrate/cmd.bat" % self_path
				cmdfile = file(cmdfile, 'w')
				cmdfile.writelines(cmd)
				cmdfile.close()
				if platform.system() != "Windows" and platform.system() != "Microsoft":
					cmdScript= """ chmod a+x '%s/bayesrate/cmd.bat' """ % self_path
					os.system(cmdScript)
				else:pass
				if platform.system() == "Darwin" and ' ' not in self_path: 
					cmdScript= """osascript -e 'tell application "Terminal" to do script "cd '%s/bayesrate' && ./cmd.bat" '""" % self_path
				elif platform.system() == "Windows" or platform.system() == "Microsoft": cmdScript= """ start "" "%s/bayesrate/cmd.bat" """ % self_path
				else: cmdScript= """xterm -e '%s' & """ % cmd
				#print cmdScript
				os.system(cmdScript)
			except: print "This option is not yet implemented"

	cmd=0
#################################### PLOT RESULTS #####################################
	if 'sub_trees' == model:
		while True:
			infile1= raw_input("> tree-file (NEXUS): ")
			if infile1 in ('h','help'): help('nexus_subsample')
			else:
				infile, input_file, name_file, path_dir=fix_path(infile1)
				burnin= cmd_1("> burnin fraction", .25, "tree_burnin")
				ntrees= cmd_1("> number of trees: ", 100, "no_tree_subsample")
				#models=zeros(2)
				#models[0]= cmd_1("> model of diversification (state 0): \n\t1) pure-birth\n\t2) birth death", 1, "pb_bd")
				#models[1]= cmd_1("> model of diversification (state 1): \n\t1) pure-birth\n\t2) birth death", 1, "pb_bd")


				rand=int(cmd_1("> Resample options: \n\t0) random sample\n\t1) evenly spaced trees", "", "rand_subsample"))
				if rand==0: rand="TRUE"
				else: rand="FALSE"
				break
		cmd= """cd "%s/bayesrate/r_functions" &&Rscript select_trees1.4.r "%s" %s %s %s""" % (self_path, infile, burnin, ntrees, rand)

	elif 'remove_outgroup' == model: 
		tree_list, input_file, name_file, path_dir=read_tree_nexus()
		for i, tree in enumerate(tree_list):
			if i<1: tree.print_plot(plot_metric='age')
		outgroup_names=raw_input("> outgroups (taxa names space separated): ")
		outgroup_names=outgroup_names.split()
		print outgroup_names
		t2=list()
		for i, tree in enumerate(tree_list):
			tree.prune_taxa_with_labels(outgroup_names)
		output="%s%s_ingroup.trees" % (path_dir, name_file)
		tree_list.write_to_path(output, schema='nexus')
		for i, tree in enumerate(tree_list):
			if i<1: tree.print_plot(plot_metric='age')
		print "The pruned trees were saved as:\n\t", output

	elif 'BMA' == model: 
		import model_averaging
		model_averaging.runBMA()
		cmd=""

	elif 'plot_rates' == model:
		while True:
			infile1= raw_input("> log-file: ")
			if infile1 in ('h','help'): help('plot_posterior')
			else: 
				infile, input_file, name_file, path_dir=fix_path(infile1)
				break
		cmd= """cd "%s/bayesrate/r_functions" &&Rscript rateHPD1.3.44.r "%s" """ % (self_path, infile)
				
	elif 'plot_rtt' == model:
		while True:
			infile1= raw_input("> log-file (*_rates): ")
			if infile1 in ('h','help'): help('plot_rtt')
			else: 
				infile, input_file, name_file, path_dir=fix_path(infile1)
				break
		cmd= """cd "%s/bayesrate/r_functions" &&Rscript rttplot1.2.r "%s" """ % (self_path, infile)

	elif 'plot_rate_shift' == model:
		while True:
			infile1= raw_input("> log-file (*_t_shift): ")
			if infile1 in ('h','help'): help('plot_shift')
			else: 
				infile, input_file, name_file, path_dir=fix_path(infile1)
				break
		cmd= """cd "%s/bayesrate/r_functions" &&Rscript rateshift_1.2.r "%s" """ % (self_path, infile)
	try: 
		if cmd!=0: os.system(cmd)
		#if cmd!=0: print cmd

	except: pass


#################################### MAIN MENU #####################################
a=0
while True: 
	try:
		def show_menu(): 
			print """> Select option:
	Parameter Estimation:
	  1) Speciation/extinction rates through time
	  2) Clade-specific rates and key innovation test
	  3) BiSSE-BMA (requires 'diversitree' R library)
	  4) MuSSE/GeoSSE/ClaSSE (requires 'diversitree', 'picante', 'optparse' R libraries)

	Marginal Likelihood (thermodynamic integration):
	  5) Speciation/extinction rates through time
	  6) Clade-specific speciation/extinction rates
	  7) BiSSE-BMA (requires 'diversitree' R library)

	Utilities:
	  8) Sub-sample trees
	  9) Remove outgroups
	 10) Joint posterior - Bayesian Model Averaging 
	 11) Plot posterior rates
	 12) Plot marginal rates through time (RTT)
	 13) Plot rate shifts\n"""
	  
		if a==0: show_menu()
		cmd = raw_input("bayesrate> ")
		#    0      1      2     3      4     5      6      7
		Q="menu", "quit", "q", "help", "h", "cite", "cit", "key"
		try: 
			cmd=int(cmd)
			if cmd<1 or cmd>13: raise()
		except: 
			if cmd not in Q: print "\tType 'help' for information on the commands that are available.\n"		
		
		if cmd in Q[1:3]: 
			if platform.system() == "Darwin": sys.stdout.write("\x1b]2;\x07")
			print "\n"
			quit()
		if cmd==Q[0]: show_menu()
		if cmd in Q[3:5]: help("main")
		if cmd in Q[5:7]: help("cite")
		if cmd==Q[7]: menu_opt(6)
		elif cmd>=1 and cmd<=12: menu_opt(cmd)
		a=1
	except(KeyboardInterrupt): 
		if platform.system() == "Darwin": sys.stdout.write("\x1b]2;BayesRate\x07")
		a=0
		print "\n\tType 'quit' to exit the program, 'help' for information on the commands that are available.\n"		

