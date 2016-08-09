#!/usr/bin/env python
# Created by Daniele Silvestro on 27/08/2010. => dsilvestro@senckenberg.de 
import sys
import os
import os.path
import platform


def help(arg):
	if arg=='nexus': print """\tInput a tree file in NEXUS format containing dated phylogenies (e.g. from BEAST output). 
	The trees can be sub-sampled and the burnin excluded using the Utility 'Sub-sample trees'.
	
	To load a tree file drag-and-drop the file onto this window.\n"""
	elif arg=='nexus_subsample': print """\tInput a tree file in NEXUS format containing dated phylogenies (e.g. from BEAST output). 
	To load a tree file drag-and-drop the file onto this window.
	Use 'Ctrl-C' to return to BayesRate main menu.\n"""
	elif arg=='no_tree_subsample': print """\tSpecify the number of trees you want to sample.\n"""
	elif arg=='rand_subsample': print """\tTrees can be sampled randomly (0) or evenly spaced (1).\n"""
	
	elif arg=='model': print """> Enter a number between 1 and 4 to select among models of diversification: 
	0: Yule process with constant rate
	1: Birth-Death process with constant rates
	2: Yule or Birth-Death with continuously varying rates: BOTHVAR, SPVAR, EXVAR (Rabosky & Lovette 2008)
	3: n-rate pure-birth: Yule process with rate shifts (Rabosky 2006)
	4: Birth-Death with incomplete taxon sampling (Yang & Rannala 1997)\n"""

	elif arg=='part_file': print """\tYou can define the data partitions interactively or via loading a partition file.
	The latter option allows a wider selection of settings, please refer to the sample file 
	for the format of this file.\n"""

	elif arg=='part_file_f': print """\tTo load a partition file drag-and-drop the file onto this window.
	Please refer to the example file and the manual for details on the format required.\n"""
	
	elif arg=='mcmc_bisse_tbl': print """\tTo load a trait file drag-and-drop the file onto this window.
	The trait file should be a table (e.g. csv) including taxa names in the first column and 
	the values of the binary trait in the second column (with headers).
	Please refer to the example file for further details on the format required.\n"""
	
	elif arg=='tax_sampling_musse': print """\tSpace-separated sampling fraction for each character state.\n"""

	elif arg=='rep': print """\tNumber of trees analyzed by the MCMC. Run on multiple trees to incorporate
	phylogenetic uncertainty in the rate estimation.\n"""

	elif arg=='runs': print """\tSpecify the number of 'heated' MCMC chains that will be run to estimate
	the marginal likelihood and the distribution of the 'temperatures' applied.
	The temperatures can be uniformly distribution as in Lartillot and Philippe (2006) or Beta distributed
	as in Xie et al. (2011). Independently of the distribution, higher number of scaling classes 
	will yield a more accurate estimation of the marginal likelihood.\n"""
	elif arg=='iterations': print """\tSpecify the number of generations the MCMC will spend sampling each tree.\n"""
	elif arg=='tree_burnin': print """\tSpecify the proportion or total number of trees from you input file
	that should be discarded as burnin.\n"""
	elif arg=='MCMC_burnin': print """\tSpecify the number of MCMC generations to be discarded as burnin.\n"""
	elif arg=='print_freq': print """\tSpecify how frequently the state of the MCMC should be logged on screen (in generations)\n"""
	elif arg=='sample_freq': print "\tSpecify how often parameter values are written in the log file (e.g. every 100 generations)\n"
	elif arg=='mod_rates': print """\tSpecify how frequently the diversification parameters are updated during the MCMC (e.g. 0.5).
	It might only be important if acceptance probability is very low.\n"""
	elif arg=='mod_times': print """\tSpecify how frequently the shift time are updated during the MCMC (e.g. 0.5). 
	Set this value to zero to fix the time of rate shift. \n"""
	elif arg=='n_rates': print "\tSpecify the number of different rates you want assume (minimum 2).\n"
	elif arg=='pb_bd': print "\tSelect between Pure-Birth and Birth-Death models\n"
	elif arg=='sample_fraction': print """\tEnter the proportion of species that is included in the phylogeny 
	(e.g. 0.85;	can be clade-specific)\n"""
	elif arg=='adv_settings': print  "\tSelect an option to customize MCMC settings or priors. Type 'run' to start the MCMC.\n"
	elif arg=='tree_selection':  print """\tTrees can be sampled:
	0: consecutively
	1: evenly spaced
	2: at random
	Note that the latter option should not be used for marginal likelihood estimation.\n"""
	elif arg=='prior_r':  print """\tSelect the appropriate prior distribution for the net diversification rate.\n"""
	elif arg=='lam_r':  print """\tShape parameter (lambda) of the exponential distribution. The expectation is
	1/lambda, but the prior allows up to infinite rates\n"""
	elif arg=='prior_m':  print """\tSelect the appropriate prior distribution for the extinction fraction (uniform or beta two-parameters).\n"""
	elif arg=='lam_m_a':  print """\tShape parameter (alpha) of the beta distribution. The expectation is alpha/(alpha+beta).\n"""
	elif arg=='lam_m_b':  print """\tShape parameter (beta) of the beta distribution. The expectation is alpha/(alpha+beta).\n"""

	elif arg=='shift_1':  print """\tMinimum age for the rate shift\n"""
	elif arg=='shift_2':  print """\tMaximum age for the rate shift\n"""

	elif arg=='bothvar':  print """> Select model of diversification: \n\t
	0: decreasing speciation, increasing extinction (BOTHVAR)
	1: constant speciation, increasing extinction (SPVAR)
	2: decreasing speciation, constant extinction (EXVAR)
	3: decreasing speciation, absent extinction (Pure-Birth SPVAR)\n"""

	elif arg=='threads': print """\tThermodynamic integration analyses require to run several MCMCs at different 'temperatures'.
	This can be carried out sequentially on a single processor (opt. 0), or in parallel either 
	on a single processor (each MCMC on a different thread; opt. 1) or on different processors
	(opt. 2). The latter option will use up to one processor for each MCMC, if available.\n"""
	
	elif arg=='n_partitions':  print """\tSpecify how many clades you would like to define. You need to have at least two clades, 
	otherwise use the spec/ext rates through time option. The rates can be constrained to be equal
	across clades by linking the partitions.\n """
	elif arg=='unlink_part':  print """Specify whether the estimation of diversification parameters of this clade should be linked 
	to those of the previous clade (i.e. have the same rates), or whether the parameters should be estimated individually."""
	elif arg=='members_clade':  print """\tTo define the members of the clade, either enter the number of all species 
	from the list above (e.g. 1-9 15 19-24), or enter their names (e.g. taxon1 taxon2 taxon4). 
	Note that numbers or names should be space separated.\n"""
	elif arg=='clade_key_inn': print """\tPerforms clade-specific analysis based on predefined clades (partition file is required).
	Option 'key innovation test' runs an estimation of the predicted species richness as 
	described by Moore and Donoghue (2009). Refer to the manual for more details 
	regarding the models implemented and the format of the partition file.\n"""
	
	elif arg=='log_BMA': print """\tInput a log file (from 'cold' MCMC chains); type 'BMA' to stop adding log-files.
	To load a log file drag-and-drop the file onto this window.
	Use 'Ctrl-C' to return to BayesRate main menu.\n"""
	
	elif arg=='BMA_ml': print """\tEnter the model marginal likelihood obtained through thermodynamic
	integration (use log units e.g. -124.65).\n"""
	
	elif arg=='plot_posterior': print """\tInput a log file to plot posterior estimates of the diversification parameters
	and calculate 95% HPDs. The results will be saved in a pdf file.
	To load a log file drag-and-drop the file onto this window.
	Use 'Ctrl-C' to return to BayesRate main menu.\n"""

	elif arg=='plot_rtt': print """\tInput a log file to plot the variation of diversification rates through time.
	The plot will be saved in a pdf file, the mean rate and 95% HPDs per Myr
	will be saved in a text file.
	To load a log file drag-and-drop the file onto this window. 
	Use 'Ctrl-C' to return to BayesRate main menu. \n"""

	elif arg=='plot_shift': print """\tInput a log file to plot the posterior distribution of the timing of each rate-shift
	(if present). The plot will be saved in a pdf file.
	To load a log file drag-and-drop the file onto this window. 
	Use 'Ctrl-C' to return to BayesRate main menu\n"""

	elif arg=='main': print """	Choose among the available options by entering a number between 1 and 12. 

	The following options are available in the main menu:

	Parameter Estimation: estimating diversificaton parameters under different 
	models (1-3). These options include fast estimation of the parameters, 
	but do not provide the model's marginal likelihood. Thus, they can not be 
	used to compare the fit of different evolutionary models. 

	Marginal Likelihood (thermodynamic integration): calculating a model's marginal
	likelihood via thermodynamic integration (TDI;  4-6). Please note that 
	these options open a new Terminal window; the TDI is performed via parallel 
	MCMCs that run over multiple processors, if available. 

	Utilities: These options (7-12) provide additional tools to prepare the input 
	data and plot the results of the analyses (pdf files will be generated).
	
	Please report any bugs or suggestions to the authors.

	Type 'help' or 'h' at any point while setting up an analysis to get information
	specific to the settings/parameters you are asked to define. 
	
	Type 'quit' to exit the program.
	Type 'cite' for the appropriate citation of the program
	Use 'Ctrl-C' to interrupt an analysis and return to BayesRate main menu.
	"""
	
	elif arg=='cite': print """	The program BayesRate is decribed in the paper:
	
	"Silvestro D., Schnitzler J. & Zizka G. (2011) A Bayesian framework to estimate
	 diversification rates and their variation through time and space. 
	 BMC Evolutionary Biology, 11, 311"
	
	
	The BiSSE-BMA method is described in:
	
	"Silvestro D., Zizka G. & Schulte K. (2013) Disentangling the effects of key innovations
	 on the diversification of Bromelioideae (Bromeliaceae). 
	 Evolution, DOI: 10.1111/evo.12236."
	
	and based on the BiSSE method and the Diversitree R library:

	"Fitzjohn R. (2012) Diversitree: Comparative phylogenetic analyses of diversification in R.
	 Methods Ecology Evolution, 3, 1084-1092",
	
	"Maddison W.P., Midford P.E. & Otto S.P. (2007) Estimating a binary character's effect on
	 speciation and extinction. 
	 Systematic Biology 56, 701-710."
	
	 
	Acknowledgments:
	
	Thanks to Ingo Michalak, Juriaan de Vos, Martha Serrano, and Anna Kostikova
	for contributions to the implementation and testing.
	
	"""
	
	else: print '\t no help available yet.'