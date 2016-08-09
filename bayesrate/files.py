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
import dendropy

###################################	
def prune_tree(retained_taxa, tree_list_t, verbose): # PYTHON 
	i=1
	retained_taxa=retained_taxa.split()
	pt=retained_taxa #[z.replace("_"," ") for z in retained_taxa]
	#print pt
	K=list()
	r=0
	#print tree_list
	for tree in tree_list_t:	
		tree.resolve_polytomies(update_bipartitions=True)
		tree.calc_node_ages(False)
		#print tree.as_string('newick')
		tree.retain_taxa_with_labels(pt)
		#print "ok"
		if i==1 and verbose==1: print tree.print_plot(plot_metric='age') #.as_string('newick')
		i+=1
		nd= tree.ageorder_node_iter(include_leaves=False, filter_fn=None, descending=True)
		ages=list()
		#print tree
		#print nd
		#print enumerate(nd)
		for n, node in enumerate(nd): 
	#		print n, node
			ages.append(node.age)
		K.append(ages)
		r += max(ages)
		#print "still ok"
	#print K[0]
	return array(K), len(ages)+1, r/len(tree_list_t)

def read_tree_nexus(): # PYTHON
	while True:
		infile= raw_input("> input tree-file (NEXUS): ")
		#try:
		if 2>1:
			if infile == "help" or infile=="h": raise(ValueError)
			infile, input_file, name_file, path_dir=fix_path(infile)
			print "\n reading trees..."
			tree_list = dendropy.TreeList.get_from_path(infile, schema="nexus", preserve_underscores=True) # , translate_dict=True
			print len(tree_list), "trees found"			
			for i, tree in enumerate(tree_list):
				if i<1: 
					tree.print_plot(plot_metric='age')
					#print len(tree.taxon_set)
			return tree_list, input_file, name_file, path_dir
			#if verbose==1: tree.print_plot()
				
			break
		#except(ValueError): help("nexus")
		#except: print "Tree format not recognized."

def get_node_ages(tree_list): # PYTHON
	K=list()
	#tree_resolved_polytomies=list()
	#for idx, tree in enumerate(tree_list):
		
	for idx, tree in enumerate(tree_list):
		#tree.calc_node_ages(False) # Adds an attribute called "age" to each node
		tree.resolve_polytomies(update_bipartitions=True)
		tree.calc_node_ages(ultrametricity_precision=0.001) # 
		nd= tree.ageorder_node_iter(include_leaves=False, filter_fn=None, descending=True)
		ages=list()
		for n, node in enumerate(nd): ages.append(node.age)
		K.append(ages)
	return array(K)

def get_root_ages(tree_infile): # PYTHON
	tree_infile, input_file, name_file, path_dir=fix_path(tree_infile)
	tree_list = dendropy.TreeList.get_from_path(tree_infile, schema="nexus")
	#tree_list = read_tree_nexus(tree_infile)
	roots=list()
	for idx, tree in enumerate(tree_list):
		#tree.calc_node_ages(False) # Adds an attribute called "age" to each node
		tree.calc_node_ages(ultrametricity_precision=0.001)
		nd= tree.ageorder_node_iter(include_leaves=False, filter_fn=None, descending=True)
		ages=list()
		for n, node in enumerate(nd): ages.append(node.age)
		roots.append(ages[0])
	return array(roots)

def load_treefile(): # PYTHON
	tree_list, input_file, name_file, path_dir=read_tree_nexus()		
	K=get_node_ages(tree_list)
	return K, input_file, name_file, path_dir
		
def parse_file(part_file, trees_ingroup, input_file, name_file, path_dir, run):
	if part_file==0: 
		while True:
			infile= raw_input("> partition file: ")
		#	infile="partition.txt"
			if infile in ('h','help', ''): help("part_file_f")
			else: break
	else: infile=part_file
	infile, input_file, name_file, path_dir=fix_path(infile)
	if run==0: print "\tparsing partition settings..."
	dataset, partition, rhos, models, part_sequence=list(), list(), list(), list(), list()
	size_clades, age_clades, name_clades=list(),list(),list()

	lists= "name_partitions", "sampling", "model", "part_sequence"
	parms= "trees", "iterations", "sampling_freq", "burnin", "output_file", "special_rule", \
	"mod_rates", "print_freq", "std_updates", "prior_r",  "prior_a(alpha)", "prior_a(beta)"
	parms=list(parms)                                    # single parameters
	parameters= name_clades, rhos, models, part_sequence # lists
	parameters=list(parameters)
	def_parms=list()

	for line in open(infile):
		line = line.partition('#')[0]
		line = line.rstrip()
		line = line.split()
		if len(line)<2: pass
		elif line[0] in lists:
			indx=[lists.index(line[0])]
			parameters[indx[0]] =(line[1:len(line)])
		elif line[0] in parms:
			indx=[parms.index(line[0])]
			try: parms[indx[0]]= int(line[1])
			except: parms[indx[0]]= line[1]

	for l in range(len(parameters)): 
		for i in range(len(parameters[l])): 
			try: parameters[l][i]=float(parameters[l][i])
			except: pass

	n_datasets=len(parameters[2])
	ntrees =        parms[0]
	iteration =     parms[1]
	sampling_freq = parms[2]
	burnin =        parms[3]
	output_file =   parms[4]
	special_rule =  parms[5]
	mod_rates =     parms[6]
	print_freq =    parms[7]
	std_updates =   parms[8]
	lam_r =			parms[9]
	a =				float(parms[10]) # beta distribution
	b =				float(parms[11])
	if a>0 and b>0:	lam_m = [a,b]
	else: lam_m=0
	#print lam_m, a, b, sum(lam_m)

	name_clades = list(parameters[0])
	taxa_list =   parameters[0]
	rhos =        parameters[1]
	models =      parameters[2]
	partition =   parameters[3]

	#print " n_datasets", n_datasets, "\n ntrees", ntrees, "\n name_clades", name_clades, \
	#"\n partition sequence", partition, "\n models", models, "\n sampling", rhos, "\n iterations", iteration, \
	#"\n beta parms", lam_m

	for line in open(infile): # get taxa names defining the clades
		line = line.partition('#')[0]
		line = line.rstrip()
		line = line.split(":")
		if len(line)<2: pass
		elif line[0].strip() in taxa_list:
			indx=[parameters[0].index(line[0].strip())]
			taxa_list[indx[0]] =(line[1]) #:len(line)])
	#print taxa_list

	trees=ntrees
	temp_trees=list()
	for n in range(n_datasets): temp_trees.append(dendropy.TreeList(trees_ingroup))
	
	for n in range(n_datasets):
		if run==0: print "\tdefining clade %s..." % (n)
		#subclade, ntrees=prune_tree(taxa_list[n], trees_ingroup, taxa_names, 0)
		#l,r=list(),0
		#for i in range(0,ntrees): 
		#	print branching_times(subclade[i])
		#	l.append(list(robjects.FloatVector(branching_times(subclade[i])))) # node ages of all trees
		#	r += max(l[i])
		#K, lenClade, root_age = array(l), len(l[0])+1, (r/ntrees)
		
		K, lenClade, root_age = prune_tree(taxa_list[n], temp_trees[n], 0)
		
		
		dataset.append(K)
		size_clades.append(lenClade)
		age_clades.append(root_age)
	#print dataset, lenClade, root_age #name_clades
	
	
	return n_datasets, dataset, size_clades, name_clades, age_clades, partition, models, rhos, \
	trees, iteration, sampling_freq, burnin, output_file, special_rule, mod_rates, 	print_freq,	std_updates, lam_r,	lam_m

def out_info(file, txt, verbose):
	file.writelines(txt)
	if verbose==0: print txt
	
def output(outfile, out_mcmc, run, log_screen): 
	outfile.writelines(out_mcmc)
	outfile.writelines("\n")
	if run==0 and log_screen !=0: print out_mcmc

def log_chain(logfile, wlog, log_mcmc):
	wlog.writerow(log_mcmc)
	logfile.flush()


def parse_bisse_file(tbl_file):
	while True:
		infile= raw_input("> BiSSE setting file: ")
		if infile in ('h','help', ''): help("part_file_bisse")
		else: break

	infile, input_file, name_file, path_dir=fix_path(infile)
	print "\tparsing BiSSE settings..."

	parms= ["output_file","rho_0","rho_1","BD_model_0","BD_model_1","Trait_model","link_speciation","link_extinction","link_trait", \
	"prior_r","prior_a(alpha)","prior_a(beta)","prior_q","trees","iterations",\
	"sampling_freq","print_freq","burnin","win_1","win_2","win_3","categories","beta_shape","path_lib"]

	for line in open(infile):
		line = line.partition('#')[0]
		line = line.rstrip()
		line = line.split()
		if len(line)<2: pass
		elif line[0] in parms:
			indx=[parms.index(line[0])]
			try: parms[indx[0]]= float(line[1]) # numbers
			except: parms[indx[0]]= line[1]     # strings
		
	S=""	
	for j in parms: S+= "%s " % (j)
	return S




def fix_path(infile1):
	if platform.system() == "Darwin":
		infile1=infile1.replace("\ ",'_spazio_')
		infile1=infile1.replace(' ','')
		infile1=infile1.replace('_spazio_',' ')
	infile=infile1.replace('"', '')
	infile=infile.replace("""' """, "'")	#else: infile=infile1
	infile=infile.replace("""'""", "")
	if platform.system() == "Windows" or platform.system() == "Microsoft": infile=infile1.replace('"', '')
	input_file = os.path.basename(infile1)		 # name input file
	name_file = os.path.splitext(input_file)[0]  # file name without extension
	path_dir = "%s/" % os.path.dirname(infile)   # path directory input file
	if path_dir=="/": path_dir=""
	return infile, input_file, name_file, path_dir




