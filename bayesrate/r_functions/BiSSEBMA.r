#!/usr/bin/Rscript
# version 20130821

arg <- commandArgs(trailingOnly=TRUE)


mcmc_BiSSE <- function(
	infile1, # tree
	infile2, # trait table
	wd="default",
	TDI=0,
	out_file_stem="default",
	rho0=1,
	rho1=1,
	BD_model_0=1,
	BD_model_1=1,
	Trait_model=0,
	link_speciation=0,
	link_extinction=0,
	link_trait=0,
	prior_r=5,
	prior_a=1,
	prior_b=1,
	prior_q=1,
	trees=1,
	IT=20000,
	sampling_freq=100,
	print_freq=100,
	burnin=1000,
	win_1=0.5,
	win_2=0.25,
	win_3=0.05,
	categories=10,
	beta_shape=0.3,
	path_lib="default"
	) {	# end args	
	
	if (wd=="default"){
		wd=dirname(infile1)
		}else{setwd(wd)}
	
	if (out_file_stem=="default"){
		st=basename(infile2)
		out_file= paste(st,"_BiSSEBMA.log", sep="")
		}else{out_file= paste(out_file_stem,".log", sep="")}
			
	if (path_lib=="default") {
		library(diversitree)
		}else{library(diversitree, , lib.loc=path_lib)}

	options(warn=-1) # hide warnings
		
	rhos=c(rho0, rho1)  #sampling fractions for state 0 and 1
	#out_file="boh.log"   #!!change in each analysis
	#TDI=0     # 0 parameter estimation, 1 TD integration (not real estimates)
	#trees=5  # one tree for TDI is advised

	#burnin=0   # run test to check the proportion
	#sampling_freq=25
	#print_freq=10
	#IT=200  # mcmc iterations (per tree)
	#categories=20
	

	use_exp=1 # 0: uniform priors; 1: exponential priors (extinction fraction has always uniform prior 0-1)
	M=c(5,5,1,1,5,5) # max parameter values (if uniform priors)
	# if exp no limit (max =1000)
	if (prior_r>0){M[1:2]=1000}	
	if (prior_q>0){M[5:6]=1000}	
	#shape_prior=c(2,2,1,1) # shape parameters of exp prior (net diversification and q rates)
	#constraints=c() # constant sp. [2], ex. [4], q. [6] const sp + ex is 2,4
	#PB=c() # 3: pb clade 1, 4: pb clade 2, 5: q0->1 =0; 6: q1->0 =0  # pure birth / irreversible rate models

	PB=NULL
	if (BD_model_0==0){PB=append(PB,3)}
	if (BD_model_1==0){PB=append(PB,4)}
	if (Trait_model==1){PB=append(PB,6)}
	if (Trait_model==2){PB=append(PB,5)}

	constraints=NULL
	if (link_speciation==1){constraints=append(constraints,2)}
	if (link_extinction==1){constraints=append(constraints,4)}
	if (link_trait==1)     {constraints=append(constraints,6)}


	shape_prior=c(prior_r,prior_a,prior_b,prior_q)
	win_size=c(win_1,win_1,win_2,win_2,win_3,win_3)
	
	Tree=read.nexus(file=infile1)
	traits=read.csv(file=infile2, header=TRUE, sep="\t")  #one single trait or create a vector
	states<-traits[,2]
	names(states)<-as.character(traits[,1])

	update_parameter <- function(i,d,M) {
		if (length(i)==1) {
			ii=abs(i+(runif(length(i),0,1)-.5)*d)
			if (ii>M) {ii=abs((M-(ii-M)))}
			if (ii>M) {ii=i}
			}
		ii}
	if (TDI>0)	{
		K=categories-1.
		k=0:K # K+1 categories
		beta=k/K
		alpha=beta_shape # categories are beta distributed
		temps=rev(beta^(1./alpha))
	} else{temps=c(1)} 

	exp_prior = function(value,l,Mv) {
		if (l>0){
			l=1./l
			log(l)-l*(value)
			}else{log(1/Mv)}
	}

	beta_prior = function(a,b,value){
		dbeta(value, a, b, log = TRUE)
	}


	true_iteration=1
	cat(sprintf("it\tlikelihood\tprior\tacceptance\tl0\tl1\tm0\tm1\tr0\tr1\ta0\ta1\tq0\tq1\ttemp\ttree\n"), append=FALSE, file=out_file)
	cat(sprintf("it\tlikelihood\tprior\tacc\tl0\tl1\tm0\tm1\ttemp\ttree\n"))
	for (J in 1:length(temps)) { # Loop temperatures
		temperature=temps[J]
		d= win_size*(4 - 3*temperature)

		for  (t_index in 1:trees) { # loop trees
			
			if (class(Tree)=="multiPhylo"){current_tree=Tree[[J]]}
			else{current_tree=Tree}
			
			BISSE=make.bisse(current_tree, states, strict=F, sampling.f=rhos)
			pars=runif(min=.1,max=.5, 6) # initialize parameters

			if (temperature<1) {burnin=0}
			LIK=0  
			acc=0
			for (iteration in 1:(IT+burnin)) { # MCMC loop
				if (iteration==1){ 
					likA = BISSE(pars)[1]
					parsA=pars}
				pars=parsA
				ind=sample(1:6,1) 
				pars[ind]=update_parameter(pars[ind],d[ind],M[ind])
				pars[PB]=0
				if (length(constraints)>0) {
					if (is.element(2, constraints) & is.element(4, constraints)) {pars[constraints] = pars[constraints-1]}
					else if (is.element(2, constraints)) {pars[2]=(pars[1]/(1-pars[3]))*(1-pars[4])}
					else if (is.element(4, constraints)) {pars[4]=pars[3]*pars[1]/(pars[2]-pars[3]*pars[2]+pars[3]*pars[1])}			
					if (is.element(6, constraints)) {pars[6] = pars[5]}
					}
				#cat("pars:", pars)
				l=pars[1:2]/(1-pars[3:4]) 
				m=-pars[3:4]*pars[1:2]/(pars[3:4]-1) 
				#lik=BISSE(c(l,m,pars[5:6]))[1]
		
				# CALC PRIORS # exponential
				if (use_exp==1) { 
					prior=sum(exp_prior(pars[1:2],shape_prior[1],M[1:2])) + sum(exp_prior(pars[5:6],shape_prior[4],M[5:6])) + sum(beta_prior(shape_prior[2],shape_prior[3],pars[3:4]))
				} else{
					prior=sum(log(1/M))
					if (min(M-temp)<0) {prior = -Inf}
					}  # uniform

				# CALC LIKELIHOOD
				if (prior> -Inf) {
					lik=try(BISSE(c(l,m,pars[5:6])))
					if (is.na(lik) | (class(lik) == "try-error" )) {lik=-Inf}
					} else{lik= -Inf}

				if (iteration==1){priorA=prior}

				tryCatch(
				{
				if ( (lik-likA)*temperature + (prior-priorA) >= log(runif(1,0,1)) ){
					likA=lik
					priorA=prior
					parsA=pars
					acc =acc + 1.					
				}
				}
				,error = function(e) NULL
				)
				

				if (iteration %% print_freq ==0) {cat(sprintf("%s\t", round(c(true_iteration, likA, priorA, parsA, temperature, t_index),2)), "\n")}
			
				if (true_iteration %% sampling_freq ==0 & iteration>=burnin) {
					l=parsA[1:2]/(1-parsA[3:4])
					m=-parsA[3:4]*parsA[1:2]/(parsA[3:4]-1)
					LIK= LIK+likA
					cat(sprintf("%s\t", c(true_iteration, likA, priorA, acc/iteration, l, m, parsA, temperature, t_index)), "\n", append=TRUE, file=out_file)
				}

				if (iteration>burnin) {true_iteration = true_iteration+1}
				}
			}
		}
	
	
	if (TDI>0)	{
		
		
		out_marg= paste(out_file_stem,"_marginal.txt",sep="")
		
		d<-read.table(out_file, header=T) #, stringsAsFactors=F,sep="\t") #, row.names = NULL)
		trapezia<-aggregate(d$lik, list(d$temp),FUN="mean")
		stds<-aggregate(d$lik, list(d$temp),FUN="sd")

		ml=0
		for (i in 1:(dim(trapezia)[[1]]-1)){
			ml=ml+( (trapezia[i,2] + trapezia[(i+1),2])/2 * (trapezia[(i+1),1] - trapezia[(i),1]))
		}
		cat("\nmean log likelihoods:", trapezia[,2], file=out_marg, append=TRUE)
		cat("\ntemperatures:", unique(d$temp), file=out_marg, append=TRUE)
		cat("\nstd(lnL):", stds[,2], file=out_marg, append=TRUE)
		cat("\n\nLog Marginal Likelihood:", ml, file=out_marg, append=TRUE)
	

		cat("\n\n Log Marginal Likelihood:", ml)
		cat("
 The marginal likelihood can be used to compare different analyses and to perform model selection
 and hypothesis testing via Bayes factors. 

 The marginal likelihood was saved to file:", out_marg,"\n\n")
		}
	
		cat("\n")
	}	
	
mcmc_BiSSE( 
	as.character(arg[1]),       # infile1, # tree
	as.character(arg[2]),       # infile2, # trait table
	as.character(arg[3]),       # wd,
	as.integer(arg[4]),         # TDI,
        as.character(arg[5]),       # out_file,
	as.double(arg[6]),	    # rho0,
	as.double(arg[7]),          # rho1,
        as.integer(arg[8]),         # BD_model_0,
        as.integer(arg[9]),         # BD_model_1,
        as.integer(arg[10]),        # Trait_model,
        as.integer(arg[11]),        # link_speciation,
        as.integer(arg[12]),        # link_extinction,
        as.double(arg[13]),         # link_trait,
        as.double(arg[14]),         # prior_r,
        as.double(arg[15]),         # prior_a,
        as.double(arg[16]),         # prior_b,
        as.double(arg[17]),         # prior_q,
        as.integer(arg[18]),        # trees,
        as.integer(arg[19]),        # iterations,
        as.integer(arg[20]),        # sampling_freq,
        as.integer(arg[21]),        # print_freq,
        as.integer(arg[22]),        # burnin,
        as.double(arg[23]),         # win_1,
        as.double(arg[24]),         # win_2,
        as.double(arg[25]),         # win_3,
        as.integer(arg[26]),        # categories,
        as.double(arg[27]),         # beta_shape
	as.character(arg[28])	    # path to diversitree library
	)

