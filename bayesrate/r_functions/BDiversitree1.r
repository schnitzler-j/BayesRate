#!/usr/bin/env Rscript
# Written by:
# Jose Rodrigo Flores Espinosa <joserodrigo.floresespinosa@unil.ch>
# Daniele Silvestro <daniele.silvestro@unil.ch>
# In colaboration with:
# Martha Serrano <martha.serrano@unil.ch>


#####################################
#####################################
### B-Diversitree is a module that runs a bayesian implementation of the Musse, Geosse and Classe models present in the package Diversitree.
###  
###
###
###
#####################################
#####################################
# Loading Libraries Required

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("diversitree"))
suppressPackageStartupMessages(library("picante"))
  
#####################################
#####################################

#####################################
# Parsing Command Line Options
        
option_list <- list(  
                  
    make_option("--i", type="integer", default=1000000,
        help=("Number of iterations to run [default %default]."),
        metavar="Iterations"),
    
    make_option("--s", type="integer", default=100,
        help="Frequency of sampling [default %default]",
        metavar="Freq_of_Sampling"),
        
    make_option("--b", type="integer", default=1000,
        help=("Skip results of the first # of iterations [default %default]."),
        metavar="Burnin"),

    make_option("--p", type="integer", default=100,
        help=("Print frequency [default %default]."),
        metavar="Prnt_freq"),
    
    make_option("--d", type="double", default=0.05,
        help=("Window size update parameters [default %default]."),
        metavar="D-Range"),
    
    make_option("--r", type="integer", default=5,
        help=("Rate to use for the calculation of the prior [default %default]."),
        metavar="Rate"),

    make_option("--t", type="integer", default=1,
        help=("Number of trees [default %default]."),
        metavar="Tree"),

        make_option("--rho",  default="1 1 1",
            help=("Taxon sampling (in quotations, space separated) [default %default]."),
            metavar="Rate"),
    
    make_option("--c", default="NULL",
        help="Parameters to be constrained. Introduce \"lamdas\", \"mus\" or \"qs\". \n\t\t Alternatively a string representing the desired contraints can be introduced. \n\t\t E.g. For three states under a Musse model: \"1,2,2,4,5,6,7,8,9,10,8,12\" \n\t\t This indicates lamda3 ~ lamda2 and q31 ~ q13 [default %default]",
        metavar="Constraints")
                
    )

# number of traits is pending
# number of trees is pending
parser_object <- OptionParser(usage = "Usage: %prog [Options] [TREE] [STATES] [MODEL]\n\n[TREE] Path to a file containig one or more phylogenetic tree(s).\n[STATES] Path to a file containig the name of each species in the pylogenetic file, along with their respective state.\n[MODEL] Model to be used for the bayesian calculations. The model can be \"musse\", \"geosse\" or \"classe\". \n", option_list=option_list, description="B-Diversitree is a module that runs a bayesian implementation of Musse, Geosse and Classe models present in the R package Diversitree.")
opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE), positional_arguments=TRUE)

#####################################
# Handling User Input and/or Option Errors

if (length(opt$args) < 3){
   cat("Error message: At least one of the input files is missing and/or the model has not been specified.\n\n"); print_help(parser_object); quit(status=1)
} 


#####################################
# Defining Variables and Objects

#set.seed(2)
file_tree <- opt$args[1]
file_states <- opt$args[2]
model <- opt$args[3]

######################################
# Reading Input Files
tree=read.nexus(file=file_tree)



if (class(tree)=="multiPhylo"){current_tree=tree[[1]]
}else{current_tree=tree}

states_frame <- read.table(file=file_states, h=T, row.names=1)
matched_tree <- match.phylo.data(current_tree,states_frame) # Drop tips with no info
states_vector <- matched_tree$data[,1] # vector of states
names(states_vector) <- row.names(matched_tree$data) # adding names to each state
no_states <- length(unique(states_vector))
current_tree <- matched_tree$phy

#####################################
# Defining Variables and Objects
print(model)
if (model == "musse"){   
  initial_parameters <- starting.point.musse(current_tree, no_states, q.div=5, yule=FALSE)   # Starting conditions
} else if (model == "geosse") {
  initial_parameters <- starting.point.geosse(current_tree, eps=0.5) 
} else if (model == "classe") {
  initial_parameters <- starting.point.classe(current_tree, no_states, eps=0.5) 
}


d <- rep(opt$options$d, length(initial_parameters))
sampling_freq <- opt$options$s
print_freq<- opt$options$p
burnin <- opt$options$b
rate <- opt$options$r
iterations <- opt$options$i
constraint <- opt$options$c
nTREES <- opt$options$t  

rhos_str <- opt$options$rho
rhos <- as.numeric(strsplit(rhos_str," ")[[1]])
print(rhos)
acc=0


######################################
# Defining Functions 
#

# Function to calculate the prior
#
exp_prior = function(values, rate){
  new_values = dexp(values, rate=1/rate,log=T) # lamda set to 5 default
  return(new_values)
}

# Function to change all parameters of one kind
update_parameters <- function(accepted_state_values, indexes, d, M) {
  i = accepted_state_values[indexes]
  ii = NULL
  ii= abs(i+(runif(length(i),0,1)-.5)*d[indexes])
  #for (j in 1:length(i)){
  #  ii[j] = abs(i[j]+(runif(1,0,1)-.5)*d[indexes[j]])
  #}
  if (max(ii) > M) {ii=abs((M-(ii-M)))}
  accepted_state_values[indexes]=ii
  return(accepted_state_values)
}  

# Function to construct the constraint vector
#
  
get_constraint <- function(temp_state, constraint, model, no_states) {

  vec = rep(1:length(temp_state))

  if (constraint == "lamdas"){
  
    if ( (model == "musse") | (model == "geosse") ) { 
      vec_temp = rep(1, no_states)                      # lamdas == to number of states (in geosse this is always 3)
      vec[1:no_states] = vec_temp
    }
    if (model == "classe") {
      x = no_states*no_states*(no_states+1)/2           # lamdas == to n*n*(n+1)/2  # length(parameters) = (k + 3) * k * k/2)
      vec_temp = rep(1, x)
      vec[1:x] = vec_temp 
    }    
    
  } else if (constraint == "mus"){
  
    if (model == "musse") {
      vec_temp = rep(no_states+1, no_states)            # mus == to number of states
      vec[(no_states+1):(no_states*2)] = vec_temp       # length(parameters) = k(k+1)
    }
    if (model == "geosse") {                            # number of mus is always 2
      vec_temp = rep(no_states+1, 2)                    # length(parameters) is always 7
      vec[(no_states+1):(no_states+2)] = vec_temp
    }
    if (model == "classe") {                            # mus = to number of states  # length(parameters) = (k + 3) * k * k/2)
      x = no_states*no_states*(no_states+1)/2           # number of lamdas
      vec_temp = rep((x+1),(no_states))                 # from #_lamdas to #_lamdas + #_states
      vec[(x+1):(x+no_states)] = vec_temp
    }
        
  } else if (constraint == "qs"){
  
    if (model == "musse") {                                                  # qs == k(k+1) - (no_states*2)
      vec_temp = rep((no_states*2+1),(length(temp_state)-no_states*2))       # length(parameters) = k(k+1)
      vec[(no_states*2+1):(length(temp_state))] = vec_temp
    }
    if (model == "geosse") {                            # number of mus is always 2
      vec_temp = rep(3+2+1, 2)                          # length(parameters) is always 7
      vec[(3+2+1):(7)] = vec_temp                       # number of lamdas is always 3
    }                                                   
    if (model == "classe") {                                            # lamdas == to n*n*(n+1)/2  # length(parameters) = (k + 3) * k * k/2)
      x = no_states*no_states*(no_states+1)/2                           # mus = to number of states
      vec_temp = rep((x+no_states+1),(length(temp_state)-x-no_states))
      vec[(x+no_states+1):(length(temp_state))] = vec_temp
    }
    
  } else {
  
    vec = as.numeric(unlist(strsplit(constraint, ",")))
    
  }
  return(vec)
}



##################################
# Main Program

cat(paste("Writing to File: ", opt$args[1], "_", opt$args[3], "_mcmc.log\n", sep=""))

real_it=0
for (J in 1:nTREES){
	if (class(tree)=="multiPhylo"){current_tree=tree[[J]]}
	else{current_tree=tree}

	states_frame <- read.table(file=file_states, h=T, row.names=1)
	matched_tree <- match.phylo.data(current_tree,states_frame) # Drop tips with no info
	states_vector <- matched_tree$data[,1] # vector of states
	names(states_vector) <- row.names(matched_tree$data) # adding names to each state
	current_tree <- matched_tree$phy
	

	for (it in 1:iterations){
   
	  # First Iteration / Initial Conditions
	  #
	  if (it==1) { 
	     if (model == "musse"){   
	       initial_parameters <- starting.point.musse(current_tree, no_states, q.div=5, yule=FALSE)   # Starting conditions
	     } else if (model == "geosse") {
	       initial_parameters <- starting.point.geosse(current_tree, eps=0.5) 
	     } else if (model == "classe") {
	       initial_parameters <- starting.point.classe(current_tree, no_states, eps=0.5) 
	     }  
 	    current_state=initial_parameters
 	    accepted_state=current_state
 	    accepted_lik=-Inf
 	    accepted_prior=-Inf
    }

   	  if (real_it==0) { 		    
	    cat(c("Iteration", "likelihood", "prior", "acceptance", "tree", names(accepted_state)), "\n", append=FALSE, file=paste(opt$args[1], "_", opt$args[3], "_mcmc.log", sep=""), sep="\t")
	  }
  
	  # Likelihood Calculation
	  #
  
	  if (model == "musse"){
		  LIKELIHOOD <- make.musse(current_tree, states_vector, no_states, strict=T, sampling.f=rhos)
	  } else if (model == "geosse") {
		  LIKELIHOOD <- make.geosse(current_tree, states_vector, strict=T, sampling.f=rhos)          
	  } else if (model == "classe") {
	    LIKELIHOOD <- make.classe(current_tree, states_vector, no_states, strict=T, sampling.f=rhos)     
	  }  else if (model == "quasse") {
	    LIKELIHOOD <- make.quasse(current_tree, states_vector, states.sd,lambda=constant.x, mu=constant.x, sampling.f=rhos)     
          }
		# More models to implement .. 
  
	  # Updating Parameters 
	  # 
  
	  flag <- sample(c("lamdas","mus","qs"), 1)

	  if (model == "musse"){
  
	    if (flag == "lamdas") { indexes = c(1:no_states)
	    } else if (flag == "mus") { indexes = c((no_states+1):(no_states*2))
	    } else if (flag == "qs")  { indexes = c((no_states*2+1):(length(initial_parameters)) ) }

	  } else if (model == "geosse") {

	    if (flag == "lamdas") { indexes = c(1:no_states)
	    } else if (flag == "mus") { indexes = c((no_states+1):(no_states+2))
	    } else if (flag == "qs")  { indexes = c((no_states+3):(length(initial_parameters)) ) }

	  } else if (model == "classe") {
	    x = no_states*no_states*(no_states+1)/2
	    if (flag == "lamdas") { indexes = c(1:x)
	    } else if (flag == "mus") { indexes = c((x+1):(x+no_states))
	    } else if (flag == "qs")  { indexes = c((x+no_states+1):(length(initial_parameters)) )}
    
  	  } else if (model == "quasse") {
  	    x = no_states*no_states*(no_states+1)/2
  	    if (flag == "lamdas") { indexes = c(1:x)
  	    } else if (flag == "mus") { indexes = c((x+1):(x+no_states))
  	    } else if (flag == "qs")  { indexes = c((x+no_states+1):(length(initial_parameters)) )}
    
  	  }

	  # Updating the parameters
	  # 
  
	  temporal_state <- update_parameters(accepted_state, indexes, d, 1000)
 
	  # Updatig the parameters according to the constraints specified, if any. 
	  #
  
	  if (constraint == "NULL") {
	    current_state = temporal_state  
	  } else { 
	    constraint_vector <- get_constraint(temporal_state, constraint, model, no_states) 
	    current_state = temporal_state[constraint_vector]  
	  }

	  # Prior calculation
	  #
		current_prior <- sum(exp_prior(current_state, rate))
  
		# Lik Calculation
	  #
		current_lik <- try(LIKELIHOOD(c(current_state)))
		if (is.na(current_lik) | (class(current_lik) == "try-error" )) {current_lik = -Inf}
  
		# Acceptance Conditions
	  #
		if ( (current_lik-accepted_lik) + (current_prior-accepted_prior) >= log(runif(1,0,1)) | it==0 ){
			accepted_lik=current_lik
			accepted_prior=current_prior
			accepted_state=current_state
			acc = acc + 1
		}
		
		if (it==0) {acc=1}

	  # Appending % of accepted states
		if ( (it %% sampling_freq == 0) & (it >= burnin) ) {    
			cat(sprintf("%s\t", c(real_it, accepted_lik, accepted_prior, acc/it, J, accepted_state)), "\n", append=TRUE, file=paste(opt$args[1], "_", opt$args[3], "_mcmc.log", sep=""))
			real_it=real_it+sampling_freq
		}
		if ( (it %% print_freq == 0)) {    
			cat(sprintf("%s\t", round(c(real_it,it, accepted_lik, accepted_prior, acc/it, J, accepted_state),3)), "\n", sep="")
		}
      	  
	
	}	

}