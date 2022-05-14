## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file

code_root="C:\\Users\\Jonathan\\OneDrive - Harvard University\\School\\Harvard\\BST249\\Project\\BST249 Project Code\\SAPHIRE-master\\"

setwd(paste0(code_root, "scripts_main"))
library(BayesianTools)

##
source(paste0(code_root, "R\\fun_SEIRpred.R"))
source(paste0(code_root, "R\\fun_SEIRsimu.R"))
source(paste0(code_root, "R\\fun_SEIRfitting_BST249.R"))
source(paste0(code_root, "R\\fun_BTSetup_BST249.R"))
source(paste0(code_root, "R\\init_cond.R"))
source(paste0(code_root, "R\\fun_R0estimate.R"))
source(paste0(code_root, "R\\fun_SEIRplot.R"))
source(paste0(code_root, "R\\fun_Findzero.R"))
##

init_sets_list=get_init_sets_list(r0 = 0.23)
randomize_startValue=T
startValue=NA 
output_ret=T 
run_id="main_analysis"
skip_MCMC=F 
panel_B_R_ylim=4
plot_combined_fig=T
pars_density=default_pars_density
pars_sampler=default_pars_sampler
pars_name=c("b12","b3","b4","b5","r12","delta3","delta4","delta5")
calc_clearance=T
n_burn_in=3001
n_iterations=100000


require(profvis)
profvis({
  ## Initialize parameters of the model
  # Number of daily new cases in the dataset
  onset_obs <- init_sets_list$daily_new_case
  
  # Initialize values of SAPHIRE
  init_states <- init_sets_list$init_states
  
  # Number of parameters in the model - default is 8
  n_pars = length(pars_name)
  
  # Number of time intervals we are studying - default is 5
  n_stage = length(init_sets_list$stage_intervals)
  
  # Create the prior for the parameters (uses beta, uniform, and normal distributions)
  # Sets a lower, upper, and best value for each of the 8 parameters
  pars_prior <- createPrior(density = pars_density, sampler = pars_sampler, 
                            lower = init_sets_list$par_lower, upper = init_sets_list$par_upper)
  
  ## Run the MCMC sampler
  if (!skip_MCMC) {
    # Initialize the Bayesian object
    # Contains the likelihood, prior, and posterior set before
    bayesSEIR <- createBayesianSetup(loglh_func, prior = pars_prior)
    
    # Run the sampler until a suitable startValue with finite likelihood is found
    if (randomize_startValue) {  
      startValue=pars_sampler()
      while (is.infinite(loglh_func(startValue))) {
        startValue=pars_sampler()
      }
    }
    
    ## Delayed Rejection Adaptive Metropolis Hastings
    # Initialize settings from BT
    mh_settings = list(startValue = startValue, optimize = T,
                       adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10)
    mcmcSampler <- Metropolis(bayesianSetup = bayesSEIR, settings = mh_settings)
    
    # Retrieve sampling parameters 
    # Burn in default is 0, thin is set to 10
    iterations <- mcmcSampler$settings$iterations
    burnin <- mcmcSampler$settings$burnin
    thin <- mcmcSampler$settings$thin
    
    # Counters for the algorithm - both set to 0
    CounterFunEvals = 0
    CounterAccept = 0
    
    # Initialize the chain (history of coefficient estimates + log probability, log likelihood, and log posterior)
    lastvalue = counter = 1
    mcmcSampler$chain = rbind(mcmcSampler$chain, array(dim = c(floor((iterations-burnin)/thin),mcmcSampler$setup$numPars+3)))
    
    # Initialize storage for the algorithm based on delayed rejection levels and number of parameters
    alpha = rep(NA, mcmcSampler$settings$DRlevels)
    proposalEval = matrix( nrow = mcmcSampler$settings$DRlevels, ncol = 3)
    proposal = matrix( nrow = mcmcSampler$settings$DRlevels, ncol = mcmcSampler$setup$numPars)
    
    # Function to calculate the MH ratio
    metropolisRatio <- function(LP2, LP1){
      if( is.na(LP2 - LP1)) out = -Inf
      else out =   min(exp( (LP2 - LP1) ), 1)
      return(out)
    } 
    
    # Run the MCMC algorithm
    for (i in 1:iterations){
      # Flag to terminate the algorithm
      accepted = F
      
      # Apply the delayed rejection
      for (j in 1:mcmcSampler$settings$DRlevels){
        # Pass in current coefficient estimates + scaling factor (default is 1 and 0.5 for delayed rejection) to get a new proposal
        # Utilizes the sampling function, number of parameters, and scaling factor to generate the move
        # move = move * scalingFactor / sqrt(num_param)
        proposal[j,] = mcmcSampler$proposalGenerator$returnProposal(x = mcmcSampler$current, scale = mcmcSampler$settings$proposalScaling[j])
        
        # Using the new proposed value, calculate the posterior density
        proposalEval[j,] <- mcmcSampler$setup$posterior$density(proposal[j,], returnAll = T)
        CounterFunEvals <- CounterFunEvals+1
        
        # For the first iteration, do normal Metropolis-Hastings
        # Alpha is the minimum of exp(par1 - par2) or 1
        if (j == 1){
          alpha[j] = metropolisRatio(proposalEval[j,1], mcmcSampler$currentLP)
          jumpProbab = alpha[1]
          # For the second iteration, do delayed rejection
        } else if (j == 2 & alpha[j-1] > 0 ){
          alpha[j] = metropolisRatio(proposalEval[j,1], proposalEval[j-1,1])
          temp <- metropolisRatio(mcmcSampler$proposalGenerator$returnDensity(proposal[1,], proposal[2,]), mcmcSampler$proposalGenerator$returnDensity(mcmcSampler$current, proposal[1,]))
          jumpProbab = metropolisRatio(proposalEval[j,1], mcmcSampler$currentLP) * temp * (1.0-alpha[j]) / (1.0-alpha[j-1]) 
        }
        
        # Given the selected ratio, generate a random Unif(0,1) to decide whether to accept
        if (runif(1) < jumpProbab){
          # Set flag
          accepted = T
          
          # Store proposed move and its log posterior density value
          mcmcSampler$current = proposal[j,]
          mcmcSampler$currentLP = proposalEval[j,1]
          
          # Based on burn-in and thinning value, update the chain
          if((i > (lastvalue+burnin)) && (i %% thin == 0) ){
            counter <- counter+1
            mcmcSampler$chain[counter,] = c(proposal[j,], proposalEval[j,])
          }
          break
        }
      }
      
      # If we fail to reject, update the chain with the previous value
      if((accepted == F) && (i > (lastvalue+burnin)) && (i %% thin == 0)){
        counter <- counter +1
        mcmcSampler$chain[counter,] = mcmcSampler$chain[counter-1,]
      } 
      
      # Increment counter if accept
      if(accepted == T) CounterAccept <- CounterAccept+1
      
      # Adaptive proposal update
      # Default adapt setting is TRUE
      # Default adaptationNotBefore = 3000 and adaptationInterval = 500
      # History of the chain is used to adapt the covariance of the proposal distribution
      if(mcmcSampler$settings$adapt  == T & i > mcmcSampler$settings$adaptationNotBefore & i %% mcmcSampler$settings$adaptationInterval == 0 ){
        start = max(1, counter - mcmcSampler$settings$adaptationDepth)
        mcmcSampler$proposalGenerator = updateProposalGenerator(proposal = mcmcSampler$proposalGenerator, chain = mcmcSampler$chain[start:counter,1:mcmcSampler$setup$numPars], message = F)
      }
    }
    
    # Make sure chain has right size
    mcmcSampler$chain <- mcmcSampler$chain[1:counter,]
    
    mcmcSampler$codaChain = coda::mcmc(mcmcSampler$chain)
    mcmcSampler$funEvals <- CounterFunEvals
    mcmcSampler$acceptanceRate <- CounterAccept/CounterFunEvals
    mh_out <- mcmcSampler
    
    # Output the results
    mcmc_pars_estimate <- BayesianTools::getSample(mh_out, start = n_burn_in+2, thin = 1)  ## set start = 2002 as the burn in period
    mcmc_pars_estimate <- round(mcmc_pars_estimate, 3)
    
    colnames(mcmc_pars_estimate) <- pars_name
    
    if (output_ret) {
      write.table(mcmc_pars_estimate, paste0("../output/pars_est_run_",run_id,".txt"), quote = F, row.names = F, sep = "\t")
    }
  }
  # Skip the MCMC sampler if parameter set to T
  else {
    mcmc_pars_estimate = read.table(paste0("../output/pars_est_run_",run_id,".txt"), header = T)
    pars_name = names(mcmc_pars_estimate)
  }
})