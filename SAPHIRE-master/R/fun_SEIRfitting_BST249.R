###########################
# BST249 Final Project Code
###########################

## Parameter description
# init_set_list: initial settings produced by init_cond.R
# randomize_startValue: this function will randomly generate an initial condition if this argument is set to T
# startValue: if randomize_startValue is set to T, you can set your own initial condition using this argument
# output_ret: Whether to output parameter estimates output by MCMC
# skip_MCMC: This is meant for redrawing all results without rerunning MCMC
SEIRfitting=function(init_sets_list, 
                     randomize_startValue=F, 
                     startValue=NA, 
                     output_ret=T, 
                     run_id=0, 
                     skip_MCMC=F, 
                     panel_B_R_ylim=4,
                     plot_combined_fig=T,
                     pars_density=default_pars_density,
                     pars_sampler=default_pars_sampler,
                     pars_name=c("b12", "b3", "b4", "b5", "r12", "delta3", "delta4", "delta5"),
                     calc_clearance=T,
                     n_burn_in=4000,
                     n_iterations=180000,
                     DRlevels = 2, 
                     adaptationInterval=500) {
   
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

   # Start timer
   ptc <- proc.time()
   
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
                         adapt = T, DRlevels = DRlevels, iterations = n_iterations, thin = 10)
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
            } else if (j == 3 & alpha[j-1] > 0 ){
               alpha[j] = metropolisRatio(proposalEval[j,1], proposalEval[j-1,1])
               temp <- metropolisRatio(mcmcSampler$proposalGenerator$returnDensity(proposal[2,], proposal[3,]), mcmcSampler$proposalGenerator$returnDensity(mcmcSampler$current, proposal[1,]))
               jumpProbab = metropolisRatio(proposalEval[j,1], mcmcSampler$currentLP) * temp * (1.0-alpha[j]) / (1.0-alpha[1]) 
            } else if (j == 4 & alpha[j-1] > 0 ){
               alpha[j] = metropolisRatio(proposalEval[j,1], proposalEval[j-1,1])
               temp <- metropolisRatio(mcmcSampler$proposalGenerator$returnDensity(proposal[3,], proposal[4,]), mcmcSampler$proposalGenerator$returnDensity(mcmcSampler$current, proposal[1,]))
               jumpProbab = metropolisRatio(proposalEval[j,1], mcmcSampler$currentLP) * temp * (1.0-alpha[j]) / (1.0-alpha[1]) 
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
         # Default adaptationDepth = 3000 and adaptationInterval = 500
         # History of the chain is used to adapt the covariance of the proposal distribution
         mcmcSampler$settings$adaptationInterval <- adaptationInterval
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
   
   
   
   
   
   
   
   ## Output results from original code
   summary_string=paste0(paste(pars_name, collapse = ","), "\n")
   
   par_str=list()
   for (i_par in 1:n_pars) {
      par_str[[i_par]]=paste0(round(mean(mcmc_pars_estimate[,i_par]),2), " (",
                              round(quantile(mcmc_pars_estimate[,i_par],0.025),2)," - " , 
                              round(quantile(mcmc_pars_estimate[,i_par],0.975),2), ")")
   }
   
   summary_string = paste0(summary_string, paste(par_str,collapse = ", "),"\n\n")
   
   estRt_mat <- apply(mcmc_pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_sets_list))
   
   summary_string = paste0(summary_string, paste0("stage",1:n_stage,collapse=","), "\n")
   
   r_str=list()
   
   if (n_stage>1) {
      for (i_stage in 1:n_stage) {
         r_str[[i_stage]]=paste0(round(mean(estRt_mat[i_stage,]),2), " (",
                                 round(quantile(estRt_mat[i_stage,],0.025),2)," - " , 
                                 round(quantile(estRt_mat[i_stage,],0.975),2), ")")
      }
   } else {
      r_str[[1]]=paste0(round(mean(estRt_mat),2), " (",
                        round(quantile(estRt_mat,0.025),2)," - " , 
                        round(quantile(estRt_mat,0.975),2), ")")
   }
   
   summary_string = paste0(summary_string, paste(r_str,collapse = ", "),"\n\n")
   
   if (calc_clearance) {
      clearance_date = Findzero(mcmc_pars_estimate, init_sets_list)
      
      summary_string = paste0(summary_string, paste(names(clearance_date), collapse = ", "))
      
      summary_string = paste0(summary_string, "\n", paste(clearance_date, collapse = ", "), "\n")
   }
   
   summary_string = paste0(summary_string, "\n", "Time elapsed: ", (proc.time()-ptc)[1] )
   write_file(summary_string, paste0("../output/summary_run_",run_id,".txt"))
   
   png(paste0("../output/par_hist_run_",run_id,".png"))
   par(mfrow = c(2, 4))
   for(i in 1:n_pars) {
      hist(mcmc_pars_estimate[, i], xlab = pars_name[i], main = "", col = "red")
      rm(i)
   }
   dev.off()
   
   png(paste0("../output/par_traj_run_",run_id,".png"), width=1000, height=500)
   par(mfrow = c(2, 4))
   for(i in 1:n_pars) {
      plot(1:nrow(mcmc_pars_estimate), mcmc_pars_estimate[, i], ylab = pars_name[i], xlab = "iter", main = "", type = "l")
      rm(i)
   }
   dev.off()
   
   if (plot_combined_fig) {
      SEIRplot(pars_estimate = mcmc_pars_estimate, file_name = run_id, init_settings = init_sets_list, panel_B_R_ylim = panel_B_R_ylim)
   }
   
   par(mfrow = c(1, 1))
}
