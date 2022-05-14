############## SEIR functions ############## 
## Density and sampler functions used to create the prior
default_pars_density <- function(pars) {
   d_vec <- rep(NA, 8)
   ##b12, b3, b4, b5
   for(i in c(1:4)) {
      d_vec[i] <- dunif(pars[i], 0, 2, log = T)
   }
   ## r12
   r12 <- pars[5]
   d_vec[5] <- dbeta(r12, beta_shape1, beta_shape2, log = T)
   ## r3
   delta_3 <- pars[6]
   d_vec[6] <- dnorm(delta_3, delta_mean, delta_sd, log = T)
   ## r4
   delta_4 <- pars[7]
   d_vec[7] <- dnorm(delta_4, delta_mean, delta_sd, log = T)
   ## r5
   delta_5 <- pars[8]
   d_vec[8] <- dnorm(delta_5, delta_mean, delta_sd, log = T)
   ##
   return(sum(d_vec))
}

default_pars_sampler <- function(n = 1) {
   s_vec <- matrix(NA, n, 8)
   ## b12, b3, b4, b5
   for(i in c(1:4)) {
      s_vec[, i] <- runif(n, 0, 2) 
   }
   ## r12 
   r12 <- rbeta(n, beta_shape1, beta_shape2)
   s_vec[, 5] <- r12
   ## r3
   s_vec[, 6] <- rnorm(n, delta_mean, delta_sd)
   ## r4
   s_vec[, 7] <- rnorm(n, delta_mean, delta_sd)
   ## r5
   s_vec[, 8] <- rnorm(n, delta_mean, delta_sd)
   return(s_vec)
}

## Log-likelihood function
loglh_func <- function(pars){
   ypred <- SEIRpred(pars, init_settings = init_sets_list)
   ypred <- ypred[, "Onset_expect"]
   onset_obs <- init_sets_list$daily_new_case
   
   # meant to suppress warnings when ypred is negative
   suppressWarnings(p <- dpois(onset_obs, ypred))
   
   if(any(p == 0) || any(is.nan(p))){
      logL <- -Inf
   }else{
      logL <- sum(log10(p))
   }
   return(logL)
}










############## BT Functions ############## 
## Create prior distribution function from BT
createPrior <- function(density = NULL, sampler = NULL, lower = NULL, upper = NULL, best = NULL){
   if(is.null(best) & ! is.null(lower) & ! is.null(upper)) best = (upper + lower) / 2
   
   parallelDensity<- function(x){
      if (is.vector(x)) return(density(x))
      else if(is.matrix(x)) return(apply(x, 1, density))
   }
   
   if(!is.null(sampler)){
      npar <- length(sampler())
      parallelSampler <- function(n=NULL){
         if(is.null(n)) out = sampler()
         else{
            if (npar == 1) out = matrix(replicate(n, sampler()))
            else if (npar >1) out = t(replicate(n, sampler(), simplify = T))
         } 
         return(out)
      }
   } 
   
   checkPrior <- function(x = NULL, z = FALSE){
      if(is.null(x)) x <- parallelSampler(1000)
      if(is.function(x)) x <- x()
      if(!is.matrix(x)) x <- parallelSampler(1000)
      check <- parallelDensity(x)
   }
   
   out<- list(density = parallelDensity, sampler = parallelSampler, lower = lower, upper = upper, best = best, originalDensity = density, checkStart = checkPrior)
   class(out) <- "prior"
   return(out)
}

## Create likelihood from BT
createLikelihood <- function(likelihood, names = NULL, parallel = F, catchDuplicates=T, 
                             sampler = NULL, parallelOptions = NULL){
   
   # check if point-wise likelihood available
   pwLikelihood = if ("sum" %in% names(as.list(args(likelihood)))) TRUE else FALSE
   
   catchingLikelihood <- function(x, ...){
      out <- tryCatch(
         {
            y = likelihood(x, ...)
            if (any(y == Inf | is.nan(y) | is.na(y) | !is.numeric(y))){
               message(paste("BayesianTools warning: positive Inf or NA / nan values, or non-numeric values occured in the likelihood. Setting likelihood to -Inf.\n Original value was", y, "for parameters", x, "\n\n "))
               y[is.infinite(y) | is.nan(y) | is.na(y) | !is.numeric(y)] = -Inf
            }
            y 
         },
         error=function(cond){
            cat(c("Parameter values ", x, "\n"))
            message("Problem encountered in the calculation of the likelihood with parameter ", x, "\n Error message was", cond, "\n set result of the parameter evaluation to -Inf ", "ParameterValues ")
            return(-Inf)
         }
      )
      return(out)
   }
   
   # initalize cl 
   cl <- NULL
   
   if (parallel == T | parallel == "auto" | is.numeric(parallel)) {
      tmp <- generateParallelExecuter(likelihood, parallel, parallelOptions) 
      parallelLikelihood <- tmp$parallelFun
      cl <- tmp$cl
      parallel = T
   }
   
   
   parallelDensity<- function(x, ...){
      if (is.vector(x)) return(catchingLikelihood(x, ...))
      else if(is.matrix(x)){
         if(catchDuplicates == TRUE){
            # Check for the rows that are not duplicated
            wn <- which(!duplicated(x))
            if(length(wn) <2) {
               return(parallelLikelihood(x, ...)) }
            else {
               # Define a output vector 
               out1 <- rep(0,length=nrow(x))
               
               # Run the likelihood function for unique values
               if (parallel == "external"){ 
                  out1[wn]<-likelihood(x[wn,], ...)
               }
               else{
                  if (parallel == T){ 
                     out1[wn]<-parallelLikelihood(x[wn,], ...)
                  }
                  else{
                     out1[wn]<-apply(x[wn,], 1, likelihood, ...)   
                  }
               }
               # Copy the values for the duplicates
               for(i in 1:length(out1)){
                  if(out1[i] != 0) next
                  else{
                     same <- numeric()
                     for(k in 1:length(out1)){
                        if(all(x[k,]== x[i,])){
                           same <- c(same,k)
                        }
                     }
                     out1[same[-1]] <- out1[same[1]]
                  }
               }
               
               return(out1)
            }}
         
         else{
            if (parallel == "external") return(likelihood(x, ...))
            else if (parallel == T){
               return(parallelLikelihood(x, ...))}
            else return(apply(x, 1, likelihood, ...))   
            
         }
      }
      else stop("parameter must be vector or matrix")
   }
   out<- list(density = parallelDensity, sampler = sampler, cl = cl, pwLikelihood = pwLikelihood, parNames = names)
   class(out) <- "likelihood"
   return(out)
}

## Create posterior from BT
createPosterior <- function(prior, likelihood){
   
   posterior <- function(x, returnAll = F){
      
      if (is.vector(x)){
         priorResult = prior$density(x) # Checking if outside prior to save calculation time
         if (! (priorResult == -Inf)) ll = likelihood$density(x)
         else ll = -Inf
         if (returnAll == F) return(ll + priorResult)    
         else return(c(ll + priorResult, ll, priorResult)) 
         
      } else if(is.matrix(x)){
         
         priorResult = prior$density(x) # Checking first if outside the prior to save calculation time
         feasible <- (! priorResult == -Inf)
         if (dim(x)[2] == 1) llResult <- likelihood$density(matrix(x[feasible, ], ncol = 1))
         else{
            if(TRUE %in% feasible) llResult <- likelihood$density(x[feasible, ])
            else llResult <- -Inf 
         }
         post = priorResult
         ll = priorResult
         ll[!feasible] = NA
         ll[feasible] = llResult
         post[feasible] = post[feasible] + llResult
         post[!feasible] = -Inf
         if (returnAll == F) return(post)    
         else{
            out <- cbind(post, ll, priorResult)
            colnames(out) = c("posterior", "likelihood", "prior")
            return(out)    
         } 
      }
      else stop("parameter must be vector or matrix")
   }
   
   out<- list(density = posterior)
   class(out) <- "posterior"
   return(out)
}

## Create setup for MCMC runs from BT
createBayesianSetup <- function(likelihood, 
                                prior = NULL,
                                priorSampler = NULL,
                                parallel = FALSE,
                                lower= NULL,
                                upper = NULL,
                                best = NULL,
                                names = NULL,
                                parallelOptions = list(variables = "all", packages = "all", dlls = NULL), 
                                catchDuplicates = FALSE,
                                plotLower = NULL,
                                plotUpper = NULL,
                                plotBest = NULL
){
   
   model <- NULL
   
   # INPUTS CHECKS
   if(is.null(parallelOptions)) parallelOptions <- list(variables = "all", packages = "all", dlls = "all")
   
   # PRIOR CHECKS
   priorClass = NULL
   if ("prior" %in% class(prior)) {
      priorClass = prior
      
   } else if (inherits(prior,"bayesianOutput")) {
      priorClass = createPriorDensity(prior)
      
   } else if ("function" %in% class(prior)) {
      if ("function" %in% class(priorSampler)) priorClass = createPrior(prior, priorSampler)
      else if (!is.null(lower) && !is.null(upper)) priorClass = createPrior(prior, lower=lower, upper=upper, best=best)
      else stop("If prior is a function, priorSampler or lower/upper is required")
      
   } else if (is.null(prior)) {
      # TODO: deprecate this
      # checks for NULL for lower/upper are already done at begin of function
      priorClass = createUniformPrior(lower = lower, upper = upper, best = best)
      
   } else stop("wrong input for prior")
   
   
   # LIKELIHOOD CHECKS
   if ("likelihood" %in% class(likelihood)) {
      likelihoodClass = likelihood
   } else if ("function" %in% class(likelihood)) {
      likelihoodClass = createLikelihood(likelihood, parallel = parallel, parallelOptions = parallelOptions, catchDuplicates = catchDuplicates)
   } else {
      stop("likelihood must be an object of class likelihood or a function")
   }
   pwLikelihood = likelihoodClass$pwLikelihood
   
   # GET NUMBER OF PARAMETERS
   numPars = length(priorClass$sampler())
   
   # CREATE POSTERIOR
   posteriorClass = createPosterior(priorClass,likelihoodClass)
   
   # CHECK FOR PLOTTING PARAMETERS
   if (is.null(plotLower)) plotLower <- priorClass$lower
   if (is.null(plotUpper)) plotUpper <- priorClass$upper
   if (is.null(plotBest)) plotBest <- priorClass$best
   
   if (is.null(plotLower) | is.null(plotUpper) | is.null(plotBest))
      print("Info is missing upper/lower/best. This can cause plotting and sensitivity analysis functions to fail. If you want to use those functions provide (plot)upper/lower/best either in createBayesianSetup or prior")
   
   # CHECK NAMES
   if (is.null(names)) {
      if (!is.null(priorClass$parNames)) names = priorClass$parNames
      else if (!is.null(likelihoodClass$parNames)) names = likelihoodClass$parNames
      else if (numPars > 0) names = paste("par", 1:numPars)
   }
   
   # CONSTRUCT OUTPUT
   info <- list(priorLower = priorClass$lower, priorUpper = priorClass$upper, priorBest = priorClass$best,
                plotLower = plotLower, plotUpper = plotUpper, plotBest = plotBest,
                parNames = names, numPars = numPars)
   out <- list(prior = priorClass, likelihood = likelihoodClass, posterior = posteriorClass,
               names = names, numPars = numPars, model = model, parallel = parallel, pwLikelihood = pwLikelihood, info = info)
   class(out) <- "BayesianSetup"
   
   return(out)
}

## Proposal generating function setup
setupStartProposal <- function(proposalGenerator = NULL, bayesianSetup, settings){
   
   # Proposal
   range = (bayesianSetup$prior$upper - bayesianSetup$prior$lower) / 50
   
   if(is.null(settings$startValue)) settings$startValue = (bayesianSetup$prior$upper + bayesianSetup$prior$lower) / 2
   
   if (length(range) != bayesianSetup$numPars) range = rep(1,bayesianSetup$numPars)
   
   if(is.null(proposalGenerator)){
      proposalGenerator = createProposalGenerator(range, gibbsProbabilities = settings$gibbsProbabilities)
   }
   
   ####### OPTIMIZATION
   
   if (settings$optimize == T){
      if(is.null(settings$message) || settings$message == TRUE){
         cat("BT runMCMC: trying to find optimal start and covariance values", "\b")
      }
      
      target <- function(x){
         out <- bayesianSetup$posterior$density(x)
         if (out == -Inf) out =  -1e20 # rnorm(1, mean = -1e20, sd = 1e-20)
         return(out)
      }
      
      try( {
         if(bayesianSetup$numPars > 1) optresul <- optim(par=settings$startValue,fn=target, method="Nelder-Mead", hessian=F, control=list("fnscale"=-1, "maxit" = 10000))
         else optresul <- optim(par=settings$startValue,fn=target, method="Brent", hessian=F, control=list("fnscale"=-1, "maxit" = 10000), lower = bayesianSetup$prior$lower, upper = bayesianSetup$prior$upper)      
         settings$startValue = optresul$par
         hessian = numDeriv::hessian(target, optresul$par)
         
         
         proposalGenerator$covariance = as.matrix(Matrix::nearPD(MASS::ginv(-hessian))$mat)
         #proposalGenerator$covariance = MASS::ginv(-optresul$hessian)
         
         # Create objects for startValues and covariance to add space between values
         startV <-covV <- character()
         
         for(i in 1:length(settings$startValue)){
            startV[i] <- paste(settings$startValue[i], "")
         } 
         for(i in 1:length( proposalGenerator$covariance)){
            covV[i] <- paste( proposalGenerator$covariance[i], "")
         } 
         
         if(is.null(settings$message) || settings$message == TRUE){
            message("BT runMCMC: Optimization finished, setting startValues to " , 
                    startV, " - Setting covariance to " , covV)
         }
         
         proposalGenerator = updateProposalGenerator(proposalGenerator)
         
      }
      , silent = FALSE)
   }  
   out = list(proposalGenerator = proposalGenerator, settings = settings)
   return(out)
}

## Create the proposal covariance matrix
createProposalGenerator <- function(
   covariance, # covariance matrix for the multivariate proposal
   gibbsProbabilities = NULL, #  changes
   gibbsWeights = NULL,
   otherDistribution = NULL,
   otherDistributionLocation = NULL, 
   otherDistributionScaled = F,
   message = F,
   method = "chol",
   scalingFactor = 2.38
) {
   
   # To provide the option of defining via sd of individual normal
   if (is.vector(covariance)) {
      covariance = diag(covariance^2) 
   }
   
   if(ncol(covariance) == 0 && nrow(covariance) == 0) covariance = 1
   
   if(is.null(otherDistribution)) numberOfParameters = max(1,nrow(covariance)) 
   else numberOfParameters = length(otherDistributionLocation)
   
   if(is.null(method) | numberOfParameters < 2){
      covarianceDecomp = NULL
      if(numberOfParameters > 1) samplingFunction = function() as.vector(mvtnorm::rmvnorm(n = 1, sigma = covariance))
      else samplingFunction = function() rnorm(n = 1, sd = sqrt(covariance)) 
   } else {
      covarianceDecomp = factorMatrice(covariance, method = method)
      samplingFunction = function() as.vector(getRmvnorm(n = 1, R = covarianceDecomp))
   }
   
   ## Assertions
   
   if(!is.null(otherDistribution)){
      stopifnot(class(otherDistribution) == "function")
      stopifnot(!is.null(otherDistributionLocation))
      if(is.numeric(otherDistributionLocation)) otherDistributionLocation = as.logical(otherDistributionLocation) 
      stopifnot(is.logical(otherDistributionLocation))      
      stopifnot(is.logical(otherDistributionScaled))  
   } 
   
   #scalingFactor = 2.38/sqrt(numberOfParameters) # CHECK ???
   #scalingFactorN = (2.38^2)/numberOfParameters 
   # note - scaling is 2.38 * sqrt because it is applied on the change, not directly on the sigma
   
   ##########################
   # Definition of proposal function
   
   returnProposal <- function(x, scale = 1){     
      
      # Possibility to mix with other distribution
      if(!is.null(otherDistribution)){
         move = rep(NA, numberOfParameters)
         move[otherDistributionLocation] = otherDistribution(x[otherDistributionLocation])
         move[!otherDistributionLocation] = samplingFunction()
      }else{
         move = samplingFunction()
      }
      
      ## Gibbs updates
      if (!is.null(gibbsProbabilities)) {
         
         nGibbs <- sample.int(length(x), size = 1, replace = F, prob = gibbsProbabilities)
         whichParametersLoc <- sample.int(length(x), nGibbs, replace = F, prob = gibbsWeights)
         move[! (1:numberOfParameters %in% whichParametersLoc)] = 0
      } else {
         nGibbs = numberOfParameters
      }
      
      if(!is.null(otherDistribution) & otherDistributionScaled == F){
         nGibbs = nrow(covariance) 
         move[!otherDistributionLocation] = move[!otherDistributionLocation] * scalingFactor / sqrt(nGibbs)
      }else{
         move = move * scalingFactor / sqrt(nGibbs)
      }
      
      newParam = x + move * scale
      
      return(newParam)
   }
   
   returnProposalMatrix <- function(x, scale = 1){
      
      numPar <- ncol(x)
      
      if (numPar == 1){
         out = matrix(apply(x, 1, returnProposal, scale = scale), ncol = 1) 
      } else {
         out = t(apply(x, 1, returnProposal, scale = scale)) 
      }
      return(out)
   }
   
   
   returnDensity <- function(x, y, scale = 1){
      if (!is.null(gibbsProbabilities) & !(is.null(otherDistribution)))stop("proposal density not implemented if Gibbs or other distribution is activated in the proposal. This error may appear if you have chosen both gibbs and delayed rejection in an MCMC algorith. This option is currently not implemented")
      
      sigmaDensity = scalingFactor^2 / numberOfParameters * covariance * scalingFactor^2
      if(length(sigmaDensity) > 1) dens = mvtnorm::dmvnorm(x, mean = y, sigma = sigmaDensity, log = T)
      else dens = dnorm(x, mean = y, sd = sqrt(sigmaDensity), log = T) 
      return(dens)
      
   }
   
   
   ##########################
   # Wrap up class fields
   
   classFields = list(
      covariance = covariance, 
      covarianceDecomp = covarianceDecomp,
      gibbsProbabilities = gibbsProbabilities, 
      gibbsWeights = gibbsWeights,
      otherDistribution = otherDistribution, 
      otherDistributionLocation = otherDistributionLocation, 
      otherDistributionScaled = otherDistributionScaled,
      returnProposal = returnProposal, 
      returnProposalMatrix = returnProposalMatrix, 
      returnDensity = returnDensity,
      updateProposalGenerator = updateProposalGenerator ,
      samplingFunction = samplingFunction
   )
   class(classFields) <- c("proposalGenerator")
   
   if(message == T){
      cat("Proposalgenerator created")
      print(classFields)
   }
   
   return(classFields)
}


## Helper function for proposal generator
factorMatrice <- function(sigma, method){
   if(method == "eigen") {
      ev <- eigen(sigma, symmetric = TRUE)
      if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
         warning("sigma is numerically not positive definite")
      }
      ## ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
      ## faster for large  nrow(sigma):
      t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
   }
   else if(method == "svd"){
      s. <- svd(sigma)
      if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
         warning("sigma is numerically not positive definite")
      }
      t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
   }
   else if(method == "chol"){
      R <- chol(sigma, pivot = TRUE)
      R[, order(attr(R, "pivot"))]
   }
}

## Helper function for proposal generator
getRmvnorm <- function(n=1, R){
   X <- matrix(rnorm(n * ncol(R)), nrow=n )%*%  R
   return(X)
}

## Helper function for proposal generator
updateProposalGenerator <- function(proposal,chain = NULL,  message = F, eps = 1e-10, manualScaleAdjustment = 1){
   
   if(!is.null(chain)){
      npar = ncol(chain)
      if(is.null(npar)) npar = 1
      
      if (npar > 1){
         covar = cov(chain) * manualScaleAdjustment
         covar = as.matrix(Matrix::nearPD(covar + diag(eps, npar))$mat)    
      }else{
         covar = var(chain) * manualScaleAdjustment
      }
      if(!any(is.na(covar))) proposal$covariance = covar
   }
   
   out <- createProposalGenerator(
      covariance = proposal$covariance, 
      gibbsProbabilities = proposal$gibbsProbabilities, 
      gibbsWeights = proposal$gibbsWeights,
      otherDistribution = proposal$otherDistribution, 
      otherDistributionLocation = proposal$otherDistributionLocation, 
      otherDistributionScaled = proposal$otherDistributionScaled
   )
   
   if(message == T){
      cat("Proposalgenerator settings changed")
      print(out)
   }
   return(out)
}

## Check and apply settings to the bayesianSetup object
checkBayesianSetup <- function(bayesianSetup, parallel = F){
   
   if(class(bayesianSetup) == "function"){
      if(is.null(parallel)) parallel = F
      bayesianSetup = createBayesianSetup(bayesianSetup, parallel = parallel)
   } 
   else if(class(bayesianSetup) == "BayesianSetup"){
      if(!is.null(parallel)) if(parallel == T & bayesianSetup$parallel == F) stop("parallel = T requested in sampler but BayesianSetup does not support parallelization. See help of BayesianSetup on how to enable parallelization")
   } 
   else stop("bayesianSetup must be class BayesianSetup or a function")
   
   return(bayesianSetup)
}

## Initialize default setting for MH MCMC
Metropolis <- function(bayesianSetup, settings){
   setup <- checkBayesianSetup(bayesianSetup, parallel = settings$parallel) # calling parallel will check if requested parallelization in settings is provided by the BayesianSetup
   settings$parallel = bayesianSetup$parallel # checking back - if no parallelization is provided, we use the parallelization in the BayesianSetup. We could also set parallel = F, but I feel it makes more sense to use the Bayesiansetup as default
   
   defaultSettings = list(proposalGenerator = NULL, 
                          consoleUpdates=100, 
                          burnin = 0,
                          adaptationInterval= 500, 
                          adaptationNotBefore = 3000,
                          proposalScaling = NULL, 
                          adaptationDepth = NULL, 
                          temperingFunction = NULL, 
                          proposalGenerator = NULL, 
                          gibbsProbabilities = NULL, 
                          currentChain = 1,
                          message = FALSE,
                          sampler = "Metropolis")
   settings = c(settings, defaultSettings)
   
   ## Parameter consistency checks 
   if(is.null(settings$adaptationDepth)){
      settings$adaptationDepth = settings$adaptationNotBefore
   } 
   
   # Decreasing scaling for DRAM by default
   if (is.null(settings$proposalScaling)) settings$proposalScaling = 0.5^(- 0:(settings$DRlevels -1))
   
   tmp <- setupStartProposal(proposalGenerator = settings$proposalGenerator, bayesianSetup = bayesianSetup, settings = settings)
   settings = tmp$settings
   proposalGenerator = tmp$proposalGenerator
   
   ####### CREATE CHAIN
   
   chain = array(dim = c(1,bayesianSetup$numPars+3))
   chain[1,1:bayesianSetup$numPars] = settings$startValue
   colnames(chain) = c(1:bayesianSetup$numPars, "LP", "LL", "LPr")
   chain[1, (bayesianSetup$numPars+1):(bayesianSetup$numPars+3)] = setup$posterior$density(settings$startValue, returnAll = T)
   
   current = settings$startValue
   currentLP = as.numeric(chain[1, (bayesianSetup$numPars+1)])
   
   ##### Sampling
   
   classFields = list(
      setup = setup,
      settings = settings,
      current = current,
      currentLP = currentLP,
      chain = chain, 
      proposalGenerator = proposalGenerator,
      funEvals = 0,
      acceptanceRate = 0
   )
   
   class(classFields) <- c("mcmcSampler", "bayesianOutput")
   return(classFields)
}
