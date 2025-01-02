sampler_RW <- nimbleFunction(
  name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    logScale            <- extractControlElement(control, 'log',                 FALSE)
    reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodes
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory      <- c(0, 0)   ## scaleHistory
    acceptanceHistory <- c(0, 0)   ## scaleHistory
    saveMCMChistory <- getNimbleOption('MCMCsaveHistory')
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
    if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
    if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
    if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run = function() {
    currentValue <- model[[target]]
    propLogScale <- 0
    if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
    propValue <- currentValue * exp(propLogScale)
    } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
    if(reflective) {
      lower <- model$getBound(target, 'lower')
      upper <- model$getBound(target, 'upper')
      while(propValue < lower | propValue > upper) {
        if(propValue < lower) propValue <- 2*lower - propValue
        if(propValue > upper) propValue <- 2*upper - propValue
      }
    }
    model[[target]] <<- propValue
    logMHR <- model$calculateDiff(target)
    if(logMHR == -Inf) {
      jump <- FALSE
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
    } else {
      logMHR <- logMHR + model$calculateDiff(calcNodesNoSelf) + propLogScale
      jump <- decide(logMHR)
      if(jump) {
        ##model$calculate(calcNodesPPomitted)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      } else {
        nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      }
    }
    if(adaptive)     adaptiveProcedure(jump)
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        if(reflective) {
          lower <- model$getBound(target, 'lower')
          upper <- model$getBound(target, 'upper')
          if(scale >= 0.5*(upper-lower)) {
            scale <<- 0.5*(upper-lower)
          }
        }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    setScale = function(newScale = double()) {
      scale         <<- newScale
      scaleOriginal <<- newScale
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
        return(numeric(1, 0))
      }
    },
    ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)