oneStepPredict <- function(obj,
                           ## Names of data objects (not all are optional)
                           observation.name = NULL,
                           data.term.indicator = NULL,
                           lower.cdf.indicator = NULL,
                           upper.cdf.indicator = NULL,
                           method=c(
                               "fullGaussian",
                               "oneStepGaussian",
                               "oneStepGeneric",
                               "cdf"),
                           subset = NULL,
                           conditional = NULL,
                           discrete = FALSE,
                           discreteSupport = NULL,
                           seed = 123,
                           trace = TRUE
                           ){
    ## Need control over selected data entries
    allDataVariables <- c(observation.name,
                          data.term.indicator,
                          lower.cdf.indicator,
                          upper.cdf.indicator
                          )
    if(missing(observation.name)) stop("'observation.name' must define a data component")
    if(!(observation.name %in% names(obj$env$data))) stop("'observation.name' must be in data component")
    if (is.null(obj$env$random)) stop("oneStepPredict is only for random effect models")
    varLengths <- sapply(obj$env$data[allDataVariables], length)
    if(min(varLengths) != max(varLengths)){
        print(varLengths)
        stop("Variable lengths differ.")
    }
    method <- match.arg(method)
    par <- obj$env$last.par.best
    par <- par[-obj$env$random]
    obs <- as.vector(obj$env$data[[observation.name]])
    ## Check
    if(!discrete){
        ndup <- sum(duplicated(obs))
        if(ndup > 0){
            warning("Observations do not look continuous. Number of duplicates = ", ndup)
        }
    }
    ## Default subset/permutation:
    if(is.null(subset)) subset <- 1:length(obs)
    ## Check
    if(!is.null(conditional)){
        if(length(intersect(subset, conditional)) > 0){
            stop("'subset' and 'conditional' have non-empty intersection")
        }
    }
    unconditional <- setdiff(1:length(obs), union(subset, conditional))

    ## Args to construct copy of 'obj'
    args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
    ## Fix all non-random components of parameter list
    tmp <- as.logical(obj$env$par * 0)
    tmp[-obj$env$random] <- TRUE
    li <- lapply(obj$env$parList(par = tmp), function(x) any(x!=0))
    fix <- names(li)[unlist(li)]
    parameters <- obj$env$parameters
    map <- lapply(parameters[fix], function(x)factor(x*NA))
    args$map <- map ## Overwrite map
    ## Find randomeffects character
    args$random <- names(li[!unlist(li)])
    ## Move data$name to parameter$name
    args$parameters[allDataVariables] <- args$data[allDataVariables]
    args$data[allDataVariables] <- NULL
    ## Pretend these are *not observed*:
    if(length(unconditional)>0){
        if(is.null(data.term.indicator))
            stop("Failed to disable some data terms (because 'data.term.indicator' missing)")
        args$parameters[[data.term.indicator]][unconditional] <- 0
    }
    ## Pretend these are *observed*:
    if(length(conditional)>0){
        if(is.null(data.term.indicator))
            stop("Failed to enable some data terms (because 'data.term.indicator' missing)")
        args$parameters[[data.term.indicator]][conditional] <- 1
    }
    ## Make map for observations and indicator variables:
    fac <- 1:length(obs)
    fac[conditional] <- NA
    fac[unconditional] <- NA
    fac[subset] <- 1:length(subset) ## Permutation
    fac <- factor(fac)
    map <- rep(list(fac), length(allDataVariables))
    names(map) <- allDataVariables
    args$map <- c(args$map, map)
    ## New object be silent
    args$silent <- TRUE
    ## Create new object
    newobj <- do.call("MakeADFun", args)

    ## Helper function to loop through observations:
    nm <- names(newobj$par)
    obs.pointer       <- which(nm == observation.name)
    data.term.pointer <- which(nm == data.term.indicator)
    lower.cdf.pointer <- which(nm == lower.cdf.indicator)
    upper.cdf.pointer <- which(nm == upper.cdf.indicator)
    observation <- local({
        obs.local <- newobj$par
        i <- 1:length(subset)
        function(k, y=NULL, lower.cdf=FALSE, upper.cdf=FALSE){
            ## Disable all observations later than k:
            obs.local[data.term.pointer[k<i]] <- 0
            ## On request, overwrite k'th observation:
            if(!is.null(y)) obs.local[obs.pointer[k]] <- y
            ## On request, get tail probs rather than point probs:
            if(lower.cdf | upper.cdf){
                obs.local[data.term.pointer[k]] <- 0
                if(lower.cdf) obs.local[lower.cdf.pointer[k]] <- 1
                if(upper.cdf) obs.local[upper.cdf.pointer[k]] <- 1
            }
            obs.local
        }
    })

    ## Trace one-step functions
    tracefun <- function(k)if(trace)print(k)
    ## Apply a one-step method and generate common output assuming
    ## the method generates at least:
    ##   * nll
    ##   * nlcdf.lower
    ##   * nlcdf.upper
    applyMethod <- function(oneStepMethod){
        pred <- as.data.frame(t(sapply(1:length(subset), oneStepMethod)))
        pred$Fx <- 1 / ( 1 + exp(pred$nlcdf.lower - pred$nlcdf.upper) )
        pred$px <- 1 / ( exp(-pred$nlcdf.lower + pred$nll) +
                         exp(-pred$nlcdf.upper + pred$nll) )
        if(discrete){
            if(!is.null(seed)) set.seed(seed)
            U <- runif(nrow(pred))
        } else {
            U <- 0
        }
        pred$residual <- qnorm(pred$Fx - U * pred$px)
        pred
    }

    ## ######################### CASE: oneStepGaussian
    if(method == "oneStepGaussian"){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        oneStepGaussian <- function(k){
            tracefun(k)
            index <- subset[k]
            f <- function(y){
                newobj$fn(observation(k, y))
            }
            g <- function(y){
                newobj$gr(observation(k, y))[obs.pointer[k]]
            }
            opt <- nlminb(obs[index], f, g)
            H <- optimHess(opt$par, f, g)
            c(observation=obs[index], mean=opt$par, sd=sqrt(1/H))
        }
        pred <- as.data.frame(t(sapply(1:length(subset), oneStepGaussian)))
        pred$residual <- (pred$observation-pred$mean)/pred$sd
    }

    ## ######################### CASE: oneStepGeneric
    if((method == "oneStepGeneric") & (!discrete)){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        oneStepGeneric <- function(k){
            tracefun(k)
            ans <- try({
                index <- subset[k]
                f <- function(y){
                    newobj$fn(observation(k, y))
                }
                nll <- f(obs[index]) ## Marginal negative log-likelihood
                F1 <- integrate(Vectorize(function(x)exp(-(f(x) - nll))), -Inf, obs[index])$value
                F2 <- integrate(Vectorize(function(x)exp(-(f(x) - nll))), obs[index], Inf)$value
                nlcdf.lower = nll - log(F1)
                nlcdf.upper = nll - log(F2)
                c(nll=nll, nlcdf.lower=nlcdf.lower, nlcdf.upper=nlcdf.upper)
            })
            if(is(ans, "try-error")) ans <- NaN
            ans
        }
        pred <- applyMethod(oneStepGeneric)
    }

    ## ######################### CASE: oneStepDiscrete
    if((method == "oneStepGeneric") & (discrete)){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        obs <- as.integer(round(obs))
        if(!is.null(seed)) set.seed(seed)
        if(is.null(discreteSupport)){
            warning("Setting 'discreteSupport' to ",min(obs),":",max(obs))
            discreteSupport <- min(obs):max(obs)
        }
        oneStepDiscrete <- function(k){
            tracefun(k)
            ans <- try({
                index <- subset[k]
                f <- function(y){
                    newobj$fn(observation(k, y))
                }
                nll <- f(obs[index]) ## Marginal negative log-likelihood
                F <- Vectorize(function(x)exp(-(f(x) - nll))) (discreteSupport)
                F1 <- sum( F[discreteSupport <= obs[index]] )
                F2 <- sum( F[discreteSupport >  obs[index]] )
                nlcdf.lower = nll - log(F1)
                nlcdf.upper = nll - log(F2)
                c(nll=nll, nlcdf.lower=nlcdf.lower, nlcdf.upper=nlcdf.upper)
            })
            if(is(ans, "try-error")) ans <- NaN
            ans
        }
        pred <- applyMethod(oneStepDiscrete)
    }

    ## ######################### CASE: fullGaussian
    if(method == "fullGaussian"){
        ## Same object with y random:
        args2 <- args
        args2$random <- c(args2$random, observation.name)
        ## Change map: Fix everything except observations
        fix <- setdiff(allDataVariables, observation.name)
        args2$map[fix] <- lapply(args2$map[fix],
                                 function(x)factor(NA*unclass(x)))
        newobj2 <- do.call("MakeADFun", args2)
        newobj2$fn() ## Test-eval to find mode
        mode <- newobj2$env$last.par
        GMRFmarginal <- function (Q, i, ...) {
            ind <- 1:nrow(Q)
            i1 <- (ind)[i]
            i0 <- setdiff(ind, i1)
            if (length(i0) == 0)
                return(Q)
            Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
            L0 <- Cholesky(Q0, ...)
            ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
                solve(Q0, Q[i0, i1, drop = FALSE])
            ans
        }
        h <- newobj2$env$spHess(mode, random=TRUE)
        i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
        Sigma <- solve( as.matrix( GMRFmarginal(h, i) ) )
        res <- obs[subset] - mode[i]
        L <- t(chol(Sigma))
        pred <- data.frame(residual = as.vector(solve(L, res)))
    }

    ## ######################### CASE: cdf
    if(method == "cdf"){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        cdf <- function(k){
            print(k)
            nll <- newobj$fn(observation(k))
            nlcdf.lower <- newobj$fn(observation(k, lower.cdf = TRUE))
            nlcdf.upper <- newobj$fn(observation(k, upper.cdf = TRUE))
            c(nll=nll, nlcdf.lower=nlcdf.lower, nlcdf.upper=nlcdf.upper)
        }
        pred <- applyMethod(cdf)
    }

    pred
}
