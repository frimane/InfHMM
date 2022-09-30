
### Code written by ``Azeddine Frimane'' (Azeddine.frimane@angstrom.uu.se)

#============================= MCMC BEAM SAMPLER ===============================
infHmm <- function(meas = NULL, n.iter = NULL, n.burn = NULL, ns.init = NULL, 
                     pri.gam = NULL, pri.alf = NULL){
#...............................................................................
  
  #######--------------------------- Initial/prior values
  N <- length(meas)
  # Hyper-parameters
  mu.0 <- mean(meas); kappa.0 <- 1
  alpha.0 <- .5/var(meas); beta.0 <- 1
  # Save Max-likelihood and some statistics
  likely_max <- -Inf
  likely_sequence <- NA
  sveZm <- NA; sveGam <- NA; sveAlf <- NA; sveLikely <- NA
  # MCMC initialization 
  alf0 <- pri.alf[1]; gam0 <- pri.gam[1] 
  alf1 <- pri.alf[2]; gam1 <- pri.gam[2]
  alf <- 1; gam <- 1 # Initial concentration parameters
  staSeq <- sample(ns.init, N, replace = T) # sample initial states
  zm <- max(staSeq) # number of states
  # sample initial emission parameters
  p. <- rgamma(zm, shape = alpha.0, rate = beta.0)
  p.. <- rnorm(zm, mean = mu.0, sd = sqrt(p.*kappa.0))
  p <- rbind(p.., 1/sqrt(p.))
  # sample the initial base distribution
  baseDist <- dirichletDistcpp(c(tabulate(staSeq), gam))
  # initial count transition-numbers
  tmpTM <- sapply(1:zm, function(t){tabulate(staSeq[which(staSeq==t)-1], zm)})
  if(!is.matrix(tmpTM)) tmpTM <- as.matrix(tmpTM)
  tMat <- t(apply(tmpTM, 1, function(x){
    x. <- c(x,0) + alf*baseDist
    dirichletDistcpp(x.)
  }))
  
  ## Update parameters
  updateParameters <- function(x){
    n <- length(x); xm <- mean(x)
    kappa.n <- kappa.0 + n
    mu.n <- (mu.0*kappa.0 + n*xm)/kappa.n
    alpha.n <- alpha.0 + .5*n
    beta.n <- beta.0 + .5*sum((x - xm)^2) + .5*kappa.0*n*((xm-mu.0)^2)/kappa.n
    pre.j <- rgamma(1, shape = alpha.n, rate  = beta.n)
    mean.j <- rnorm(1, mean = mu.n, sd = 1/sqrt(pre.j*kappa.n))
    return(c(mean.j, 1/sqrt(pre.j)))
  }
  
  #####---------------------------------- MCMC iterations
  for(k in 1:n.iter){
    # -----------------------------------------------
    # -----------------------------------------------
    # Sample auxiliary variables
    auxVar <- 1
    for (i in 2:N) auxVar[i] <- runif(1, 0, tMat[staSeq[i-1], staSeq[i]])
    # -----------------------------------------------
    # # Expand the transition matrix and # of States if necessary
    tMatlist <- eXpand(tMat, auxVar, gam, alf, baseDist)
    zm <- tMatlist[[3]]
    supZm <- tMatlist[[4]]
    tMat <- tMatlist[[1]][ ,-(zm+1)]
    baseDist <- tMatlist[[2]]
    
    if(!is.matrix(tMat)) tMat <- as.matrix(tMat) # case where there is only one state
    tMat <- tMat/rowSums(tMat)
    # -----------------------------------------------
    # Sample emission parameters
    p <- sapply(seq(zm-supZm), function(y) updateParameters(meas[staSeq==y]), 
                simplify = TRUE)
    if(is.null(dim(p))){ # If there is just one state, this will overcome the 
      p <- matrix(p) }   # conversion to numeric class
    # sample -if needed- new emission parameters
    if(supZm != 0){
      p. <- rgamma(supZm, shape = alpha.0, rate = beta.0)
      p.. <- rnorm(supZm, mean = mu.0, sd = 1/sqrt(p.*kappa.0))
      p <- cbind(p, rbind(p..,1/sqrt(p.)))
    } 
    # -----------------------------------------------
    # Sample new trajectories, FF-BS algorithm
    lData <- sapply(meas, function(x) dnorm(x, p[1,], p[2,]), simplify = TRUE)
    if(is.vector(lData)) lData <- matrix(lData, nrow = 1)
    staSeq0 <- as.numeric(FF_BS(lData, tMat, rep(1, zm), auxVar))
    # -----------------------------------------------
    # Remove empty states, and Relabel them
    sveStat <- which(tabulate(staSeq0)!=0)
    fstates <- factor(staSeq0)
    levels(fstates) <- 1:length(sveStat)
    staSeq <- as.numeric(fstates) # 
    # -----------------------------------------------
    # Compute transition counts
    zm <- max(staSeq) 
    tmpTM <- sapply(1:zm, function(t){tabulate(staSeq[which(staSeq==t)-1], zm) })
    if(!is.matrix(tmpTM)) tmpTM <- as.matrix(tmpTM)
    # -----------------------------------------------
    # Sample the auxiliary counts
    tmpTMlist <- split(tmpTM, col(tmpTM, as.factor = TRUE)) 
    cts <- mapply(function(x, y){
      sum(countsCpp(x, y, alf))
    }, tmpTMlist, baseDist[sveStat]) 
    # -----------------------------------------------
    # Update Gamma-concentration parameter:
    tc <- sum(cts)
    pAux1 <- tc/(tc+gam)
    aux1 <- rbinom(1,1,pAux1)
    aux2 <- rbeta(1, gam+1, tc)
    gam <- rgamma(1, gam0+zm-aux1, gam1-log(aux2))
    # Update Alpha-concentration parameter:
    auxVars <- apply(tmpTM, 1, function(x){
      v1 <- sum(x)/(sum(x)+alf)
      c(rbinom(1,1,v1), log(rbeta(1, alf+1, sum(x))) )
    }) 
    alf <- rgamma(1, alf0+tc-sum(auxVars[1,]), alf1-sum(auxVars[2,]))
    # And lastly Update the shared base distribution 
    baseDist <- dirichletDistcpp(c(cts, gam))
    # -----------------------------------------------
    # Sample the transition matrix
    tMat <- t(apply(tmpTM, 1, function(x){
      x. <- c(x,0) + alf*baseDist
      dirichletDistcpp(x.)
    }))
    # -----------------------------------------------
    # -----------------------------------------------
    # Compute state-likelihoods after burn-in iterations
    if(k > n.burn){
    p <- as.matrix(p[,sveStat])
      likely <- log(dnorm(meas[1], p[1,staSeq[1]], p[2,staSeq[1]])) +
        sum(sapply(2:N, function(x)
          log(dnorm(meas[x], p[1, staSeq[x]], p[2, staSeq[x]])) +
            log(tMat[staSeq[x-1], staSeq[x]])))
      sveLikely <- c(sveLikely, likely)
      if (likely > likely_max) {
        likely_max <- likely
        likely_sequence <- staSeq
        likely_transition_matrix <- tMat
        likely_parameters <- p
        likely_iter <- k
        likely_alf <- alf
        likely_gam <- gam
        likely_baseDist <- baseDist
      }
    }
    # some statistics
    sveZm <- c(sveZm,zm)
    sveGam <- c(sveGam, gam)
    sveAlf <- c(sveAlf, alf)
    
    # print outputs
    cat("\nIteration",k,"","Alpha",alf,"","Gamma",gam,"",
        "Number of states is",zm,"\n")
  }
  return(
    list(states.sequence = likely_sequence,
         likely.num.elements.by.state = table(likely_sequence),
         likely.transition.matrix = likely_transition_matrix,
         likely.emission.parameters = likely_parameters,
         likely.concentration.alpha = likely_alf,
         likely.concentration.gamma = likely_gam,
         likely.iteration = likely_iter,
         likely.base.Dist = likely_baseDist,
         saved.num.states = sveZm[-1],
         saved.Gam = sveGam[-1],
         saved.Alf = sveAlf[-1],
         saved.Likelihood = sveLikely[-1]))
}
#===============================================================================
