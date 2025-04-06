Rcpp::sourceCpp("/GitHub/CONGA/Modified sampling steps/pocal.cpp")

#lambdashrk is the shrinkage parameter, similar to the one in Wang (2012).

CONGAfitNewer <- function(X, Total_itr = 5000, lambdashrk=1, burn = 2500){
  library(mvtnorm)
  library(bayesm)
  library(armspp)
  
  Ti <- nrow(X)
  c  <- ncol(X)
  po <- max(X)
  # #Fitting the model starts here
  
  # tes <- rep(0, po)
  # 
  # #Tune power paramter (theta) from the paper.
  # for(i in 1:(2*po)){
  #   tes[i] <- mean((cov(atan(X)^(i/2))-cov(X))^2)
  # }
  
  po <- posel(po, X)#which.min(tes)/10
  
  library(combinat)
  
  apnum <- 100
  
  numerlike = 0
  
  index <- as.matrix(combinat::combn(1:c, 2))
  
  #Define useful functions
  ll <- function(ind, lambda, beta){
    return(sum(dpois(ind, lambda, log = T) + lambda)+beta*atan(ind[index[1, ]])*atan(ind[index[2, ]]))
  }
  
  numerlike <- function(lambda, beta){
    numll <- sum(exp(apply(mat, 1, FUN = ll, lambda, beta)))
    return(numll)
  }
  
  numerlike1 <- function(i, j, lambda, beta, apnum = 100){
    Beta <- matrix(0, c, c)
    Beta[row(Beta)>col(Beta)] <- beta
    Beta <- Beta + t(Beta)
    poisd <- dpois(0:apnum, lambda[j]) * exp(lambda[j])
    betapart <- exp(sum(-Beta[j, -j]*(atan(X[i, -j])^po))*(atan(0:apnum))^po)
    return(sum(poisd * betapart))
  }
  
  llhoodl <- function(j, X, lambda, beta, apnum = 100){
    #ret <- rep(0, Ti)
    #for(k in 1:Ti){
    Beta <- matrix(0, c, c)
    Beta[row(Beta)>col(Beta)] <- beta
    Beta <- Beta + t(Beta)
    
    ret <- -lambda[j]+X[,j]*log(lambda[j]) + (lambda[j])
    ret <- (ret) - unlist(sapply(1:Ti, FUN=function(i){return(log(numerlike1(i, j, lambda, beta)))})) 
    
    return(sum(ret))
  } 
  
  llhoodb <- function(j, X, lambda, beta, apnum = 100){
    Beta <- matrix(0, c, c)
    Beta[row(Beta)>col(Beta)] <- beta
    Beta <- Beta + t(Beta)
    sum <- 0
    for(i in 1:Ti){
      sum = sum + sum(-Beta[j, -j]*(atan(X[i, -j])^po)*(atan(X[i, j])^po))- log(numerlike1(i, j, lambda, beta))
    }
    return(sum)
  }
  
  atanmean <- function(theta, apnum = 100){
    return(sum(atan(0:apnum)^po * dpois(0:100, theta)))
  }
  
  #lambda0   <- lambda
  
  #MCMC starts here
  
  lambda_p  <- list()
  beta_p    <- list()
  itr       <- 0
  newcount  <- 0
  M         <- rep(2, c)
  alpha     <- 1
  betalam   <- 1
  lambda    <- colMeans(X) + 1e-10
  matan     <- colMeans(atan(X)^po) 
  Me        <- matrix(matan, Ti, c, byrow = T)
  pdx       <- crossprod(atan(X)^po-Me) 
  me        <- pdx[t(index)]
  
  betalen   <- c*(c-1)/2
  beta      <- rnorm(betalen)
  ap        <- 1/betalen
  psi       <- rexp(betalen, 0.5)
  phi       <- rdirichlet(rep(ap, betalen))
  tau       <- rgamma(1, betalen*ap, 0.5)
  sigma     <- rep(1, betalen)
  newc      <- rep(0, Total_itr)
  newcount  <- 0
  Beta <- matrix(0, c, c)
  Beta[row(Beta)>col(Beta)] <- beta
  Beta <- Beta + t(Beta)
  
  s1    <- 100  
  s0    <- 100
  Z     <- rep(1, betalen)
  
  prq   <- 1/betalen
  consb <- 1
  
  aclam <- 0
  acbeta <- 0
  
  tryout <- try(pdxid <- diag(solve(cov(atan(X)^(po)))), silent = T)#pdx/(Ti)))
  const <- 1e-20
  while(class(tryout)=="try-error"){
    tryout <- try(pdxid <- diag(solve(cov(atan(X)^(po)) + const*diag(c))), silent = T)
    const <- 10*const
  }
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  gamout <- dgamma(lambda,alpha+X, betalam + 1)
  sdlam <- 1e-3
  
  uplam <- function(k){
    lambdakc <- rgamma(1, alpha+sum(X[, k]), betalam + Ti)#lambda[k] + rnorm(1, sd=sdlam)
    
    lambdac <- lambda
    lambdac[k] <- lambdakc
    
    R <- llhoodl(k, X, lambdac, beta) - llhoodl(k, X, lambda, beta)
    R <- R + (dgamma(lambdakc, alpha, betalam, log = T) - dgamma(lambda[k], alpha, betalam, log = T)) 
    
    #R <- R - (dgamma(lambdakc, alpha+sum(X[, k]), betalam + Ti, log = T) + dgamma(lambda[k], alpha+ sum(X[, k]), betalam + Ti, log = T)) 
    
    Q <- dgamma(lambda[k], alpha+sum(X[, k]), betalam + Ti, log = T) - dgamma(lambdakc, alpha+sum(X[, k]), betalam + Ti, log = T)
    R <- R + Q
    u <- runif(1)
    
    ret <- lambda[k]
    if(is.na(R) == T || is.nan(R) == T){R = 1}
    if(log(u) < R)
    {ret <- lambdac[k]}
   
    return(ret) 
  }
  
  while(itr < Total_itr){
    itr <- itr + 1
    
    out <- lapply(1:c, uplam)
    lambda <- unlist(out)
    
    lambda_p[[itr]] <- lambda
    
    bsigma <- matrix(0, c, c)
    bsigma[row(bsigma)>col(bsigma)] <- Z * rep(s1, betalen) + (1-Z) * rep(s0, betalen)
    bsigma <- bsigma + t(bsigma)
    Betac  <- Beta
    
    #t1 <- proc.time()
    for(i in 1:c){
      mean <- - pdx[i, -i] 
      varc <- bsigma[i,-i] 
      varctemp <- matrix(0, c-1, c-1)
      #Beta[-i, -i] is the Omega_11
      varctemp <- Beta[-i, -i]
      
      #Unlike Wang the diagonal will be precomputed as above
      diag(varctemp) <- pdxid[-i]
      
      #Calculating inverse of Omega_11 using eigen first get eigen
      
      varctempei <- eigen(varctemp)
      
      #Get the eigen values of Cinv
      
      varcei <- (var(atan(X[, i])^po)*Ti+lambdashrk) /abs(varctempei$values) + 1 / varc
      
      #Preparing to get the inverse using UD^{-1}t(U) = crossprod(t(U)/sqrt(diag(D)))
      varceiU <- t(varctempei$vectors)/sqrt(abs(varcei))
      
      ##Calculate the C matrix
      varc <- crossprod(varceiU)#varceiU %*% diag(1/abs(varceiD)) %*% t(varceiU)
      
      #Generate the Candidate beta
      betac <- array(rmvnorm(1, varc %*% mean, varc))
      
      ###Next do MH step whether to accept or reject unlike Wang
      if(length(is.na(betac))) betac[which(is.na(betac))] <- 0
      #const <- as.numeric(sdb / sqrt(crossprod(betac - Beta[i, -i])))
      #if(const>1){const = 1}
      #betac <- Beta[i, -i] + const*(betac - Beta[i, -i])
      Betac[i, -i] <- betac 
      Betac[- i, i] <- betac 
      #betac   <-  Betac[i, -i]
      R <- llhoodb(i, X, lambda, Betac[t(index)]) - llhoodb(i, X, lambda, beta)
      R <- R + sum((dnorm(Betac[i, -i], 0, bsigma[i, -i], log = T) - dnorm(Beta[i, -i], 0, bsigma[i, -i], log = T)))
      
      Q <- dmvnorm(Beta[- i, i], varc %*% mean, varc, log = T) - dmvnorm(betac, varc %*% mean, varc, log = T)
      R <- R + Q
      u <- runif(1)
      
      if(is.na(R) == T || is.nan(R) == T){R = 1}
      if(log(u) < R)
      {beta <- Betac[t(index)]
      Beta <- Betac
      acbeta <- acbeta + 1}
    }
    #t2 <- proc.time()
    #t2-t1
    beta_p[[itr]] <- beta
    
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
    
  }
  close(pb)
  out <- list(BetaMCMC = beta_p[(burn+1):Total_itr], LambdaMCMC = lambda_p[(burn+1):Total_itr])
  return(out)
}

