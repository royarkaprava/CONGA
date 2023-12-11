
CONGAfitNew <- function(X, Total_itr = 5000, lambdashrk=1, burn = 2500){
  library(mvtnorm)
  library(bayesm)
  
  Ti <- nrow(X)
  c  <- ncol(X)
  po <- max(X)
  #Fitting the model starts here
  tes <- rep(0, po)
  
  #Tune power paramter (theta) from the paper.
  for(i in 1:(2*po)){
    tes[i] <- mean((cov(atan(X)^(i/2))-cov(X))^2)
  }
  
  po <- which.min(tes)/10
  
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
    poisd <- dpois(0:apnum, lambda[i, j]) * exp(lambda[i, j])
    betapart <- exp(sum(-Beta[j, -j]*(atan(X[i, -j])^po))*(atan(0:apnum))^po)
    return(sum(poisd * betapart))
  }
  
  llhoodl <- function(i, j, X, lambda, beta, apnum = 100){
    #ret <- rep(0, Ti)
    #for(k in 1:Ti){
    Beta <- matrix(0, c, c)
    Beta[row(Beta)>col(Beta)] <- beta
    Beta <- Beta + t(Beta)
    meant <- rep(0, c)
    for(k in 1:c){
      meant[k] <- atanmean(lambda[i, k])
    }
    
    ret <- -lambda[i, j]+X[i, j]*log(lambda[i, j]) + (lambda[i, j])
    ret <- (ret) - log(numerlike1(i, j, lambda, beta))
    
    return((ret))
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
  lambda    <- X 
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
  
  pdxid <- diag(solve(cov(atan(X)^(po))))#pdx/(Ti)))
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  while(itr < Total_itr){
    itr <- itr + 1
    
    for(k in 1:c){
      
      lambdac       <- lambda
      Q <- matrix(0, Ti, Ti)
      for(i in 1:Ti){
        #for(j in 1:Ti){
        Q[i, ] <- exp(-lambda[,k]+log(lambda[, k])*X[i, k]-lgamma(X[i, k]+1))+1e-100
        #}
        Q[i, i] <- M[k] * dgamma(lambda[i, k],alpha+X[i, k], betalam + 1)
        Q[i, ] <- Q[i, ] / sum(Q[i, ])
        new_ind <- sample(x = 1:Ti, 1, replace = T, prob = Q[i,])
        if(new_ind == i){
          lamc          <- rgamma(1, alpha+X[i, k], betalam + 1)
          lambdac[i, k] <- lamc
          newcount <- newcount + 1
          
          R <- llhoodl(i, k, X, lambdac, beta) - llhoodl(i, k, X, lambda, beta)
          R <- R + (dgamma(lamc, alpha, betalam, log = T) - dgamma(lambda[i, k], alpha, betalam, log = T)) 
          
          R <- R - (dgamma(lamc, alpha+X[i, k], betalam + 1, log = T) + dgamma(lambda[i, k], alpha+X[i, k], betalam + 1, log = T)) 
          
          u <- runif(1)
          
          if(is.na(R) == T || is.nan(R) == T){R = 1}
          if(log(u) < R)
          {lambda[, k] <- lambdac[, k]
          aclam <- aclam + 1}
        }
        else
        {lambda[i, k]= lambda[new_ind, k]
        lambdac[i, k] = lambda[new_ind, k]}
      }
      if(sum(lambdac[, k]!=lambda[, k])>0){}
      ka <- length(unique(lambda[, k]))
      delta2 <- rbeta(1, M[k], Ti)
      M[k] <- rgamma(1, 10 + ka, 10 - log(delta2))
    }
    
    newc[itr] <- newcount
    lambda_p[[itr]] <- lambda
    
    bsigma <- matrix(0, c, c)
    bsigma[row(bsigma)>col(bsigma)] <- Z * rep(s1, betalen) + (1-Z) * rep(s0, betalen)
    bsigma <- bsigma + t(bsigma)
    Betac  <- Beta
    
    for(i in 1:c){
      mean <- - pdx[i, -i] 
      varc <- bsigma[i,-i] 
      varctemp <- matrix(0, c-1, c-1)
      varctemp <- Beta[-i, -i]
      diag(varctemp) <- pdxid[-i]
      varctempei <- eigen(varctemp)
      varctempeiU <- varctempei$vectors
      varctempeiD <- varctempei$values
      varctemp <- varctempeiU %*% diag(1/abs(varctempeiD)) %*% t(varctempeiU)
      varcei <- eigen((var(atan(X[, i])^po+lambdashrk) * Ti)* varctemp + diag(1 / varc))
      varceiU <- varcei$vectors
      varceiD <- varcei$values
      varc <- varceiU %*% diag(1/abs(varceiD)) %*% t(varceiU)
      #varc <- (varc+t(varc))/2
      betac <- array(rmvnorm(1, varc %*% mean, varc))
      if(length(is.na(betac))) betac[which(is.na(betac))] <- 0
      #const <- as.numeric(sdb / sqrt(crossprod(betac - Beta[i, -i])))
      #if(const>1){const = 1}
      #betac <- Beta[i, -i] + const*(betac - Beta[i, -i])
      Betac[i, -i] <- betac 
      Betac[- i, i] <- betac 
      betac   <-  Betac[i, -i]
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
    
    beta_p[[itr]] <- beta
  
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
  }
  close(pb)
  out <- list(BetaMCMC = beta_p[(burn+1):Total_itr], LambdaMCMC = lambda_p[(burn+1):Total_itr])
  return(out)
}

