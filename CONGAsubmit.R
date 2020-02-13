library(mvtnorm)
library(statmod)
library(GIGrvg)
library(MCMCpack)
library(tmvtnorm)
library(XMRF)
#install.packages("lcmix", repos="http://R-Forge.R-project.org")

Ti <- 100

c <- 30

Posdef <- function (n, ev = runif(n, 0, 3))
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

pdmat0 <- Posdef(n=c, ev=runif(c, 0, 10))

pdmat0[which(abs(pdmat0) < 1)] <- 0 #Generate sparse precision matrix

data_p <- list()
beta_pp <- list()
#for(rep in 1:100){
#Generate count data from copula type model
  rawvars <- rmvnorm(Ti, mean=rep(0, c), sigma=solve(pdmat))
  pvars <- pnorm(rawvars)
  X <- qpois(pvars, 5)
  
#Fitting the model starts here
  tes <- rep(0, po)
  
  #Tune power paramter (theta) from the paper.
  for(i in 1:(10*po)){
    tes[i] <- mean((cov(atan(X)^(i/10))-cov(X))^2)
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
    
    ret <- dpois(X[i, j], lambda[i, j], log = T) + (lambda[i, j])
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
  
  lambda0   <- lambda
  
  #MCMC starts here
  
  Total_itr <- 5000
  lambda_p  <- list()
  beta_p    <- list()
  itr       <- 0
  aclam     <- 0
  newcount  <- 0
  M         <- rep(2, c)
  alpha     <- 1
  betalam   <- 1
  lambda    <- X 
  matan     <- colMeans(atan(X)^po) 
  Me        <- matrix(matan, Ti, c, byrow = T)
  pdx       <- t(atan(X)^po-Me)%*%(atan(X)^po-Me) 
  me        <- pdx[t(index)]
  
  betalen   <- c*(c-1)/2
  beta      <- rnorm(betalen)
  sdl       <- 0.5
  sdb       <- 1
  acbeta    <- 0
  ap        <- 1/betalen
  psi       <- rexp(betalen, 0.5)
  phi       <- rdirichlet(1, rep(ap, betalen))
  tau       <- rgamma(1, betalen*ap, 0.5)
  sigma     <- rep(1, betalen)
  newc      <- rep(0, Total_itr)
  newcount  <- 0
  Beta <- matrix(0, c, c)
  Beta[row(Beta)>col(Beta)] <- beta
  Beta <- Beta + t(Beta)

  s1    <- 1  
  s0    <- 1
  Z     <- rep(1, betalen)
  
  prq   <- 1/betalen
  consb <- 1
  
  pdxid <- diag(solve(pdx/(Ti)))
  while(itr < Total_itr){
    itr <- itr + 1
    
    for(k in 1:c){
      
      lambdac       <- lambda
      Q <- matrix(0, Ti, Ti)
      for(i in 1:Ti){
        #for(j in 1:Ti){
        Q[i, ] <- dpois(X[i, k], lambda[, k])
        #}
        Q[i, i] <- M[k] * dgamma(lambda[i, k],alpha+X[i, k], betalam + 1)
        Q[i, ] <- Q[i, ] / sum(Q[i, ])
        new_ind <- sample(x = 1:Ti, 1, replace = T, prob = Q[i,])
        if(new_ind == i){
          lamc          <- rgamma(1, alpha+X[i, k], betalam + 1)
          lambdac[i, k] <- lamc
          newcount <- newcount + 1
          
          lambdac[, k]         <- lambda[, k] + sdl*(lambdac[, k] - lambda[, k])
          
          lambdac[which(lambdac[, k] < 0), k] <- 0
          
          R <- llhoodl(i, k, X, lambdac, beta) - llhoodl(i, k, X, lambda, beta)
          R <- R + (dgamma(lamc, alpha, betalam, log = T) - dgamma(lambda[i, k], alpha, betalam, log = T)) 
          
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
    if((itr%%100)==0){
      ar <- aclam / (newcount)
      if(newcount==0){ar = 0.35}
      
      if(ar<.25){sdl <- sdl/2}
      if(ar>.45){sdl <- sdl*2}
      if(sdl>1){sdl = 1}
    }
    
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
      varctemp <- solve(varctemp)
      varc <- solve((var(atan(X[, i])^po+4) * Ti)* varctemp + diag(1 / varc))
      varc <- (varc+t(varc))/2
      betac <- array(rmvnorm(1, varc %*% mean, varc))
      betac[which(is.na(betac))] <- 0
      const <- as.numeric(sdb / sqrt(crossprod(betac - Beta[i, -i])))
      if(const>1){const = 1}
      betac <- Beta[i, -i] + const*(betac - Beta[i, -i])
      Betac[i, -i] <- betac 
      Betac[- i, i] <- betac 
      betac   <-  Betac[i, -i]
      R <- llhoodb(i, X, lambda, Betac[t(index)]) - llhoodb(i, X, lambda, beta)
      R <- R + sum((dnorm(Betac[t(index)], 0, sigma, log = T) - dnorm(beta, 0, sigma, log = T)))
      u <- runif(1)
      
      if(is.na(R) == T || is.nan(R) == T){R = 1}
      if(log(u) < R)
      {beta <- Betac[t(index)]
      Beta <- Betac
      acbeta <- acbeta + 1}
    }
    
    beta_p[[itr]] <- beta
    if((itr%%100)==0){
      ar <- acbeta / (c*itr)
      
      if(ar<.25){sdb <- sdb/2}
      if(ar>.45){sdb <- sdb*2}
      
    }
    print(itr)
  }
  #print(rep)
  
  #Generated MCMC samples of the coefficient vector (beta)
  Betamat <- matrix(unlist(beta_p), betalen, 5000)
#  beta_pp[[rep]] <- retu
  
#}

#The range of postburn samples
M1 <- 2501
M2 <- 5000

#Post process to construct the graph. This is a different appraoch from the paper.
#This works better to generate the ROC curve. This will be updated in the next revision of the paper

retu <- matrix(unlist(beta_p[M1:M2]), betalen, M2-M1+1)

posmat <- function(x){
  mean(x >= 0)
}

Cutoff <- 0.7 #To construct ROC, we need to vary this cutoff
Qvec   <- apply(retu, 1, posmat)
Qvec   <- (abs(Qvec - 0.5)/0.5>Cutoff)

pdmatind <- (pdmat0 != 0)
pdmatt <- pdmatind[t(index)]
ind1 <- which(pdmatt==1) #Find the indices where there is an edge
ind0 <- which(pdmatt==0) #Find the indices where there is no edge

mean(Qvec[ind0]==1) #False positive
mean(Qvec[ind1]==1) #True positive