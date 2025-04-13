Rcpp::sourceCpp("/GitHub/CONGA/Modified sampling steps/pocal.cpp")

#lambdashrk is the shrinkage parameter, similar to the one in Wang (2012).

CONGAfitNewerHS <- function(X, Total_itr = 5000, lambdashrk=1, burn = 2500){
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
  
  numerlike1 <- function(i, j, lambda, Beta, apnum = 100){
    poisd <- dpois(0:apnum, lambda[j]) * exp(lambda[j])
    betapart <- exp(sum(-Beta[j, -j]*(atan(X[i, -j])^po))*(atan(0:apnum))^po)
    return(sum(poisd * betapart))
  }
  
  llhoodl <- function(j, X, lambda, Beta, apnum = 100){
    #ret <- rep(0, Ti)
    #for(k in 1:Ti){
   
    ret <- -lambda[j]+X[,j]*log(lambda[j]) + (lambda[j])
    ret <- (ret) - unlist(sapply(1:Ti, FUN=function(i){return(log(numerlike1(i, j, lambda, Beta)))})) 
    
    return(sum(ret))
  } 
  
  llhoodb <- function(j, X, lambda, Beta, apnum = 100){
    #sum <- 0
    ret <- sum(-Beta[j, -j]*(atan(X[, -j])^po)*(atan(X[, j])^po))- unlist(sapply(1:Ti, FUN=function(i){return(log(numerlike1(i, j, lambda, Beta)))})) 
    
    return(sum(ret))
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
  select_all_but_diag <- function(x) matrix(x[lower.tri(x, diag = F) | upper.tri(x, diag = F)], nrow = nrow(x) - 1, ncol = ncol(x))
  
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  gamout <- dgamma(lambda,alpha+X, betalam + 1)
  sdlam <- 1e-3
  
  lamhorse <- matrix(1, c, c-1)
  tauhorse <- 10#rep(10, m)
  xihorse <- 1#tauhorse/10
  nuhorse  <- matrix(1, c, c-1)
  uplam <- function(k){
    lambdakc <- rgamma(1, alpha+sum(X[, k]), betalam + Ti)#lambda[k] + rnorm(1, sd=sdlam)
    
    lambdac <- lambda
    lambdac[k] <- lambdakc
    
    R <- llhoodl(k, X, lambdac, Beta) - llhoodl(k, X, lambda, Beta)
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
    
    lambdac <- rgamma(c, alpha+colSums(X), betalam + Ti)
    R <- lapply(1:c, FUN=function(k){llhoodl(k, X, lambdac, Beta) - llhoodl(k, X, lambda, Beta)})
    
    R <- unlist(R) + dgamma(lambdac, alpha, betalam, log = T) - dgamma(lambda, alpha, betalam, log = T) 
    
    #R <- R - (dgamma(lambdakc, alpha+sum(X[, k]), betalam + Ti, log = T) + dgamma(lambda[k], alpha+ sum(X[, k]), betalam + Ti, log = T)) 
    
    Q <- dgamma(lambda, alpha+colSums(X), betalam + Ti, log = T) - dgamma(lambdac, alpha+colSums(X), betalam + Ti, log = T)
    R <- R + Q
    u <- runif(c)
    
    lambda <- (log(u) > R)*lambda + (log(u) < R)*lambdac
    
    
    lambda_p[[itr]] <- lambda
    
    bsigma <- matrix(0, c, c)
    bsigma[row(bsigma)>col(bsigma)] <- Z * rep(s1, betalen) + (1-Z) * rep(s0, betalen)
    bsigma <- bsigma + t(bsigma)
    Betac  <- Beta
    
    B1 <- t(select_all_but_diag(t(Beta)))
    tauhorse   <- 10#sqrt(1/rgamma(1, c*(c-1)/2+1/2, 1/xihorse + sum(B1^2/2/lamhorse^2)))
    #xihorse   <- 1/rgamma(1, 1, 1+1/tauhorse^2)
    
    #t1 <- proc.time()
    for(i in 1:c){
      #Update the Horse shoe parameters
      lamhorse[i, ] <- sqrt(1/rgamma(c-1, 1, 1/nuhorse[i, ] + Beta[i, -i]^2/2/tauhorse^2))
      nuhorse[i, ] <- 1/rgamma(c-1, 1, 1+1/lamhorse[i,]^2)
      
      priorbeta <- 1/(lamhorse[i, ]^2 * tauhorse^2)
      
      mean <- - pdx[i, -i] 
      varc <- bsigma[i,-i] 
      varctemp <- matrix(0, c-1, c-1)
      #Beta[-i, -i] is the Omega_11
      varctemp <- Beta[-i, -i]
      
      #Unlike Wang the diagonal will be precomputed as above
      diag(varctemp) <- pdxid[-i]
      
      #Calculating inverse of Omega_11 using eigen; first get eigen
      
      varctempei <- eigen(varctemp)
      if(prod(varctempei$values > -0.01)){
        #Get the inverse as UD^{-1}t(U) = crossprod(t(U)/sqrt(diag(D)))
        varctempeiU <- t(varctempei$vectors)/sqrt(abs(varctempei$values))
        #varctempeiD <- varctempei$values
        varctemp <- crossprod(varctempeiU) #varctempeiU %*% diag(1/abs(varctempeiD)) %*% t(varctempeiU)
        
        ##Now calculate Cinv and its eigen following 1(b)
        varinv <- (var(atan(X[, i])^po)* Ti)* varctemp + diag(priorbeta)
        varcei <- eigen(varinv)
        
        #Get the inverse of Cinv using again UD^{-1}t(U) = crossprod(t(U)/sqrt(diag(D)))
        varceiU <- t(varcei$vectors)/sqrt(abs(varcei$values))
        #varceiD <- varcei$values
        varc <- crossprod(varceiU)
        
        ###Generate the Candidate beta
        betac <- array(mvtnorm::rmvnorm(1, varc %*% mean, varc))
        #rmvnorm.canonical(1, mean, varcei)#
        
        ###Next do MH step whether to accept or reject unlike Wang
        if(length(is.na(betac))) betac[which(is.na(betac))] <- 0
        
        Betac[i, -i] <- betac 
        Betac[- i, i] <- betac 
        betac   <-  Betac[i, -i]
        R <- llhoodb(i, X, lambda, Betac) - llhoodb(i, X, lambda, Beta)
        R <- R + sum((dnorm(Betac[i, -i], 0, 1/priorbeta, log = T) - dnorm(Beta[i, -i], 0, 1/priorbeta, log = T)))
        
        Q <- Beta[- i, i]*(-varinv %*% Beta[- i, i]/2+mean)-betac*(-varinv %*% betac/2 +mean) #dmvnorm(Beta[- i, i], varc %*% mean, varc, log = T) - dmvnorm(betac, varc %*% mean, varc, log = T)
        R <- R + sum(Q)
        u <- runif(1)
        
        if(is.na(R) == T || is.nan(R) == T){R = 1}
        if(log(u) < R)
        {beta <- Betac[t(index)]
        Beta <- Betac
        acbeta <- acbeta + 1}
      }
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

