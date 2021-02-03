library(combinat)

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

#for(rep in 1:100){
#Generate count data from copula type model
  rawvars <- rmvnorm(Ti, mean=rep(0, c), sigma=solve(pdmat0))
  pvars <- pnorm(rawvars)
  X <- qpois(pvars, 5)

  fit <- CONGAfit(X)
  beta_p <- fit$BetaMCMC
#Post process to construct the graph. This is a different appraoch from the paper.
#This works better to generate the ROC curve. This will be updated in the next revision of the paper

retu <- matrix(unlist(beta_p), lenght(beta_p[[1]]), length(beta_p))

posmat <- function(x){
  mean(x >= 0)
}

Cutoff <- 0.7 #To construct ROC, we need to vary this cutoff
Qvec   <- apply(retu, 1, posmat)
Qvec   <- (abs(Qvec - 0.5)/0.5>Cutoff)

pdmatind <- (pdmat0 != 0)
index <- as.matrix(combinat::combn(1:c, 2))
pdmatt <- pdmatind[t(index)]
ind1 <- which(pdmatt==1) #Find the indices where there is an edge
ind0 <- which(pdmatt==0) #Find the indices where there is no edge

mean(Qvec[ind0]==1) #False positive
mean(Qvec[ind1]==1) #True positive