rocvec <- function(fitN, pdmat0){
  cut <- c(1:18)/20
  k <- 0
  FP = 0
  TP = 0
  for(Cutoff in cut){
    beta_p <- fitN$BetaMCMC[1:1000]
    #Post process to construct the graph. This is a different appraoch from the paper.
    #This works better to generate the ROC curve. This will be updated in the next revision of the paper
    
    
    #The following line builds 'number of parameters' X 'number of MCMC samples' matrix. 
    #Hence each row of this matrix corresponds to the MCMC samples of each parameter in beta. 
    #If there are c many variables as above, the length of beta is c\choose 2.
    
    retu <- matrix(unlist(beta_p), length(beta_p[[1]]), length(beta_p))
    
    #From above matrix, we can compute the posterior probability of being on the positive part of the real line. 
    #This is computed in the following few lines.
    
    
    posmat <- function(x){
      mean(x >= 0)
    }
    
    Qvec   <- apply(retu, 1, posmat)
    
    #Now there are different ways to use these probabilities.
    #We set a cutoff say 0.7 and compute the following,
    #Increase in the cutoff would lead to sparser graph and decrease in it
    #would lower the sparsity.
    
    #Cutoff <- 0.5 #To construct ROC, we need to vary this cutoff
    Qvec   <- (abs(Qvec - 0.5)/0.5>Cutoff)
    
    
    #All possible edges are the columns in the following index matrix
    index <- as.matrix(combinat::combn(1:c, 2))
    
    #Estimated adjacency under a cutoff 0.7
    estiAdj <- matrix(0, c, c)
    estiAdj[t(index)] <- Qvec
    
    estiAdj <- estiAdj + t(estiAdj)
    
    #Now under the cutoff 0.7 index[,which(Qvec==1)] gives you the detected edges. 
    
    #If we have a true precision matrix encoding the graphical dependence. 
    #Then,
    
    pdmatind <- (pdmat0 != 0) #Binarize the entries either zero or non-zero.
    
    pdmatt <- pdmatind[t(index)] #extracts the entries correspond to different edges
    ind1 <- which(pdmatt==1) #Find the indices where there is an edge based on the true precision matrix
    ind0 <- which(pdmatt==0) #Find the indices where there is no edge based on the true precision matrix
    
    k <- k+1
    FP[k] <- mean(Qvec[ind0]==1) #False positive
    TP[k] <- mean(Qvec[ind1]==1) #True positive
  }
  
  return(list(FP=FP, TP=TP))
}

rc <- rocvec(fitN, pdmat0)

plot(rc$FP, rc$TP, type="l", xlim=c(0,1), ylim=c(0,1), xlab="False Positive", ylab="True Positive")
