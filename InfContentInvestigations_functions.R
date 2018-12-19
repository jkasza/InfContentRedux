###########################################
# Information content investigations 
# Helper functions
# J Kasza, 2018-08-03
###########################################

#This file contains functions to calculate
#information content of cluster-period cells 
#when there is treatment effect heterogeneity
#and when some cells of a multiple-period CRT
#are missing by design (e.g. transition periods in a SW)

library(Matrix)


#######################
# Functions for generating design matrices
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}

###SWdesmat with implementation periods
SWdesmatIMP <- function(T, nIMP) {
  Xsw <- SWdesmat(T+nIMP)
  for(i in 1:T) Xsw[i, (i+1):(i+nIMP)] <- NA
 
  return(Xsw[1:(T-1),])
}



SWdesmat2sandwich <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1) + 2)
  for(i in 1:(T-1)) {
    Xsw[i+1,(i+1):T] <- 1
  }
  Xsw[1, ] <- 1
  return(Xsw)
}



pllelbasedesmat2 <- function(T, K) {
  if(K%%2 == 0) {
    Xpllelbase <- matrix(data=0, ncol = T, nrow = K)
    Xpllelbase[1:K/2,2:T] <- 1
    return(Xpllelbase)
  }
  if(K%%2 == 1) {
    Xpllelbase <- matrix(data=0, ncol = T, nrow = K)
    Xpllelbase[1:(K+1)/2,2:T] <- 1
    return(Xpllelbase)
  }
}

crxodesmat2 <- function(T, K) {
  if(K%%2 == 0) {
    Xcrxo <- matrix(data=0, ncol = T, nrow = K)
    Xcrxo[1:K/2, seq(1,T,2)] <- 1
    Xcrxo[(K/2 + 1):K, seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  if(K%%2 == 1) {
    Xcrxo <- matrix(data=0, ncol = T, nrow = K)
    Xcrxo[1:(K+1)/2, seq(1,T,2)] <- 1
    Xcrxo[((K+1)/2 + 1):(K), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  
}
###############################

###############################
# Functions to incorporate treatment effect heterogeneity
# Complete designs, or designs missing a cell or a cluster:
CRTVar_TEHet <- function(Xmat, m, sigu2, sigv2, siguv =0, sig2E) {
  #Xmat: the design matrix
  #m: number if subjects in each cluster-period
  #sigu2: variance of cluster random intercepts
  #sigv2: variance of treatment effect slopes
  #siguv: covariance of random intercept and slope
  #sig2E: variance of subject-level errors
  
  sig2 <- sig2E/m
  
  T <- ncol(Xmat)
  K <- nrow(Xmat)
  
  
  Xvec <- as.vector(t(Xmat))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  
  #If type = 0 (HH or HG)
  #Variance matrix
  Vmat <- matrix(data=c(sigu2, siguv, siguv, sigv2), nrow=2, byrow=TRUE)
  
  Dtemp <- cbind(matrix(data=1, nrow=K, ncol=T), Xmat)
  #Generate the D matrix for all clusters
  # D <- kronecker(diag(1,2), Dvec)
  
  D<- blkDiag(z=array(as.vector(t(Dtemp)), c(T,2,K)))
  #blkDiag assumes square matrices so puffs up the matrices with zero columns.
  #need to remove these zero columns
  #To avoid removing those columns corresponding to clusters that are never
  #exposed, add a header row, which is later deleted.
  D<-rbind(rep(c(1,1,rep(0, T-2)), K),D)
  D<- D[,c(colSums(D, na.rm = TRUE)!=0)]
  D<-D[-1,]
  
  #Remove those rows of D which correspond to clusters not observed in 
  #particular periods
  D<- D[!is.na(Xvec),]
  
  varYbar <- D%*%kronecker(diag(1,K), Vmat)%*%t(D) + sig2*diag(1, nrow(D))
  
  #Variance of the treatment effect estimator is then given by:
  return(solve(t(Zmat)%*%solve(varYbar)%*%Zmat)[ncol(Zmat),ncol(Zmat)] )
  
  #Think about the problems when there are multiple removed observations later:
  ##there will be problems if Zmat is not of full column rank
  ##if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  ##else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
  
  #return(vartheta)
  
}

#Designs missing an entire period:
CRTVar_TEHet_ocol <- function(Xmat, ocol,  m, sigu2, sigv2, siguv =0, sig2E) {
  #Xmat: the design matrix
  #m: number if subjects in each cluster-period
  #sigu2: variance of cluster random intercepts
  #sigv2: variance of treatment effect slopes
  #siguv: covariance of random intercept and slope
  #sig2E: variance of subject-level errors
  
  sig2 <- sig2E/m
  
  T <- ncol(Xmat[,-ocol])
  K <- nrow(Xmat)
  Xvec <- as.vector(t(Xmat[,-ocol]))
  
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  
  #Variance matrix
  Vmat <- matrix(data=c(sigu2, siguv, siguv, sigv2), nrow=2, byrow=TRUE)
  
  Dtemp <- cbind(matrix(data=1, nrow=K, ncol=T), Xmat[,-ocol])
  #Generate the D matrix for all clusters
  # D <- kronecker(diag(1,2), Dvec)
  
  D<- blkDiag(z=array(as.vector(t(Dtemp)), c(T,2,K)))
  #blkDiag assumes square matrices so puffs up the matrices with zero columns.
  #need to remove these zero columns
  #To avoid removing those columns corresponding to clusters that are never
  #exposed, add a header row, which is later deleted.
  D<-rbind(rep(c(1,1,rep(0, T-2)), K),D)
  D<- D[,c(colSums(D, na.rm = TRUE)!=0)]
  D<-D[-1,]
  
  #Remove those rows of D which correspond to clusters not observed in 
  #particular periods
  D<- D[!is.na(Xvec),]
  
  varYbar <- D%*%kronecker(diag(1,K), Vmat)%*%t(D) + sig2*diag(1, nrow(D))
  
  #Variance of the treatment effect estimator is then given by:
  return(solve(t(Zmat)%*%solve(varYbar)%*%Zmat)[ncol(Zmat),ncol(Zmat)] )
  
  #Think about the problems when there are multiple removed observations later:
  ##there will be problems if Zmat is not of full column rank
  ##if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  ##else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
  
  #return(vartheta)
  
}
###############################

#Incorporate discrete-time decay
CRTVar_TEHet_decay <- function(Xmat, m, sigu2, sigv2, siguv =0, sig2E, type = 0, r) {
  #Xmat: the design matrix
  #m: number if subjects in each cluster-period
  #sigu2: variance of cluster random intercepts
  #sigv2: variance of treatment effect slopes
  #siguv: covariance of random intercept and slope
  #sig2E: variance of subject-level errors
  #type: the type of correlation structure. 
  #   if type = 0, Hussey and Hughes assumed if r=1, Hooper/Girling if r< 1
  #   if type = 1, discrete-time decay model assumed.
  
  sig2 <- sig2E/m
  
  T <- ncol(Xmat)
  K <- nrow(Xmat)
  
  
  Xvec <- as.vector(t(Xmat))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  

  #Constant decay var if type==0
  if(type==0) { 
    Vi <-diag((1-r)*sigu2, T) + matrix(data=sigu2*r, nrow=T, ncol=T)
  }
  #exponential decay structure
  if(type==1) { 
    Vi <- sigu2*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  }
  
  #Need to add on the treatment effect components:
  Vi <- rbind(Vi, rep(siguv, T))
  Vi <- cbind(Vi, c(rep(siguv, T), sigv2))
  
  #cbind(matrix(data=rep(diag(1,T), K), nrow=K*T, ncol=T, byrow=TRUE), Xvec)
  #above as below
  #Dtemp<- cbind(kronecker(rep(1, K), diag(1,T)), Xvec)
  
  Ttemp<-matrix(data=rep(as.vector(diag(1,T)), K), nrow=K, ncol=T*T, byrow=TRUE)
  Dtemp<- cbind(Ttemp, Xmat)


  #Generate the D matrix for all clusters
  # D <- kronecker(diag(1,2), Dvec)
  
  Darray <- array(as.vector(t(Dtemp)), c(T,T+1,K))
  
  D <- Darray[,,1]
  for(i in 2:K) D <- bdiag(D, Darray[,,i])
  
  #Remove those rows of D which correspond to clusters not observed in 
  #particular periods
  D<- D[!is.na(Xvec),]
  
  varYbar <- D%*%kronecker(diag(1,K), Vi)%*%t(D) + sig2*diag(1, nrow(D))
  
  #Variance of the treatment effect estimator is then given by:
  return(solve(t(Zmat)%*%solve(varYbar)%*%Zmat)[ncol(Zmat),ncol(Zmat)] )

  
}

#Designs missing an entire period incorporating discrete-time decay

CRTVar_TEHet_ocol_decay <- function(Xmat, ocol,  m, sigu2, sigv2, siguv =0, sig2E, type = 0, r) {
  #Xmat: the design matrix
  #m: number if subjects in each cluster-period
  #sigu2: variance of cluster random intercepts
  #sigv2: variance of treatment effect slopes
  #siguv: covariance of random intercept and slope
  #sig2E: variance of subject-level errors
  
  sig2 <- sig2E/m
  
  T <- ncol(Xmat[,-ocol])
  K <- nrow(Xmat)
  Xvec <- as.vector(t(Xmat[,-ocol]))
  
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  
  
  #Need to create the Vi for the *original* design matrix
  #and then chop out the relevant row and column
  #T+1 is to allow for this.
  
  #Constant decay var if type==0
  if(type==0) { 
    Vi <-diag((1-r)*sigu2, T+1) + matrix(data=sigu2*r, nrow=T+1, ncol=T+1)
  }
  #exponential decay structure
  if(type==1) { 
    Vi <- sigu2*(r^abs(matrix(1:(T+1),nrow=T+1, ncol=T+1, byrow=FALSE) - matrix(1:(T+1),nrow=T+1, ncol=T+1, byrow=TRUE)))
  }
  Vi <- Vi[-ocol, -ocol]
  
  #Need to add on the treatment effect components:
  Vi <- rbind(Vi, rep(siguv, T))
  Vi <- cbind(Vi, c(rep(siguv, T), sigv2))
  
  #cbind(matrix(data=rep(diag(1,T), K), nrow=K*T, ncol=T, byrow=TRUE), Xvec)
  #above as below
  #Dtemp<- cbind(kronecker(rep(1, K), diag(1,T)), Xvec)
  
  Ttemp<-matrix(data=rep(as.vector(diag(1,T)), K), nrow=K, ncol=T*T, byrow=TRUE)
  Dtemp<- cbind(Ttemp, Xmat[,-ocol])
  
  
  #Generate the D matrix for all clusters
  # D <- kronecker(diag(1,2), Dvec)
  
  Darray <- array(as.vector(t(Dtemp)), c(T,T+1,K))
  
  D <- Darray[,,1]
  for(i in 2:K) D <- bdiag(D, Darray[,,i])
  
  #Remove those rows of D which correspond to clusters not observed in 
  #particular periods
  D<- D[!is.na(Xvec),]
  
  varYbar <- D%*%kronecker(diag(1,K), Vi)%*%t(D) + sig2*diag(1, nrow(D))
  
  #Variance of the treatment effect estimator is then given by:
  return(solve(t(Zmat)%*%solve(varYbar)%*%Zmat)[ncol(Zmat),ncol(Zmat)] )
  
  
#  #Constant decay var if type==0
#  if(type==0) { 
#    Vi <-diag((1-r)*sigu2, T) + matrix(data=sigu2*r, nrow=T, ncol=T)
#  }
#  #exponential decay structure
#  if(type==1) { 
#    Vi <- sigu2*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
#  }
#  
#  #Need to add on the treatment effect components:
#  Vi <- rbind(Vi, rep(siguv, T))
#  Vi <- cbind(Vi, c(rep(siguv, T), sigv2))
#
#  Ttemp<-matrix(data=rep(as.vector(diag(1,T)), K), nrow=K, ncol=T*T, byrow=TRUE)
#  Dtemp<- cbind(Ttemp, Xmat)
#  
#  
#  #Generate the D matrix for all clusters
#  # D <- kronecker(diag(1,2), Dvec)
#  
#  Darray <- array(as.vector(t(Dtemp)), c(T,T+1,K))
#  
#  D <- Darray[,,1]
#  for(i in 2:K) D <- bdiag(D, Darray[,,i])
#  
#  #Remove those rows of D which correspond to clusters not observed in 
#  #particular periods
#  D<- D[!is.na(Xvec),]
#  
#  varYbar <- D%*%kronecker(diag(1,K), Vi)%*%t(D) + sig2*diag(1, nrow(D))
#  
#  #Variance of the treatment effect estimator is then given by:
#  return(solve(t(Zmat)%*%solve(varYbar)%*%Zmat)[ncol(Zmat),ncol(Zmat)] )
#  
#  #Think about the problems when there are multiple removed observations later:
#  ##there will be problems if Zmat is not of full column rank
#  ##if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
#  ##else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
#  
#  #return(vartheta)
  
}


#################################
#General function (taken from ExpDecayVar.r in Information content folder)
CRTVarGeneral <- function(Xmat, m, rho0, r, type) {
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  
  T <- ncol(Xmat)
  K <- nrow(Xmat)
  
  
  Xvec <- as.vector(t(Xmat))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  
  #Variance matrix for one cluster, with decay in correlation over time
  #Vi <- diag(sig2,T) + matrix(sig2CP,nrow=T, ncol=T)
  #Constant decay var if type==0
  if(type==0) { 
    Vi <-diag(sig2 +(1-r)*sig2CP, T) + matrix(data=sig2CP*r, nrow=T, ncol=T)
  }
  #exponential decay structure
  if(type==1) { 
    Vi <- diag(sig2,T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  }
  #Variance matrix for all clusters
  Vall <- kronecker(diag(1,K), Vi)
  Vall <- Vall[!is.na(Xvec),!is.na(Xvec)]
  
  
  #vartheta <- solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)]
  
  #xtry <- try(solve((t(Zmat)%*%solve(Vall)%*%Zmat))) 
  #if('try-error' %in% class(xtry)) return(NA)
  
  #there will be problems if Zmat is not of full column rank
  if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
  
  #return(vartheta)
  
}

