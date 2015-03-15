

BayesPen.MCMC<-function(y,X,keep.samples=FALSE,
                         prior=list(a1=0.1,b1=0.1,a2=0.1,b2=0.2),eps=.Machine$double.eps*100,
                         nIter=1000,burnIn=100){

    tick <- proc.time()[3]
    n <- nrow(X)
    p <- ncol(X)

    int  <- mean(y)
    beta <- rep(0,p)
    tau1 <- 1/var(y)    
    tau2 <- 100

    Xb        <- rep(0,n)
    keep.beta <- NULL
    if(keep.samples){keep.beta<-matrix(0,nIter,p)}

    B1 <- B2 <- 0

    E   <- eigen(t(X)%*%X)
    G   <- E$vectors
    D   <- E$values
    D   <- ifelse(D<0,0,D)
    XG  <- X%*%G
    M1  <- as.vector(t(XG)%*%y)
    M2  <- as.vector(t(XG)%*%rep(1,n))

    for(iter in 1:nIter){

   
      # UPDATE BETA

       VVV       <- 1/(tau1*D + tau2)
       beta      <- tau1*VVV*(M1-int*M2)+rnorm(p,0,sqrt(VVV))
       beta      <- as.vector(G%*%beta)
       Xb        <- X%*%beta

      # UPDATE THE VARIANCES

       tau1 <- rgamma(1,n/2+prior$a1,sum((y-int-Xb)^2)/2+prior$b1)
       tau2 <- rgamma(1,p/2+prior$a2,sum(beta^2)/2+prior$b2)

      # UPDATE THE INTERCEPT

       VVV  <- tau1*n + eps
       MMM  <- tau1*sum(y-Xb)
       int  <- rnorm(1,MMM/VVV,1/sqrt(VVV))

      if(keep.samples){
        keep.beta[iter,]<-beta
      }
      if(iter>burnIn){
         B1 <- B1 + beta/(nIter-burnIn)
         B2 <- B2 + outer(beta,beta)/(nIter-burnIn)
      }

    }

    MN  <- as.vector(B1)
    COV <- B2 - outer(MN,MN)
    tock <- proc.time()[3]
    out <- list(MN=MN,COV=COV,beta=keep.beta,time=tock-tick)

return(out)}



