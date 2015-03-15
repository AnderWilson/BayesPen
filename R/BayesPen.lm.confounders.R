BayesPen.lm.confounders <-
function(y,x,u,prior, nIter=500, burnIn=100, force, max.steps=NULL, max.refit, include.me=FALSE, z.score=FALSE){
    
  if(missing(prior)) prior<- list(a1=0.1,b1=0.1,a2=0.1,b2=0.2)
  if(missing(force)) force <- NULL
  if(missing(max.refit)) max.refit <- NULL
  
  if (nIter < ncol(x)){
    stop("The number of MCMC draws must be set to be larger than the dimension for cronfounder selection.")
  }
  
  x <- as.matrix(x)
  fit.out<-BayesPen.MCMC(y=y, X=cbind(x,u), prior=prior, nIter = nIter, burnIn = burnIn)
  if(is.null(dim(x)) | ncol(x)==1){
    confounder.weights <- c(0,abs(BayesPen.MCMC(y=x, X=u, prior=prior, nIter = nIter, burnIn = burnIn)$MN))
    force <- c(1,force+1)
  }else{
    confounder.weights<-rep(0,ncol(u)+ncol(x))
    for(i in 1:ncol(x)){
      if(include.me){
        fit.temp <- BayesPen.MCMC(y=x[,i], X=cbind(x[,-i],u), prior=prior, nIter = nIter, burnIn = burnIn)
        weight.temp <- abs(fit.temp$MN)
        if(z.score) weight.temp <- weight.temp/sqrt(diag(fit.temp$COV))*sqrt(diag(fit.out$COV))[-1]
        confounder.weights <- confounder.weights+c(rep(0,ncol(x)),weight.temp[-c(1:c(ncol(x)-1))])
      }else{
        fit.temp <- BayesPen.MCMC(y=x[,i], X=u, prior=prior, nIter = nIter, burnIn = burnIn)
        weight.temp <- abs(fit.temp$MN)
        if(z.score) weight.temp <- weight.temp/sqrt(diag(fit.temp$COV))*sqrt(diag(fit.out$COV)[-c(1:3)])
        confounder.weights <- confounder.weights+c(rep(0,ncol(x)),weight.temp)
      }
    }
    force <- c(1:ncol(x),force+ncol(x))
  }
  print("Model fitting complete, start post processing.")
  print("------------------------------------------------------------")

  fit <- BayesPen(beta=fit.out$MN, beta_cov=fit.out$COV,joint=TRUE,confounder.weights=confounder.weights,force=force,max.steps=max.steps)
  
  print("Refit Model")
  print("------------------------------------------------------------")
  refit <- BayesPen.refit(y-mean(y),cbind(x,u),fit=fit, max.refit=max.refit)

  print("Complete")
  print("------------------------------------------------------------")
  out <- c(fit,refit[names(refit)[-grep("joint",names(refit))]])
  out$lm <- fit.out
  out$confounder.weights <- confounder.weights
  return(out)
  
}
