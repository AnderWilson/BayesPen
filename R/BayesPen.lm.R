BayesPen.lm <-
function(y,x,prior, nIter=500, burnIn=100, joint, force=NULL, max.steps=NULL, max.refit){
  if(missing(prior)) prior<- list(a1=0.1,b1=0.1,a2=0.1,b2=0.2)
  p <- ncol(x)
  if(missing(joint) & p<=1000){
        print("============================================================")
        print("Joint not specified.")
        print("Joint credible sets will be used for variable selection.")
        print("============================================================")
        joint<-FALSE
    }else if(missing(joint)){
        print("============================================================")
        print("Joint not specified and p>1000.")
        print("Margingal credible sets will be used for variable selection.")
        print("============================================================")
        joint<-FALSE
    }
  if(missing(max.refit)) max.refit <- NULL

  
  if (nIter < ncol(x) && joint){
    stop("The number of MCMC draws must be set to be larger than the dimension to use the joint posterior covariance matrix.")
  }

  fit.lm <- BayesPen.MCMC(y=y, X=x, prior=prior, nIter = nIter, burnIn = burnIn)
  
  print("Model fitting complete.  Start post processing")
  print("------------------------------------------------------------")
  fit <- BayesPen(beta=fit.lm$MN, beta_cov=fit.lm$COV,joint=joint,force=force,max.steps=max.steps)

  print("Refit Model")
  print("------------------------------------------------------------")
  refit <- BayesPen.refit(y-mean(y),x,fit=fit, max.refit=max.refit)

  print("Complete")
  print("------------------------------------------------------------")
  out <- c(fit,refit[names(refit)[-grep("joint",names(refit))]])
  out$lm <- fit.lm
  return(out)

}
