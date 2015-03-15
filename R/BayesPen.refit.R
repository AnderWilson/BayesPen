BayesPen.refit <-
function(y,x,fit,joint,max.refit,...){
  
  n <- length(y)
  p <- ncol(x)
  
  if(missing(joint)){
    joint<-fit$joint
  }else if(joint==TRUE & fit$joint==FALSE){
    print("Joint credible regions are not provided.  Marginal will be used.")
    joint <- FALSE
  }
    

    if(missing(max.refit)){
        if(joint){
            max.refit <- min(200,n-1,nrow(fit$joint.path))
        }else{
            max.refit <- min(200,(n-1),p,max(fit$marginal.path))
        }
        print("============================================================")
        print("max.refit missing.")
        print(paste("max.refit will be set to ",max.refit,".",sep=""))
        print("============================================================")
    }else if(is.null(max.refit)){
        if(joint){
            max.refit <- min(200,n-1,nrow(fit$joint.path))
        }else{
            max.refit <- min(200,(n-1),p,max(fit$marginal.path))
        }
        print("============================================================")
        print("max.refit null.")
        print(paste("max.refit will be set to ",max.refit,".",sep=""))
        print("============================================================")
    }
  
  if(joint){
    adapt.ind<-fit$joint.path
    df.joint<-apply(adapt.ind,1,sum)
    
    coefs.joint<-matrix(0,nrow=min(max.refit,n-1,length(df.joint)),ncol=p)
    SSE.joint<-rep(0,min(max.refit,n-1,length(df.joint))) 
    dev.joint<-rep(0,min(max.refit,n-1,length(df.joint)))
    df.joint <- df.joint[1:nrow(coefs.joint)]
    
    for (dloc in 1:min(max.refit,n-1,length(df.joint)))
    {
      if (df.joint[dloc]>0)
      {
        X_curr<-x[,adapt.ind[dloc,]!=0]
        fits.joint<-glm(y~X_curr-1,...)
        coefs.joint[dloc,adapt.ind[dloc,]!=0]<-fits.joint$coef
        SSE.joint[dloc]<-t(fits.joint$res)%*%fits.joint$res
        dev.joint[dloc] <- fits.joint$deviance
      }else{
        SSE.joint[dloc]<-t(y)%*%y
        dev.joint[dloc] <- NA
      }
    }
    coefs.joint=coefs.joint[df.joint<=min((n-1),p),]
    SSE.joint=SSE.joint[df.joint<=min((n-1),p)]
    dev.joint=dev.joint[df.joint<=min((n-1),p)]
    return(list(coefs=coefs.joint,SSE=SSE.joint,dev=dev.joint, df=df.joint,joint=joint))
  }else{
    ordered.z.s <- fit$marginal.path
    coefs.marg<-matrix(0,nrow=min(max.refit,(n-1),p,max(ordered.z.s)),ncol=p)
    SSE.marg<-dev.marg<-df.marg<-rep(0,min(max.refit,(n-1),p,max(ordered.z.s)))
    
    for (dloc in 1:min(max.refit,(n-1),p,max(ordered.z.s)))
    {
      X_curr<-x[,(ordered.z.s<=dloc)]
      fits.marg<-glm(y~X_curr-1,...)
      coefs.marg.temp<-rep(0,p)
      coefs.marg[dloc,(ordered.z.s<=dloc)]<-fits.marg$coef
      SSE.marg[dloc]<-t(fits.marg$res)%*%fits.marg$res
      dev.marg[dloc]<-fits.marg$deviance
      df.marg[dloc]<-length(which(ordered.z.s<=dloc))
    }
    return(list(coefs=coefs.marg,SSE=SSE.marg, dev=dev.marg ,df=df.marg,joint=joint))
    
  }
}
