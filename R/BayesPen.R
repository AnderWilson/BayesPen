BayesPen <-
function(beta, beta_cov, joint, force=NULL, confounder.weights, max.steps=NULL){
  
  p <- length(beta)
  
  if(missing(joint)){
    if(missing(confounder.weights) & p<=1000){
        print("============================================================")
        print("Joint not specified.")
        print("Joint credible sets will be used for variable selection.")
        print("============================================================")
        joint<-FALSE
    }else if(missing(confounder.weights)){
        print("============================================================")
        print("Joint not specified and p>1000.")
        print("Margingal credible sets will be used for variable selection.")
        print("============================================================")
        joint<-FALSE
    }else{
        joint<-TRUE
    }
  }
  
  
  if(joint){
    prec_est<-solve(beta_cov)
    
    x1 <- chol(prec_est, pivot=TRUE)
    pivot <- attr(x1, "pivot")
    x1 <- x1[, order(pivot)]
    
    y1 <- x1%*%beta
    if(missing(confounder.weights)){
      w <- beta^2 
    }else{
      if(is.null(dim(confounder.weights))){
        w <- (abs(beta)+abs(confounder.weights))^2         
      }else{
        w <- (abs(beta)+rowSums(abs(confounder.weights)))^2 
      }
    }
    w[force] <- 1
    
    xs<-scale(x1,center=FALSE,scale=1/w)
    
    if(!is.null(force)){
      if(max(force)>p | min(force<1)){ print("ERROR: force improperly specified."); return;}
      IPx.force <- diag(rep(1,p)) - xs[,force]%*%solve(t(xs[,force])%*%xs[,force])%*%t(xs[,force])
      xs <- IPx.force%*%xs
      y1 <- IPx.force%*%y1
    }
    
    if(is.null(max.steps)) max.steps = 8 * p
    adapt.cred.set<-lars(xs,y1,type="lasso",normalize=FALSE,intercept=FALSE,max.steps=max.steps,use.Gram=F)
    
    adapt.ind<-1*(coef(adapt.cred.set)!=0)
    
    if(!is.null(force)){
      adapt.ind[,force] <- 1
    }
    order.joint<-unlist(adapt.cred.set$action)
  }else{
    adapt.ind <- NULL
    order.joint <- NULL
  }

  
  if(missing(confounder.weights)){
    z.stats<-beta/sqrt(diag(beta_cov))
    if(!is.null(force)){
      ordered.z.s<-rep(NA,length(z.stats))
      ordered.z.s[-force]<-rank(-abs(z.stats[-force]),ties.method="random")+1
      ordered.z.s[force]<-1
      order.marg<-rep(NA,length(z.stats))
      order.marg[-force]<-sort(ordered.z.s[-force],index.return=T)$ix+1
      order.marg[force]<-1
    }else{
      ordered.z.s<-rank(-abs(z.stats),ties.method="random")
      order.marg=sort(ordered.z.s,index.return=T)$ix
    }
    
  }else{
    ordered.z.s <- NULL
    order.marg <- NULL
    if(!joint){
      print("Confounder weights can only be applied to the joint method.")
      return;
    }
  }

  
  
  
  return(list(joint.path=adapt.ind,marginal.path=ordered.z.s,order.joint=order.joint,order.marg=order.marg,joint=joint, force=force))
}
