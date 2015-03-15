SimExample <-
function(n=100, p, model, rho){
    
    if(missing(model)){
      print("ERROR: no model specified.")
      return();
    }else if(toupper(model)=="WPD2"){
        

      
      #number of variables
      if(missing(p)) p <- 57
      if(p<14){
        print("ERROR: p ge 14 is required for model WPD2.")
        return();
      }
      n_true <- 14
      n_noise <- p - n_true 
      if(n<p) print("warnings: n < p")
      
      #covariance struction for data
      if(missing(rho)) rho = 0.7
      S <- diag(15)
      for(k in 1:8){
        for(l in 1:8){
          S[k,l] <- S[l,k] <- rho^(k+l-2)
        }
        S[k,k] <- 1
      }
      
      true <- 1:15
      beta <- rep(0,p)
      beta[true] <- rep(0.1,dim(S)[1])
      
      
      #generate data
      C <- chol(S)
      U <-matrix(rnorm(n*(n_true+1)),n,n_true+1) %*% C
      y <- rnorm(n, U%*%beta[true],1)
      X <- as.matrix(U[,1],n,1)
      U <- as.matrix(cbind(U[,-1],matrix(rnorm(n_noise*n),n,n_noise)))
      
      return(list(y=y,X=X,U=U,p=p,beta=beta, rho=rho, model=model))
    }else if(toupper(model)=="BR1"){
      
        


      if(missing(p)) p <- 50
      if(p<40){
        print("ERROR: p ge 40 is required for model BR1.")
        return();
      }
      n_true <- 10
      n_noise <- p - n_true
      if(n<p) print("warnings: n < p")
      
      
      #covariance struction for data
      if(missing(rho)) rho = 0.5
      S <- diag(p)
      for(k in 1:p){
        for(l in 1:p){
          S[k,l] <- S[l,k] <- rho^abs(k-l)
        }
        S[k,k] <- 1
      }
      
      true <- c(11:15,36:40)
      beta <- rep(0,p)
      beta[true] <- runif(length(true))*5
      
      
      #generate data
      C <- chol(S)
      X <- matrix(rnorm(n*p),n,p) %*% C
      y <- rnorm(n, X[,true]%*%beta[true],1)
      
      
      return(list(y=y,X=X,p=p,beta=beta, rho=rho, model=model))
    }else if(toupper(model)=="BR2"){
      
        

        
      if(missing(p)) p <- 50
      if(p<40){
        print("ERROR: p ge 3 is required for model BR2.")
        return();
      }
      n_true <- 3
      n_noise <- p - n_true
      if(n<p) print("warnings: n < p")
      
      
      #covariance struction for data
      S <- rwish(p,diag(p))
      true <- c(1:3)
      beta <- rep(0,p)
      beta[true] <- runif(length(true))
      
      
      #generate data
      C <- chol(S)
      X <- matrix(rnorm(n*p),n,p) %*% C
      y <- rnorm(n, X[,true]%*%beta[true],1)
      
      
      return(list(y=y,X=X,p=p,beta=beta, model=model))
    }
  }
