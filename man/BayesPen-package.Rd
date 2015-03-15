\name{BayesPen-package}
\alias{BayesPen-package}
\docType{package}
\title{
Bayesian Penalized Credible Regions
}
\description{
Performs the Bayesian penalized credible regions methods for variable selection (Bondell and Reich 2012) and confounder selection (Wilson and Reich 2014). 
The solution path can be computed from the posterior mean and covariance matrix for a variety of models. 
The function BayesPen performs confounder or variable selection from the posterior means and covariances of a Bayesian regression model. 
The wrapper functions BayesPen.lm and BayesPen.lm.confounders fit linear models with code based on the BLM package (Campos and Rodriguez 2012) and run BayesPen together in one step.
}




\details{
\tabular{ll}{
Package: \tab BayesPen\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-02-16\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Ander Wilson, Howard D. Bondell, and Brian J. Reich

Maintainer: Ander Wilson <awilson@hsph.harvard.edu>
}
\references{
Bondell, H. D. and Reich, B. J. (2012). Consistent high-dimensional Bayesian variable selection via penalized credible regions. \emph{J. Am. Statist. Assoc}. 107, 1610-1624.

Wilson A., Reich B. J. (2014). Confounder selection via penalized credible regions. \emph{Biometrics} 70: 852-861.
}
\examples{

######################
#Variable Selection
set.seed(1234)
dat <- SimExample(500,model="BR1")
X <- dat$X
y <- dat$y

#fit the full model assuming flat priors on beta
fit1 <- lm(y~X-1)
betahat <- coef(fit1)
cov <- vcov(fit1)

#find solution path
fit.BayesPen <- BayesPen(beta=betahat, beta_cov=cov)

#refit the model
refit <- BayesPen.refit(y,X,fit.BayesPen)

#plot it
BayesPen.plot(refit)


######################
#Confounder Selection
set.seed(1234)
dat <- SimExample(500,model="WPD2")
X <- dat$X
U <- dat$U
W <- cbind(X,U)
y <- dat$y

#fit the full outcome model assuming flat priors on beta
fit1 <- lm(y~W-1)
betahat <- coef(fit1)
cov <- vcov(fit1)

#fit the full exposure model assuming flat priors on beta
fit2 <- lm(X~U-1)
gammahat <- coef(fit2)

#find solution path
fit.BayesPen <- BayesPen(beta=betahat, beta_cov=cov, confounder.weights=c(0,gammahat), force=1)

#refit the model
refit <- BayesPen.refit(y,W,fit.BayesPen)

#plot it
BayesPen.plot(refit)
}

