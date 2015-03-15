\name{BayesPen}
\alias{BayesPen}
\title{ 
Bayesian Penalized Credible Regions 
}
\description{
Fits the Bayesian penalized credible regions method from the posterior mean and covariance of a Bayesian model. 
The function performs variable selection (Bondell and Reich 2012) and confounder selection (Wilson and Reich 2014). 
The default is variable selection and confounder selection is only performed if confounder weights are provided.
}
\usage{
BayesPen(beta, beta_cov, joint, force = NULL, confounder.weights, max.steps = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
p-vector of posterior means from the fitting of the full regression model. For confounder selection this is the outcome model.
}
  \item{beta_cov}{
Posterior covariance matrix corresponding to beta.
}
  \item{joint}{
Indicator if joint credible regions approach should be used. 
If joint=FALSE the marginal approach of Bondell and Reich (2012) will be used.
The marginal approach is only available for variable selection, not confounder selection.
}
  \item{force}{
A vector of columns corresponding to variables that are forced to  be included in the model. For example, this may include an intercept. 
For confounder selection the exposure(s) of interest should be forced into the model.}
  \item{confounder.weights}{
The posterior mean from the exposure model for confounder selection. 
For a single exposure this is a p-vector with exposure model regression coefficients in the same order as in beta. 
For multiple exposures this is a matrix with p rows and a column for each exposure.
The locations corresponding to exposure(s) in beta can be set to any numeric value; they are not used.
For variable selection this is omitted.
}
  \item{max.steps}{
Maximum number of steps to be performed in the LARS algorithm (Hastie and Efron 2013).
}
}
\value{

\item{joint.path}{A complete solution path for the joint credible regions approach. Each row is a model in the solution path with a 1 indicating a variable is included and a 0 indicating it is not included.}
\item{marginal.path}{A complete solution path for the marginal credible regions approach. The p-vector denotes the step at which each covariate is included in the model.}
\item{order.path}{The action returned from lars that shows when each covariate is added to the model.}
\item{order.marg}{The the covariate added at each step.}
\item{joint}{Returns a vector indicating which variables are forced into the model. }
\item{force}{Returns the logical joint.}
}
\references{
Bondell, H. D. and Reich, B. J. (2012). Consistent high-dimensional Bayesian variable selection via penalized credible regions. \emph{J. Am. Statist. Assoc}. 107, 1610-1624.

Trevor Hastie and Brad Efron (2013). lars: Least Angle Regression, Lasso and Forward Stagewise. R package version 1.2. http://CRAN.R-project.org/package=lars

Wilson A., Reich B. J. (2014). Confounder selection via penalized credible regions. \emph{Biometrics} 70: 852-861.
}
\author{
Ander Wilson, Howard D. Bondell, and Brian J. Reich
}
\seealso{
 \code{\link{BayesPen.lm}}, \code{\link{BayesPen.lm.confounders}}
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

\keyword{ regression }
