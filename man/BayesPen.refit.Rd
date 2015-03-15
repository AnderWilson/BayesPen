\name{BayesPen.refit}
\alias{BayesPen.refit}
\title{ Bayesian Penalized Credible Regions Solution Path Refit }
\description{ Refits the solution path given by BayesPen. }
\usage{
BayesPen.refit(y, x, fit, joint, max.refit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{	A n-vector of responses. If fit is a list from BayesPen.lm or BayesPen.lm.confounders then y is not required.}
  \item{x}{ A n x p design matrix that includes all potential covariates. In the confounder selection case this includes the exposures and confounders, i.e. cbind(x,u). If fit is a list from BayesPen.lm or BayesPen.lm.confounders then x is not required.}
  \item{fit}{ A list returned from BayesPen, BayesPen.lm, or BayesPen.lm.confounders. }
  \item{joint}{ For variable selection this indicates if the joint or marginal solution path should be used. Joint must be TRUE for confounder selection.}
  \item{max.refit}{ The maximum number of models to be refit.}
  \item{\dots}{ These are additional terms passed to glm to refit a glm other than linear. The default is linear.}
}
\details{
This refits each model in the solution path with the frequentist model using the glm function.    
}
\value{
\item{coefs }{A matrix of regression coefficients for each model in the solution path. The regression coefficients for parameters omitted from a model are set to 0.}
\item{SSE }{SSE of each refitted model.}
\item{dev }{Deviance of each refitted model.}
\item{df }{ Error degrees of freedom from each refitted model.}
\item{joint }{Returns the logical joint.}
}
\references{
Bondell, H. D. and Reich, B. J. (2012). Consistent high-dimensional Bayesian variable selection via penalized credible regions. \emph{J. Am. Statist. Assoc}. 107, 1610-1624.

Wilson A., Reich B. J. (2014). Confounder selection via penalized credible regions. \emph{Biometrics} 70: 852-861.
}
\author{
Ander Wilson, Howard D. Bondell, and Brian J. Reich
}

\seealso{
 \code{\link{BayesPen}}
 }
\examples{

######################
#Variable Selection
set.seed(1234)
dat <- SimExample(500,model="BR1")
X <- dat$X
y <- dat$y

#fit the full model assuming flat priors on beta
betahat <- solve(t(X)\%*\%X) \%*\% t(X) \%*\% y
cov <- solve(t(X)\%*\%X) * sum((X\%*\%betahat-y)^2)/(length(y)-length(betahat))

#find solution path
fit.BayesPen <- BayesPen(beta=betahat, beta_cov=cov)

#refit the model
refit <- BayesPen.refit(y,X,fit.BayesPen)

#plot it
BayesPen.plot(refit)

}

\keyword{ regression }
