\name{BayesPen.lm}
\alias{BayesPen.lm}
\title{
Variable Selection via Penalized Credible Regions for Linear Models
}
\description{
Fits a Bayesian linear model and performs variable selection via penalized credible regions (Bondell and Reich 2012).
}
\usage{
BayesPen.lm(y, x, prior, nIter=500, burnIn=100, joint, force = NULL,
		    max.steps = NULL, max.refit)
}
\arguments{
  \item{y}{
A n-vector of responses.
}
  \item{x}{
A n x p design matrix.  An intercept is automatically added.
}
  \item{prior}{
A list specifying the priors for the regression coefficients and the error variance. The four elements of the list are a1, b1, a2, and b2. The residuals are assumed to be iid normal mean 0 with with a gamma(a1,b1) hyperprior on the precision. The regression coeficients are iid mean 0 and have a gamma(a2,b2) hyperprior on the precision.
}
  \item{nIter}{
The number of MCMC iterations (integer).
}
  \item{burnIn}{
The number of MCMC iterations to be discarded as burnin (integer).}
  \item{joint}{
Indicator if joint credible regions approach should be used. If joint=FALSE the marginal approach of Bondell and Reich (2012) will be used.}
  \item{force}{
An optional vector indexing which covariates variables should be forced into the model.}
  \item{max.steps}{
Maximum number of steps to be performed in the LARS algorithm (Hastie and Efron 2013).
}
  \item{max.refit}{ The maximum number of models to be refit.}
}
\value{
\item{joint.path}{A complete solution path for the joint credible regions approach. Each row is a model in the solution path with a 1 indicating a variable is included and a 0 indicating it is not included.}
\item{marginal.path}{A complete solution path for the marginal credible regions approach. The p-vector denotes the step at which each covariate is included in the model.}
\item{order.path}{The action returned from lars that shows when each covariate is added to the model.}
\item{order.marg}{The the covariate added at each step.}
\item{joint}{Returns a vector indicating which variables are forced into the model. }
\item{force}{Returns the logical joint.}
\item{coefs }{A matrix of regression coefficients for each model in the solution path. The regression coefficients for parameters omitted from a model are set to 0.}
\item{SSE }{SSE of each refitted model.}
\item{dev }{Deviance of each refitted model.}
\item{df }{ Error degrees of freedom from each refitted model.}
 \item{lm}{ Full fitting of the model.}
}
\references{
Bondell, H. D. and Reich, B. J. (2012). Consistent high-dimensional Bayesian variable selection via penalized credible regions. \emph{J. Am. Statist. Assoc}. 107, 1610-1624.
}
\author{
Ander Wilson, Howard D. Bondell, and Brian J. Reich
}

\examples{
set.seed(1234)
dat <- SimExample(500,model="BR1")
fit.BRnew <- BayesPen.lm(y=dat$y,x=dat$X)
BayesPen.plot(fit.BRnew)
}

\keyword{ regression }
