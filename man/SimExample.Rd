\name{SimExample}
\alias{SimExample}
\title{
Simulate Example Data for BayesPen
}
\description{
Simulates example data used in Bondell and Reich (2012), Wang et. al (2012), and Wilson and Reich (2014).
}
\usage{
SimExample(n = 100, p, model, rho)
}
\arguments{
  \item{n}{ Number of observations. }  
\item{p}{ The total number of covariates (including the exposure of interest for WPD2). }
\item{model}{ What model to simulate data from.  WPD2 is design 2 from Wang et. al (2012).  BR1 and BR2 are designs 1 and 2 from Bondell and Reich (2012). }
\item{rho}{ This specifies the correlation between covariates in WPD2 and BR2. }
}
\value{
 \item{y }{n vector of responses.}
 \item{X }{n vector of exposures for WPD and the n x p design matrix for BR1 and BR2.}
 \item{U }{n x p matrix of potential confounders for WPD2.This is missing for BR1 and BR2.}
 \item{p }{Total number of potential confounders}
 \item{beta }{The regression coeficients.  For WPD2 the first beta corresponds to X.}
 \item{rho}{Returns rho.}
 \item{model}{Returns model.}
}
\references{
Bondell, H. D. and Reich, B. J. (2012). Consistent high-dimensional Bayesian variable selection via penalized credible regions. \emph{J. Am. Statist. Assoc}. 107, 1610-1624.

Wang, C., Parmigiani, G., and Dominici, F. (2012). Bayesian effect estimation accounting for adjustment uncertainty. \emph{Biometrics} 68, 661-671.

Wilson A., Reich B. J. (2014). Confounder selection via penalized credible regions. \emph{Biometrics} 70: 852-861.
}
\author{
Ander Wilson, Howard D. Bondell, and Brian J. Reich
}
\examples{
set.seed(1234)
dat <- SimExample(500,model="BR1")
lm.fit <- lm(dat$y~dat$X)
}
\keyword{ datagen }