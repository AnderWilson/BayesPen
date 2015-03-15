\name{BayesPen.plot}
\alias{BayesPen.plot}
\title{
Plot Bayesian Penalized Credible Region Solution Path.
}
\description{
Plots the solution path for Bayesian penalized credible regions from a an object returned by BayesPen.refit.
}
\usage{
BayesPen.plot(refit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refit}{Object returned from BayesPen.refit.}
  \item{\dots}{Additional graphics parameters from matplot.}
}
\references{
Bondell, H. D. and Reich, B. J. (2012). Consistent high-dimensional Bayesian variable selection via penalized credible regions. \emph{J. Am. Statist. Assoc}. 107, 1610-1624.

Wilson A., Reich B. J. (2014). Confounder selection via penalized credible regions. \emph{Biometrics} 70: 852-861.
}
\author{
Ander Wilson, Howard D. Bondell, and Brian J. Reich
}
\seealso{
 \code{\link{BayesPen}},\code{\link{BayesPen.refit}}
 }
\examples{
set.seed(1234)
dat <- SimExample(500,model="WPD2")
fit.BRnew <- BayesPen.lm.confounders(y=dat$y,x=dat$X, u=dat$U)
BayesPen.plot(fit.BRnew,
			  lty=1, #change line type
			  lwd=c(10,rep(2,nrow(fit.BRnew$coef)-1)), #change width
			  col=c("black","grey",rep("blue",7),rep("red",7),
			  rep("grey",nrow(fit.BRnew$coef)-15)), #change the colors
			   ylim=c(-.1,.4) #set limits
			  )
legend("topright", lty=1, lwd=c(10,2,2,2), col=c("black","blue","red","grey"), 
	   legend=c("Exposure","Confounder","Covariate","Noise"))
}
\keyword{ hplot }

