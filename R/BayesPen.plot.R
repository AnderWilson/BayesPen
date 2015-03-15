BayesPen.plot <-
function(refit, ...) matplot(refit$coefs, type="l", las=1, main="BPCR Solution Path", ylab="Beta hat",xlab="Step", ...)
