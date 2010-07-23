\name{COZIGAM.cts}
\alias{COZIGAM.cts}
\title{Fitting Component-Specific-Proportionally Constrained 
Zero-Inflated Generalized Additive Models with Continuous Responses}
\description{
\code{COZIGAM.cts} fits a COZIGAM with component-specific-proportional constraint 
for continuous exponential family responses. 

\code{\link{COZIGAM.dis}} fits a COZIGAM with discrete exponential family 
responses. One of the two funcions is called from \code{\link{cozigam}}
automatically depending on the distribution family if
the argument \code{constraint="component"}. 
}
\author{
Hai Liu and Kung-Sik Chan
}

\references{
Liu, H and Chan, K.S. (2010)  Introducing COZIGAM: An R Package for Unconstrained and Constrained Zero-Inflated Generalized Additive Model Analysis. Journal of Statistical Software, 35(11), 1-26.
\url{http://www.jstatsoft.org/v35/i11/}
}
\seealso{
   \code{\link{cozigam}}
}
\keyword{smooth} \keyword{models} \keyword{regression}


