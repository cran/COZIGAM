\name{PCOZIGAM.cts}
\alias{PCOZIGAM.cts}
\title{Fitting Proportionally Constrained Zero-Inflated Generalized Additive Models with
Continuous Responses}
\description{
\code{PCOZIGAM.cts} fits a COZIGAM with proportional constraint for continuous exponential family 
responses. 

\code{\link{PCOZIGAM.dis}} fits a COZIGAM with proportional constraint for discrete exponential family 
responses. One of the two funcions is called from \code{\link{cozigam}}
automatically depending on the distribution family if
the argument \code{constraint="proportional"}. 
}
\author{
Hai Liu and Kung-Sik Chan
}

\references{
Liu, H and Chan, K.S. (2008) Constrained Generalized Additive Model with Zero-Inflated Data. 
Technical Report 388, Department of Statisics and Actuarial Science, The Unversity of Iowa. 
\url{http://www.stat.uiowa.edu/techrep/tr388.pdf} 
}
\seealso{
   \code{\link{cozigam}}
}
\keyword{smooth} \keyword{models} \keyword{regression}


