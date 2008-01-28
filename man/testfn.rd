\name{testfn}
\alias{testfn}
\alias{f0}
\alias{f1}
\title{Some Test Functions}
\description{
\code{testfn} defines a 2-D test function.

\code{f0} and \code{f1} define two 1-D test functions. These test
functions may be used in simulation study of COZIGAM fits.
}
\details{
\deqn{f_0(x)=0.2x^{11}(10(1-x))^{6}+10(10x)^{3}(1-x)^{10}}
\deqn{f_1(x)=\sin(\pi x)}
\deqn{testfn(x,z)=0.3\times 0.4\pi\{1.2e^{-(x-0.2)^2/0.3^2
-(z-0.3)^2}+0.8e^{-(x-0.7)^2/0.3^2-(z-0.8)^2/0.4^2}\}}.
}

\author{
Hai Liu and Kung-Sik Chan
}

\seealso{
   \code{\link{cozigam}}
}
\keyword{smooth} \keyword{models} \keyword{regression}

