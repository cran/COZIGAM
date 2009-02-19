\name{summary.cozigam}
\alias{summary.cozigam}
\alias{print.summary.cozigam}
\title{Summary for a COZIGAM fit}
\description{
   Takes a fitted \code{cozigam} object produced by \code{cozigam()}
   and produces various useful summaries from it.
}
\usage{
\method{summary}{cozigam}(object, dispersion = NULL, ...)

\method{print}{summary.cozigam}(x,digits = max(3, getOption("digits") - 3),
                  signif.stars = getOption("show.signif.stars"), ...)
}
\arguments{
   \item{object}{a fitted \code{cozigam} object as produced 
   by \code{cozigam()}.}
   \item{dispersion}{a known dispersion parameter. \code{NULL} to use 
   estimate or default (e.g. 1 for Poisson).}
   \item{x}{a \code{summary.cozigam} object produced by 
   \code{summary.cozigam()}.}
   \item{digits}{the number of significant digits to use when printing.}
   \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" 
   are printed for each coefficient.}
   \item{...}{other arguments.}
}
\details{
\code{print.summary.cozigam} tries to be smart about formatting the coefficients, standard errors, etc.
and additionally gives ``significance stars" if signif.stars is \code{TRUE}.
}
\value{
\code{summary.cozigam} produces a list of summary information 
for a fitted \code{cozigam} object.
   \item{p.coeff}{an array of estimates of the strictly parametric model 
   coefficients, including the linear constraints parameters.}
   \item{p.t}{an array of the \code{p.coeff}'s divided by their standard errors.}
   \item{p.pv}{an array of p-values for the null hypothesis that the corresponding parameter is zero.
   Calculated with reference to the t distribution with the estimated residual degrees of freedom for
   the model fit if the dispersion parameter has been estimated, and the standard normal if not.}
   \item{m}{the number of smooth terms in the model.}
   \item{chi.sq}{an array of test statistics for assessing the significance of model smooth terms.
   If \eqn{b_i} is the parameter vector for the \code{i}th smooth term, and this term has 
   estimated covariance matrix \eqn{V_i} then the statistic is \eqn{b_i'V_i^{k-}b_i}, 
   where \eqn{V_i^{k-}} is the rank \code{k} pseudo-inverse of \eqn{V_i}, and \code{k} is 
   estimated rank of \eqn{V_i}.}
   \item{s.pv}{an array of approximate p-values for the null hypotheses that each smooth term is zero.
   Be warned, these are only approximate. In the case of known dispersion parameter, 
   they are obtained by comparing the chi.sq statistic given above to the chi-squared
   distribution with \code{k} degrees of freedom, where \code{k} is the estimated rank of 
   \eqn{V_i}. If the dispersion parameter is unknown (in which case it will have been estimated)
   the statistic is compared to an F distribution with \code{k} upper d.f. and lower d.f. 
   given by the residual degrees of freedom for the model.Typically the p-values will be 
   somewhat too low, because they are conditional on the smoothing parameters,
   which are usually uncertain, but note that the statistic can also have low power if the 
   rank, \code{k}, is too high relative to the EDF of the term.}
   \item{se}{array of standard error estimates for all parameter estimates.}
   \item{edf}{array of estimated degrees of freedom for the model terms.}
   \item{residual.df}{estimated residual degrees of freedom.}
   \item{n}{number of data.}
   \item{family}{the family used.}
   \item{formula}{the original GAM formula.}
   \item{dispersion}{estimated (or given) scale parameter.}
   \item{cov.unscaled}{the estimated covariance matrix of the parameters, divided by scale parameter.}
   \item{cov.scaled}{the estimated covariance matrix of the parameters.}
   \item{p.table}{significance table for parameters.}
   \item{s.table}{significance table for smooths.}
}
\author{
Hai Liu and Kung-Sik Chan
}
\seealso{
   \code{\link{cozigam}}, \code{\link{predict.cozigam}}
}
\examples{
set.seed(1)
n <- 600
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- test(x1, x2)*2 + f0(x3)/5
eta0 <- f/1.1
mu0 <- exp(eta0)  

eta.p10 <- (test(x1,x2) - mean(test(x1,x2)))*2/1.1
eta.p20 <- (f0(x3) - mean(f0(x3)))/5/1.1

alpha0 <- 0.5
delta10 <- 1
delta20 <- 0
eta.p0 <- delta10*eta.p10 + delta20*eta.p20 
p0 <- .Call("logit_linkinv", alpha0 + eta.p0, PACKAGE = "stats")

z <- rbinom(rep(1,n), 1, p0)
y <- rpois(rep(1,n), mu0)
y[z==0] <- 0; rm(z)

res <- cozigam(y~s(x1,x2)+s(x3), constraint="component", zero.delta=c(NA, 0), family=poisson)
summary(res)

}
\keyword{models} \keyword{smooth} \keyword{regression}
