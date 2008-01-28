\name{summary.cozigam}
\alias{summary.cozigam}
\alias{print.summary.cozigam}
\title{Summary for a COZIGAM fit}
\description{
   Produces various useful summaries of a fitted COZIGAM.
}
\usage{
\method{summary}{cozigam}(object, dispersion = NULL, ...)

\method{print}{summary.cozigam}(x,digits = max(3, getOption("digits") - 3),
                  signif.stars = getOption("show.signif.stars"), ...)
}
\arguments{
   \item{object}{A fitted \code{cozigam} object produced by the \code{cozigam} function.}
   \item{dispersion}{A known dispersion parameter value. Set it to \code{NULL},
   the dispersion parameter is set to an 
   estimate or default value from the COZIGAM (e.g. 1 for Poisson).}
   \item{x}{A \code{summary.cozigam} object produced by the \code{summary.cozigam} function.}
   \item{digits}{The number of significant digits to use when printing.}
   \item{signif.stars}{Logical. If \code{TRUE}, "significance stars" 
   are printed for each coefficient.}
   \item{...}{Other arguments.}
}
\details{
\code{print.summary.cozigam} smartly formats the coefficients, standard errors, etc.
and additionally labels coefficients with "significance stars" if \code{signif.stars} is \code{TRUE}.
}
\value{
\code{summary.cozigam} produces a list of summary information 
for a fitted \code{cozigam} object.
   \item{p.coeff}{An array of the strictly parametric estimates, including those of the linear constraints parameters.}
   \item{p.t}{An array of the t-ratios, i.e. estimates in the \code{p.coeff} divided by their standard errors.}
   \item{p.pv}{An array of p-values for the null hypothesis that the corresponding parameter is zero.
   Calculated with reference to the t distribution with the estimated residual degrees of freedom for
   the model fit if the dispersion parameter has been estimated, but otherwise with reference to the standard normal distribution.}
   \item{m}{The number of smooth terms in the model.}
   \item{chi.sq}{An array of test statistics for assessing the significance of model smooth terms.
   If \eqn{b_i} is the parameter vector for the \code{i}th smooth term, and this term has 
   estimated covariance matrix \eqn{V_i} then the statistic is \eqn{b_i'V_i^{k-}b_i}, 
   where \eqn{V_i^{k-}} is the rank \code{k} pseudo-inverse of \eqn{V_i}, and \code{k} is 
   estimated rank of \eqn{V_i}.}
   \item{s.pv}{An array of approximate p-values for the null hypotheses that each smooth term is zero.
   Be warned, these are only approximate values. In the case of known dispersion parameter, 
   they are obtained by comparing the above \code{chi.sq} statistic to the chi-squared
   distribution with \code{k} degrees of freedom, where \code{k} is the estimated rank of 
   \eqn{V_i}. If the dispersion parameter is unknown (in which case it will have been estimated)
   the statistic is compared to an F distribution with \code{k} upper d.f. and lower d.f. 
   given by the residual degrees of freedom for the model.Typically the p-values will be 
   somewhat too low, because they are conditional on the smoothing parameters,
   which are usually uncertain, but note that the statistic can also have low power if the 
   rank, \code{k}, is too high relative to the EDF (Estimated Degrees of Freedom) of the term.}
   \item{se}{An array of standard error estimates for all parameter estimates.}
   \item{edf}{An array of estimated degrees of freedom for the model terms.}
   \item{residual.df}{The estimated residual degrees of freedom.}
   \item{n}{The number of data.}
   \item{family}{The family used.}
   \item{formula}{The original GAM formula.}
   \item{dispersion}{The estimated (or given) scale parameter.}
   \item{cov.unscaled}{The estimated covariance matrix of the parameters, divided by the scale parameter.}
   \item{cov.scaled}{The estimated covariance matrix of the parameters.}
   \item{p.table}{The significance table for the parameter estimates.}
   \item{s.table}{The significance table for smooths.}
}
\author{
Hai Liu and Kung-Sik Chan
}
\seealso{
   \code{\link{cozigam}},\code{\link{predict.cozigam}}
}
\examples{
set.seed(11)
n <- 600
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- testfn(x1, x2)*4-mean(testfn(x1, x2)*4) + f0(x3)/2-mean(f0(x3)/2)
sig <- 0.5
mu0 <- f + 3
y <- mu0 + rnorm(n, 0, sig)

alpha0 <- -2.2
delta0 <- 1.2
p0 <- .Call("logit_linkinv", alpha0 + delta0 * mu0, PACKAGE = "stats")
z <- rbinom(rep(1,n), 1, p0)
y[z==0] <- 0

res <- cozigam(y~s(x1,x2)+s(x3), conv.crit.out = 1e-4, family = gaussian)

summary(res)
}
\keyword{models} \keyword{smooth} \keyword{regression}
