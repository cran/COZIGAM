\name{disgam}
\alias{disgam}
\title{Fitting Discrete Generalized Additive Models with Model Selection Criterion}
\description{
   Fit a discrete Generalized Additive Model (GAM) to data and calculate the logarithmic marginal likelihood.
}
\usage{
disgam (formula, size=NULL, family = poisson(), ...)
}
\arguments{
   \item{formula}{A GAM formula. This is exactly like the formula for a GLM except that smooth terms
   can be added to the right hand side of the formula (and a formula of the form \code{y ~ .} is not allowed).
   Smooth terms are specified by expressions of the form: \code{s(var1,var2,...,k=12,fx=FALSE,bs="tp",by=a.var)}
   where \code{var1}, \code{var2}, etc. are the covariates which the smooth is a function of and
   \code{k} is the dimension of the basis used to represent the smooth term. \code{by} can be used to
   specify a variable by which the smooth should be multiplied.}
   \item{size}{Optional. Must be specified when \code{family} is \code{binomial}.}
   \item{family}{This is a family object specifying the distribution and link to use in fitting etc.
   See \code{glm} and \code{family} for more details. Currently support Poisson and binomial distributions.}
   \item{...}{Additional arguments to be passed to the low level regression fitting functions.}
}

\details{
It is necessary to assess whether there is zero-inflation in count data, e.g., Poisson or binomial data. 
The model selection approach
can be used to determine whether a zero-inflated model is needed. To do that, we can fit both a ZIGAM and a
regular GAM to the data, and then compare the logarithmic marginal likelihoods from these two models. Higher
logarithmic marginal likelihood from the ZIGAM would indicate that there is zero-inflation in the count data.
Otherwise, we can simply fit a regular GAM instead of a ZIGAM. See Liu and Chan (2008) for more detail.
}

\value{
A list containing the following components:
   \item{fit.gam}{A fitted GAM assuming there is no zero-inflation in the data.}
   \item{V.beta}{The estimated covariance matrix of the GAM.}
   \item{mu}{The fitted mean values.}
   \item{formula}{Model formula.}
   \item{family}{The family used.}
   \item{loglik, ploglik}{The (penalized) log-likelihood of the fitted model.}
   \item{logE}{Approximated logarithmic marginal likelihood by Laplace method used for model selection.}
}

\author{
Hai Liu and Kung-Sik Chan
}

\references{
Liu, H and Chan, K.S. (2008) Constrained Generalized Additive Model with Zero-Inflated Data.
Technical Report 388, Department of Statistics and Actuarial Science, The University of Iowa.
\url{http://www.stat.uiowa.edu/techrep/tr388.pdf}}
\seealso{
   \code{\link{zigam}}
}
\examples{
## Poisson Response
set.seed(11)
n <- 200
x1 <- runif(n, 0, 1)

eta0 <- f0(x1)/4 - 0.5
mu0 <- exp(eta0)

y <- rpois(rep(1,n), mu0) # generating non-zero-inflated data

res.gam <- disgam(y~s(x1), family=poisson) # fit a regular GAM
res.zigam <- zigam(y~s(x1), maxiter=10, family=poisson) # fit a ZIGAM

res.gam$logE > res.zigam$logE # compare the model selction criterion

# Another example
set.seed(11)
n <- 200
x1 <- runif(n, 0, 1)
eta0 <- f0(x1)/4 - 0.5
mu0 <- exp(eta0)

alpha0 <- 0.4
delta0 <- 0.8
p0 <- .Call("logit_linkinv", alpha0 + delta0 * eta0, PACKAGE = "stats")

# Generating zero-inflated Poisson count data
z <- rbinom(rep(1,n), 1, p0)
y <- rpois(rep(1,n), mu0)
y[z==0] <- 0

res.gam <- disgam(y~s(x1), family=poisson) # fit a regular GAM
res.zigam <- zigam(y~s(x1), family=poisson) # fit a ZIGAM

res.gam$logE < res.zigam$logE # compare the model selction criterion

}

\keyword{smooth} \keyword{models} \keyword{regression}
