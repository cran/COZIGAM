\name{cozigam}
\alias{cozigam}
\title{Fitting Constrained Zero-Inflated Generalized Additive Models}
\description{
   Fits a COnstrained Zero-Inflated Generalized Additive Model (COZIGAM) to data.
}
\usage{
cozigam (formula, maxiter = 30, conv.crit.in = 1e-5, 
    conv.crit.out = 1e-4, size = NULL, log.tran = FALSE, family, ...)
}
\arguments{
   \item{formula}{a GAM formula. This is exactly like the formula for a GLM except that smooth terms 
   can be added to the right hand side of the formula (and a formula of the form y ~ . is not allowed). 
   Smooth terms are specified by expressions of the form: \code{s(var1,var2,...,k=12,fx=FALSE,bs="tp",by=a.var)} 
   where \code{var1}, \code{var2}, etc. are the covariates which the smooth is a function of and 
   \code{k} is the dimension of the basis used to represent the smooth term. \code{by} can be used to 
   specify a variable by which the smooth should be multiplied.}
   \item{maxiter}{the maximum number of iterations allowed in the estimation procedure.}
   \item{conv.crit.in}{convergence criterion in the inner loop.}
   \item{conv.crit.out}{convergence criterion in the outer loop.}
   \item{size}{optional. Must be specified when \code{family} is \code{binomial}.}
   \item{log.tran}{logical. \code{TRUE} if log-transformation is needed for the response.}
   \item{family}{This is a family object specifying the distribution and link to use in fitting etc. 
   See \code{glm} and \code{family} for more details. Currently support Gaussian/lognormal, Gamma, 
   Poisson and binomial distributions.}
   \item{...}{additional arguments to be passed to the low level regression fitting functions.}
}

\details{
A COnstrained Zero-Inflated Generalized Additive Model (COZIGAM) assumes the response variable \eqn{Y_i} 
is distributed from a zero-inflated 1-parameter exponential family with covariate \eqn{x_i} (could be high dimensional). 
More specifically, \eqn{Y_i} comes from a regular exponential family distribution \eqn{f(x_i)} with 
probability \eqn{p_i} and equals zero with probability \eqn{1-p_i}, with the further assumption that the 
probability of zero-inflation \eqn{p_i} is some monotone function of the mean response function on the link scale, 
e.g., if a logit link is used, \eqn{p_i} is linearly constrained by: \deqn{logit(p_i) = \alpha + \delta  g(\mu_i)},
where \eqn{g()} is the link function of the assumed exponential family and \eqn{\alpha} and \eqn{delta} are 
two unknown parameters to be estimated. 
The mean response is estimated nonparametrically by a Generalized Additive Model (GAM), 
i.e., \eqn{g(\mu_i) = s(x_i)} for some smooth function. 
This bypasses the problems of two popular methods for analyzing zero-inflated data that either focus only 
on the non-zero data or model the presence-absence data and the non-zero data separately.

The estimation approach is a modified penalized iteratively reweighted least squares algorithm (P-IRLS) which requires 
the \code{magic()} function in the \code{mgcv} library. EM algorithm is used if the underlying 
regular exponential family has strictly positive probability mass at zero. 
The covariance matrix of the estimated coefficients is obtained by inverting the 
observed information from Louis' method. See Liu and Chan (2008) for more detail.
}

\value{
An object of class \code{"cozigam"} is a list containing the following components: 
   \item{coefficients}{a named vector of coefficients including the linear constraint parameters.}
   \item{V.theta}{The covariance matrix of the estimated parameters.}
   \item{beta, V.beta}{the estimated parameters of smooth terms with associated covariance matrix V.beta.}
   \item{mu}{the fitted mean values.}
   \item{linear.predictor}{the fitted values on the link scale.}
   \item{dispersion}{(estimated) dispersion parameter.}
   \item{formula}{model formula.}
   \item{p}{the fitted zero-inflation rates.}
   \item{G}{an objected returned by \code{gam()} containing information used in the estimation procedure.}
   \item{psi}{conditional expectation of the zero-inflation indicator (only for the discrete case which involves EM algorithm).}
   \item{family}{the family used.}
   \item{loglik, ploglik}{the (penalized) log-likelihood of the fitted model.}
   \item{y}{the response used.}
   \item{converged}{logical. \code{TRUE} if the iterative procedure is converged.}
   \item{est.disp}{logical. \code{TRUE} if the dispersion parameter is estimated.}
   \item{fit.nonzero}{fitted GAM using only the nonzero data as in the presence/absence analysis}
   \item{score}{the GCV or UBRE score returned by \code{magic()} in the iterative procedure.}
   \item{edf}{estimated degrees of freedom for each model parameter. 
   Penalization means that many of these are less than 1.}
   \item{edf.smooth}{estimated degrees of freedom of each smooth term} 
}

\author{
Hai Liu and Kung-Sik Chan
}

\references{
Liu, H and Chan, K.S. (2008) Constrained Generalized Additive Model with Zero-Inflated Data. 
}
\seealso{
   \code{\link{plot.cozigam}},\code{\link{predict.cozigam}},
   \code{\link{summary.cozigam}}
}
\examples{
## Normal/Log-Normal Response

set.seed(11)
n <- 600
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- test(x1, x2)*4-mean(test(x1, x2)*4) + f0(x3)/2-mean(f0(x3)/2)
sig <- 0.5
mu0 <- f + 3
y <- mu0 + rnorm(n, 0, sig)

alpha0 <- -2.2
delta0 <- 1.2
p0 <- .Call("logit_linkinv", alpha0 + delta0 * mu0, PACKAGE = "stats")
z <- rbinom(rep(1,n), 1, p0)
y[z==0] <- 0

res <- cozigam(y~s(x1,x2)+s(x3), conv.crit.out = 1e-4, log.tran = FALSE, family = gaussian)

plot(res, too.far = 0.1)

## Poisson Response

n <- 500
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- test(x1, x2)*2 + f0(x3)/5
eta0 <- f/1.1
mu0 <- exp(eta0)  

alpha0 <- -1
delta0 <- 2
p0 <- .Call("logit_linkinv", alpha0 + delta0 * eta0, PACKAGE = "stats")

z <- rbinom(rep(1,n), 1, p0)
y <- rpois(rep(1,n), mu0)
y[z==0] <- 0

res <- cozigam(y~s(x1,x2)+s(x3), maxiter = 30, conv.crit.out = 1e-3, family = poisson)

## A Tensor Product Smooth Example
test.ten <- function(x,z,sx=0.3,sz=0.4)  
{ x<-x*20
  (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
  0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}

set.seed(11)
n <- 400
x1 <- runif(n, 0, 1)/20
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- test.ten(x1, x2)*4-mean(test.ten(x1, x2)*4) + f0(x3)/2-mean(f0(x3)/2)
sig <- 0.5
mu0 <- f + 3
y <- mu0 + rnorm(n, 0, sig)

alpha0 <- -2.2
delta0 <- 1.2
p0 <- .Call("logit_linkinv", alpha0 + delta0 * mu0, PACKAGE = "stats")

z <- rbinom(rep(1,n), 1, p0)
y[z==0] <- 0

# If use thin plate spline ...
res.tps <- cozigam(y~s(x1,x2)+s(x3), conv.crit.out = 1e-4, family=gaussian)
par(mfrow=c(1,2))
plot(res.tps, select=1, too.far=0.1)
# Compare with tensor product spline
res.ten <- cozigam(y~te(x1,x2)+s(x3), conv.crit.out = 1e-4, family=gaussian)
plot(res.ten, select=1, too.far=0.1)
}

\keyword{smooth} \keyword{models} \keyword{regression}


