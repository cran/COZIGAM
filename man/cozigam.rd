\name{cozigam}
\alias{cozigam}
\title{Fitting Constrained Zero-Inflated Generalized Additive Models}
\description{
   Fit a COnstrained Zero-Inflated Generalized Additive Model (COZIGAM) to data.
}
\usage{
cozigam (formula, constraint = "proportional", zero.delta = NULL, maxiter = 20, 
    conv.crit.in = 1e-5, conv.crit.out = 1e-3, size = NULL, log.tran = FALSE, family, ...)
}
\arguments{
   \item{formula}{A GAM formula. This is exactly like the formula for a GLM except that smooth terms 
   can be added to the right hand side of the formula (and a formula of the form \code{y ~ .} is not allowed). 
   Smooth terms are specified by expressions of the form: \code{s(var1,var2,...,k=12,fx=FALSE,bs="tp",by=a.var)} 
   where \code{var1}, \code{var2}, etc. are the covariates which the smooth is a function of and 
   \code{k} is the dimension of the basis used to represent the smooth term. \code{by} can be used to 
   specify a variable by which the smooth should be multiplied.}
   \item{constraint}{Type of constraint on the zero-inflation probability, can be either ``proportional" or ``component". 
   See details for more information.}
   \item{zero.delta}{A vector specifying which subset of constraint parameter delta set to be zero when the
   \code{constraint} is set to ``component". For instance, \code{zero.delta=c(NA, 0)}
   means only the first smooth function is included in the zero-inflation constraint. See details for the model formulation.}
   \item{maxiter}{The maximum number of iterations allowed in the estimation procedure.}
   \item{conv.crit.in}{Convergence criterion in the inner loop.}
   \item{conv.crit.out}{Convergence criterion in the outer loop.}
   \item{size}{Optional. Must be specified when \code{family} is \code{binomial}.}
   \item{log.tran}{Logical. \code{TRUE} if log-transformation is needed for the response.}
   \item{family}{This is a family object specifying the distribution and link to use in fitting etc. 
   See \code{glm} and \code{family} for more details. Currently support Gaussian/lognormal, Gamma, 
   Poisson and binomial distributions.}
   \item{...}{Additional arguments to be passed to the low level regression fitting functions.}
}

\details{
A COnstrained Zero-Inflated Generalized Additive Model (COZIGAM) assumes the response variable \eqn{Y_i} 
is distributed from a zero-inflated 1-parameter exponential family with covariate \eqn{x_i} (could be high dimensional). 
More specifically, \eqn{Y_i} comes from a regular exponential family distribution \eqn{f(x_i)} with 
probability \eqn{p_i} and equals zero with probability \eqn{1-p_i}, with the further assumption that the 
probability of non-zero-inflation \eqn{p_i} is related to the regular mean response or its smooth components, depending on
the type of the \code{constraint}. If \code{constraint} is ``proportional", \eqn{p_i} is assumed to be 
some monotone function of the regular mean response function on the link scale,
e.g., if a logit link is used, \eqn{p_i} is linearly constrained by: \deqn{logit(p_i) = \alpha + \delta  g(\mu_i),}
where \eqn{g()} is the link function of the regular exponential family with mean \eqn{\mu_i}
 and \eqn{\alpha} and \eqn{\delta} are 
two unknown parameters to be estimated; If \code{constraint} is ``component", \eqn{p_i} is assumed to be
linearly related to the smooth components of the regular mean response on the link scale, i.e.,
\deqn{logit(p_i) = \alpha + \delta_1 s(x_{1i},x_{2i}) + \delta_2 s(x_{3i}) + ... ,}
and the regular mean response has the structure
\deqn{g(\mu_i) =  b_0 + s(x_{1i},x_{2i}) + s(x_{3i}) + ... ,}
where \eqn{s()}'s are some centered smooth functions to be estimated nonparametrically by a Generalized Additive Model (GAM) and \eqn{b_0} is the overall intercept.
In the component-specific-proportional constraint case, a vector \code{zero.delta} must be specified. 
For instance, in the above model setting, \code{zero.delta=c(NA, 0)} would estimate \eqn{\delta_1} 
as an unknown parameter and set \eqn{\delta_2} to be
fixed at 0. That is, the model would assume that \eqn{p_i} is not dependent on the covariate \eqn{x_3}.  

The estimation approach is a modified penalized iteratively reweighted least squares algorithm (P-IRLS) which requires 
the \code{magic()} function in the \code{mgcv} library. EM algorithm is used if the underlying 
regular exponential family has strictly positive probability mass at zero. 
The covariance matrix of the estimated coefficients is obtained by inverting the 
observed information matrix. See Liu and Chan (2008) for more detail.
}

\value{
An object of class \code{"cozigam"} is a list containing the following components: 
   \item{coefficients}{A named vector of coefficients including the linear constraint parameters.}
   \item{V.theta}{The covariance matrix of the estimated parameters.}
   \item{V.beta}{The covariance matrix of the estimated parameters associated with the smooth functions.}
   \item{se.alpha, se.delta}{Estimated standard errors of parameters alpha and delta.}
   \item{mu}{The fitted regular mean values.}
   \item{linear.predictor}{The fitted values on the link scale.}
   \item{dispersion}{(Estimated) dispersion parameter.}
   \item{formula}{Model formula.}
   \item{p}{The fitted non-zero-inflation probabilities.}
   \item{G}{An objected returned by \code{gam()} containing information used in the estimation procedure.}
   \item{psi}{Conditional expectation of the zero-inflation indicator (only for the discrete case which involves EM algorithm).}
   \item{family}{The family used.}
   \item{loglik, ploglik}{The (penalized) log-likelihood of the fitted model.}
   \item{y}{The response used.}
   \item{converged}{Logical. \code{TRUE} if the iterative procedure is converged.}
   \item{est.disp}{Logical. \code{TRUE} if the dispersion parameter is estimated.}
   \item{fit.nonzero}{Fitted GAM using only the nonzero data as in the presence/absence analysis}
   \item{score}{The GCV or UBRE score returned by \code{magic()} in the iterative procedure.}
   \item{edf}{Estimated degrees of freedom for each model parameter. 
   Penalization means that many of these are less than 1.}
   \item{edf.smooth}{Estimated degrees of freedom of each smooth term.}
   \item{sp}{Estimated smoothing parameters.}
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
   \code{\link{plot.cozigam}}, \code{\link{predict.cozigam}},
   \code{\link{summary.cozigam}}, \code{\link{zigam}}
}
\examples{
## Normal/Log-Normal Response with proportional constraint
set.seed(11)
n <- 400
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- test(x1,x2)*4-mean(test(x1,x2)*4) + f0(x3)/2-mean(f0(x3)/2)
sig <- 0.5
mu0 <- f + 3
y <- mu0 + rnorm(n, 0, sig)

alpha0 <- -2.2
delta0 <- 1.2
p0 <- .Call("logit_linkinv", alpha0 + delta0 * mu0, PACKAGE = "stats")
z <- rbinom(rep(1,n), 1, p0)
y[z==0] <- 0

res <- cozigam(y~s(x1,x2)+s(x3), constraint = "proportional", family = gaussian)

plot(res)

## Poisson Response with component-specific constraint
set.seed(11)
n <- 600
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- test(x1, x2)*2 + f0(x3)/5
eta0 <- f/1.1
mu0 <- exp(eta0)  

eta.p10 <- (test(x1,x2) - mean(test(x1,x2)))*2/1.1
eta.p20 <- (f0(x3) - mean(f0(x3)))/5/1.1

alpha0 <- 0.2
delta10 <- 1.2
delta20 <- 0
eta.p0 <- delta10*eta.p10 + delta20*eta.p20 
p0 <- .Call("logit_linkinv", alpha0 + eta.p0, PACKAGE = "stats")

z <- rbinom(rep(1,n), 1, p0)
y <- rpois(rep(1,n), mu0)
y[z==0] <- 0; rm(z)

res <- cozigam(y~s(x1,x2)+s(x3), constraint="component", zero.delta=c(NA, 0),
  conv.crit.out = 1e-3, family=poisson)


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
res.tps <- cozigam(y~s(x1,x2)+s(x3), constraint = "proportional", family=gaussian)
par(mfrow=c(1,2))
plot(res.tps, select=1)
# Compare with tensor product spline
res.ten <- cozigam(y~te(x1,x2)+s(x3), constraint = "proportional", family=gaussian)
plot(res.ten, select=1)
}

\keyword{smooth} \keyword{models} \keyword{regression}


