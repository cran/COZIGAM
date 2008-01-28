\name{cozigam}
\alias{cozigam}
\title{Fit a Constrained Zero-Inflated Generalized Additive Model}
\description{
   Fit a COnstrained Zero-Inflated Generalized Additive Model (COZIGAM).
}
\usage{
cozigam (formula, maxiter = 30, conv.crit.in = 1e-5, 
    conv.crit.out = 1e-4, size = NULL, log.tran = FALSE, family, ...)
}
\arguments{
   \item{formula}{A GAM formula. This is exactly like the formula for a GLM except that smooth terms 
   can be added to the right hand side of the formula (and a formula of the form \code{y~.} is not allowed); 
   c.f. the formula argument of the \code{gam} function of the \code{mgcv} library.  
   Smooth terms are specified by expressions of the form: \code{s(var1,var2,...,k=12,fx=FALSE,bs="tp",by=a.var)} 
   where \code{var1}, \code{var2}, etc. are the covariates of which the smooth is a function and 
   \code{k} is the dimension of the basis used to represent the smooth term. The argument \code{by} can be used to 
   specify a variable by which the smooth should be multiplied.}
   \item{maxiter}{The maximum number of iterations allowed in the estimation procedure.}
   \item{conv.crit.in}{The convergence criterion in the inner loop.}
   \item{conv.crit.out}{The convergence criterion in the outer loop.}
   \item{size}{Optional. Must be specified when \code{family} is specified as \code{binomial}.}
   \item{log.tran}{Logical. \code{TRUE} if log-transformation is needed for the response.}
   \item{family}{This is a family object for specifying the distribution and link to use in model fitting. 
   See \code{glm} and \code{family} for more details. Currently supported families include Gaussian/lognormal, Gamma, 
   Poisson and binomial distributions.}
   \item{...}{Additional arguments to be passed to the lower level regression fitting functions.}
}

\details{
A COnstrained Zero-Inflated Generalized Additive Model (COZIGAM) assumes that the response variable \eqn{Y_i} 
is distributed from a zero-inflated 1-parameter exponential family with covariate \eqn{x_i} (could be high dimensional). 
More specifically, \eqn{Y_i} comes from a regular exponential family distribution \eqn{f(x_i)} (for example, a Poisson distribution) with 
probability \eqn{p_i} and equals zero with probability \eqn{1-p_i}, with the further assumption that the 
probability of zero-inflation \eqn{p_i} is some monotone function of the mean response function on the link scale, 
e.g., if a logit link is used, \eqn{p_i} is linearly constrained by: \deqn{logit(p_i) = \alpha + \delta g(\mu_i),}
where \eqn{g()} is the link function of the assumed exponential family and \eqn{\alpha} and \eqn{\delta} are 
two unknown parameters to be estimated. 
The regular mean response (for example, the Poisson mean) is specified as \eqn{g(\mu_i) = s(x_i)} for some smooth function of the covariate \eqn{x}. 

The estimation approach is a modified penalized iteratively reweighted least squares algorithm (P-IRLS) which requires 
the \code{magic} function in the \code{mgcv} library. The EM algorithm is used if the underlying 
regular exponential family has strictly positive probability mass at zero. 
The covariance matrix of the estimated coefficients is obtained by inverting the 
observed information from Louis' method. See Liu and Chan (2008) for more detail.
}

\value{
An object of class \code{"cozigam"} is a list containing the following components: 
   \item{coefficients}{A named vector of coefficients including the linear constraint parameters.}
   \item{V.theta}{The covariance matrix of the estimated parameters.}
   \item{beta, V.beta}{The estimated parameters of smooth terms with associated covariance matrix V.beta.}
   \item{mu}{The fitted mean values.}
   \item{linear.predictor}{The fitted values on the link scale.}
   \item{dispersion}{(Estimated) dispersion parameter.}
   \item{formula}{Model formula.}
   \item{p}{The fitted zero-inflation rates.}
   \item{G}{An objected returned by \code{gam()} containing information used in the estimation procedure.}
   \item{psi}{Conditional expectation of the zero-inflation indicator (only for the discrete case which involves EM algorithm).}
   \item{family}{The family used.}
   \item{loglik, ploglik}{The (penalized) log-likelihood of the fitted model.}
   \item{y}{The response used.}
   \item{converged}{Logical. \code{TRUE} if the iterative procedure has converged.}
   \item{est.disp}{Logical. \code{TRUE} if the dispersion parameter is estimated.}
   \item{fit.nonzero}{Fitted GAM using only the nonzero data as in the presence/absence analysis}
   \item{score}{The GCV or UBRE score returned by \code{magic()} in the iterative procedure.}
   \item{edf}{Estimated degrees of freedom for each model parameter. 
   Because of the use of the penalized likelihood, the \code{edf} may be less than 1.}
   \item{edf.smooth}{Estimated degrees of freedom of each smooth term.} 
}

\author{
Hai Liu and Kung-Sik Chan
}

\references{
Liu, H and Chan, K.S. (2008) Constrained Generalized Additive Model with Zero-Inflated Data. Technical Report, 
Department of Statisics and Actuarial Science, The Unversity of Iowa. 
\url{http://www.stat.uiowa.edu/techrep/tr388.pdf} 
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

f <- testfn(x1, x2)*4-mean(testfn(x1, x2)*4) + f0(x3)/2-mean(f0(x3)/2)
sig <- 0.5
mu0 <- f + 3
y <- mu0 + rnorm(n, 0, sig)

alpha0 <- -2.2
delta0 <- 1.2
p0 <- .Call("logit_linkinv", alpha0 + delta0 * mu0, PACKAGE = "stats")
z <- rbinom(rep(1,n), 1, p0)
y[z==0] <- 0

res <- cozigam(y~s(x1,x2)+s(x3), conv.crit.out = 1e-4, log.tran = FALSE, family = gaussian)

plot(res)

## Poisson Response

n <- 500
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)

f <- testfn(x1, x2)*2 + f0(x3)/5
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
plot(res.tps, select=1)
# Compare with tensor product spline
res.ten <- cozigam(y~te(x1,x2)+s(x3), conv.crit.out = 1e-4, family=gaussian)
plot(res.ten, select=1)


## A Real Application of Pollock Egg Density Estimation
data(eggdata) # load Pollock egg dataset
y <- eggdata$catch # response, egg catch
lon <- eggdata$lon; lat <- eggdata$lat # longitude & latitude
day <- eggdata$j.day # Julian day
bot <- log(eggdata$bottom) # log-transformed bottom depth
# Fit COZIGAM
egg.res <- cozigam(y~s(lon,lat)+s(day)+s(bot), log.tran=TRUE, family=gaussian)
summary(egg.res)
plot(egg.res, too.far=0.05)
# Compare with the results from the presence/absense analysis
pa <- as.numeric(y>0)
egg.bin <- gam(pa~s(lon,lat)+s(day)+s(bot), family=binomial)
plot(egg.bin)

}

\keyword{smooth} \keyword{models} \keyword{regression}


