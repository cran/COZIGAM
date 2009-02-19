\name{zigam}
\alias{zigam}
\title{Fitting (Unconstrained) Zero-Inflated Generalized Additive Models}
\description{
   Fit a Zero-Inflated Generalized Additive Model (ZIGAM) to data.
}
\usage{
zigam (formula, maxiter = 20, conv.crit = 1e-3,
    size = NULL, log.tran = FALSE, family, data=list(), ...)
}
\arguments{
   \item{formula}{A GAM formula. This is exactly like the formula for a GLM except that smooth terms 
   can be added to the right hand side of the formula (and a formula of the form \code{y ~ .} is not allowed). 
   Smooth terms are specified by expressions of the form: \code{s(var1,var2,...,k=12,fx=FALSE,bs="tp",by=a.var)} 
   where \code{var1}, \code{var2}, etc. are the covariates which the smooth is a function of and 
   \code{k} is the dimension of the basis used to represent the smooth term. \code{by} can be used to 
   specify a variable by which the smooth should be multiplied.}
   \item{maxiter}{The maximum number of iterations allowed in the EM algorithm in estimating discrete ZIGAMs.}
   \item{conv.crit}{Convergence criterion in the iterative estimation algorithm.}
   \item{size}{Number of trials. Must be specified when \code{family} is \code{binomial}.}
   \item{log.tran}{Logical. \code{TRUE} if log-transformation is needed for the response.}
   \item{family}{This is a family object specifying the distribution and link to use in fitting etc. 
   See \code{glm} and \code{family} for more details. Currently support Gaussian/lognormal, Gamma, 
   Poisson and binomial distributions.}
   \item{data}{A data frame or list containing the model response variable and covariates required by the formula.}
   \item{...}{Additional arguments to be passed to the low level regression fitting functions.}
}

\details{
A Zero-Inflated Generalized Additive Model (ZIGAM) assumes the response variable \eqn{Y_i} 
is distributed from a zero-inflated 1-parameter exponential family with covariate \eqn{x_i} (could be high dimensional). 
More specifically, \eqn{Y_i} comes from a non-zero-inflated exponential family distribution \eqn{f(x_i)} 
(regular component) with probability \eqn{p_i} and equals zero with probability \eqn{1-p_i}. 
The probability of non-zero-inflation \eqn{p_i} also depends on the covariates through some unknown smooth funtions.
Different from the COnstrained Zero-Inflated Generalized Additive Model (COZIGAM), the process of generating the 
non-zero-inflated responses and the zero-inflation process are assumed to be independent.
The mean of the non-zero-inflated exponential family distribution is assumed to be
 \deqn{g(\mu_i) = s_1(x_i),} and the non-zero-inflation probability is linked to the covariates by
 \deqn{logit(p_i) = s_2(x_i),} where \eqn{s_1} and \eqn{s_2} are two possibly distinct smooth functions to
be estimated nonparametrically by Generalized Additive Models. See Liu and Chan (2008) for more detail.
}

\value{
A list containing the following components: 
   \item{fit.gam}{A fitted GAM of the regular component, i.e., non-zero-inflated exponential family regression model.}
   \item{fit.lr}{A logistic regression model on the zero-inflation process.}
   \item{V.beta}{The covariance matrix of the estimated parameters associated with the smooth functions in the 
       non-zero-inflated data generating process.}
   \item{V.gamma}{The covariance matrix of the estimated parameters associated with the smooth functions in the 
       zero-inflation process.}
   \item{mu}{The fitted regular mean values.}
   \item{dispersion}{(Estimated) dispersion parameter.}
   \item{formula}{Model formula.}
   \item{p}{The fitted non-zero-inflation probabilities.}
   \item{psi}{Conditional expectation of the zero-inflation indicator (only for the discrete case which involves EM algorithm).}
   \item{family}{The family used.}
   \item{loglik, ploglik}{The (penalized) log-likelihood of the fitted model.}
   \item{logE}{Approximated logarithmic marginal likelihood by Laplace method used for model selection.} 
   \item{X1}{Design matrix in the non-zero-inflated exponential family regression model.}
   \item{X2}{Design matrix in the logistic regression model.}
}

\author{
Hai Liu and Kung-Sik Chan
}

\references{
Liu, H and Chan, K.S. (2008) Constrained Generalized Additive Model with Zero-Inflated Data. 
Technical Report 388, Department of Statistics and Actuarial Science, The University of Iowa.
\url{http://www.stat.uiowa.edu/techrep/tr388.pdf}}
\seealso{
   \code{\link{cozigam}}
}
\examples{
## Gaussian Response 
set.seed(11)
n <- 200
x1 <- runif(n, 0, 1)
f <- (f0(x1)-mean(f0(x1)))/2
sig <- 0.3
mu0 <- f + 1.5
y <- mu0 + rnorm(n, 0, sig)

eta.p0 <- f1(x1)*2 - 1 # true function used in zero-infation process
p0 <- .Call("logit_linkinv", eta.p0, PACKAGE = "stats")

z <- rbinom(rep(1,n), 1, p0)
y[z==0] <- 0

# Fit a ZIGAM
res.un <- zigam(y~s(x1), family=gaussian)

# Compare with a COZIGAM
res <- cozigam(y~s(x1), family=gaussian)
res.un$logE > res$logE

## Poisson Response
set.seed(11)
n <- 400
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)

eta0 <- test(x1,x2)*3
mu0 <- exp(eta0)  

alpha0 <- -0.2
delta0 <- 1.0
p0 <- .Call("logit_linkinv", alpha0 + delta0 * eta0, PACKAGE = "stats")


z <- rbinom(rep(1,n), 1, p0)
y <- rpois(rep(1,n), mu0)
y[z==0] <- 0

res.un <- zigam(y~s(x1,x2), maxiter=30, family=poisson) # fit a ZIGAM
res <- cozigam(y~s(x1,x2), maxiter=30, family=poisson) # fit a COZIGAM
res.un$logE < res$logE # compare the model selction criterion

}

\keyword{smooth} \keyword{models} \keyword{regression}


