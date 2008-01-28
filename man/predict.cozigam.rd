\name{predict.cozigam}
\alias{predict.cozigam}
\title{Prediction from a fitted COZIGAM}
\description{
   Predict from a COZIGAM  fitted by the \code{cozigam} function, 
   given a new set of covariate values, if supplied, or 
   the original set of covariate values used for the model fit. 
   The standard errors of the point predictors are also computed. 
}
\usage{
\method{predict}{cozigam}(object, newdata, type="link", se.fit=FALSE, ...)
}
\arguments{
   \item{object}{A fitted \code{cozigam} object as produced by the \code{cozigam} function.}
   \item{newdata}{A data frame containing the values of the model covariates at which 
   predictions are required. If this is not provided then predictions corresponding to 
   the original data are returned. 
   If \code{newdata} is provided then it must contain all the variables needed for prediction.}
   \item{type}{When this has the value \code{"link"} (default) the linear predictor (possibly 
   with associated standard errors) is returned. 
   When \code{type="terms"} each component of the linear predictor is returned seperately (possibly 
   with standard errors): this excludes any intercept. 
   When \code{type="response"} predictions on the scale of the response are returned (possibly 
   with approximate standard errors).}
   \item{se.fit}{Logical. If \code{TRUE} (not default), the standard errors are returned.}
   \item{...}{Other arguments.}
}
\details{
The standard errors produced by \code{predict.cozigam()} are based on the covariance matrix of 
the parameters obtained by Louis' method in the fitted \code{gam} object.
}
\value{
If \code{se.fit} is \code{TRUE} then a three-item list is returned with items (both arrays) \code{fit}, 
\code{se.fit} containing the point predictors and associated standard error estimates and \code{p} 
containing the corresponding zero-inflation rates, 
otherwise a two-item list without the array of standard error estimates is returned. 
The dimensions of the returned arrays depends on whether \code{type} is \code{"terms"} or not: 
if it is then the array is 2-dimensional with each term in the linear predictor included separately, 
otherwise the array is 1-dimensional and contains the linear predictor/predicted values (or 
corresponding standard errors). The linear predictor returned termwise excludes the intercept.
}
\author{
Hai Liu and Kung-Sik Chan
}
\references{
Liu, H and Chan, K.S. (2008) Constrained Generalized Additive Model with Zero-Inflated Data. Technical Report, 
Department of Statisics and Actuarial Science, The Unversity of Iowa. 
\url{http://www.stat.uiowa.edu/techrep/tr388.pdf}  

Louis, T. A. (1982) Finding the Observed Information Matrix When Using EM Algorithm. J. R. Statist. Soc. B, 44, 226-233
}
\seealso{
   \code{\link{cozigam}},\code{\link{plot.cozigam}}
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

newdata <- data.frame(x1=c(0.5,0.8), x2=c(0.2,0.1), x3=c(0.3,0.7))
pred <- predict(res, newdata=newdata, se.fit=TRUE, type="response")
}
\keyword{smooth} \keyword{models} \keyword{regression}
 
