\name{plot.cozigam}
\alias{plot.cozigam}
\title{Default COZIGAM plotting}
\description{
   Takes a fitted \code{cozigam} object produced by \code{cozigam()} and plots 
   the component smooth functions that make it up, on the scale of the linear predictor.
}
\usage{
\method{plot}{cozigam}(x, plot.2d = "contour", too.far = 0.05, 
    n.1d = 100, n.2d = 30, theta = 30, phi = 30, select = NULL, image.col = "topo", 
    persp.col = "lightblue", contour.col = "red", n.Col = 100, shade.ci = FALSE,
    shade.col = "gray80", Rug = TRUE, ...)
}
\arguments{
   \item{x}{a fitted \code{cozigam} object produced by \code{cozigam()}.}
   \item{plot.2d}{one of "contour" (default) or "persp".}
   \item{select}{allows the plot for a single model term to be selected for printing. 
   e.g. if you just want the plot for the second smooth term set \code{select=2.}}
   \item{n.1d}{number of points used for each 1-D plot. Default value 100.}
   \item{n.2d}{square root of number of points used to grid estimates of 2-D functions 
   for contouring.}
   \item{theta}{one of the perspective plot angles.}
   \item{phi}{the other perspective plot angle.}
   \item{too.far}{if greater than 0 then this is used to determine when a location 
   is too far from data to be plotted when plotting 2-D smooths. This is useful 
   since smooths tend to go wild away from data. The data are scaled into the unit 
   square before deciding what to exclude, and \code{too.far} is a distance within 
   the unit square.}
   \item{shade.ci}{logical. If \code{TRUE}, produce shaded regions as confidence bands for smooths.}
   \item{shade.col}{define the color used for shading confidence bands.}
   \item{image.col}{define the color used for 2-D image plots.}
   \item{persp.col}{define the color used for 2-D perspective plots.}
   \item{contour.col}{define the color used for the 2-D contour lines.}
   \item{n.Col}{control the number of colors in 2-D image plots.}
   \item{Rug}{logical, if \code{TRUE} (default) then the covariate to which the 
   plot applies is displayed as a rug plot at the foot of each plot of a 1-D smooth, 
   and the locations of the covariates are plotted as points on the contour plot
   representing a 2-D smooth.}
   \item{...}{other graphics parameters to pass on to plotting commands.}
}
\details{
Produces default plot showing the smooth components of a fitted \code{cozigam}.

Smooths of more than 2 variables are not currently dealt with, 
but simply generate a warning.
}
\value{
The function simply generates plots.
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

f <- test(x1, x2)*4-mean(test(x1, x2)*4) + f0(x3)/2-mean(f0(x3)/2)
sig <- 0.5
mu0 <- f + 3
y <- mu0 + rnorm(n, 0, sig)

alpha0 <- -2.2
delta0 <- 1.2
p0 <- .Call("logit_linkinv", alpha0 + delta0 * mu0, PACKAGE = "stats")
z <- rbinom(rep(1,n), 1, p0)
y[z==0] <- 0

res <- cozigam(y~s(x1,x2)+s(x3), conv.crit.out = 1e-4, family = gaussian)

plot(res, plot.2d = "contour", too.far = 0.1, image.col="topo")
plot(res, plot.2d = "persp", too.far = 0.1, select=1)
}
\keyword{hplot} \keyword{models} \keyword{regression} \keyword{smooth}
