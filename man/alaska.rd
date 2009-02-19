\name{alaska}
\docType{data}
\alias{alaska}
\title{Alaska Peninsula Coastline}
\description{
The coastline of the Alaska Peninsula.
}
\author{
Hai Liu and Kung-Sik Chan
}
\usage{data(alaska)}
\format{A data frame with the following 2 variables.
   \item{lat}{A numeric vector: latitude}
   \item{lon}{A numeric vector: longitude}
}
\details{
This data frame can be used to draw the coastline of the Alaska peninsula. See example below.
}
\examples{
data(alaska)
plot(alaska$lon,alaska$lat,type="l",col="green")
}
\keyword{datasets}