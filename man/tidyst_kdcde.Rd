\name{tidyst_kdcde}
\alias{tidy_kdcde}
\alias{st_kdcde}
\title{Tidy and geospatial kernel deconvolved density estimates}

\description{
    Tidy and geospatial versions of kernel deconvolved density estimates for 1- and 2-dimensional data.
}

\usage{
tidy_kdcde(data, ...)
st_kdcde(x,...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{...}{other parameters in \code{ks::kdcde} function}
}
  
\value{
    The output from \code{*_kdcde} have the same structure as the standard kernel density estimate from \code{\link[=st_kde]{*_kde}}. 
}

\details{
    A deconvolved kernel density estimate is a modification of the standard density estimate for data observed with error. This version is based on a weighted kernel density estimate. For details of the computation of the kernel deconvolved density estimate and the bandwidth selector procedure, see \code{?ks::kdcde}. 
}

\examples{\donttest{## tidy 2-d deconvolved density estimate
library(ggplot2)
data(air, package="ks")
air <- na.omit(air[, c("time","co2","pm10")])
air <- dplyr::filter(air, time=="20:00")
air <- dplyr::select(air, co2, pm10)
## for details on computation of Sigma.air, see ?ks::kdcde
Sigma.air <- diag(c(6705.765, 957.664)) 

t1 <- tidy_kde(air)
t2 <- tidy_kdcde(air, Sigma=Sigma.air, reg=0.00021)
t3 <- rbind(t1, t2, labels=c("Standard KDE","Deconvolved KDE"))
tb <- contour_breaks(t3, cont=seq(10,90,by=20), group=FALSE)

## deconvolved estimate is more clearly bimodal than standard KDE 
gt <- ggplot(t3, aes(x=co2, y=pm10)) 
gt + geom_contour_filled_ks(colour=1) + facet_wrap(~group)
gt + geom_contour_filled_ks(colour=1, breaks=tb) + facet_wrap(~group)} 
}

\keyword{smooth}
