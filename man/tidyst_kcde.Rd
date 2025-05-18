\name{tidyst_kcde}
\alias{tidy_kcde}
\alias{tidy_kcopula}
\alias{st_kcde}
\title{Tidy and geospatial kernel cumulative distribution and copula estimates}

\description{
    Tidy and geospatial versions of kernel cumulative distribution estimates for 1- and 2-dimensional data, and kernel copula estimates for 2-dimensional data.
}

\usage{
tidy_kcde(data, ...)
tidy_kcopula(data, ...)
st_kcde(x, ...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{...}{other parameters in \code{ks::kcde}, \code{ks::kcopula} functions}
}

\value{
    The outputs from \code{*_kcde} have the same structure as the kernel density estimate from \code{\link[=st_kde]{*_kde}}, except that \code{estimate} indicates the cumulative distribution rather than the density values. Likewise for \code{tidy_kcopula}.
}

\details{
 	For details of the computation of the kernel distribution and copula estimates, and of the bandwidth selector procedures, see \code{?ks::kcde}, \code{?ks::kcopula}.
}

\examples{library(ggplot2)
data(crabs, package="MASS")
## tidy 1-d distribution estimate 
crabs1 <- dplyr::select(crabs, FL)
t1 <- tidy_kcde(crabs1)
gt2 <- ggplot(t1, aes(x=FL)) 
gt2 + geom_line(colour=1) + geom_rug_ks(colour=4) 

## tidy 2-d copula estimate
crabs2 <- dplyr::select(crabs, FL, CW)
t2 <- tidy_kcopula(crabs2)
gt2 <- ggplot(t2, aes(x=FL, y=CW)) 
gt2 + geom_contour_filled_ks(colour=1, cont=seq(10,90,by=10))

## geospatial distribution estimate
\donttest{data(wa)
data(grevilleasf)
paradoxa <- dplyr::filter(grevilleasf, species=="paradoxa")
paradoxa_bbox <- sf::st_bbox(c(xmin=4e5, xmax=8e5, ymin=6.4e6, ymax=6.65e6), 
    crs=sf::st_crs(paradoxa))
xminb <- paradoxa_bbox[1:2]; xmaxb <- paradoxa_bbox[3:4]
paradoxa_bbox <- sf::st_as_sfc(paradoxa_bbox)
paradoxa_crop <- sf::st_filter(paradoxa, paradoxa_bbox)
s1 <- st_kcde(paradoxa_crop, xmin=xminb, xmax=xmaxb)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(paradoxa_bbox, add=TRUE, lty=3, lwd=2)
plot(s1, add=TRUE)

## geom_sf plot
gs1 <- ggplot(s1) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map() 
gs1 + geom_sf(data=paradoxa_bbox, linewidth=1.2, linetype="dotted", fill=NA) +
    geom_sf(data=st_get_contour(s1), aes(fill=contperc)) +
    coord_sf(xlim=xlim, ylim=ylim)}
}

\keyword{smooth}
