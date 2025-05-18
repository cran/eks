\name{tidyst_kms}
\alias{tidy_kms}
\alias{st_kms}
\title{Tidy and geospatial kernel mean shift clustering}

\description{
    Tidy and geospatial versions of a kernel mean shift clustering for 1- and 2-dimensional data.
}

\usage{
tidy_kms(data, ...)
st_kms(x, ...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{...}{other parameters in \code{ks::kms} function}
}
  
\value{
    The output from \code{*_kms} have the same structure as the kernel density estimate from \code{\link[=st_kde]{*_kde}}, except that \code{x,y} indicate the data points rather than the grid points, and \code{estimate} indicates the mean shift cluster label of the data points, rather than the density values. 
}

\details{
    Mean shift clustering is a generalisation of \eqn{k}{k}-means clustering (aka unsupervised learning) which allows for non-ellipsoidal clusters and does not require the number of clusters to be pre-specified. The mean shift clusters are determined by following the initial observations along the density gradient ascent paths to the cluster centre.

    For details of the computation and the bandwidth selection procedure of the kernel mean shift clustering, see \code{?ks::kms}. The bandwidth matrix of smoothing parameters is computed as in \code{ks::kdde(deriv_order=1)}. 
}

\examples{
## tidy 2-d mean shift clustering 
library(ggplot2) 
data(crabs, package="MASS")
crabs2 <- dplyr::select(crabs, FL, CW)
t1 <- tidy_kms(crabs2)
## convex hulls of clusters
t2 <- dplyr::group_by(t1, label)
t2 <- dplyr::slice(t2, chull(FL,CW))

gt <- ggplot(t1, aes(x=FL, y=CW)) 
gt + geom_point(aes(colour=label)) +
    geom_polygon(data=t2, aes(fill=label), alpha=0.1, col=1, linetype="dotted")

## geospatial mean shift clustering 
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
s1 <- st_kms(hakeoides)
## convex hulls of clusters
s2 <- dplyr::group_by(s1$sf, label)
s2 <- dplyr::summarise(s2, geometry=sf::st_combine(geometry))
s2 <- sf::st_convex_hull(s2)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, add=TRUE, pch=16)
plot(s2, add=TRUE, lty=3, pal=function(.){
    colorspace::qualitative_hcl(n=., palette="Set2", alpha=0.15)})

## geom_sf plot
gs <- ggplot(s1) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map()
gs + geom_sf(data=s1$sf, aes(colour=label), alpha=0.5) + 
    geom_sf(data=s2, aes(fill=label), linetype="dotted", alpha=0.15) + 
    coord_sf(xlim=xlim, ylim=ylim) 
}


\keyword{smooth}
