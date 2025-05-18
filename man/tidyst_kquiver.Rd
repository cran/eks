\name{tidyst_kquiver}
\alias{tidy_kquiver}
\alias{st_kquiver}
\title{Tidy and geospatial kernel density quiver estimate}

\description{
  Tidy and geospatial versions of a kernel density quiver estimate for 2-dimensional data.
}

\usage{
tidy_kquiver(data, thin=5, transf=1/4, neg.grad=FALSE)
st_kquiver(x, thin=5, transf=1/4, neg.grad=FALSE, scale=1)
}

\arguments{
  \item{data}{tidy kernel density gradient estimate (output from \code{\link{tidy_kdde}(deriv_order=1)})}
  \item{x}{geospatial kernel density gradient estimate (output from \code{\link{st_kdde}(deriv_order=1)})}
  \item{thin}{number to thin out estimation grid. Default is 5.}
  \item{transf}{power index in transformation. Default is 1/4.} 
  \item{neg.grad}{flag to compute arrows in negative gradient direction. Default is FALSE.}
  \item{scale}{scale factor to normalise length of arrows. Default is 1.}
}
  
\value{
    The output from \code{tidy_kquiver} has the same structure as the input kernel density gradient estimate, with the added columns \code{u,v} for the density gradient value in the \eqn{x}{x}-, \eqn{y}{y}-axis. This structure is compatible with the \code{ggquiver::geom_quiver} layer function for quiver plots. 

    Since \code{ggquiver::geom_quiver} is not compatible with \code{geom_sf} layers, the output from \code{st_kquiver} has added columns \code{lon}, \code{lat}, \code{lon_end}, \code{lat_end}, \code{len} which are required in \code{geom_segment}. 
}
 
\details{
    A kernel quiver estimate is a modification of the standard kernel density gradient estimate in \code{\link[=st_kdde]{*_kdde}} where the density derivatives are not given in the separate groups as indexed in \code{deriv_group}, but as extra columns \code{u} (for \code{deriv_group=(1,0)}) and \code{v} (for \code{deriv_group=(0,1)}). 

    The bandwidth matrix of smoothing parameters is computed as in \code{ks::kdde(deriv_order=1)}.
}

\examples{
## tidy kernel quiver estimate
\donttest{library(ggplot2)
data(crabs, package="MASS")
crabs2 <- dplyr::select(crabs, FL, CW)
t1 <- tidy_kde(crabs2)
t2 <- tidy_kdde(crabs2, deriv_order=1)
t3 <- tidy_kquiver(t2, thin=5) 
gt <- ggplot(t1, aes(x=FL, y=CW)) 
gt + geom_contour_filled_ks(alpha=0.5) + 
    ggquiver::geom_quiver(data=t3, aes(u=u, v=v), colour=1)

## geospatial kernel `quiver' estimate
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
s1 <- st_kde(hakeoides)
s2 <- st_kdde(hakeoides, deriv_order=1)
s3 <- st_kquiver(s2, thin=5)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, add=TRUE, alpha=0.5, border="grey50")
plot(s3$tidy_ks$ks[[1]], thin=5, add=TRUE, display="quiver")

## geom_sf plot 
gs <- ggplot(s1) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map()
gs + geom_sf(data=st_get_contour(s1), aes(fill=contperc), alpha=0.5) +
   ggquiver::geom_quiver(data=s3$sf, aes(x=lon, y=lat, u=u, v=v), colour=1) + 
    coord_sf(xlim=xlim, ylim=ylim)}
}

\keyword{smooth}
