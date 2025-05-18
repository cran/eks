\name{tidyst_ksupp}
\alias{tidy_ksupp}
\alias{st_ksupp}
\title{Tidy and geospatial kernel support estimate}

\description{
  Tidy and geospatial versions of a kernel support estimate for 2-dimensional data.
}

\usage{
tidy_ksupp(data, cont=95, convex_hull=TRUE, ...)
st_ksupp(x, cont=95, convex_hull=TRUE, ...)
}

\arguments{
    \item{data}{tidy kernel density estimate (output from \code{\link{tidy_kde}})}
    \item{x}{spatial kernel density estimate (output from \code{\link{st_kde}})}
    \item{cont}{scalar contour level. Default is 95.}
    \item{convex_hull}{flag to compute convex hull of contour region. Default is TRUE.}
    \item{...}{other parameters in \code{ks::ksupp} function}
}

\value{
    The output from \code{*_ksupp} have the same structure as the kernel density estimate from \code{\link[=st_kde]{*_kde}}, except that \code{x,y} indicate the boundary of the density support estimate (if \code{convex_hull=TRUE}) or the grid points inside the density support (if \code{convex_hull=FALSE}), rather than the complete grid points themselves.

    For \code{st_ksupp}, the density support estimate is stored as a (multi)polygon \code{sf} object.
}

\details{
    The kernel support estimate is considered to be the \code{cont}\% probability contour of the kernel density estimate, with an additional convex hull calculation if \code{convex_hull=TRUE}. For details of the computation of the kernel support estimate, see \code{?ks::ksupp}.
}

\examples{
## tidy density support estimate
library(ggplot2)
data(crabs, package="MASS")
crabs2 <- dplyr::select(crabs, FL, CW)
t1 <- tidy_kde(crabs2)
t2 <- tidy_ksupp(t1)
ggplot(t1, aes(x=FL, y=CW)) + 
    geom_contour_filled_ks(cont=c(25,50,75,95), colour="grey50") +
    geom_polygon(data=t2, aes(linetype=label), fill=NA, colour=1) +
    scale_linetype_manual(values="dashed", name="Support\nconvex hull")

## geospatial density support estimate
\donttest{data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
s1 <- st_kde(hakeoides)
s2 <- st_ksupp(s1, cont=97.5)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, cont=c(25,50,75,97.5), add=TRUE, border="grey50")
plot(s2, add=TRUE, lty=2, pos="bottom")

## geom_sf plot
ggplot(s1) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map() + 
    scale_linetype_manual(values=c("dashed"), name="Support\nconvex hull") + 
    geom_sf(data=st_get_contour(s1, cont=c(25,50,75,97.5)), aes(fill=contperc)) +
    geom_sf(data=st_get_contour(s2), aes(linetype=contperc), fill=NA, col=1) + 
    coord_sf(xlim=xlim, ylim=ylim)}
}

\keyword{smooth}
