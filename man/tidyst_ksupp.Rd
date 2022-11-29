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
    The output from \code{*_ksupp} have the same structure as the kernel density estimate from \code{\link[=st_kde]{*_kde}}, except that \code{x,y} indicate the boundary of the density support estimate (if \code{convex.hull=TRUE}) or the grid points inside the density support (if \code{convex.hull=FALSE}), rather than the complete grid points themselves.

    For \code{st_kdr}, the density support estimate is stored as a (multi)polygon \code{sf} object.
}

\details{
    The kernel support estimate is considered to be the \code{cont}\% probability contour of the kernel density estimate, with an additional convex hull calculation if \code{convex_hull=TRUE}. For details of the computation of the kernel support estimate, see \code{?ks::ksupp}.
}

\examples{
## tidy density support estimate
data(crabs, package="MASS")
crabs2 <- dplyr::select(crabs, FL, CW)
t1 <- tidy_kde(crabs2)
t2 <- tidy_ksupp(t1)
gt <- ggplot2::ggplot(t1, ggplot2::aes(x=FL, y=CW)) 
gt + geom_contour_filled_ks(cont=c(25,50,75,95), colour=1) +
    ggplot2::geom_polygon(data=t2, ggplot2::aes(linetype=label), fill=NA, colour=1,
    linewidth=1.2) + colorspace::scale_fill_discrete_sequential() +
    ggplot2::scale_linetype_manual(values="dashed", name="Support\nconvex hull")

gt + ggplot2::geom_polygon(data=t2, ggplot2::aes(colour=label), fill=NA,
    linetype="dashed") + geom_point_ks(data=t1, colour=3) +
    ggplot2::scale_colour_manual(values=1)

## geospatial density support estimate
\donttest{data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
s1 <- st_kde(hakeoides)
s2 <- st_ksupp(s1)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, cont=c(25,50,75,95), add=TRUE)
plot(sf::st_geometry(hakeoides), add=TRUE, pch=16, col=1)
plot(s2, add=TRUE, lty=2, lwd=2)

## geom_sf plot
gs1 <- ggplot2::ggplot(s1) + ggplot2::geom_sf(data=wa, fill=NA)  + 
    ggthemes::theme_map()
gs1 + ggplot2::geom_sf(data=st_get_contour(s1, cont=c(25,50,75,95)),
    ggplot2::aes(fill=label_percent(contlabel))) +
    ggplot2::geom_sf(data=st_get_contour(s2,cont=95), 
    ggplot2::aes(linetype=contlabel), fill=NA, colour=1) + 
    ggplot2::geom_sf(data=hakeoides, ggplot2::aes(colour=species)) +
    colorspace::scale_fill_discrete_sequential(palette="Heat2") +
    ggplot2::scale_colour_manual(values=1) + 
    ggplot2::guides(colour=ggplot2::guide_legend(title="Locations")) +
    ggplot2::scale_linetype_manual(values="dashed", name="Support\nconvex hull") +
    ggplot2::coord_sf(xlim=xlim, ylim=ylim)}
}

\keyword{smooth}
