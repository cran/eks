\name{tidyst_kcurv}
\alias{tidy_kcurv}
\alias{st_kcurv}

\title{Tidy and geospatial kernel summary density curvature estimates}

\description{
  Tidy and geospatial versions of kernel summary density curvature estimates for 2-dimensional data.
}

\usage{
tidy_kcurv(data, ...)
st_kcurv(x, ...)
}

\arguments{
  \item{data}{tidy kernel density curvature estimate (output from \code{\link{tidy_kdde}(deriv_order=2)})}
  \item{x}{geospatial density curvature estimate (output from \code{\link{st_kdde}(deriv_order=2)})}
  \item{...}{other parameters in \code{ks::kcurv} function}
}


\value{
    The output from \code{*_kcurv} have the same structure as the input kernel density curvature estimate from \code{\link[=st_kdde]{*_kdde}}, except that \code{estimate} indicates the summary curvature values rather than the density curvature values, and that \code{deriv_group} for each of the partial derivatives is collapsed into a single grouping.
}


\details{
    A kernel density summary curvature estimate is a modification of a kernel density curvature estimate where the matrix of second order partial derivative values is summarised as a scalar value. For details of the computation of the kernel density summary curvature estimate, see \code{?ks::kcurv}. The bandwidth matrix of smoothing parameters is computed as in \code{ks::kdde(deriv_order=2)}. 
}

\examples{
\donttest{## tidy kernel summary density curvature estimate
data(crabs, package="MASS")
crabs2 <- dplyr::select(crabs, FL, CW)
t1 <- tidy_kdde(crabs2, deriv_order=2)
t2 <- tidy_kcurv(t1)
gt1 <- ggplot2::ggplot(t2, ggplot2::aes(x=FL, y=CW)) 
gt1 + geom_contour_filled_ks(colour=1) + 
    ggplot2::scale_fill_brewer(palette="Oranges")

## geospatial kernel summary density curvature estimate
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
s1 <- st_kdde(hakeoides, deriv_order=2)
s2 <- st_kcurv(s1)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s2, add=TRUE)

## geom_sf plot
gs1 <- ggplot2::ggplot(s2) + ggplot2::geom_sf(data=wa, fill=NA) 
gs1 + ggplot2::geom_sf(data=st_get_contour(s2), ggplot2::aes(fill=label_percent(contlabel))) +
    colorspace::scale_fill_discrete_sequential(h1=30,c1=360,c2=30) +
    ggplot2::coord_sf(xlim=xlim, ylim=ylim) + 
    ggthemes::theme_map()}
}

\keyword{smooth}
