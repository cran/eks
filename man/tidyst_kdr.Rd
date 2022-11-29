\name{tidyst_kdr}
\alias{tidy_kdr}
\alias{st_kdr}
\title{Tidy and geospatial kernel density ridge estimates}

\description{
  Tidy and geospatial versions of kernel density ridge estimates for 2-dimensional data.
}

\usage{
tidy_kdr(data, dTolerance, ...)
st_kdr(x, dTolerance, ...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{dTolerance}{tolerance parameter in \code{sf::st_simplify} for reducing complexity of density ridge}
    \item{...}{other parameters in \code{ks::kdr} function}
}

\value{
    The output from \code{*_kdr} have the same structure as the kernel density estimate from \code{\link[=st_kde]{*_kde}}, except that \code{x,y} indicate the points on the density ridge, rather than the grid points themselves, and \code{estimate} becomes NA. For \code{st_kdr}, the density ridge is stored as a multipoints \code{sf} object.
}

\details{
    A density ridge can be interpreted as the line connecting the peaks in the kernel density estimate, like for a mountain range. It can also be interpreted as the filament generalisation of 2-d principal components.  For details of the computation and the bandwidth selection procedure of the kernel density ridge estimate, see \code{?ks::kdr}. The bandwidth matrix of smoothing parameters is computed as in \code{ks::kdde(deriv_order=2)}.

    To reduce the complexity of the density ridge, a call to \code{sf::st_simplify(,dTolerance)} is made. If \code{dTolerance} is missing, then it defaults to approximately the mean distance between each pair of consecutive points in each segment of the density ridge. If \code{dTolerance=0} then this step of Ramer-Douglas-Peucker simplification is not carried out.
}

\examples{
## tidy density ridge estimate
\donttest{data(cardio, package="ks")
cardio <- dplyr::as_tibble(cardio[,c("ASTV","Mean")])
set.seed(8192)
cardio <- cardio[sample(1:nrow(cardio), round(nrow(cardio)/4,0)),]
## gridsize=c(21,21) is for illustrative purposes only
## remove for more complete KDR
t1 <- tidy_kdr(cardio, gridsize=c(21,21))
gt <- ggplot2::ggplot(t1, ggplot2::aes(x=ASTV, y=Mean)) 
gt + geom_point_ks(colour=3, alpha=0.8) + 
    ggplot2::geom_path(ggplot2::aes(colour=label, group=segment), size=1.2, alpha=0.8) +
    ggplot2::scale_colour_manual(values=6)}

## geospatial density ridge estimate
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
## gridsize=c(21,21) is for illustrative purposes only 
## remove for more complete KDR
s1 <- st_kdr(hakeoides, gridsize=c(21,21))

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(sf::st_geometry(hakeoides), add=TRUE, col=3, pch=16)
plot(s1, add=TRUE, col=6, lwd=3, alpha=0.8)

## geom_sf plot
gs <- ggplot2::ggplot(s1) + ggplot2::geom_sf(data=wa, fill=NA) +
    ggthemes::theme_map()
gs + ggplot2::geom_sf(data=hakeoides, colour=3, alpha=0.5) +
    ggplot2::geom_sf(data=s1$sf, ggplot2::aes(colour=label), linewidth=1.2, alpha=0.8) +
    ggplot2::scale_colour_manual(values=6) +
    ggplot2::coord_sf(xlim=xlim, ylim=ylim) 
}

\keyword{smooth}
