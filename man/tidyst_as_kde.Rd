\name{tidyst_as_kde}
\alias{tidy_as_kde}
\alias{st_as_kde}
\title{Tidy and geospatial quasi kernel density estimates from gridded data}

\description{
    Tidy and geospatial versions of quasi kernel density estimates for 1- and 2-dimensional gridded data.
}

\usage{
tidy_as_kde(data, density, ...)
st_as_kde(x, attrib=1, density, ...)
}

\arguments{
    \item{data}{data frame/tibble, where first d columns are evaluation points of an estimation grid, and last column is estimate value}
    \item{x}{sf object, where \code{attrib} is estimate value and estimation grid is rectangle polygon geometry}
    \item{attrib}{name or position of estimate variable. Default is 1.}
    \item{density}{flag to compute bona fide density. Default is TRUE.}
    \item{...}{other parameters}
}

\value{
    The input gridded data should be a tidy form of a complete Cartesian product that forms a regular rectangular grid with no missing values or grid points/polygons. The input is converted to a quasi density estimate so that it can be subsequently treated like an output from \code{tidy_kde} or \code{st_kde}.  
}

\details{
    If \code{density=TRUE} then all the estimate values are non-negative. If \code{density=FALSE} then the estimate values can be positive or negative.
}

\examples{
## tidy quasi density estimate
library(ggplot2)
data(faithfuld, package="ggplot2")
t1 <- tidy_as_kde(faithfuld)
## probabilistic density contour levels
ggplot(t1, aes(x=eruptions, y=waiting)) + 
    geom_contour_filled_ks(colour=1, aes(fill=after_stat(contperc)))

## non-probabilistic contour levels in ggplot2::geom_contour_filled
ggplot(faithfuld, aes(x=eruptions, y=waiting, z=density)) +
    geom_contour_filled(bins=4) + 
    colorspace::scale_fill_discrete_sequential(palette="Heat")

## geospatial quasi density estimate
## see example in ? st_intergrid 
}

\keyword{smooth}
