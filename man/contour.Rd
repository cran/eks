\name{contour}
\alias{contour_breaks}
\alias{contourLevels.tidy_ks}
\alias{contourLevels.sf_ks}
\alias{st_get_contour}
\title{Contour functions for tidy and geospatial kernel estimates}

\description{
    Contour functions for tidy and geospatial kernel estimates.
}

\usage{
\method{contourLevels}{tidy_ks}(x, ..., cont=c(25,50,75))
\method{contourLevels}{sf_ks}(x, ..., cont=c(25,50,75))
contour_breaks(data, cont=c(25,50,75), n=3, group=FALSE, type="density")
st_get_contour(x, cont=c(25,50,75), breaks, disjoint=TRUE, digits)
}

\arguments{
    \item{x,data}{tidy kernel estimate (output from \code{tidy_k*}) or geospatial kernel estimate (output from \code{st_k*})}
    \item{cont}{vector of contour levels. Default is c(25,50,75). Used only for \code{level="density", "quantile"}.}
    \item{n}{number of contour levels. Default is 3. Used only for \code{level="length", "natural"}.}
    \item{type}{type of contour levels: one of \code{"density", "length", "quantile", "natural"}. Default is "density".}
    \item{group}{flag to compute contour levels per group. Default is FALSE.}
    \item{breaks}{tibble or vector of contour levels (e.g. output from \code{contour_breaks})}
    \item{disjoint}{flag to compute disjoint contours. Default is TRUE.}
    \item{digits}{number of significant digits in output labels. If missing, default is 4.} 
    \item{...}{other parameters (not implemented)}    
}

\value{
    --The output from \code{contour_breaks} is a tibble of the values of the contour breaks. If \code{group=FALSE} then a single set of contour breaks is returned. If \code{group=TRUE} then a set of contour breaks for each group is returned. 

    --The output from \code{st_get_contour} is an sf object of the contours as multipolygons, with added attributes
    \item{estimate}{density value or factor with \code{digits} s.f.}
    \item{contlabel}{factor of contour probability level}
    \item{contperc}{factor of contour probability level with added \%}
    \item{contregion}{factor of contour probability level with added \eqn{\leq}{<=} or \eqn{\geq}{>=}}
    These are used to render the contours in the base \proglang{R} base graphics, and via the aesthetics in \code{ggplot}. 
}

\details{
    By default, the 5\%, 10\%, ..., 95\% contours are computed for an \code{st_k*} output, though a plot of 19 of them would be too crowded. \code{st_get_contour} selects a subset of these, as specified by \code{cont}. If a contour level in \code{cont} does not already exist or if absolute contour levels are specified in \code{breaks}, then the corresponding contours are computed. If \code{disjoint=TRUE} (default) then the contours are computed as a set of disjoint multipolygons: this allows for plotting without overlapping transparent colours. If \code{disjoint=FALSE} then the contours are overlapping and so their colours alpha-mixed, but they strictly satisfy the probabilistic definition, e.g. a 25\% contour region is the smallest region that contains 25\% of the probability mass defined by the kernel estimate, see \code{\link{geom_contour_ks}}.

    Since these default probability contours are relative contour levels, they aren't suitable for producing a contour plot with fixed contour levels across all groups. It may require trial and error to obtain a single set of contour levels which is appropriate for all groups: one possible choice is provided by \code{contour_breaks}. 
}

\seealso{\code{\link{geom_contour_ks}}}

\examples{
library(ggplot2)
data(crabs, package="MASS")
## tidy data 
crabs2 <- dplyr::select(crabs, FL, CW)
t1 <- tidy_kde(crabs2)
ggplot(t1, aes(x=FL, y=CW)) + geom_contour_filled_ks()

\donttest{## geospatial data
crabs2s <- sf::st_as_sf(crabs2, coords=c("FL","CW"))
s1 <- st_kde(crabs2s)
s1_cont <- st_get_contour(s1, cont=seq(10,90, by=20))

## base R plot
vc <- function(.) colorspace::sequential_hcl(palette="viridis", n=.)
plot(s1_cont[,"contperc"], pal=vc)

## ggplot
ggplot() + scale_fill_viridis_d() + 
    geom_sf(data=s1_cont, aes(fill=contperc))}
}

\keyword{smooth}
