\name{contour}
\alias{contour_breaks}
\alias{contourLevels.tidy_ks}
\alias{contourLevels.sf_ks}
\alias{st_get_contour}
\alias{label_percent}
\title{Contour functions for tidy and geospatial kernel estimates}

\description{
    Contour functions for tidy and geospatial kernel estimates.
}

\usage{
\method{contourLevels}{tidy_ks}(x, cont=c(25,50,75), group=FALSE, ...)
\method{contourLevels}{sf_ks}(x, cont=c(25,50,75), group=FALSE, ...)
contour_breaks(data, cont=c(25,50,75), group=FALSE)
label_percent(y)

st_get_contour(x, cont=c(25,50,75), breaks, which_deriv_ind, disjoint=TRUE,
    as_point=FALSE)
}

\arguments{
    \item{x,data}{tidy kernel estimate (output from \code{tidy_k*}) or geospatial kernel estimate (output from \code{st_k*})}
    \item{cont}{vector of contour levels. Default is c(25,50,75).}
    \item{group}{flag to compute contour levels per group. Default is FALSE.}
    \item{breaks}{tibble or vector of contour levels (e.g. output from \code{contour_breaks})}
    \item{which_deriv_ind}{derivative index (only required for \code{st_kdde} objects)}
    \item{disjoint}{flag to compute disjoint contours. Default is TRUE.}
    \item{as_point}{flag to return polygons as point coordinates in tidy format. Default is TRUE.}  
    \item{y}{factor variable}
    \item{...}{other parameters (not implemented)}    
}

\value{
    The output from \code{contour_breaks} is a tibble of the values of the contour breaks. The output from \code{st_get_contour} is an sf object of the contours as multipolygons.
}

\details{
    By default, the 1\% to 99\% contours are computed for an \code{st_k*} output, though a plot of all 99 of them would be too crowded. \code{st_get_contour} selects a subset of these, as specified by \code{cont}. If a contour level in \code{cont} does not already exist or if absolute contour levels are specified in \code{breaks}, then the corresponding contours are computed. If \code{disjoint=TRUE} (default) then the contours are computed as a set of disjoint multipolygons: this allows for plotting without overlapping transparent colours. If \code{disjoint=FALSE} then the contours are overlapping and so their colours alpha-mixed, but they strictly satisfy the probabilistic definition, e.g. a 25\% contour region is the smallest region that contains 25\% of the probability mass defined by the kernel estimate, see \code{\link{geom_contour_ks}}.

    Since these default probability contours are relative contour levels, they aren't suitable for producing a contour plot with fixed contour levels across all groups. It may require trial and error to obtain a single set of contour levels which is appropriate for all groups: one possible choice is provided by \code{contour_breaks}. 
}

\seealso{\code{\link{geom_contour_ks}}}

\examples{
data(crabs, package="MASS")
crabs2 <- dplyr::select(crabs, FL, CW, sex)
crabs2 <- dplyr::group_by(crabs2, sex)
t1 <- tidy_kde(crabs2)
b <- contour_breaks(t1)
ggplot2::ggplot(t1, ggplot2::aes(x=FL, y=CW)) + 
    geom_contour_filled_ks(colour=1, breaks=b) + 
    ggplot2::facet_wrap(~sex)

## extract contour polygons
\donttest{crabs2s <- sf::st_as_sf(crabs2, coords=c("FL","CW"))
t2 <- st_kde(crabs2s)
t2 <- st_get_contour(t2, breaks=b, as_point=TRUE)
t2 <- dplyr::rename(t2, FL=X, CW=Y)
ggplot2::ggplot(t2, ggplot2::aes(x=FL, y=CW)) + 
    ggplot2::geom_polygon(ggplot2::aes(fill=contlabel, 
    subgroup=contlabel_group), col=1) + 
    ggplot2::scale_fill_viridis_d() + 
    ggplot2::guides(fill=ggplot2::guide_legend(reverse=TRUE)) + 
    ggplot2::facet_wrap(~sex)} 
}

\keyword{smooth}
