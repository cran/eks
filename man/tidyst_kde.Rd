\name{tidyst_kde}
\alias{tidy_kde}
\alias{st_kde}
\title{Tidy and geospatial kernel density estimates}

\description{
    Tidy and geospatial versions of kernel density estimates for 1- and 2-dimensional data.
}

\usage{
tidy_kde(data, ...)
st_kde(x, ...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{...}{other parameters in \code{ks::kde} function}
}

\value{
    --For \code{tidy_kde}, the output is an object of class \code{tidy_ks}, which is a tibble with columns:
    \item{x}{evaluation points in x-axis (name is taken from 1st input variable in \code{data})}
    \item{y}{evaluation points in y-axis (2-d) (name is taken from 2nd input variable in \code{data})}
    \item{estimate}{kernel estimate values}
    \item{ks}{first row (within each \code{group}) contains the untidy kernel estimate from \code{ks::kde}}
    \item{tks}{short object class label derived from the \code{ks} object class}
    \item{label}{long object class label}
    \item{group}{grouping variable (if grouped input) (name is taken from grouping variable in \code{data}).}

    --For \code{st_kde}, the output is an object of class \code{st_ks}, which is a list with fields:
    \item{tidy_ks}{tibble of simplified output (\code{ks}, \code{tks}, \code{label}, \code{group}) from \code{tidy_kde}}
    \item{grid}{sf object of grid of kernel density estimate values, as polygons, with attributes \code{estimate}, \code{group} copied from the \code{tidy_ks} object}
    \item{sf}{sf object of 5\%, 10\%, ..., 95\% 
    contour regions of kernel density estimate, as multipolygons, with attributes \code{estimate}, \code{contlabel}, \code{contperc}, \code{contregion} derived from the contour level, and \code{group} copied from the \code{tidy_ks} object.}
}

\details{
    For \code{tidy_kde}, the first columns of the output tibble are copied from \code{aes(x)} (1-d) or \code{aes(x,y)} (2-d). These columns are the evaluation grid points. The \code{estimate} column is the kernel density values at these grid points. The \code{group} column is a copy of the grouping variable of the input data. The \code{ks} column is a copy of the untidy kernel estimate from \code{ks::kde}, since the calculations for the layer functions \code{\link{geom_contour_ks}}, \code{\link{geom_contour_filled_ks}} require both the observations \code{data} and the kernel estimate as a \code{kde} object. For this reason, it is advised to compute a tidy kernel estimate first and then to create a \code{ggplot} with this tidy kernel estimate as the default \code{data} in the layer. 

    For \code{st_kde}, the output list contains the field \code{tidy_ks} which is the output from \code{tidy_ks}. The field \code{grid} is the kernel estimate values, with rectangular polygons. The field \code{sf} is the 5\%, 10\%, ..., 95\% probability contour regions as multipolygons, with the derived attribute \code{contlabel = 100\%-cont}.

    The structure of the \code{tidy_kde} output is inherited from the input, i.e. if the input is a data frame/ (grouped) tibble then the output is a data frame/(grouped) tibble. Likewise for the \code{sf} object outputs for \code{st_kde}.

    The default bandwidth matrix is the unconstrained plug-in selector \code{ks::Hpi}, which is suitable for a wide range of data sets, since it is not restrained to smoothing along the coordinate axes. This produces a kernel estimate which is more representative of the data than with the default bandwidth in \code{geom_density_2d} and \code{geom_density_2d_filled}. For further details of the computation of the kernel density estimate and the bandwidth selector procedure, see \code{?ks::kde}.

    If \code{breaks} is missing, then for \code{geom_contour_ks}, \code{geom_filled_contour_ks}, the density contours at \code{cont} percent are displayed, and \code{contperc} is the default aesthetic, which shows legend keys with percentages. 
    For \pkg{eks} \eqn{\geq}{>=} 1.1.0, the aesthetic \code{label_percent(contlabel)} is deprecated in favour of \code{contperc}. The other derived quantities \code{contregion}, \code{estimate} computed via \code{stat_filled_contour_ks}, \code{stat_contour_ks}, are also suitable aesthetics and these display the contour levels and \eqn{\geq}{>=} contour levels in the legend keys respectively. For \code{geom_sf}, the default aesthetic is not suitable, and should be chosen from \code{contperc}, \code{contregion}, \code{estimate}, to produce analogous plots.  

    If \code{breaks} is not missing, then the density contours at \code{breaks} are displayed, and \code{estimate} is the default aesthetic. The other derived quantity \code{contregion} is also suitable, whereas \code{contperc} is NA so is not suitable. 

    For facetted plots, it is recommended to set the same \code{breaks} for all facets, and to use \code{estimate} or \code{contregion} as aesthetics. The same density contour percentages \code{cont} can lead to different contour levels for different densities in different facets which is not suitable for a facetted plot.
}

\examples{
## tidy density estimates
library(ggplot2)
data(crabs, package="MASS")
## tidy 1-d density estimate per species
crabs1 <- dplyr::select(crabs, FL, sp)
crabs1 <- dplyr::group_by(crabs1, sp)
t1 <- tidy_kde(crabs1)
gt1 <- ggplot(t1, aes(x=FL)) 
gt1 + geom_line(colour=1) + geom_rug_ks(colour=4) + facet_wrap(~sp)

## tidy 2-d density estimate
## suitable smoothing matrix gives bimodal estimate
crabs2 <- dplyr::select(crabs, FL, CW)
t2 <- tidy_kde(crabs2)
gt2 <- ggplot(t2, aes(x=FL, y=CW)) 

## filled contour regions with density percentage labels
gt2 + geom_contour_filled_ks(colour=1) 
## filled contour regions with density level labels
gt2 + geom_contour_filled_ks(colour=1, aes(fill=after_stat(estimate)))
gt2 + geom_contour_filled_ks(colour=1, aes(fill=after_stat(contregion)))

## contour lines only with density percentage labels
gt2 + geom_contour_ks(aes(colour=after_stat(contperc)))
## contour lines only with density level labels
gt2 + geom_contour_ks(aes(colour=after_stat(estimate)))
gt2 + geom_contour_ks(aes(colour=after_stat(contregion)))

## default smoothing in geom_density_2d_filled() gives unimodal estimate
gt3 <- ggplot(crabs2, aes(x=FL, y=CW)) 
gt3 + geom_density_2d_filled(bins=4, colour=1) +
    colorspace::scale_fill_discrete_sequential(palette="Heat") +
    guides(fill=guide_legend(title="Density", reverse=TRUE))

\donttest{## facet wrapped geom_sf plot with fixed contour levels for all facets
crabs3 <- dplyr::select(crabs, FL, CW, sex)
t3 <- tidy_kde(dplyr::group_by(crabs3, sex))
tb <- contour_breaks(t3, group=FALSE)
gt4 <- ggplot(t3, aes(x=FL, y=CW)) + facet_wrap(~sex)
## filled contours
gt4 + geom_contour_filled_ks(colour=1, breaks=tb, aes(fill=after_stat(estimate))) 
gt4 + geom_contour_filled_ks(colour=1, breaks=tb, aes(fill=after_stat(contregion))) 
## contour lines 
gt4 + geom_contour_ks(breaks=tb, aes(colour=after_stat(estimate)))
gt4 + geom_contour_ks(breaks=tb, aes(colour=after_stat(contregion))) 

## geospatial density estimate
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
hakeoides_coord <- data.frame(sf::st_coordinates(hakeoides))
s1 <- st_kde(hakeoides)

## base R plot
## filled contour regions with density percentage labels
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, add=TRUE)

## geom_sf plot
## suitable smoothing matrix gives optimally smoothed contours
gs1 <- ggplot(s1) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map() 
glim <- coord_sf(xlim=xlim, ylim=ylim)
## filled contour regions with density percentage labels
gs1 + geom_sf(data=st_get_contour(s1), aes(fill=contperc)) + glim 
## filled contour regions with density level labels
gs1 + geom_sf(data=st_get_contour(s1), aes(fill=estimate)) + glim
gs1 + geom_sf(data=st_get_contour(s1), aes(fill=contregion)) + glim

## default smoothing in geom_density_2d_filled() is oversmoothed
ggplot(hakeoides_coord) + geom_sf(data=wa, fill=NA) + 
    ggthemes::theme_map() + 
    colorspace::scale_fill_discrete_sequential(palette="Heat") + 
    geom_density_2d_filled(aes(x=X, y=Y), bins=4, colour=1) +
    guides(fill=guide_legend(title="Density", reverse=TRUE)) + glim
    
## facet wrapped geom_sf plot with fixed contour levels for all facets
grevillea_gr <- dplyr::filter(grevilleasf, species=="hakeoides" |
    species=="paradoxa")
grevillea_gr <- dplyr::mutate(grevillea_gr, species=factor(species))  
grevillea_gr <- dplyr::group_by(grevillea_gr, species)
s2 <- st_kde(grevillea_gr)
sb <- contour_breaks(s2, group=TRUE)

gs3 <- ggplot(s2) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map() +
    facet_wrap(~species)
## filled contour regions with density percentage labels w/o breaks
gs3 + geom_sf(data=st_get_contour(s2), aes(fill=contperc)) +
    glim 
## filled contour regions with density level labels with breaks
gs3 + geom_sf(data=st_get_contour(s2, breaks=sb), aes(fill=estimate)) +
    glim
gs3 + geom_sf(data=st_get_contour(s2, breaks=sb), aes(fill=contregion)) +
    glim} 

\dontrun{## export as geopackage for external GIS software
sf::write_sf(wa, dsn="grevillea.gpkg", layer="wa")
sf::write_sf(hakeoides, dsn="grevillea.gpkg", layer="hakeoides")
sf::write_sf(gs1_cont, dsn="grevillea.gpkg", layer="hakeoides_cont")
sf::write_sf(s1$grid, dsn="grevillea.gpkg", layer="hakeoides_grid")}
}

\keyword{smooth}
