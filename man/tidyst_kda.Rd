\name{tidyst_kda}
\alias{tidy_kda}
\alias{st_kda}
\title{Tidy and geospatial kernel discrimination analysis (classification)}

\description{
  Tidy and geospatial versions of kernel discrimination analysis (classification) for 1- and 2-dimensional data.
}

\usage{
tidy_kda(data, ...)
st_kda(x, ...)
}

\arguments{
    \item{data}{grouped tibble of data values}
    \item{x}{sf object with grouping attribute and with point geometry}
    \item{...}{other parameters in \code{ks::kda} function}
}
  
\value{
    --For \code{tidy_kda}, the output is an object of class \code{tidy_ks}, which is a tibble with columns: 
    \item{x}{evaluation points in x-axis (name is taken from 1st input variable in \code{data})}
    \item{y}{evaluation points in y-axis (2-d) (name is taken from 2nd input variable in \code{data})}
    \item{estimate}{weighted kernel density estimate values}
    \item{prior_prob}{prior probabilities for each group}
    \item{ks}{first row (within each \code{group}) contains the untidy kernel estimate from \code{ks::kda}}
    \item{tks}{short object class label derived from the \code{ks} object class}
    \item{label}{estimated KDA group label at \code{(x,y)}}
    \item{group}{grouping variable (same as input).}

    --For \code{st_kda}, the output is an object of class \code{st_ks}, which is a list with fields:
    \item{tidy_ks}{tibble of simplified output (\code{ks}, \code{tks}, \code{label}, \code{group}) from \code{tidy_kda}}
    \item{grid}{sf object of grid of weighted kernel density estimate values, as polygons, with attributes \code{estimate}, \code{label}, \code{group} copied from the \code{tidy_ks} object}
    \item{sf}{sf object of 1\% to 99\% contour regions of weighted kernel density estimate, as multipolygons, with attributes \code{contlabel} derived from the contour level; and \code{estimate}, \code{group} copied from the \code{tidy_ks} object.}
}

\details{
    A kernel discriminant analysis (aka classification or supervised learning) assigns each grid point to the group with the highest density value, weighted by the prior probabilities. 

    The output from \code{*_kda} have the same structure as the kernel density estimate from \code{\link[=st_kde]{*_kde}}, except that \code{estimate} is the weighted kernel density values at the grid points (weighted by \code{prior_prob}), and \code{label} becomes the KDA grouping variable that indicates to which of the groups the grid points belong. The output is a grouped tibble, grouped by the input grouping variable.

    For details of the computation of the kernel discriminant analysis and the bandwidth selector procedure, see \code{?ks::kda}. The bandwidth matrix of smoothing parameters is computed as in \code{ks::kde} per group. 
}

\examples{
## tidy discriminant analysis (classification)
library(ggplot2)
data(cardio, package="ks")
cardio <- dplyr::as_tibble(cardio[,c("ASTV","Mean","NSP")])
cardio <- dplyr::mutate(cardio, NSP=ordered(NSP))
cardio <- dplyr::group_by(cardio, NSP)
set.seed(8192)
cardio.train.ind <- sample(1:nrow(cardio), round(nrow(cardio)/4,0))
cardio.train <- cardio[cardio.train.ind,]
cardio.train1 <- dplyr::select(cardio.train, ASTV, NSP)
cardio.train2 <- dplyr::select(cardio.train, ASTV, Mean, NSP)

## tidy 1-d classification
t1 <- tidy_kda(cardio.train1) 
gt1 <- ggplot(t1, aes(x=ASTV)) 
gt1 + geom_line(aes(colour=NSP)) + 
    geom_rug(aes(colour=label), sides="b", linewidth=1.5) +
    scale_colour_brewer(palette="Dark2", na.translate=FALSE) 

\donttest{## tidy 2-d classification
t2 <- tidy_kda(cardio.train2)
gt2 <- ggplot(t2, aes(x=ASTV, y=Mean)) + theme_bw()
gt2 + geom_contour_ks(aes(colour=NSP)) + 
    geom_tile(aes(fill=label), alpha=0.2) +
    scale_fill_brewer(palette="Dark2", na.translate=FALSE) +
    scale_colour_brewer(palette="Dark2")}

## geospatial classification
\donttest{data(wa)
data(grevilleasf)
grevillea_gr <- dplyr::filter(grevilleasf, species=="hakeoides" |
    species=="paradoxa")
grevillea_gr <- dplyr::mutate(grevillea_gr, species=factor(species))  
grevillea_gr <- dplyr::group_by(grevillea_gr, species)
s1 <- st_kda(grevillea_gr)
s2 <- st_ksupp(st_kde(grevillea_gr))
s1$grid <- sf::st_filter(s1$grid, sf::st_convex_hull(sf::st_union(s2$sf)))

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, which_geometry="grid", add=TRUE, border=NA, legend=FALSE)
plot(s1, add=TRUE, lwd=2, border=rep(colorspace::qualitative_hcl(
    palette="Dark2", n=2, alpha=0.5), each=3))

## geom_sf plot
gs1 <- ggplot(s1) + geom_sf(data=wa, fill=NA) + 
    geom_sf(data=dplyr::mutate(s1$grid, species=label), aes(fill=species), 
    alpha=0.1, colour=NA) + ggthemes::theme_map()
gs1 + geom_sf(data=st_get_contour(s1), aes(colour=species), fill=NA) +
    colorspace::scale_colour_discrete_qualitative(palette="Dark2") +
    colorspace::scale_fill_discrete_qualitative(palette="Dark2") +
    facet_wrap(~species) + coord_sf(xlim=xlim, ylim=ylim)}
}

\keyword{smooth}
