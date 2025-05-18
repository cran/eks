\name{tidyst_kde_balloon}
\alias{tidy_kde_balloon}
\alias{tidy_kde_sp}
\alias{st_kde_balloon}
\alias{st_kde_sp}
\title{Tidy and geospatial kernel density estimates with variable kernels}

\description{
    Tidy and geospatial versions of kernel density estimates with variable kernels for 2-dimensional data.
}

\usage{
tidy_kde_balloon(data, ...)
tidy_kde_sp(data, ...)
st_kde_balloon(x, ...)
st_kde_sp(x, ...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{...}{other parameters in \code{ks::kde.balloon}, \code{ks::kde.sp} functions}
}
  
\value{
    The outputs from \code{*_kde_balloon}, \code{*_kde_sp} have the same structure as the standard kernel density estimate from \code{\link[=st_kde]{*_kde}}. 
}

\details{
    A variable kernel density estimate is a modification of the standard density estimate where the bandwidth matrix is variable. There are two main types: balloon kernel estimates (\code{*_kde_balloon}) where the bandwidth varies with the grid point, and sample point kernel estimates (\code{*_kde_sp}) where the bandwidth varies with the data points. For details of the computation of the variable kernel estimates and of the bandwidth selector procedure, see \code{ks::kde.balloon}, \code{ks::kde.sp}.
}

\examples{
## tidy variable density estimates
\donttest{library(ggplot2)
data(worldbank, package="ks")
worldbank <- dplyr::as_tibble(worldbank)
wb2 <- na.omit(worldbank[,c("GDP.growth", "inflation")])
xmin <- c(-70,-25); xmax <- c(25,70)

## standard density estimate
t1 <- tidy_kde(wb2, xmin=xmin, xmax=xmax)
## sample point variable density estimate
t2 <- tidy_kde_sp(wb2, xmin=xmin, xmax=xmax)
## balloon variable density estimate
## gridsize=c(21,21) only for illustrative purposes
t3 <- tidy_kde_balloon(wb2, xmin=xmin, xmax=xmax, gridsize=c(21,21))
t4 <- rbind(t1, t2, t3, labels=c("Standard KDE","Sample point KDE","Balloon KDE"))
tb <- contour_breaks(t4, group=FALSE)
gt <- ggplot(t4, aes(x=GDP.growth, y=inflation))
gt + geom_contour_filled_ks(breaks=tb, colour=1) + facet_wrap(~group)

## geospatial variable density estimates
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")

## standard density estimate
s1 <- st_kde(hakeoides)
## sample point variable density estimate
s2 <- st_kde_sp(hakeoides)   
s3 <- c(s1, s2, labels=c("Standard KDE","Sample point KDE"))
sb <- contour_breaks(s3, group=FALSE)
bcols <- colorspace::sequential_hcl(nrow(sb), palette="Heat", rev=TRUE)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
layout(matrix(1:2, ncol=2))
par(mar=(c(0,0,0,0)))
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, add=TRUE, col=bcols[1:2], breaks=sb, legend=FALSE)
par(mar=(c(0,0,0,0)))
plot(wa, xlim=xlim, ylim=ylim)
plot(s2, add=TRUE, col=bcols, breaks=sb)
layout(1)

## geom_sf plot
gs <- ggplot(s3) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map()
gs + geom_sf(data=st_get_contour(s3, breaks=sb), aes(fill=estimate)) +
    coord_sf(xlim=xlim, ylim=ylim) + facet_wrap(~group)}
}

\keyword{smooth}
