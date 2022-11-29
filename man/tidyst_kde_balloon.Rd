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
\donttest{data(worldbank, package="ks")
worldbank <- dplyr::as_tibble(worldbank)
wb2 <- na.omit(worldbank[,c("GDP.growth", "inflation")])
xmin <- c(-70,-25); xmax <- c(25,70)

## standard density estimate
t1 <- tidy_kde(wb2, xmin=xmin, xmax=xmax)
## sample point variable density estimate
t2 <- tidy_kde_sp(wb2, xmin=xmin, xmax=xmax)
tt <- c(t1, t2, labels=c("Standard KDE","Sample point KDE"))

## fixed contour levels for all three plots
b <- contour_breaks(tt)
gt <- ggplot2::ggplot(tt, ggplot2::aes(x=GDP.growth, y=inflation)) 
gt + geom_contour_filled_ks(breaks=b, colour=1) + 
    colorspace::scale_fill_discrete_sequential() +
    ggplot2::facet_wrap(~group)

## balloon variable density estimate
## gridsize=c(21,21) only for illustrative purposes
t3 <- tidy_kde_balloon(wb2, xmin=xmin, xmax=xmax, gridsize=c(21,21))
tt <- c(t1, t2, t3, labels=c("Standard KDE","Sample point KDE","Balloon KDE"))
b <- contour_breaks(tt, cont=seq(10,90,by=10))
gt + geom_contour_filled_ks(data=tt, breaks=b, colour=1) + 
    colorspace::scale_fill_discrete_sequential() +
    ggplot2::facet_wrap(~group)

## geospatial variable density estimates
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")

## standard density estimate
s1 <- st_kde(hakeoides)
## sample point variable density estimate
s2 <- st_kde_sp(hakeoides)   
s3 <- c(s1, s2, labels=c("Standard KDE","Sample point KDE"))
b <- contour_breaks(s3)
bcols <- colorspace::sequential_hcl(nrow(b), palette="Heat2", rev=TRUE)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, add=TRUE, col=bcols[1:2], breaks=b)
plot(wa, xlim=xlim, ylim=ylim)
plot(s2, add=TRUE, col=bcols, breaks=b)

## geom_sf plot
gs <- ggplot2::ggplot(s3) + ggplot2::geom_sf(data=wa, fill=NA) + 
    ggthemes::theme_map()
gs + ggplot2::geom_sf(data=st_get_contour(s3, breaks=b),  
    ggplot2::aes(fill=contlabel)) + 
    colorspace::scale_fill_discrete_sequential(palette="Heat2") +
    ggplot2::coord_sf(xlim=xlim, ylim=ylim) +
    ggplot2::facet_wrap(~group)}
}

\keyword{smooth}
