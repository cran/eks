\name{tidyst_kde_boundary}
\alias{tidy_kde_boundary}
\alias{tidy_kde_truncate}
\alias{st_kde_boundary}
\alias{st_kde_truncate}
\title{Tidy and geospatial kernel density estimates with boundary and truncated kernels}

\description{
    Tidy and geospatial versions of kernel density estimates with boundary and truncated kernels for 1- and 2-dimensional data. 
}

\usage{
tidy_kde_boundary(data, ...)
tidy_kde_truncate(data, boundary, ...)
st_kde_boundary(x, ...)
st_kde_truncate(x, boundary, ...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{boundary}{data frame/sf point geometry of boundary} 
    \item{...}{other parameters in \code{ks::kde.boundary} function}
}
  
\value{
   The outputs from \code{*_kde_boundary}, \code{*_kde_truncate} have the same structure as the standard kernel density estimate from \code{\link[=st_kde]{*_kde}}. 
}

\details{
    A boundary kernel density estimate is a modification of the standard density estimate for bounded data. There are two main types: beta kernels (\code{boundary.kernel="beta"}) and linear kernels (\code{boundary.kernel="linear"}). For details of the computation of the boundary kernel estimates and of the bandwidth selector procedure, see \code{ks::kde.boundary}. Currently only rectangular boundaries are supported, as defined by \code{xmin} x \code{xmax}. 
    
    A truncated kernel density estimate is a modification of the standard density estimate for bounded data. All the probability mass outside of \code{boundary} is set to zero and re-assigned over the regions inside in the boundary. The boundary can be any polygon. For further details of the computation of the truncated kernel estimate, see \code{ks::kde.truncate}. 

    For details of the computation of the boundary kernel estimates and the truncated kernel density estimates, and of the bandwidth selector procedure, see \code{ks::kde.boundary}, \code{ks::kde.truncate}.
}

\examples{
## tidy boundary density estimates
library(ggplot2)
data(worldbank, package="ks")
worldbank <- dplyr::as_tibble(worldbank)
## percentage data is bounded on [0,100] x [0,100]
wb2 <- na.omit(worldbank[,c("internet", "ag.value")])
xmin <- c(0,0); xmax <- c(100,100)
rectb <- data.frame(x=c(0,100,100,0,0), y=c(0,0,100,100,0))

## standard density estimate
t1 <- tidy_kde(wb2)
## beta boundary density estimate
t2 <- tidy_kde_boundary(wb2, boundary.kernel="beta", xmin=xmin, xmax=xmax)
## linear boundary density estimate
t3 <- tidy_kde_boundary(wb2, boundary.kernel="linear", xmin=xmin, xmax=xmax)
## tidy truncated density estimate
t4 <- tidy_kde_truncate(wb2, boundary=rectb)

t5 <- c(t1, t2, t3, t4, labels=c("Standard KDE","Beta bd KDE", "Linear bd KDE",
    "Truncated KDE"))

## standard estimate exceeds boundary but not boundary or truncated estimates
gr <- geom_polygon(data=rectb, inherit.aes=FALSE, aes(x=x,y=y), 
    fill=NA, colour=1, linetype="dashed")
gt <- ggplot(t5, aes(x=internet,y=ag.value)) 
gt + geom_contour_ks(aes(colour=group)) + gr + facet_wrap(~group)

## geospatial boundary density estimates
\donttest{data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
hakeoides_crop <- sf::st_crop(hakeoides, xmin=4e5, xmax=5.7e5, ymin=6.47e6, ymax=7e6)
hakeoides_bbox <- sf::st_as_sfc(sf::st_bbox(hakeoides_crop))
xminb <- sf::st_bbox(hakeoides_crop)[1:2]
xmaxb <- sf::st_bbox(hakeoides_crop)[3:4]
s1 <- st_kde(hakeoides_crop)
s2 <- st_kde_boundary(hakeoides_crop, boundary.kernel="beta", 
    xmin=xminb, xmax=xmaxb)
s3 <- st_kde_boundary(hakeoides_crop, boundary.kernel="linear", 
    xmin=xminb, xmax=xmaxb)
## geospatial truncated density estimate    
s4 <- st_kde_truncate(x=hakeoides_crop, boundary=hakeoides_bbox)
s5 <- c(s1, s2, s3, s4, labels=c("Standard KDE","Beta bd KDE", "Linear bd KDE",
    "Truncated KDE"))

## base R plots
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(hakeoides_bbox, add=TRUE, lty=3, lwd=2)
plot(s1, add=TRUE, border=1, col="transparent", legend=FALSE)
plot(s2, add=TRUE, border=2, col="transparent", legend=FALSE)
mapsf::mf_legend(type="symb", val=c("Standard KDE","Beta bd KDE"), pal=c(1,2), 
    cex=rep(3,2), pch=rep("-",2), title="Density", pos="bottomleft")

plot(wa, xlim=xlim, ylim=ylim)
plot(hakeoides_bbox, add=TRUE, lty=3, lwd=2)
plot(s1, add=TRUE, border=1, col="transparent", legend=FALSE)
plot(s3, add=TRUE, border=3, col="transparent", legend=FALSE)
mapsf::mf_legend(type="symb", val=c("Standard KDE","Linear bd KDE"), pal=c(1,3), 
    cex=rep(3,2), pch=rep("-",2), title="Density", pos="bottomleft")

plot(wa, xlim=xlim, ylim=ylim)
plot(hakeoides_bbox, add=TRUE, lty=3, lwd=2)
plot(s1, add=TRUE, border=1, col="transparent", legend=FALSE)
plot(s4, add=TRUE, border=4, col="transparent", legend=FALSE)
mapsf::mf_legend(type="symb", val=c("Standard KDE","Truncated KDE"), pal=c(1,4), 
    cex=rep(3,2), pch=rep("-",2), title="Density", pos="bottomleft")

## geom_sf plots
gs <- ggplot(s1) + geom_sf(data=wa, fill=NA) + 
    geom_sf(data=hakeoides_bbox, linewidth=1.2, linetype="dotted", fill=NA) + 
    geom_sf(data=dplyr::mutate(st_get_contour(s1), gr="Standard KDE"), 
        aes(colour=gr), fill=NA) + ggthemes::theme_map()
gs + geom_sf(data=st_get_contour(s5), aes(colour=group), fill=NA) +
    coord_sf(xlim=xlim, ylim=ylim) + facet_wrap(~group)}
}

\keyword{smooth}
