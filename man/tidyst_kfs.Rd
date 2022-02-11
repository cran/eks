\name{tidyst_kfs}
\alias{tidy_kfs}
\alias{st_kfs}

\title{Tidy and geospatial kernel feature significance}

\description{
  Tidy and geospatial versions of kernel feature significance for 1- and 2-dimensional data.
}

\usage{
tidy_kfs(data, ...)
st_kfs(x, ...)
}

\arguments{
    \item{data}{data frame/tibble of data values}
    \item{x}{sf object with point geometry}
    \item{...}{other parameters in \code{ks::kfs} function}
}
  

\value{
   The output from \code{tidy_kfs} has the same structure as the kernel density estimate from \code{\link{tidy_kde}}, except that all values of \code{estimate} outside of the significant curvature regions are set to zero, and the \code{label} indicates whether the corresponding \code{x,y} point is inside a significant curvature region. 
   
   The output from \code{st_kfs} has a single contour (with \code{contlabel=50}) as a multipolygon which delimits significant curvature regions. 
}
 
\details{
    A significant kernel curvature region consist of all points whose density curvature value is significantly different less than zero (i.e. forms a bump surrounding a local maximum). A Hochberg procedure is employed to control the significance level for multiple significance tests. 

    For details of the computation of the significant kernel curvature regions, see \code{?ks::kfs}. The bandwidth matrix of smoothing parameters is computed as in \code{ks::kdde(deriv_order=2)}.  
}

\examples{
\donttest{## tidy significant curvature regions
data(hsct, package="ks")
hsct <- dplyr::as_tibble(hsct)
hsct <- dplyr::filter(hsct, PE.Ly65Mac1>0 & APC.CD45.2>0)
hsct12 <- dplyr::filter(hsct, subject==12)   
hsct12 <- dplyr::select(hsct12, PE.Ly65Mac1, APC.CD45.2)
t1 <- tidy_kfs(hsct12)
gt <- ggplot2::ggplot(t1, ggplot2::aes(x=PE.Ly65Mac1, y=APC.CD45.2)) 
gt + geom_contour_filled_ks(ggplot2::aes(colour=label), colour=1) +
    ggplot2::scale_fill_manual(values=7) 

## geospatial significant curvature regions
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
s1 <- st_kfs(hakeoides)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, add=TRUE)

## geom_sf plot
gs <- ggplot2::ggplot(s1) + ggplot2::geom_sf(data=wa, fill=NA) 
gs + ggplot2::geom_sf(data=st_get_contour(s1, cont=50), ggplot2::aes(fill=label)) +
    ggplot2::scale_fill_manual(values=7) +
    ggplot2::coord_sf(xlim=xlim, ylim=ylim) + 
    ggthemes::theme_map()}
}

\keyword{smooth}
