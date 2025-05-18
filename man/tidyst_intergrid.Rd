\name{tidyst_intergrid}
\alias{tidy_intergrid}
\alias{st_intergrid}
\title{Tidy and geospatial interpolated grids}

\description{
    Tidy and geospatial interpolated grids for 2-dimensional gridded data.
}

\usage{
tidy_intergrid(data, attrib, cellsize, verbose=FALSE)
st_intergrid(x, attrib, cellsize, verbose=FALSE)
}

\arguments{
    \item{data,x}{sf object, where attribute is estimate value and incomplete/irregular estimation grid is polygon geometry}
    \item{attrib}{name or position of estimate variable. Default is 1.}
    \item{cellsize}{cell size. If missing then automatically calculated from data.}
    \item{verbose}{flag for verbose output. Default is FALSE.}
}

\value{
    The input gridded data is interpolated to a complete regular rectangular grid. 
}

\details{
   Any missing grid cells are inferred so there are no gaps in the output, and the attribute value to set to 0. For any other grid cells with missing attribute values, the attributate is also set to 0. \code{tidy_overgrid}/\code{st_overgrid} is usually deployed on gridded data from third parties, where geometries are excluded/varying to reduce storage requirements, but \code{tidy_as_kde}/\code{st_as_kde} require complete regular rectangular grids.  
}

\examples{
## geospatial quasi density estimate
library(ggplot2)
data(ales_grid, package="eks")

## incomplete 1 km x 1 km grid
## ind = #individuals in grid cells 
gs <- ggplot() + ggthemes::theme_map() + 
    colorspace::scale_fill_continuous_sequential(palette="Heat", breaks=seq(0,6000,by=1000))  
gs + geom_sf(data=ales_grid, aes(fill=ind))

## complete regular interpolated 1 km x 1 km grid
ales_sgrid <- st_intergrid(ales_grid, attrib="ind", cellsize=c(1000,1000))
gs + geom_sf(data=ales_sgrid, aes(fill=ind)) 

## geom_sf KDE plot
ales_skde <- st_as_kde(ales_sgrid)
ggplot(ales_skde) + ggthemes::theme_map() + 
    geom_sf(data=st_get_contour(ales_skde), aes(fill=contperc))
}

\keyword{smooth}
