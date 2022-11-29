\name{st_add_coordinates}
\alias{st_add_coordinates}
\title{Add coordinates as attributes to geospatial data}

\description{
    Add coordinates as attributes to geospatial data.
}

\usage{
st_add_coordinates(x, as_sf=FALSE, as_tibble=FALSE, rename=TRUE)
}

\arguments{
    \item{x}{\code{sf} object with point geometry}
    \item{as_sf}{flag for output as \code{sf} object. Default is TRUE.}
    \item{as_tibble}{flag for output as tibble. Default is TRUE.}
    \item{rename}{flag to rename output from  \code{X,Y} to \code{lon,lat}. Default is TRUE.}
}

\value{
    The longitude and latitude of the point geometry are added as attributes.  
}

\details{
    The \code{sf::st_coordinates} is applied to obtain the longitude and latitude coordinates. 
}

\examples{
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
hakeoides_coord <- st_add_coordinates(hakeoides)
}

\keyword{smooth}
