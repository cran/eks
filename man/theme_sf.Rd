\name{theme_sf}
\alias{theme_sf}

\title{Theme for geospatial maps}

\description{
    Theme for geospatial maps. 
}

\usage{
theme_sf()
}

\details{
    A partial copy of \code{ggthemes::theme_map()}.
}

\examples{
library(ggplot2)
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
ggplot() + theme_sf() + geom_sf(data=wa, fill=NA) + 
    geom_sf(data=hakeoides, alpha=0.4) 
}

\keyword{hplot}
