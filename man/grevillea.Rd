\name{grevilleasf}
\docType{data}
\alias{grevilleasf}
\alias{wa}
\title{Geographical locations of Grevillea plants in Western Australia}

\description{
The \code{wa} data set contains the polygon of the administrative boundary of Western Australia (excluding islands). The \code{grevillea} data set contains the locations of 22303 grevillea plants in Western Australia. 
}

\usage{data(wa)
data(grevilleasf)
}

\format{\code{wa} is an sf object, whose geometry is the polygon in the 
  EPSG:7850 (GDA2020/MGA zone 50) projection.

  \code{grevilleasf} is an sf object with 22303 rows and 2 attributes.
  Each row corresponds to an observed plant. 
  The first column is the full scientific name,
  the second is the species name.
  The geometry is the point location of the plant in the EPSG:7850 (GDA2020/MGA zone 50) projection. This is a superset of the \code{grevillea} dataset in the \pkg{ks} package. 
} 

\source{
Atlas of Living Australia (2021). Grevillea occurrence. \url{https://biocache.ala.org.au}. Accessed on 2021-08-18. 

Geoscape Australia (2021). WA State Boundary -- Geoscape Administrative Boundaries. \url{https://data.gov.au/data/dataset/wa-state-boundary-geoscape-administrative-boundaries}.  Accessed on 2021-08-18. 
}

\keyword{datasets}

