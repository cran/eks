\name{ales_grid}
\docType{data}
\alias{ales_grid}

\title{Gridded census data in Alès, France}

\description{
The \code{ales_grid} data set contains census data of 806 1 km x 1 km grid cells which cover the city of Alès and its surrounding regions in southern France. The grid cell size was chosen to be a compromise between ensuring privacy and providing detailed public statistics from the 2019 census.
}

\usage{data(ales_grid)
}

\format{\code{ales_grid} is an sf object, whose geometry are 1 km x 1 km square polygons in the European-wide standard EPSG:3035 (ETRS89-extended / LAEA Europe) projection. There are 31 attributes collected by the French national statistical agency (Insee), e.g. \code{ind} is the number of individuals. See the URL below for more details.
} 

\source{
Insee (2023). Revenus, pauvreté et niveau de vie en 2019 - Données carroyées [In French]. \url{https://www.insee.fr/fr/statistiques/7655464?sommaire=7655515}.  Accessed on 2025-01-16. 
}

\keyword{datasets}

