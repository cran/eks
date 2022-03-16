\name{eks-package}
\alias{eks}
\docType{package}
\title{eks}

\description{
    This package extends the functionality of the kernel smoothing functions from the \pkg{ks} package in base \R to the tidyverse and to GIS (Geographical Information Systems) ecosystems.     
}

\details{
    As the kernel smoothers from the \pkg{ks} package are prefixed as \code{k*}, their equivalents in \pkg{eks} are systematically named as follows:   
    \itemize{
        \item \code{tidy_k*} for 1- and 2-d tidy data 
        \item \code{st_k*} for 2-d geospatial data.
    }
      
    The output data tibbles (tidy data frames provided by the \pkg{tibble} package) from \code{tidy_k*} can be visualised within the \pkg{ggplot2} graphical interface, using the usual layer functions and the custom ones supplied in this package. These \code{tidy_k*} functions are analogous to those in the \pkg{broom} and related packages, though the latter tend to focus on tidying the summary diagnostic output from model fitting (and not on tidying the underlying estimates themselves), whereas \code{tidy_k*} are more substantive since they do compute tidy estimates.   
  
    The output simple feature geometries (provided by the \pkg{sf} package) from \code{st_k*} can be visualised in the (i) \pkg{ggplot2} graphical interface using primarily the \code{geom_sf} layer function, or (ii) in the base \R graphical interface using the \code{plot} method supplied in this package. These simple feature geometries can also be exported as standard geospatial formats (e.g. shapefile, GEOS geometry) for use in external GIS software such as ArcGIS and QGIS. 
}   
  
  
\author{ Tarn Duong }

\references{
  Chacon, J.E. & Duong, T. (2018) \emph{Multivariate Kernel Smoothing
  and Its Applications}. Chapman & Hall/CRC, Boca Raton. 

  Duong, T. (2022) \emph{Statistical visualisation for tidy and geospatial data in R via kernel smoothing methods in the eks package}. Submitted. DOI:10.48550/arXiv.2203.01686
}

\keyword{package}



