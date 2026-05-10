## Introduction 

This package extends the functionality of the kernel smoothing functions from the `ks` package in base `R` to the tidyverse and to GIS (Geographical Information Systems) ecosystems. 

As the kernel smoothers from the `ks` package are prefixed as `k*`, their equivalents in `eks` are systematically named as follows:   

* `tidy_k*` for 1- and 2-d tidy data 
* `st_k*` for 2-d geospatial data.
    
The output data tibbles (tidy data frames provided by the `tibble` package) from `tidy_k*` can be visualised within the `ggplot2` graphical interface, using the usual layer functions and the custom ones supplied in this package. These `tidy_k*` functions are analogous to those in the `broom` and related packages, though the latter tend to focus on tidying the summary diagnostic output from model fitting (and not on tidying the underlying estimates themselves), whereas `tidy_k*` are more substantive since they do compute tidy estimates.   
  
The output simple feature geometries (provided by the `sf` package) from `st_k*` can be visualised in the (i) `ggplot2` graphical interface using primarily the `ggplot2::geom_sf` layer function, or (ii) in the base `R` graphical interface using the `plot` method supplied in this package. These simple feature geometries can also be exported as standard geospatial formats (e.g. shapefile, GEOS geometry) for use in external GIS software such as [ArcGIS](https://www.arcgis.com/index.html) and [QGIS](https://qgis.org/). 

For an illustration of kernel density estimates for tidy and geospatial data, see `vignette("eks")`.

## Installation

Install from CRAN:

```r
install.packages("eks")
```

## Further reading

Chacon, J. E. and Duong, T. (2018) [_Multivariate Kernel Smoothing and Its Applications_](https://mvstat.net/mvksa/) Chapman & Hall/CRC Press, Boca Raton.

Duong, T. (2025) [Statistical visualisation for tidy and geospatial data in R via kernel smoothing methods in the eks package](https://mvstat.net/tduong/research/publications/duong-2025-compstat.pdf) _Computational Statistics_ **40**, 2825--2847.