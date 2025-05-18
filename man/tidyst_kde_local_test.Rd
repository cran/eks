\name{tidyst_kde_local_test}
\alias{tidy_kde_local_test}
\alias{st_kde_local_test}
\title{Tidy and geospatial kernel density based local two-sample comparison tests}

\description{
  Tidy and geospatial versions of kernel density based local two-sample comparison tests for 1- and 2-dimensional data.
}

\usage{
tidy_kde_local_test(data1, data2, labels, ...)
st_kde_local_test(x1, x2, labels, ...)
}

\arguments{
    \item{data1, data2}{data frames/tibbles of data values}
    \item{x1,x2}{sf objects with point geometry}
    \item{labels}{flag or vector of strings for legend labels}
    \item{...}{other parameters in \code{ks::kde.local.test} function}
}
  
\value{
    The output has the same structure as the kernel density estimate from \code{\link[=st_kde]{*_kde}}, except that \code{estimate} is the difference between the density values \code{data1-data2} rather than the density values, and \code{label} becomes an indicator factor of the local comparison test result: "f1<f2" = \code{data1 < data2}, 0 = \code{data1 = data2}, "f2>f1" = \code{data1 > data2}.

    The output from \code{st_kde_local_test} has two contours, with \code{contlabel=-50} (for f1<f2) and \code{contlabel=50} (for f1>f2), as multipolygons which delimit the significant difference regions.
}

\details{
    A kernel local density based two-sample comparison is a modification of the standard kernel density estimate where the two data samples are compared. A Hochberg procedure is employed to control the significance level for multiple comparison tests. 

    For details of the computation of the kernel local density based two-sample comparison test and the bandwidth selector procedure, see \code{?ks::kde.local.test}. The bandwidth matrix of smoothing parameters is computed as in \code{ks::kde} per data sample. 

    If \code{labels} is missing, then the first sample label is taken from \code{x1}, and the second sample label from \code{x2}. If \code{labels="default"} then these are "f1" and "f2". Otherwise, they are assigned to the values of the input vector of strings.
}

\examples{
## tidy local test between unsuccessful and successful grafts
library(ggplot2)
data(hsct, package="ks")
hsct <- dplyr::as_tibble(hsct)
hsct <- dplyr::filter(hsct, PE.Ly65Mac1 >0 & APC.CD45.2>0)
hsct6 <- dplyr::filter(hsct, subject==6)   ## unsuccessful graft 
hsct6 <- dplyr::select(hsct6, PE.Ly65Mac1, APC.CD45.2)
hsct12 <- dplyr::filter(hsct, subject==12) ## successful graft 
hsct12 <- dplyr::select(hsct12, PE.Ly65Mac1, APC.CD45.2)
t1 <- tidy_kde_local_test(data1=hsct6, data2=hsct12)
gt <- ggplot(t1, aes(x=PE.Ly65Mac1, y=APC.CD45.2)) + coord_fixed() 
gt + geom_contour_filled_ks(aes(fill=after_stat(contregion)))

\donttest{t2 <- tidy_kfs(hsct12)
gt + geom_contour_ks(data=t1, aes(colour=after_stat(contregion))) +
    geom_contour_filled_ks(data=t2, aes(fill=after_stat(contregion))) + 
    colorspace::scale_fill_discrete_qualitative(palette="Dark2", 
        rev=TRUE, nmax=2) + 
    guides(fill=guide_legend(title="Signif curv"))}

## geospatial local test between Grevillea species
\donttest{data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
paradoxa <- dplyr::filter(grevilleasf, species=="paradoxa")
s1 <- st_kde_local_test(x1=hakeoides, x2=paradoxa)
s2 <- st_kfs(hakeoides)

## base R plot
xlim <- c(1.2e5, 1.1e6); ylim <- c(6.1e6, 7.2e6)
plot(wa, xlim=xlim, ylim=ylim)
plot(s1, add=TRUE)
plot(s2, add=TRUE, pos="bottom")

## geom_sf plot
gs <- ggplot(s1) + geom_sf(data=wa, fill=NA) + ggthemes::theme_map() 
gs + geom_sf(data=st_get_contour(s1), aes(colour=contregion), fill=NA) +
    geom_sf(data=st_get_contour(s2), aes(fill=contregion)) + 
    guides(fill=guide_legend(title="Signif curv")) + 
    scale_fill_manual(values=7) + coord_sf(xlim=xlim, ylim=ylim)}
}

\keyword{smooth}
