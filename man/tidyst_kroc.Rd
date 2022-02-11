\name{tidyst_kroc}
\alias{tidy_kroc}
\alias{st_kroc}
\title{Tidy and geospatial kernel receiver operating characteristic (ROC) curve}

\description{
  Tidy and geospatial versions of kernel receiver operating characteristic (ROC) curve for 1- and 2-dimensional data.
}

\usage{
tidy_kroc(data1, data2, ...)
st_kroc(x1, x2, ...)
}

\arguments{
    \item{data1, data2}{data frames/tibbles of data values}
    \item{x1,x2}{sf objects with point geometry}
    \item{...}{other parameters in \code{ks::kroc} function}
}
  

\value{
    The output has the same structure as the 1-d kernel distribution estimate from \code{\link[=st_kcde]{*_kcde}}, except that \code{fpr} (\eqn{x}{x}-variable) is the false positive rate (complement of specificity) and \code{estimate} is the true positive rate (sensitivity), rather than the usual estimation grid points and cdf values. 
}

\details{
    A kernel ROC curve is a modification of the standard kernel distribution estimate where the two data samples are compared. For details of the computation and the bandwidth selection procedure of the kernel density ROC curve, see \code{?ks::kroc}. The bandwidth matrix of smoothing parameters is computed as in \code{ks::kcde} per data sample.
}

\examples{
## 2-d kernel ROC curve between unsuccessful and successful grafts
data(hsct, package="ks")
hsct <- dplyr::as_tibble(hsct)
hsct <- dplyr::filter(hsct, PE.Ly65Mac1 >0 & APC.CD45.2>0)
hsct6 <- dplyr::filter(hsct, subject==6)   ## unsuccessful graft 
hsct6 <- dplyr::select(hsct6, PE.Ly65Mac1, APC.CD45.2)
hsct12 <- dplyr::filter(hsct, subject==12) ## successful graft 
hsct12 <- dplyr::select(hsct12, PE.Ly65Mac1, APC.CD45.2)
t1 <- tidy_kroc(data1=hsct6, data2=hsct12)
gt <- ggplot2::ggplot(t1, ggplot2::aes(x=fpr)) 
gt + ggplot2::geom_line(colour=1) 

## geospatial local test between Grevillea species
data(wa)
data(grevilleasf)
hakeoides <- dplyr::filter(grevilleasf, species=="hakeoides")
paradoxa <- dplyr::filter(grevilleasf, species=="paradoxa")
s1 <- st_kroc(x1=hakeoides, x2=paradoxa)
gs <- ggplot2::ggplot(s1, ggplot2::aes(x=fpr)) 
gs + ggplot2::geom_line(colour=1) 
}

\keyword{smooth}
