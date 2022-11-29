\name{geom_contour_ks}
\alias{geom_contour_ks}
\alias{stat_contour_ks}
\alias{geom_contour_filled_ks}
\alias{stat_contour_filled_ks}
\alias{GeomContourKs}
\alias{StatContourKs}
\alias{GeomContourFilledKs}
\alias{StatContourFilledKs}
\title{Contour and filled contour plots for tidy kernel estimates}

\description{
    Contour and filled contour plots for tidy kernel estimates for 2-dimensional data.
}

\usage{
geom_contour_ks(mapping=NULL, data=NULL, stat="contour_ks",
    position="identity", ..., cont=c(25,50,75), label_percent=NULL,
    breaks=NULL, show.legend=NA, inherit.aes=TRUE)
stat_contour_ks(mapping=NULL, data=NULL, geom="contour_ks",
    position="identity", ..., cont=c(25,50,75), label_percent=NULL,
    breaks=NULL, show.legend=NA, inherit.aes=TRUE)
geom_contour_filled_ks(mapping=NULL, data=NULL, stat="contour_filled_ks",
    position="identity", ..., cont=c(25,50,75), label_percent=NULL,
    breaks=NULL, show.legend=NA, inherit.aes=TRUE)
stat_contour_filled_ks(mapping=NULL, data=NULL, geom="contour_filled_ks",
    position="identity", ..., cont=c(25,50,75), label_percent=NULL,
    breaks=NULL, show.legend=NA, inherit.aes=TRUE)
}

\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{aes()}. If specified and \code{inherit.aes = TRUE} (the
default), it is combined with the default mapping at the top level of the
plot. You must supply \code{mapping} if there is no plot mapping.}

\item{data}{The data to be displayed in this layer. There are three
options:

If \code{NULL}, the default, the data is inherited from the plot
data as specified in the call to \code{ggplot()}.

A \code{data.frame}, or other object, will override the plot
data. All objects will be fortified to produce a data frame. See
\code{fortify()} for which variables will be created.

A \code{function} will be called with a single argument,
the plot data. The return value must be a \code{data.frame}, and
will be used as the layer data. A \code{function} can be created
from a \code{formula} (e.g. \code{~ head(.x, 10)}).}

\item{stat}{The statistical transformation to use on the data for this
layer, as a string.}

\item{position}{Position adjustment, either as a string, or the result of
a call to a position adjustment function.}

\item{...}{Other arguments passed on to \code{layer()}. These are
often aesthetics, used to set an aesthetic to a fixed value, like
\code{colour="red"} or \code{size=3}. They may also be parameters
to the paired geom/stat.}

\item{cont}{Vector of contour probabilities. Default value is \code{cont=c(25,50,75}).}

\item{label_percent}{Flag for legend label as percentage. Default is TRUE.}

\item{breaks}{Numeric vector to set the contour breaks e.g. output from \code{contour_breaks}. Overrides \code{cont}.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{borders()}.}

\item{geom}{The geometric object to use display the data.}
}

\value{
    Similar output as the standard layer functions \code{ggplot2::geom_contour}, \code{geom_contour_filled} and  \code{ggplot2::stat_contour}, \code{stat_contour_filled}. 
}

\details{
    These layer functions are modifications of the standard layer functions \code{ggplot2::geom_contour}, \code{geom_contour_filled} and  \code{ggplot2::stat_contour}, \code{stat_contour_filled}. Their usage and output are similar, except that they require a tidy kernel estimate as the input, rather than the observations themselves, and that the underlying choice of the contour levels is different. For most cases, \code{geom_contour_ks} is equivalent to \code{geom_contour(stat="contour_ks")}, and likewise for \code{geom_contour_filled_ks}.

    The choice of the contour levels are based on probability contours. A 25\% contour region is the smallest region that contains 25\% of the probability mass defined by the kernel estimate. Probability contours offer a more intuitive approach to selecting the contour levels that reveal the pertinent characteristics of the kernel estimates. See Chacon & Duong (2018, Chapter 2.2). They are specified by the \code{cont} parameter: the default value is \code{cont=c(25,50,75)}, which computes the upper quartile, median and lower quartile probability contours. If \code{percent_label=TRUE}, then the legend labels are given as these percentage in \code{cont}. Otherwise, the labels are the contour levels themselves.

    Since these probability contours are computed for each group of the grouping variable in \code{data}, then these relative contour levels are different for each group. To produce a contour plot with fixed contour levels across all groups, then these can be supplied in \code{breaks}: a possible choice is provided by \code{\link{contour_breaks}}. 
}

\examples{
data(crabs, package="MASS")
crabs2 <- dplyr::select(crabs, FL, CW, sp)
crabs2 <- dplyr::group_by(crabs2, sp)
tt <- tidy_kde(crabs2)
gt <- ggplot2::ggplot(tt, ggplot2::aes(x=FL, y=CW))
gt + geom_contour_ks() + ggplot2::facet_wrap(~sp)
gt + ggplot2::geom_contour(stat="contour_ks") +
    ggplot2::facet_wrap(~sp) ## same output
gt + geom_contour_filled_ks(colour=1) + ggplot2::facet_wrap(~sp)
gt + ggplot2::geom_contour_filled(stat="contour_filled_ks", colour=1) +
    ggplot2::facet_wrap(~sp) ## same output
}

\seealso{\code{\link{contour_breaks}}}

\references{
  Chacon, J.E. & Duong, T. (2018) \emph{Multivariate Kernel Smoothing
  and Its Applications}. Chapman & Hall/CRC, Boca Raton.
}

\keyword{hplot}
