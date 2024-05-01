\name{geom_point_ks}
\alias{geom_point_ks}
\alias{stat_point_ks}
\alias{geom_rug_ks}
\alias{stat_rug_ks}
\alias{GeomPointKs}
\alias{GeomRugKs}
\alias{StatPointKs}
\alias{StatRugKs}
\title{Rug and scatter plots for tidy kernel estimates}

\description{
    Rug and scatter plots for tidy kernel estimates for 1- and 2-dimensional data. 
}

\usage{
geom_point_ks(mapping=NULL, data=NULL, stat="point_ks", position="identity", 
    ..., na.rm=FALSE, jitter=FALSE, show.legend=NA, inherit.aes=TRUE) 
stat_point_ks(mapping=NULL, data=NULL, geom="point_ks", position="identity", 
    ..., na.rm=FALSE, show.legend=NA, inherit.aes=TRUE)
geom_rug_ks(mapping=NULL, data=NULL, stat="rug_ks", position="identity", 
    ..., outside=FALSE, sides="bl", length=unit(0.03, "npc"), na.rm=FALSE, 
    jitter=FALSE, show.legend=NA, inherit.aes=TRUE) 
stat_rug_ks(mapping=NULL, data=NULL, geom="rug_ks", position="identity", 
    ..., na.rm=FALSE, show.legend=NA, inherit.aes=TRUE)   
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

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with
a warning. If \code{TRUE}, missing values are silently removed.}

\item{jitter}{Flag to jitter data before plot. Default value is FALSE.}

\item{outside}{logical that controls whether to move the rug tassels outside of the plot area. Default is off (FALSE). You will also need to use \code{coord_cartesian(clip = "off")}. When set to TRUE, also consider changing the sides argument to "tr". See examples.}

\item{sides}{A string that controls which sides of the plot the rugs appear on.
It can be set to a string containing any of \code{"trbl"}, for top, right,
bottom, and left.}

\item{length}{A \code{grid::unit()} object that sets the length of the rug lines. Use scale expansion to avoid overplotting of data.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{borders()}.}

\item{geom}{The geometric object to use display the data}
}

\value{
    Similar output as the standard layer functions \code{ggplot2::geom_point}, \code{ggplot2::geom_rug} and \code{ggplot2::stat_point}.
}

\details{
    These layer functions are modifications of the standard layer functions \code{ggplot2::geom_point}, \code{ggplot2::geom_rug} and \code{ggplot2::stat_point}. Their usage and output are similar, except that they require a tidy kernel estimate as the input, rather than the observations themselves. For most cases, \code{geom_rug_ks} is equivalent to \code{geom_rug(stat="rug_ks")}, and likewise for \code{geom_point_ks}. 
}

\examples{library(ggplot2)
data(crabs, package="MASS")

## rug plot for tidy 1-d kernel density estimate
crabs1 <- dplyr::select(crabs, FL)
t1 <- tidy_kde(crabs1)
g1 <- ggplot(t1, aes(x=FL)) + geom_line()
g1 + geom_rug_ks(colour=4) 
g1 + geom_rug(stat="rug_ks", colour=4) ## same output

## scatter plot for tidy 2-d kernel density estimate
crabs2 <- dplyr::select(crabs, FL, CW)
t2 <- tidy_kde(crabs2)
g2 <- ggplot(t2, aes(x=FL, y=CW)) 
g2 + geom_contour_ks(colour=1) + geom_point_ks(colour=4) 
g2 + geom_contour(stat="contour_ks", colour=1) +
    geom_point(stat="point_ks", colour=4) ## same output
}

\keyword{hplot}   
