\name{scale_transparent}
\alias{scale_transparent}
\title{Change individual colours in discrete colour scale to transparent}

\description{
    Change individual colours in discrete colour scale to transparent.
}

\usage{
scale_transparent(x, ind=NULL)
}

\arguments{
    \item{x}{discrete colour scale}
    \item{ind}{index of colour scale to change to transparent. Default is median.}
}

\value{
    The output is the same colour scale, except that the colours at the indices enumerated by \code{ind} are changed to transparent.  
}

\examples{
## see ? tidy_kdde
}

\keyword{graphics}
