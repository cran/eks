\name{tidyst_plot}
\alias{ggplot.tidy_ks}
\alias{ggplot.sf_ks}
\alias{plot.sf_ks}
\alias{aes_ks}
\alias{guides_ks}
\alias{labs_ks}

\title{Plots for tidy and geospatial kernel estimates}
\description{
  Plots for tidy and geospatial kernel estimates.
}

\usage{
\method{ggplot}{tidy_ks}(data=NULL, mapping=aes(), ...)
\method{ggplot}{sf_ks}(data=NULL, mapping=aes(), ..., which_geometry="sf")
\method{plot}{sf_ks}(x, ...)
}

\arguments{
  \item{data,x}{object of class \code{tidy_ks} (output from \code{tidy_k*}) or object of class \code{sf_ks} (output from \code{st_k*})}
  \item{mapping}{default list of aesthetic mappings to use for plot.}
   \item{which_geometry}{type of geometry to display: one of \code{c("sf", "grid")}. Default is \code{"sf"}.}
  \item{...}{other graphics parameters. See below.
    }
}
  

\value{
  \code{ggplot} plot object is created. Base \R plot is sent to graphics window. 
}


\details{
    For \code{tidy_ks} objects, the \code{ggplot} method adds some default aesthetics based on derived variables in the computed kernel estimate. These are \code{aes(y=estimate, weight=ks)} (1-d) and are \code{aes(z=estimate, weight=ks)} (2-d). These derived variables computed in the tibble output from \code{tidy_k*} are: \code{estimate} is the kernel estimate value and \code{ks} is the untidy version of the kernel estimate, which is required to compute contour levels. The \code{ggplot} method also adds some default labels for the axes and grouping variable (via \code{labs_ks}), and some default formatting for the legends. These defaults replicate the appearance of the corresponding plots from the \pkg{ks} package.  


    For \code{sf_ks} objects, the \code{ggplot} method is similar to the above method, except no default aesthetics are added. The function header for the \code{plot} method is 
    \preformatted{
    plot(x, which_geometry="sf", percent_label=TRUE, cont=c(25,50,75), 
      abs_cont, which_deriv_ind=1, main="", legend=FALSE, ...)}
    where
    \describe{  
      \item{\code{which_geometry}}{type of geometry to display: one of \code{c("sf", "grid")}. Default is \code{"sf"}.}
      \item{\code{cont}}{vector of percentages for contour heights}
      \item{\code{abs_cont}}{vector of values for contour heights}
      \item{\code{which_deriv_ind}}{index for partial derivative for density derivative estimate. Default is 1.}
      \item{\code{main}}{main plot label. Default is "".}
      \item{\code{legend}}{flag to add legend, created by \code{mapfs::mf_legend}. Default is FALSE. ***Experimental: The output from \code{mapfs::mf_legend} in base \R is not as robust as the legend output in \code{ggplot2}.***  }
      \item{...}{other graphics parameters in the \code{plot} method for \code{sf} objects.}
      }    
      
}

\seealso{\code{\link[=st_kde]{tidy_kde,st_kde}}}


\keyword{hplot}
