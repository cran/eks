#####################################################################
## Auxiliary functions for tidy_* functions
#####################################################################

ggplot.tidy_ks <- function(data=NULL, mapping=aes(), ...)
{
    data1 <- data
    if (!is.null(data)) class(data1) <- class(data1)[class(data1)!="tidy_ks"]
    
    p <- ggplot2::ggplot(data=data1, mapping=mapping, ...)
    p <- p + aes_ks(data) + labs_ks(data) + guides_ks(data)

    return(p)
}

## default aesthetics
## object = ggplot or tidy ks object

aes_ks <- function(object)
{
    if (inherits(object, "ggplot")) object <- object$data
    d <- dim_ks(object)
    if (d==1) aesd <- aes(y=!!sym("estimate"), weight=!!sym("ks"))
    else if (d==2)
    {
   		if (all(c("u","v") %in% names(object))) aesd <- aes(u=!!sym("u"), v=!!sym("v"))
        else aesd <- aes(z=!!sym("estimate"), weight=!!sym("ks"))
    }

    return(aesd)
}

## default labels
## object = ggplot or tidy ks object

labs_ks <- function(object, x, y, ...)
{
    if (inherits(object, "ggplot")) object <- object$data
    d <- dim_ks(object)
; 
    if (inherits(object, "sf_ks"))
    {
        object <- object$tidy_ks[1,]
        oc <- head(object$tks,1)
        if (missing(x)) x <- "Longitude"
        if (missing(y)) y <- "Latitude"
    }
    else
    {
        oc <- head(object$tks,1)
        if (oc %in% "ksupp") { xnames <- names(untidy_ks(object)) }
        else if (oc %in% "kdr") 
        { 
            gv <- dplyr::group_vars(object)
            gv <- setdiff(gv, "segment")
            if (length(gv)>0) object <- dplyr::group_by(object, dplyr::across(dplyr::all_of(gv)))
            else object <- dplyr::ungroup(object)
            xnames <- untidy_ks(object, "names") 
        }
        else xnames <- untidy_ks(object, "names")
        while (is.list(xnames)) { xnames <- xnames[[1]] }

        if (is.null(xnames)) xnames <- c("x", "y")
        if (missing(x)) x <- xnames[1]
        if (missing(y)) y <- xnames[2]
    }

    if (d==1)
    {
        if (oc %in% c("kde", "kfs", "kde.boundary")) y <- "Density function"
        else if (oc %in% "kda") y <- "Weighted density"
        else if (oc %in% "kdde") y <- "Density derivative"
        else if (oc %in% "kcde") y <- "Distribution function"
        else if (oc %in% "kde.loctest") y <- expression("Density difference "*f1-f2)
    }
    if (oc %in% "kroc")
    {
        x <- expression("False positive rate"~~group("(", list(bar(specificity)), ")"))
        y <- "True positive rate (sensitivity)"
    }

    labsd <- ggplot2::labs(x=x, y=y, ...)

    return(labsd)
}

## default guides parameters

guides_ks <- function(object, kda_part=TRUE, colour_contour=TRUE)
{
    if (inherits(object, "ggplot")) object <- object$data
    d <- dim_ks(object)
    if (inherits(object, "sf_ks")) object <- object$tidy_ks[1,]
    oc <- head(object$tks,1)

    if (d==1)
    {
        titlec <- NULL
        if (oc %in% "kda") { titlec <- names(dplyr::select(object, dplyr::group_vars(object)))[1] }
        else if (oc %in% "kde.loctest") { titlec <- "Signif difference" }
        else if (oc %in% "kfs") { titlec <- "" }
        if (!is.null(titlec)) gu <- ggplot2::guides(colour=ggplot2::guide_legend(title=titlec))
        else gu <- ggplot2::guides()
    }
    else if (d==2)
    {
        titlef <- NULL; titlec <- NULL; rev <- TRUE
        if (oc %in% c("kde", "kdcde", "kde.boundary", "kde.truncate", "kde.sp", "kde.balloon")) titlef <- "Density"
        else if (oc %in% "kcurv") titlef <- "Density curvature"
        else if (oc %in% "kda")
        {
            if (kda_part)
            {
                titlef <- "Classif\nlabel"; titlec <- "Density"; rev <- FALSE
            }
            else titlef <- "Density"
        }
        else if (oc %in% "kdde") titlef <- "Density derivative"
        else if (oc %in% c("kcde", "kcopula")) titlef <- "Distribution"
        else if (oc %in% "kde.loctest") { titlef <- "Signif difference" }
        else if (oc %in% "kms") { titlef <- "Cluster"; rev <- FALSE }
        else if (oc %in% c("kdr", "kfs")) { titlef <- ""; titlec <- ""; rev <- FALSE }
        else if (oc %in% "ksupp") titlef <- expression(atop("Support", "convex hull"))
        if (colour_contour) titlec <- titlef

        if ((titlef=="" & titlec=="") ) gu <- ggplot2::guides(fill=ggplot2::guide_legend(title=NULL, reverse=TRUE & rev), colour=ggplot2::guide_legend(title=NULL, reverse=colour_contour & rev))
        else if (!is.null(titlef) & !is.null(titlec)) gu <- ggplot2::guides(fill=ggplot2::guide_legend(title=titlef, reverse=TRUE & rev), colour=ggplot2::guide_legend(title=titlec, reverse=colour_contour & rev))
        else if (!is.null(titlef) & is.null(titlec)) gu <- ggplot2::guides(fill=ggplot2::guide_legend(title=titlef, reverse=TRUE & rev))
        else if (is.null(titlef) & !is.null(titlec)) gu <- ggplot2::guides(colour=ggplot2::guide_legend(title=titlec, reverse=colour_contour & rev))
        else gu <- ggplot2::guides()
    }

    return(gu)
}

## make individual colours in colour scale to be transparent 

scale_transparent <- function(x, ind=NULL)
{  
    Scale <- ggplot2::ggproto("Scale", x, 
        palette = function (self, n)
        { 
            pal <- x$palette(n)
            if (is.null(ind)) ind <- round(stats::median(seq_len(length(pal))),0)
            pal[ind] <- "transparent"
            return(pal)
        }
    )

    return(Scale)
}

## dimension of tidy_ks/sf_ks objects

dim_ks <- function(object)
{
    if (inherits(object, "sf_ks"))
    {
        object <- object$tidy_ks[1,]
        if (object$tks %in% "kda") d <- ncol(untidy_ks(object)$H[[1]])
        else if (object$tks %in% "ksupp") d <- ncol(untidy_ks(object))
        else d <- ncol(get_data_ks(object))
    }
    else
        d <- unlist(tail(getElement(dplyr::ungroup(object), "ks"),n=1))

    return(as.integer(d))
}

## move group_vars in x to last column
## group_vars can be those in y  

move_group_vars <- function(x, y)
{
    if (missing(y)) y <- x
    if (dplyr::is_grouped_df(y))
    {
        gvy <- dplyr::group_vars(y)
        x <- dplyr::group_by(x, dplyr::across(dplyr::all_of(gvy)))
        x <- dplyr::relocate(x, dplyr::all_of(gvy), .after=dplyr::last_col())
    }

    return(x)
}

## extract ks object from tidy_ks object
## object = tidy_ks object

untidy_ks <- function(object, name, collapse=FALSE, var_ks="ks")
{
    ge <- dplyr::pull(dplyr::slice_head(object), var_ks)
    n <- length(ge)

    if (n==1)
    {
        ge <- ge[[1]]
        if (!missing(name)) ge <- getElement(ge, name)
    }
    else
    {
        if (!missing(name)) ge <- lapply(ge, getElement, name)
    }

    return(ge)
}

## extract data matrix from ks object from tidy_ks object
## object = tidy_ks object

get_data_ks <- function(object, var_ks="ks")
{
    ge <- untidy_ks(object, "x", var_ks=var_ks)
    if (is.null(ge)) ge <- rbind(untidy_ks(object, "fhat1", var_ks=var_ks)$x, untidy_ks(object, "fhat2", var_ks=var_ks)$x)
    oc <- class(getElement(object, var_ks)[[1]])
    
    if (oc %in% "kda")
    {
        gr <- unique(object$group)
        ge <- data.frame(ge[[gr]])
        colnames(ge) <- c("x","y","z")[1:ncol(ge)]
        ge <- data.frame(ge, group=gr)
    }
    else
    {
        if (!is.list(ge))
        {
            ge <- data.frame(ge)
            colnames(ge) <- c("x","y","z")[1:ncol(ge)]
        }
        else
        {
            gr <- rep(unique(object$group), lengths(ge))
            d <- unlist(tail(getElement(dplyr::ungroup(object), "ks"),n=1))
            if (d==1) ge <- data.frame(do.call("c",ge))
            else ge <- do.call("rbind",ge)
            colnames(ge) <- c("x","y","z")[1:ncol(ge)]
            rownames(ge) <- NULL
            ge <- data.frame(ge, group=gr)
        }
    }

    return(ge)
}

## renaming function

rename_ks <- function(data, gg, d)
{
    gv <- dplyr::group_vars(data)
    if (length(gv)>0) nd <- names(dplyr::select(dplyr::ungroup(data), -dplyr::all_of(gv)))
    else nd <- names(data)

    if (d>=1) names(gg)[names(gg)=="x"] <- nd[1]
    if (d>=2) names(gg)[names(gg)=="y"] <- nd[2]

    return(gg)
}

## cast data as "tidy_ks" object

as_tidy_ks <- function(data)
{
    class(data) <- unique(c("tidy_ks", class(data)))

    return(data)
}

## concatenate tidy_ks objects

c.tidy_ks <- function (..., labels) { c_tidy_ks(..., labels=labels) } 

c_tidy_ks <- function(..., labels)
{
    args <- list(...)
    nargs <- length(args)
    if (nargs>=1) 
    {
        tt <- dplyr::mutate(args[[1]], group=1L)
        if (nargs>=2)
        {
            for (i in 2:nargs) 
                tt <- dplyr::bind_rows(tt, dplyr::mutate(args[[i]], group=i))
        }
        if (missing(labels)) tt$group <- factor(tt$group)
        else tt$group <- factor(tt$group, labels=labels)
        tt <- dplyr::group_by(tt, .data$group)
        tt <- as_tidy_ks(tt)
    }
    else tt <- NULL

    return(tt)
}

#####################################################################
## Probability contour regions
#####################################################################

## compute contour breaks for ks objects
## data = tidy ks object

contour_breaks <- function(data, cont=c(25,50,75), group=FALSE) 
{ 
    cont.aug <- cont[cont>=1 & cont<=99]
    cb <- ks::contourLevels(x=data, cont=cont.aug, group=TRUE)
    cb <- dplyr::slice(cb, 1:(dplyr::n()-1)) 
    
    if (!group) 
    {
        oc <- class(data)       
        cb <- dplyr::arrange(cb, .data$breaks)
        cb <- dplyr::ungroup(cb) 
        cb <- dplyr::distinct(cb, .data$breaks)
           
        cont.prob.diff <- function(cb)
        {
            if (any(oc %in% "sf_ks")) data2 <- data$tidy_ks
            else data2 <- dplyr::slice_head(data, n=1) 
            cb.mat <- matrix(0, ncol=nrow(data2), nrow=nrow(cb))
            for (i in 1:ncol(cb.mat)) 
            {
                if (inherits(untidy_ks(data2[i,]), "kda"))
                {
                    kdei <- untidy_ks(data2[i,])
                    kdei$x <- kdei$x[[i]]
                    kdei$estimate <- kdei$estimate[[i]]
                    kdei$H <- kdei$H[[i]]
                    kdei$cont <- kdei$cont[[i]]   
                    kdei$w <- kdei$w[[i]]
                    class(kdei) <- "kde"
                    data2[i,]$ks <- list(kdei)
                }
                if (inherits(untidy_ks(data2[i,]), "kdde")) cb.mat[,i] <- cb$breaks
                else cb.mat[,i] <- ks::contourProbs(untidy_ks(data2[i,]), abs.cont=as.vector(unlist(cb)))
            }
            return(as.vector(cb.mat))
        }
        cb2 <- cont.prob.diff(cb)
        cont.cutoff <- abs(head(diff(cont.aug),n=1))
        if (any(oc %in% "sf_ks")) data2 <- data$tidy_ks else data2 <- data
        if (any(data2$tks %in% "kdde")) {cb2 <- (cb2-min(cb2))/(max(cb2)-min(cb2))} 
        clust <- cutree(hclust(dist(cb2), method="complete"), h=cont.cutoff/100)
        cb <- data.frame(breaks=rep(cb$breaks, times=length(cb2)/nrow(cb)))
        cb$contgr <- clust
        cb <- dplyr::distinct(cb)
        cb <- dplyr::group_by(cb, .data$contgr)
        cb.summary <- dplyr::summarise(cb, min.breaks=min(.data$breaks), mean.breaks=mean(.data$breaks), max.breaks=max(.data$breaks))
        cb.summary <- dplyr::mutate(cb.summary, breaks=ifelse(dplyr::row_number()==1, .data$min.breaks, ifelse(dplyr::row_number()==dplyr::n(), .data$max.breaks, .data$mean.breaks)))
        cb.summary <- dplyr::select(cb.summary, dplyr::all_of(c("contgr", "breaks"))) 
        cb.summary <- dplyr::arrange(cb.summary, .data$breaks)   
    }
    else cb.summary <- cb

    return(cb.summary)    
}

## contour levels method for tidy_ks opjects

contourLevels.tidy_ks <- function(x, cont=c(25,50,75), group=FALSE)
{
    data <- x
    oc <- head(x$tks,1) 
    data$.group <- dplyr::group_indices(data)
    gg <- dplyr::slice_head(data)
    
    if (oc %in% "kdde")
        gg <- dplyr::group_modify(.data=gg, .f=~dplyr::tibble(breaks=.contourLevels(untidy_ks(.x), cont=cont, which.deriv.ind=unique(.x$deriv.ind), approx.cont=TRUE)))
    else if (oc %in% "kda")
        gg <- dplyr::group_modify(.data=gg, .f=~dplyr::tibble(breaks=.contourLevels(untidy_ks(.x), cont=cont, approx.cont=TRUE)[[.x$.group]]))
    else
        gg <- dplyr::group_modify(.data=gg, .f=~dplyr::tibble(breaks=.contourLevels(untidy_ks(.x), cont=cont, approx.cont=TRUE)))

    gg <- dplyr::mutate(gg, breaks=ifelse(dplyr::row_number()==dplyr::n(), max(gg$breaks), .data$breaks))
    if (!group) {
        gg <- dplyr::mutate(gg, group_breaks=dplyr::row_number())
        gg <- dplyr::group_by(gg, .data$group_breaks)
        ggbreaks <- dplyr::summarise(gg, min=min(.data$breaks), max=max(.data$breaks), .groups="drop_last")
        ggbreaks <- dplyr::transmute(ggbreaks, breaks=ifelse(dplyr::row_number()==1, min, ifelse(max<0 & dplyr::row_number()<dplyr::n(), min, max)))

    }
    else ggbreaks <- gg

    return(ggbreaks)

}

## fhat = untidy ks object

.contourLevels <- function(fhat, cont, which.deriv.ind=1, approx.cont=TRUE)
{
    cont <- sort(cont)
    oc <- class(fhat)
    if (oc %in% "kdde")
    {
        hts <- ks::contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont, which.deriv.ind=fhat$deriv.ind.scalar)
    }
    else if (oc %in% "kdecurv")
    {
        class(fhat) <- "kde"
        fhat.temp <- fhat
        fhat.temp$x <- fhat$x[predict(fhat, x=fhat$x) > 0, ]
        hts <- ks::contourLevels(fhat.temp, prob=(100-cont)/100, approx=approx.cont)
    }
    else if (oc %in% "kcde")
    {
        hts <- sort(cont/100)
        names(hts) <- paste0(cont,"%")
    }
    else if (oc %in% "kcopula")
    {
        hts <- sort(unique(c(cont,100))/100)
        names(hts) <- paste0(cont,"%")
    }
    else if (oc %in% "kde.loctest")
        hts <- c(-10,-0.5, 0.5, 10)
    else if (oc %in% "kfs")
        hts <- c(0.5, 10)
    else
        hts <- ks::contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)

    if (oc %in% c("kde", "kdecurv", "kcde"))
    {
        hts <- sort(hts)
        hts <- c(hts, max(c(fhat$estimate, hts)) + 0.01 * max(abs(fhat$estimate)))
    }
    else if (oc %in% "kda")
    {
        hts <- lapply(hts, sort)
        pp <- fhat$prior.prob
        for (i in 1:length(hts))
             hts[[i]] <- c(pp[i]*hts[[i]], max(c(fhat$estimate[[i]], pp[i]*hts[[i]])) + 0.01 * max(abs(fhat$estimate[[i]])))
    }
    else if (oc %in% "kdde")
    {
        hts1 <- hts[1,]; hts2 <- hts[2,]
        names(hts1) <- paste0("-", names(hts1))
        hts <- c(sort(hts1), sort(hts2))
        hts <- c(min(c(sapply(fhat$estimate, min), hts)) - 0.01 * max(abs(sapply(fhat$estimate, max))), hts, max(c(sapply(fhat$estimate, max), hts)) + 0.01 * max(abs(sapply(fhat$estimate, max))))
    }

    return (hts)
}

## local copies of unexported functions from ggplot2

.ggplot2_new_data_frame <- function(x=list(), n=NULL)
{
    if (length(x) != 0 && is.null(names(x))) { stop("Elements must be named") }
    lengths <- vapply(x, length, integer(1))
    if (is.null(n)) 
    {
        n <- if (length(x)==0 || min(lengths)==0) 0
        else max(lengths)
    }
    for (i in seq_along(x)) 
    {
        if (lengths[i] == n) next
        if (lengths[i] != 1) { stop("Elements must equal the number of rows or 1") }
        x[[i]] <- rep(x[[i]], n)
    }
    class(x) <- "data.frame"
    attr(x, "row.names") <- .set_row_names(n)
    
    return(x)
}

.ggplot2_xyz_to_isolines <- function(data, breaks)
{
    isoband::isolines(x=sort(unique(data$x)), y=sort(unique(data$y)), z=.ggplot2_isoband_z_matrix(data), levels=breaks)
}

.ggplot2_xyz_to_isobands <- function (data, breaks)
{
    isoband::isobands(x=sort(unique(data$x)), y=sort(unique(data$y)),
        z=.ggplot2_isoband_z_matrix(data), levels_low=breaks[-length(breaks)],
        levels_high=breaks[-1])
}

.ggplot2_isoband_z_matrix <- function (data)
{
    x_pos <- as.integer(factor(data$x, levels=sort(unique(data$x))))
    y_pos <- as.integer(factor(data$y, levels=sort(unique(data$y))))
    nrow <- max(y_pos)
    ncol <- max(x_pos)
    raster <- matrix(NA_real_, nrow=nrow, ncol=ncol)
    raster[cbind(y_pos, x_pos)] <- data$z
    
    return(raster)
}

.ggplot2_iso_to_path <- function (iso, group=1)
{
    lengths <- vapply(iso, function(x) length(x$x), integer(1))
    if (all(lengths == 0)) 
    {
        warning("stat_contour(): Zero contours were generated")
        return(.ggplot2_new_data_frame())
    }
    levels <- names(iso)
    xs <- unlist(lapply(iso, "[[", "x"), use.names=FALSE)
    ys <- unlist(lapply(iso, "[[", "y"), use.names=FALSE)
    ids <- unlist(lapply(iso, "[[", "id"), use.names=FALSE)
    item_id <- rep(seq_along(iso), lengths)
    groups <- paste(group, sprintf("%03d", item_id), sprintf("%03d", ids), sep="-")
    groups <- factor(groups)
    .ggplot2_new_data_frame(list(level=rep(levels, lengths), x=xs,
        y=ys, piece=as.integer(groups), group=groups),
        n=length(xs))
}

.ggplot2_iso_to_polygon <- function (iso, group=1)
{
    lengths <- vapply(iso, function(x) length(x$x), integer(1))
    if (all(lengths == 0)) {
        warning("stat_contour(): Zero contours were generated")
        return(.ggplot2_new_data_frame())
    }
    levels <- names(iso)
    xs <- unlist(lapply(iso, "[[", "x"), use.names=FALSE)
    ys <- unlist(lapply(iso, "[[", "y"), use.names=FALSE)
    ids <- unlist(lapply(iso, "[[", "id"), use.names=FALSE)
    item_id <- rep(seq_along(iso), lengths)
    groups <- paste(group, sprintf("%03d", item_id), sep="-")
    groups <- factor(groups)
    .ggplot2_new_data_frame(list(level=rep(levels, lengths), x=xs,
        y=ys, piece=as.integer(groups), group=groups, subgroup=ids),
        n=length(xs))
}

#####################################################################
## stat_* and geom_* auxiliary functions for ggplot2 layers
#####################################################################

## geom and stat function for extracting data points from tidy_ks object
## required for adding scatter plots of original data in geom_point and geom_rug

geom_point_ks <- function(mapping=NULL, data=NULL, stat="point_ks", position="identity", ..., na.rm=FALSE, jitter=FALSE, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=stat, geom=GeomPointKs, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(na.rm=na.rm, jitter=jitter, ...))
}

GeomPointKs <- ggplot2::ggproto("GeomPointKs", GeomPoint)

stat_point_ks <- function(mapping=NULL, data=NULL, geom="point_ks", position="identity", ..., na.rm=FALSE, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=StatPointKs, geom=geom, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(na.rm=na.rm, jitter=jitter, ...))
}

StatPointKs <- ggplot2::ggproto("StatPointKs", ggplot2::Stat,
compute_group=function(self, data, scales, jitter=FALSE)
{
    x <- get_data_ks(data, var_ks="weight")
    if (jitter) x <- jitter(as.matrix(x))
    return(x)
},
required_aes=c("x","y"),
dropped_aes=c("z", "weight")
)

stat_rug_ks <- function(mapping=NULL, data=NULL, geom="rug_ks", position="identity", ..., na.rm=FALSE, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=StatRugKs, geom=geom, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(na.rm=na.rm, jitter=jitter, ...))
}

StatRugKs <- ggplot2::ggproto("StatPointKs", ggplot2::Stat,
compute_group=function(self, data, scales, jitter=FALSE)
{
    x <- get_data_ks(data, var_ks="weight")
    if (jitter) x <- jitter(as.matrix(x))
    return(x)
},
required_aes=c("x"),
dropped_aes=c("y", "z", "weight") 
)

geom_rug_ks <- function(mapping=NULL, data=NULL, stat="rug_ks", position="identity", ..., outside=FALSE, sides="bl", length=unit(0.03, "npc"), na.rm=FALSE, jitter=FALSE, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=stat, geom=GeomRugKs, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(outside=outside, sides=sides, length=length, na.rm=na.rm, jitter=jitter, ...))
}

GeomRugKs <- ggplot2::ggproto("GeomRug", GeomRug)

## stat and geom functions for computing probability contours
## replace ggplot2::geom_contour and ggplot2::stat_contour

geom_contour_ks <- function(mapping=NULL, data=NULL, stat="contour_ks", position="identity", ..., cont=c(25,50,75), label_percent=NULL, breaks=NULL, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=stat, geom=GeomContourKs, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(cont=cont, label_percent=label_percent, breaks=breaks, ...))
}

GeomContourKs <- ggplot2::ggproto("GeomContourFilled", GeomContour)

stat_contour_ks <- function(mapping=NULL, data=NULL, geom="contour_ks", position="identity", ..., cont=c(25,50,75), label_percent=NULL, breaks=NULL, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=StatContourKs, geom=geom, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(cont=cont, label_percent=label_percent, breaks=breaks, ...))
}

## replaces ggplot2::StatContour
## data = tidy_ks object

StatContourKs <- ggplot2::ggproto("StatContourKs", ggplot2::StatContour,
dropped_aes="weight",
compute_group=function (data, scales, cont=c(25,50,75), label_percent=NULL, breaks, bins=NULL, binwidth=NULL)
{
    fhat <- untidy_ks(data, var_ks="weight")
    oc <- class(fhat)
    if (is.null(label_percent)) label_percent <- TRUE

    ## compute ks probability contour levels
    if (!is.null(cont) & is.null(breaks))
    {
        breaks <- .contourLevels(fhat=fhat, cont=cont)

        if (oc %in% "kcopula") names(breaks)[length(breaks)] <- "100%"
        if (oc %in% "kde.loctest")
        {
            fhat.diff.lab <- fhat$fhat1
            fhat.diff.lab$estimate <- fhat$fhat.diff.pos$estimate - fhat$fhat.diff.neg$estimate
            d <- ncol(fhat$fhat1$H)
            diff.lab <- compute.tidy(fhat.diff.lab, d=d, tidy=FALSE)$estimate
            data$z <- diff.lab
            breaks <- c(-10,-0.5,0.5,10)
        }

        if (is.list(breaks)) breaks <- breaks[[unique(data$group)]]
        isolines <- .ggplot2_xyz_to_isolines(data, breaks)
        path_df <- .ggplot2_iso_to_path(isolines, data$group[1])
        path_df$level <- as.numeric(path_df$level)

        if (oc %in% "kdde")
        {
            if (label_percent)
            {
                nc <- length(cont)
                names.isolines <- paste0(c(-100 - as.numeric(gsub("%", "",tail(head(names(breaks), n=nc+1),n=-1))), 100 - as.numeric(gsub("%", "",head(tail(names(breaks), n=nc+1),n=-1)))), "%")
            }
            else
                names.isolines <- format(head(tail(breaks, n=-1),n=-1), digits=3, trim=TRUE)
            path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
        else if (oc %in% "kde.loctest")
        {
            names.isolines <- c("-1", "1")
            path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
        else if (oc %in% "kcopula")
        {
            path_df$level <- findInterval(path_df$level, as.numeric(names(table(path_df$level))))
            names.isolines <- 100-as.numeric(gsub("%", "",names(breaks), "%"))
            path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
        else
        {
            path_df$level <- findInterval(path_df$level, as.numeric(names(table(path_df$level))))
            names.isolines <- 100-as.numeric(gsub("%", "", head(names(breaks), n=-1), "%"))
            path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
    }
    else
    {
        if (!is.vector(breaks))
        {
            if ("group" %in% names(breaks))
            {
                breaks_ind <- breaks$group %in% levels(breaks$group)[unique(data$group)]
                breaks_orig <- breaks
                breaks <- breaks_orig$breaks[breaks_ind]

                if (is.list(breaks)) breaks <- breaks[[unique(data$group)]]
            }
            else breaks <- breaks$breaks
        }

        breaks <- sort(breaks)
        isolines <- .ggplot2_xyz_to_isolines(data, breaks)
        path_df <- .ggplot2_iso_to_path(isolines, data$group[1])
        path_df$level <- signif(as.numeric(path_df$level),3)
        path_df$level <- format(path_df$level, digits=3, trim=TRUE)
    }
    return(path_df)
},
required_aes=c("x","y","z"),
dropped_aes=c("z","weight") 
)

## stat and geom functions for computing filled probability contours
## replace ggplot2::geom_contour_filled & ggplot2::stat_contour_filled

geom_contour_filled_ks <- function(mapping=NULL, data=NULL, stat="contour_filled_ks", position="identity", ..., cont=c(25,50,75), label_percent=NULL, breaks=NULL, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=stat, geom=GeomContourFilled, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(cont=cont, label_percent=label_percent, breaks=breaks, ...))
}

GeomContourFilledKs <- ggplot2::ggproto("GeomContourFilled", GeomContourFilled)

## replace ggplot2::stat_contour_filled & ggplot2::StatContourFilled

stat_contour_filled_ks <- function(mapping=NULL, data=NULL, geom="contour_filled_ks", position="identity", ..., cont=c(25,50,75), label_percent=NULL, breaks=NULL, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=StatContourFilledKs, geom=geom, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(breaks=breaks, na.rm=FALSE, ...))
}

StatContourFilledKs <- ggplot2::ggproto("StatContourFilledKs", StatContourFilled,
compute_group=function(data, scales, cont=c(25,50,75), label_percent=NULL, breaks, bins=NULL, binwidth=NULL)
{
    fhat <- untidy_ks(data, var_ks="weight")
    oc <- class(fhat)
    if (is.null(label_percent)) label_percent <- TRUE

    ## compute ks probability contour levels
    if (!is.null(cont) & is.null(breaks))
    {
        if (oc %in% "kde.loctest")
        {
            d <- ncol(fhat$fhat1$H)
            fhat.diff.lab <- fhat$fhat1
            fhat.diff.lab$estimate <- fhat$fhat.diff.pos$estimate - fhat$fhat.diff.neg$estimate
            diff.lab <- compute.tidy(fhat.diff.lab, d=d, tidy=FALSE)$estimate
            data$z <- diff.lab
        }
        else if (oc %in% "kfs")
        {
            d <- ncol(fhat$H)
            data$z <- compute.tidy(fhat, d=d, tidy=FALSE)$estimate
        }

        breaks <- .contourLevels(fhat=fhat, cont=cont)
        if (oc %in% "kda") breaks <- breaks[[unique(data$PANEL)]]
        isobands <- .ggplot2_xyz_to_isobands(data, breaks)

        if (oc %in% "kdde")
        {
            if (label_percent)
            {
                nc <- length(cont)
                names(isobands) <- paste0(c(-100 - as.numeric(gsub("%", "",tail(head(names(breaks), n=nc+1),n=-1))),0, 100 - as.numeric(gsub("%", "",head(tail(names(breaks), n=nc+1),n=-1)))), "%")
            }
            else
                 names(isobands) <- format(head(breaks, n=-1), digits=3, trim=TRUE)
        }
        else if (oc %in% "kde.loctest")
        {
            label_percent <- FALSE
            labels <- na.omit(data$weight[[1]]$label)
            names(isobands) <- c(labels[1], "", labels[2]) 
        }
        else if (oc %in% "kfs")
        {
            label_percent <- FALSE
            names(isobands) <- c("Signif curv")
        }
        else
        {
            if (label_percent) names(isobands) <- paste0(100-as.numeric(gsub("%", "", head(names(breaks), n=-1), "%")),"%")
            else names(isobands) <- format(head(breaks, n=-1), digits=3, trim=TRUE)
        }
    }
    else
    {
        if (!is.vector(breaks))
        {
            if ("group" %in% names(breaks))
            {
                breaks_ind <- breaks$group %in% levels(breaks$group)[unique(data$group)]
                breaks_orig <- breaks
                breaks <- breaks_orig$breaks[breaks_ind]

                if (is.list(breaks)) breaks <- breaks[[unique(data$group)]]
            }
            else breaks <- breaks$breaks
        }

        if ((oc %in% c("kdde","kda")))
        {
            breaks <- c(sort(breaks), max(c(sapply(fhat$estimate,max), breaks)) + 0.01 * max(abs(sapply(fhat$estimate,max))))
        }
        else
           breaks <- c(sort(breaks), max(c(fhat$estimate, breaks)) + 0.01 * max(abs(fhat$estimate)))
        isobands <- .ggplot2_xyz_to_isobands(data, breaks)
        names(isobands) <- format(signif(as.numeric(format(head(breaks, n=-1), dig.lab=3, trim=TRUE)),3), dig.lab=3, trim=TRUE)
        if (exists("breaks_orig")) breaks_isobands <- format(unique(as.vector(sort(c(breaks_orig$breaks, max(breaks))))), digits=3, trim=TRUE)
    }

    ## local copy of unexported function from scales
    .scales_rescale_max <- function (x, to=c(0, 1), from=range(x, na.rm=TRUE)) { x/from[2] * to[2] }

    path_df <- .ggplot2_iso_to_polygon(isobands, data$group[1])
    path_df$level <- ordered(path_df$level, levels=names(isobands))
    path_df$level_low <- breaks[as.numeric(path_df$level)]
    path_df$level_high <- breaks[as.numeric(path_df$level) + 1]
    path_df$level_mid <- 0.5 * (path_df$level_low + path_df$level_high)
    path_df$nlevel <- .scales_rescale_max(path_df$level_high)

    if (exists("breaks_isobands"))
        path_df$level <- ordered(path_df$level, levels=breaks_isobands)

    return(path_df)
},
required_aes=c("x","y","z"),
dropped_aes=c("z","weight") 
)
