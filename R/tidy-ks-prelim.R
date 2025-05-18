#####################################################################
## Auxiliary functions for tidy_* functions
#####################################################################

ggplot.tidy_ks <- function(data=NULL, mapping=aes(), ...)
{
    data1 <- data
    if (!is.null(data)) 
    {
        class(data1) <- class(data1)[class(data1)!="tidy_ks"]
        transp_neutral <- any(object_class(data, type="ks") %in% "kdde")
    }
    else transp_neutral <- FALSE
    p <- ggplot2::ggplot(data=data1, mapping=mapping, ...)
    p <- p + aes_ks(data) + labs_ks(data) + guides_ks(data)
    p <- p + scale_colour_ks(data) + scale_fill_ks(data, transp_neutral=transp_neutral)
    
    return(p)
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
c.tidy_ks <- function (..., add_group=TRUE, labels) 
{
    args <- list(...)
    nargs <- length(args)
    
    if (nargs>=1) 
    {
        if (add_group)
        {
        
            tt <- dplyr::mutate(args[[1]], group=1L)
            if (nargs>=2)
            {
                for (i in 2:nargs) 
                    tt <- dplyr::bind_rows(tt, dplyr::mutate(args[[i]], group=i))
            }
            if (any(names(tt) %in% "group"))
            {
                if (missing(labels)) tt$group <- factor(tt$group)
                else tt$group <- factor(tt$group, labels=labels)
                tt <- dplyr::group_by(tt, .data$group)
            }
            tt <- as_tidy_ks(tt)
        }
        else 
        {   
            tt <- dplyr::bind_rows(...)
            tt <- as_tidy_ks(tt)
        }
    }
    else tt <- NULL

    return(tt)
}

rbind.tidy_ks <- c.tidy_ks

#####################################################################
## Probability contour regions
#####################################################################
## compute contour breaks for ks objects
## data = tidy ks/sf_ks object
contour_breaks <- function(data, cont=c(25,50,75), n=3, group=FALSE, type="density") 
{ 
    cont.aug <- cont[cont>=1 & cont<=99]
    oct <- object_class(data, type="tks")
    type <- match.arg(type, c("density", "length", "quantile", "natural"))

    ## sf_ks object
    if (inherits(data, "sf_ks")) 
    {
        data$tidy_ks <- as_tidy_ks(data$tidy_ks)
        is.kdde <- oct %in% "kdde"
       
        if (type=="density") 
        {
            ## contour levels only w/o min, w/o max 
            cb <- contourLevels(x=data, cont=cont.aug)
            if (is.kdde) cb <- dplyr::group_modify(cb, .f=~head(tail(.x, n=-1), n=-1))
            else cb <- dplyr::group_modify(cb, .f=~head(.x, n=-1))
        }
        else 
        {
            ## contour levels, kde: with min, w/o max, kdde: w/o min, w/o max
            cb <- dplyr::group_modify(data$grid, .f=~data.frame(breaks=ks::.contourBreaks(.x, cont=cont.aug, n=n, type=type, is.kdde=is.kdde)))
            if (is.kdde) cb <- dplyr::group_modify(cb, .f=~head(tail(.x, n=-1), n=-1))
        }
    }
    ## tidy_ks object
    else
    {
        data <- as_tidy_ks(data) 
        is.kdde <- oct %in% "kdde"
       
        if (type=="density") 
        {
            cb <- ks::contourLevels(x=data, cont=cont.aug)
            cb <- dplyr::group_modify(cb, .f=~head(.x, n=-1))
        }
        else 
        {
            cb <- dplyr::group_modify(data, .f=~data.frame(breaks=ks::.contourBreaks(.x, cont=cont.aug, n=n, type=type, is.kdde=is.kdde)))
            if (is.kdde) cb <- dplyr::group_modify(cb, .f=~head(.x, n=-1))
        }    
    }
    
    if (group)
    {
        cb.summary <- dplyr::mutate(cb, contgr=dplyr::row_number(), .before=1)
        cb.summary <- dplyr::select(cb.summary, dplyr::all_of(c("contgr", "breaks", dplyr::group_vars(cb.summary))))
        if (inherits(data, "sf_ks")) data2 <- sf::st_drop_geometry(data$sf) else data2 <- data
        if (is.kdde & any(names(data2) %in% "deriv_ind")) 
        { 
            cb.summary <- dplyr::left_join(cb.summary, dplyr::distinct(dplyr::ungroup(data2)[,c("deriv_ind", "deriv_group")]), by="deriv_group")
            cb.summary <- dplyr::relocate(cb.summary, "contgr", "breaks", .before=1)
        }
        cb.summary <- dplyr::arrange(cb.summary, .data$breaks, .by_group=TRUE)  
    }
    else 
    {    
        cb <- dplyr::arrange(cb, .data$breaks)
        cb <- dplyr::ungroup(cb) 
        cb <- dplyr::distinct(cb, .data$breaks)

        cont.prob.diff <- function(cb)
        {
            if (inherits(data, "sf_ks")) data2 <- data$tidy_ks
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

        if (nrow(cb)>2) 
        {
            ## always keep min value
            if (inherits(data, "tidy_ks")) 
            {
                cb1 <- head(cb, n=1); cb1$contgr <- 0;
                cb <- tail(cb, n=-1)
            }
            cb2 <- cont.prob.diff(cb)
            cont.cutoff <- abs(head(diff(cont.aug),n=1))
            if (inherits(data, "sf_ks")) data2 <- data$tidy_ks else data2 <- data
            if (any(data2$tks %in% "kdde")) {cb2 <- (cb2-min(cb2))/(max(cb2)-min(cb2))} 
            clust <- cutree(hclust(dist(cb2), method="complete"), h=cont.cutoff/100)
            if (length(clust)==0) clust <- 1:length(cb2)
            cb <- data.frame(breaks=rep(cb$breaks, times=length(cb2)/nrow(cb)))
            cb$contgr <- clust
            cb <- dplyr::distinct(cb)
            if (inherits(data, "tidy_ks")) { cb <- rbind(cb1, cb); cb$contgr <- cb$contgr+1 }
        }
        else 
            cb$contgr <- 1:nrow(cb)

        cb <- dplyr::group_by(cb, .data$contgr)
        cb.summary <- dplyr::summarise(cb, min.breaks=min(.data$breaks), mean.breaks=mean(.data$breaks), max.breaks=max(.data$breaks))
        cb.summary <- dplyr::mutate(cb.summary, breaks=ifelse(dplyr::row_number()==1, .data$min.breaks, ifelse(dplyr::row_number()==dplyr::n(), .data$max.breaks, .data$mean.breaks)))
        cb.summary <- dplyr::select(cb.summary, dplyr::all_of(c("contgr", "breaks"))) 
        cb.summary <- dplyr::arrange(cb.summary, .data$breaks) 
        cb.summary <- dplyr::distinct(cb.summary, .data$breaks, .keep_all=TRUE)
        cb.summary$contgr <- 1:nrow(cb.summary)  
    }
    
    return(cb.summary)    
}

## contour levels method for tidy_ks opjects
contourLevels.tidy_ks <- function(x, ..., cont=c(25,50,75))
{
    oct <- object_class(x, type="tks")
    if (oct %in% "kcurv") 
    { 
        ftemp <- x$ks[lengths(x$ks)>1]
        ftemp <- lapply(ftemp, function(.) { class(.) <- "kdecurv"; . })
        x$ks[lengths(x$ks)>1] <- ftemp
    }
    data <- x
    data$.group <- dplyr::group_indices(data)
    gg <- dplyr::slice_head(data)
    add_minmax <- TRUE 

    if (oct %in% "kdde")
        gg <- dplyr::group_modify(.data=gg, .f=~{ breaks <- .contourLevels(untidy_ks(.x), cont=cont, which.deriv.ind=unique(.x$deriv.ind), add_minmax=add_minmax); data.frame(breaks, cont=as.numeric(names(breaks))) })
    else if (oct %in% "kda")
        gg <- dplyr::group_modify(.data=gg, .f=~{ breaks <- .contourLevels(untidy_ks(.x), cont=cont, add_minmax=add_minmax)[[.x$.group]]; data.frame(breaks, cont=as.numeric(names(breaks))) })
    else
        gg <- dplyr::group_modify(.data=gg, .f=~{ breaks <- .contourLevels(untidy_ks(.x), cont=cont, add_minmax=add_minmax); data.frame(breaks, cont=as.numeric(names(breaks))) })
    
    return(gg)
}

## fhat = untidy ks object
## add_minmax = add min and/or max values of fhat$estimate to output 
.contourLevels <- function(fhat, cont, which.deriv.ind=1, approx.cont=TRUE, add_minmax=TRUE)
{
    cont <- sort(cont)
    oc <- head(class(fhat),n=1) ## equivalent to object_class(type="ks")
    
    if (oc %in% c("kdde","kqdde"))
    {
        hts <- ks::contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont, which.deriv.ind=fhat$deriv.ind.scalar)
    }
    else if (oc %in% "kdecurv")
    {
        class(fhat) <- "kde"
        fhat.temp <- fhat
        fhat.temp$x <- fhat$x[predict(fhat, x=fhat$x)>0, ]
        hts <- ks::contourLevels(fhat.temp, prob=(100-cont)/100, approx=approx.cont)
    }
    else if (oc %in% c("kcde", "kcopula"))
    {
        hts <- sort(cont/100)
        names(hts) <- paste0(cont,"%")
    }
    else if (oc %in% "kde.loctest")
        hts <- c(-10,-0.5, 0.5, 10)
    else if (oc %in% "kfs")
        hts <- c(0.5, 10)
    else
        hts <- ks::contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)

    if (oc %in% "kda")
    {
        hts <- lapply(hts, sort)
        pp <- fhat$prior.prob
        if (add_minmax) 
        {
            for (i in 1:length(hts)) 
            { 
                hts[[i]] <- c(pp[i]*hts[[i]], max(c(fhat$estimate[[i]], pp[i]*hts[[i]])) + 0.01 * max(abs(fhat$estimate[[i]])))
                names(hts[[i]])[length(hts[[i]])] <- "101%"
            }
        }
    }
    else if (oc %in% c("kdde","kqdde"))
    {
        hts1 <- hts[1,]; hts2 <- hts[2,]
        names(hts1) <- paste0("-", names(hts1))
        hts <- c(sort(hts1), sort(hts2))
        if (add_minmax) 
        {
            hts <- c(min(c(sapply(fhat$estimate, min), hts)) - 0.01 * max(abs(sapply(fhat$estimate, max))), hts, max(c(sapply(fhat$estimate, max), hts)) + 0.01 * max(abs(sapply(fhat$estimate, max))))
            names(hts)[1] <- "-101%"
            names(hts)[length(hts)] <- "101%"
        }
    }
    else 
    {
        hts <- sort(hts)
        if (add_minmax) 
        {
            hts <- c(hts, max(c(fhat$estimate, hts)) + 0.01 * max(abs(fhat$estimate)))
            names(hts)[length(hts)] <- "101%"
        }
    }
    names(hts) <- gsub("%", "", names(hts))

    return (hts)
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
    if (jitter) x <- data.frame(do.call(cbind,lapply(x, jitter, fac=1)))
    return(x)
},
required_aes=c("x","y"),
dropped_aes=c("z", "weight")
)

geom_rug_ks <- function(mapping=NULL, data=NULL, stat="rug_ks", position="identity", ..., outside=FALSE, sides="bl", length=unit(0.03, "npc"), na.rm=FALSE, jitter=FALSE, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=stat, geom=GeomRugKs, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(outside=outside, sides=sides, length=length, na.rm=na.rm, jitter=jitter, ...))
}

GeomRugKs <- ggplot2::ggproto("GeomRug", GeomRug)

stat_rug_ks <- function(mapping=NULL, data=NULL, geom="rug_ks", position="identity", ..., na.rm=FALSE, show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=StatRugKs, geom=geom, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(na.rm=na.rm, ...))
}

StatRugKs <- ggplot2::ggproto("StatPointKs", ggplot2::Stat,
compute_group=function(self, data, scales, jitter=FALSE)
{
    x <- get_data_ks(data, var_ks="weight")
    if (jitter) x <- data.frame(do.call(cbind,lapply(x, jitter, fac=1)))
     return(x)
},
required_aes=c("x"),
dropped_aes=c("y", "z", "weight") 
)

## stat and geom functions for computing probability contours
## replace ggplot2::geom_contour and ggplot2::stat_contour
geom_contour_ks <- function(mapping=NULL, data=NULL, stat="contour_ks", position="identity", ..., cont=c(25,50,75), contperc=TRUE, breaks=NULL, digits=NULL, show.legend=NA, inherit.aes=TRUE)
{
    aest <- dplyr::as_label(mapping$colour)
    if (aest!="NULL") aest <- gsub("\\)", "", gsub("after_stat\\(", "", aest))
    
    ## breaks_panel = set of common breaks for all panels
    breaks_panel <- NULL
    glayer <- ggplot2::layer(data=data, mapping=mapping, stat=stat, geom=GeomContourKs, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(cont=cont, contperc=contperc, breaks=breaks, digits=digits, breaks_panel=breaks_panel, aest=aest, ...))
    
    ## clear previous fill scale + create new fill scale with breaks
    if (!is.null(data) & !is.null(breaks))
    {
        if (object_class(data, type="tks") %in% "kdde") 
        {
            if (any(names(breaks) %in% "breaks")) breaks1 <- breaks$breaks else breaks1 <- breaks
            if (inherits(breaks, c("matrix","data.frame","tbl_df"))) 
                breaks <- breaks[breaks1 > min(data$estimate) & breaks1 < max(data$estimate),]
            else
                breaks <- breaks[breaks1 > min(data$estimate) & breaks1 < max(data$estimate)] 
            transp_neutral <- object_class(data, type="ks") %in% "kdde" 
            glayer <- c(scale_asymmetric(data, breaks=breaks, transp_neutral=transp_neutral), glayer)
        }
    }
    
    return(glayer)
}

GeomContourKs <- ggplot2::ggproto("GeomContourFilled", GeomContour)

stat_contour_ks <- function(mapping=NULL, data=NULL, geom="contour_ks", position="identity", ..., show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=StatContourKs, geom=geom, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(...))
}

## replaces ggplot2::StatContour
## data = tidy_ks object
## default values needed here since this StatContourKs is invoked directly from geom_contour(stat="contour_ks") rather than reading in the parameters from geom_contour_ks()
StatContourKs <- ggplot2::ggproto("StatContourKs", ggplot2::StatContour,
dropped_aes="weight",
compute_group=function(self, data, scales, cont=c(25,50,75), contperc=FALSE, breaks=NULL, breaks_panel=NULL, digits=NULL, aest=NULL, bins=NULL, binwidth=NULL)
{
    fhat <- untidy_ks(data, var_ks="weight")
    oc <- head(class(fhat),n=1) ## equivalent to object_class(fhat, type="ks")
    null_breaks <- is.null(breaks)
    if (is.null(digits)) digits <- ifelse(oc %in% c("kcde", "kcopula", "kde.loctest", "kfs", "kdr"), 2, 4)
   
    ## compute ks probability contour levels
    if (is.null(breaks))
    {
        if (oc %in% "kde.loctest") 
        {
            breaks <- c(-10,-0.5,0.5,10)
            breaks.label <- fhat$label
        }
        else if (oc %in% "kfs")
        {
            breaks <- c(-10, 0.5, 10)
            breaks.label <- fhat$label
        }
        else 
        {
            if (isTRUE(fhat$type=="kcurv")) class(fhat) <- "kdecurv"
            breaks <- .contourLevels(fhat=fhat, cont=cont)
            if (is.list(breaks)) breaks <- breaks[[unique(data$group)]]
        }

        if (oc %in% "kcopula") names(breaks)[length(breaks)] <- "100%"
        else if (oc %in% "kde.loctest")
        {
            fhat.diff.lab <- fhat$fhat1
            fhat.diff.lab$estimate <- fhat$fhat.diff.pos$estimate - fhat$fhat.diff.neg$estimate
            d <- ncol(fhat$fhat1$H)
            diff.lab <- compute.tidy(fhat.diff.lab, d=d, tidy=FALSE)$estimate
            data$z <- diff.lab     
        }
        else if (oc %in% "kfs")
        {
            d <- ncol(fhat$H)
            data$z <- compute.tidy(fhat, d=d, tidy=FALSE)$estimate
        }
           
        isolines <- .ggplot2_xyz_to_isolines(data, breaks)
        path_df <- .ggplot2_iso_to_path(isolines, data$group[1])
        path_df$estimate <- path_df$level <- as.numeric(path_df$level)
        
        ## if #distinct values in path_df$level != #names.isolines    
        if (oc %in% c("kdde", "kqdde"))
        {
            nc <- length(cont)
            names.isolines <- c(-100 - as.numeric(gsub("%", "",tail(head(names(breaks), n=nc+1),n=-1))), 100 - as.numeric(gsub("%", "",head(tail(names(breaks), n=nc+1),n=-1))))
            if (length(unique(path_df$estimate)) != length(names.isolines)) 
            {
                breaks2 <- signif(breaks[names(breaks) %in% names.isolines],8)
                names.isolines <- names.isolines[!duplicated(breaks2)]
                breaks2 <- breaks2[!duplicated(breaks2)]
                levels2 <- as.character(signif(unique(path_df$level),8))
                path_df$level <- ordered(path_df$level, labels=names.isolines[breaks2 %in% levels2])
            }
            else  
                path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
        else if (oc %in% "kde.loctest")
        {
            names.isolines <- c("-0.5", "0.5")
            path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
        else if (oc %in% "kfs")
        {
            names.isolines <- c("0.5")
            path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
        else
        {
            path_df$level <- findInterval(path_df$level, as.numeric(names(table(path_df$level))))
            names.isolines <- 100-as.numeric(gsub("%", "", head(names(breaks), n=-1)))
            path_df$level <- ordered(path_df$level, labels=names.isolines)
        }
    }
    else
    {
        if (any(names(breaks) %in% "breaks")) breaks <- breaks$breaks
        else if (is.list(breaks)) breaks <- breaks[[unique(data$group)]]
        breaks <- sort(breaks)

        isolines <- .ggplot2_xyz_to_isolines(data, breaks)
        path_df <- .ggplot2_iso_to_path(isolines, data$group[1])
        path_df$estimate <- path_df$level <- as.numeric(path_df$level)
        path_df$level <- round_signif(as.numeric(path_df$level), digits=digits) 
        path_df$level <- ordered(path_df$level)
    }
    
    ## ordered factor used in original ggplot2::StatContourFilled
    ## can cause issues when combining with unordered factors in ggplots 
    class(path_df$level) <- "factor"
    path_df$contlabel <- path_df$level
     
    ## create more informative labels
    path_df <- create_label(path_df, digits=digits, add_contperc=null_breaks, is_kdde=oc %in% c("kdde", "kqdde", "kde.loctest"), is_filled_contour=FALSE)   
    path_df$estimate <- round_signif(path_df$estimate, digits=digits)

    ## breaks_panel = common breaks for all panels   
    if (is.null(breaks_panel))
        path_df$estimate <- factor(path_df$estimate, levels=unique(path_df$estimate)) 
    else
    {    
        breaks_panel_br <- sort(unique(breaks_panel$breaks))
        breaks_panel_df <- create_label(data.frame(estimate=breaks_panel_br, contlabel=breaks_panel_br), digits=digits, add_contperc=null_breaks, is_kdde=oc %in% c("kdde", "kde.loctest"), is_filled_contour=FALSE) 
        breaks_panel_df$estimate <- round_signif(breaks_panel_df$estimate, digits=digits)
        path_df$estimate <- factor(path_df$estimate, levels=breaks_panel_df$estimate)
        path_df$contregion <- factor(path_df$contregion, levels=breaks_panel_df$contregion)
    }

    if (oc %in% "kde.loctest") 
    {
        levels(path_df$contlabel) <- breaks.label
        levels(path_df$contperc) <- c("-50%", "50%")
        path_df$contregion <- path_df$contlabel
    }
    else if (oc %in% "kfs")
    {
        levels(path_df$contperc) <- "50%"
        levels(path_df$contlabel) <- breaks.label
        path_df$contregion <- path_df$contlabel
    }
    if (aest=="NULL") path_df$level <- path_df$contperc
    else if (aest %in% c("contregion", "contlabel", "contline", "estimate"))path_df$level <- path_df[[aest]]
    
    return(path_df)
},
required_aes=c("x","y","z"),
dropped_aes=c("z","weight") 
)

## stat and geom functions for computing filled probability contours
## replace ggplot2::geom_contour_filled & ggplot2::stat_contour_filled
geom_contour_filled_ks <- function(mapping=NULL, data=NULL, stat="contour_filled_ks", position="identity", ..., cont=c(25,50,75), contperc=TRUE, breaks=NULL, transp_neutral=NULL, digits=NULL, show.legend=NA, inherit.aes=TRUE)
{
    aest <- dplyr::as_label(mapping$fill)
    if (aest!="NULL") aest <- gsub("\\)", "", gsub("after_stat\\(", "", aest))
    
    breaks_panel <- NULL
    glayer <- ggplot2::layer(data=data, mapping=mapping, stat=stat, geom=GeomContourFilled, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(cont=cont, contperc=contperc, digits=digits, breaks=breaks, breaks_panel=breaks_panel, aest=aest, ...))
 
    ## clear previous fill scale + create new fill scale with breaks
    if (!is.null(data) & !is.null(breaks))
    {
        if (object_class(data, type="tks") %in% "kdde") 
        {
            if (any(names(breaks) %in% "breaks")) breaks1 <- breaks$breaks else breaks1 <- breaks
            if (inherits(breaks, c("matrix","data.frame","tbl_df"))) 
                breaks <- breaks[breaks1 > min(data$estimate) & breaks1 < max(data$estimate),]
            else
                breaks <- breaks[breaks1 > min(data$estimate) & breaks1 < max(data$estimate)] 
            transp_neutral <- object_class(data, type="ks") %in% "kdde" 
            glayer <- c(scale_asymmetric(data, breaks=breaks, transp_neutral=transp_neutral), glayer)
        }
    }
    
    ## if contline is aesthetic then don't display its legend
    if (any(names(mapping) %in% "linetype")) 
    {
        if (length(grep("contline", dplyr::as_label(mapping$linetype)))>0)
            glayer <- c(glayer, ggplot2::guides(linetype="none"))
    }

    return(glayer)
}

GeomContourFilledKs <- ggplot2::ggproto("GeomContourFilled", GeomContourFilled)

## replace ggplot2::stat_contour_filled & ggplot2::StatContourFilled
stat_contour_filled_ks <- function(mapping=NULL, data=NULL, geom="contour_filled_ks", position="identity", ..., show.legend=NA, inherit.aes=TRUE)
{
    ggplot2::layer(data=data, mapping=mapping, stat=StatContourFilledKs, geom=geom, position=position, show.legend=show.legend, inherit.aes=inherit.aes, params=list(na.rm=FALSE, ...))
}

## default values needed here since this StatContourFilledKs is invoked directly from geom_contour_filled(stat="contour_filled_ks") rather than reading in the parameters from geom_contour_filled_ks()
StatContourFilledKs <- ggplot2::ggproto("StatContourFilledKs", StatContourFilled,
compute_group=function(data, scales, cont=c(25,50,75), contperc=TRUE, breaks=NULL, breaks_panel=NULL, digits=NULL, aest=NULL, bins=NULL, binwidth=NULL)
{
    fhat <- untidy_ks(data, var_ks="weight")
    oc <- head(class(fhat),n=1) ## equivalent to object_class(type="ks")
    null_breaks <- is.null(breaks)
    if (is.null(digits)) digits <- ifelse(oc %in% c("kcde", "kcopula", "kde.loctest", "kfs", "kdr"), 2, 4)
  
    ## compute ks probability contour levels
    if (is.null(breaks))
    {
        if (oc %in% "kde.loctest")
        {
            d <- ncol(fhat$fhat1$H)
            fhat.diff.lab <- fhat$fhat1
            fhat.diff.lab$estimate <- fhat$fhat.diff.pos$estimate - fhat$fhat.diff.neg$estimate
            diff.lab <- compute.tidy(fhat.diff.lab, d=d, tidy=FALSE)$estimate
            data$z <- diff.lab
            breaks <- c(-1, 0, 1, 10)
        }
        else if (oc %in% "kfs")
        {
            d <- ncol(fhat$H)
            data$z <- compute.tidy(fhat, d=d, tidy=FALSE)$estimate
            breaks <- c(-1.0, 10, 10.1)
        }

        if (!(oc %in% "kde.loctest")) 
        {
            if (isTRUE(fhat$type=="kcurv")) class(fhat) <- "kdecurv"
            breaks <- .contourLevels(fhat=fhat, cont=cont)
        }
        if (oc %in% "kda") breaks <- breaks[[unique(data$PANEL)]]
       
        isobands <- .ggplot2_xyz_to_isobands(data, breaks)
        
        if (oc %in% c("kdde", "kqdde"))
        {
            nc <- length(cont)
            names(isobands) <- c(-100 - as.numeric(gsub("%", "",tail(head(names(breaks), n=nc+1),n=-1))), 100, 100 - as.numeric(gsub("%", "",head(tail(names(breaks), n=nc+1),n=-1))))
        }
        else if (oc %in% "kde.loctest")
        {
            contperc <- FALSE
            labels <- na.omit(data$weight[[1]]$label)
            names(isobands) <- c(labels[1], "", labels[2]) 
        }
        else if (oc %in% "kfs")
        {
            contperc <- FALSE
            labels <- na.omit(data$weight[[1]]$label)
            names(isobands) <- head(breaks, n=-1)
        }
        else
        {
            names(isobands) <- 100-as.numeric(gsub("%", "", head(names(breaks), n=-1)), "%")
        }
    }
    else
    {
        if (any(names(breaks) %in% "breaks")) 
            breaks <- breaks$breaks
        else if (is.list(breaks)) breaks <- breaks[[unique(data$group)]]
        breaks <- sort(breaks)
       
        if ((oc %in% c("kdde","kqdde","kda")))
        {
            breaks <- c(sort(breaks), max(c(sapply(fhat$estimate,max), breaks)) + 0.01 * max(abs(sapply(fhat$estimate,max))))
        }
        else
           breaks <- c(sort(breaks), max(c(fhat$estimate, breaks)) + 0.01 * max(abs(fhat$estimate)))
        isobands <- .ggplot2_xyz_to_isobands(data, breaks)
        names(isobands) <- round_signif(as.numeric(head(breaks, n=-1)), digits=digits)
    }

    ## local copy of unexported function from scales
    .scales_rescale_max <- function (x, to=c(0, 1), from=range(x, na.rm=TRUE)) { x/from[2] * to[2] }
    
    ## copied from original StatContourFilled
    path_df <- .ggplot2_iso_to_polygon(isobands, data$group[1])
    path_df$level <- ordered(path_df$level, levels=names(isobands))    
    path_df$level_low <- breaks[as.numeric(path_df$level)]
    path_df$level_high <- breaks[as.numeric(path_df$level) + 1]
    path_df$level_mid <- 0.5 * (path_df$level_low + path_df$level_high)
    path_df$nlevel <- .scales_rescale_max(path_df$level_high)
    
    ## create more informative labels
    ## contlabel = contour percentage (without %)
    ## contperc = contour percentage (with %)
    ## estimate/level = contour level
    ## contregion = contour level (with >= or <=)
    path_df$estimate <- as.numeric(path_df$level_low)
    path_df$estimate_high <- as.numeric(path_df$level_high)
    path_df$contlabel <- path_df$level
    path_df.temp <- create_label(path_df, digits=digits, add_contperc=null_breaks, is_kdde=oc %in% c("kdde", "kqdde", "kde.loctest"))
    
    if (oc %in% "kde.loctest")
    {   
        levels(path_df$contlabel)[levels(path_df$contlabel)==""] <- NA
        path_df$contperc <- factor(path_df$estimate)
        levels(path_df$contperc) <- c("-50%", NA, "50%")
        path_df$contregion <- path_df$contlabel
        path_df$estimate <- factor(path_df$estimate, labels=c(-0.5, 0, 0.5))
        levels(path_df$estimate)[2] <- NA
    }
    else if (oc %in% "kfs")
    {
        path_df$contperc <- path_df$contlabel
        levels(path_df$contperc) <- c("50%",NA)
        path_df$contlabel <- labels
        path_df$contlabel <- factor(path_df$contlabel)
        path_df$contregion <- path_df$contlabel
        path_df$estimate <- factor(path_df.temp$estimate)
    }
    else 
    {
        path_df$contperc <- path_df.temp$contperc 
        path_df$contline <- path_df.temp$contline
        path_df$contregion <- path_df.temp$contregion
        path_df$estimate <- round_signif(path_df$estimate, digits=digits)
     
        ## breaks_panel = common breaks for all panels   
        if (is.null(breaks_panel))
            path_df$estimate <- factor(path_df$estimate, levels=unique(path_df$estimate)) 
        else
        {    
            breaks_panel_br <- sort(unique(breaks_panel$breaks))
            breaks_panel_df <- create_label(data.frame(estimate=breaks_panel_br, contlabel=breaks_panel_br), digits=digits, add_contperc=null_breaks, is_kdde=oc %in% c("kdde", "kde.loctest"), is_filled_contour=FALSE) 
            breaks_panel_df$estimate <- round_signif(breaks_panel_df$estimate, digits=digits)
            path_df$estimate <- factor(path_df$estimate, levels=breaks_panel_df$estimate)
            path_df$contregion <- factor(path_df$contregion, levels=breaks_panel_df$contregion)
        }
    }

    if (aest=="NULL") 
    {
        if (null_breaks) path_df$level <- path_df$contperc
        else path_df$level <- path_df$estimate
    }
    else if (aest %in% c("contregion", "contlabel", "contline", "estimate")) path_df$level <- path_df[[aest]]

    return(path_df)
},
required_aes=c("x","y","z"),
dropped_aes=c("z","weight") 
)

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
