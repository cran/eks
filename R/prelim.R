#####################################################################
## Auxiliary functions for tidy_*/st_* functions
#####################################################################

## create more informative legend labels
## contlabel = contour percentage (without %)
## contperc = contour percentage (with %)
## estimate/level = contour level
## contregion = contour level (with >= or <=)
create_label <- function(x, digits, add_contperc, is_kdde, is_filled_contour=TRUE)
{
    est.ord <- order(x$estimate)
    
    ## create contregion label
    if (is_kdde)
    {
        if (is_filled_contour)
        {
            x$contregion_low <- paste0(ifelse(x$estimate>=0, ">", "<"), round_signif(x$estimate, digits=digits)) 
            x$contregion_high <- paste0(ifelse(x$estimate_high>=0, ">", "<"), round_signif(x$estimate_high, digits=digits)) 
            x$contregion <- ifelse(x$estimate_high<0 & x$estimate<0, x$contregion_high, x$contregion_low)
            x$contregion[x$estimate_high>=0 & x$estimate<0] <- gsub("<",">", x$contregion)[x$estimate_high>=0 & x$estimate<0]
        }
        else 
            x$contregion <- paste0(ifelse(x$estimate>=0, ">", "<"), round_signif(x$estimate, digits=digits))
    }
    else 
        x$contregion <- paste0(">", round_signif(x$estimate, digits=digits))

    x$contregion <- gsub("-", "\u2013", gsub(">", "\u2265 ", gsub("<", "\u2264 ", x$contregion)))
    x$contlabel <- factor(x$contlabel, levels=unique(x$contlabel))
    x$contregion <- factor(x$contregion, levels=unique(x$contregion))
   
    ## create contline label
    if (is_kdde & is_filled_contour)
    {
        x$contline <- ifelse(sign(x$estimate_high)==1 & sign(x$estimate)==-1, NA, 1L)
        x$contline <- factor(x$contline)
    }

    ## create cont percentage label  
    if (add_contperc) 
    {   
        x$contperc <- x$contlabel
        x$contperc <- factor(x$contperc, levels=unique(x$contperc[est.ord]))
        
        levels(x$contperc) <- paste0(levels(x$contperc), "%")
    }
    else 
    {
        x$contperc <- NA  
        x$contperc <- factor(x$contper)
    }
    
    return(x)
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

round_signif <- function(x, digits=4, format=FALSE) 
{
    x <- sapply(x, function(.) ifelse(.>=10^(digits-1), round(.), signif(.,digits))) 
    if (format) x <- format(x, nsmall=digits)

    return(x)
}

unfactor <- function(f) { as.numeric(levels(f))[f] }

## object class 
## type="tks" includes "kde", "kdde", "kda" ks object classes etc.
## type="ks" includes above + "kqde", "kqdde" quasi ks object classes
object_class <- function(x, type="tks")
{
    type <- match.arg(type, c("ks", "tks"))
    if (type=="ks") oc <- ifelse(inherits(x, "sf_ks"), class(x$tidy_ks$ks[[1]]), class(x$ks[[1]]))
    else if (type=="tks") oc <- ifelse(inherits(x, "sf_ks"), x$tidy_ks$tks[1], x$tks[1])

    return(oc)
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
    oct <- object_class(object, type="tks")

    if (inherits(object, "sf_ks"))
    {
        object <- object$tidy_ks[1,]
        if (missing(x)) x <- "Longitude"
        if (missing(y)) y <- "Latitude"
    }
    else
    {
        if (oct %in% "ksupp") { xnames <- names(untidy_ks(object)) }
        else if (oct %in% "kdr") 
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
        if (oct %in% c("kde", "kfs", "kde.boundary")) y <- "Density function"
        else if (oct %in% "kda") y <- "Weighted density"
        else if (oct %in% "kdde") y <- "Density derivative"
        else if (oct %in% "kcde") y <- "Distribution function"
        else if (oct %in% "kde.loctest") y <- expression("Density difference "*f1-f2)
    }
    if (oct %in% "kroc")
    {
        x <- expression("False positive rate"~~group("(", list(bar(specificity)), ")"))
        y <- "True positive rate (sensitivity)"
    }

    labsd <- ggplot2::labs(x=x, y=y, ...)

    return(labsd)
}

## default guides parameters
guides_ks <- function(object, kda_part=TRUE)
{
    if (inherits(object, "ggplot")) object <- object$data
    d <- dim_ks(object)
    if (inherits(object, "sf_ks")) object <- object$tidy_ks[1,]
    oct <- object_class(object, type="tks")

    if (d==1)
    {
        titlec <- NULL
        if (oct %in% "kda") { titlef <- "Classif label"; titlec <- "Classif label"; rev <- FALSE }
        else if (oct %in% "kde.loctest") { titlec <- "Signif difference" }
        else if (oct %in% "kfs") { titlec <- "" }
        if (!is.null(titlec)) gu <- ggplot2::guides(colour=ggplot2::guide_legend(title=titlec))
        else gu <- ggplot2::guides()
    }
    else if (d==2)
    {
        titlef <- NULL; titlec <- NULL; rev <- TRUE
        if (oct %in% c("kde", "kdcde", "kde.boundary", "kde.truncate", "kde.sp", "kde.balloon")) titlef <- "Density"
        else if (oct %in% "kcurv") titlef <- "Density curvature"
        else if (oct %in% "kda")
        {
            if (kda_part) { titlef <- "Classif label"; titlec <- "Density"; rev <- FALSE }
            else titlef <- "Density"
        }
        else if (oct %in% "kdde") { titlef <- ifelse(object_class(object, type="ks") %in% "kqdde", "Density", "Density derivative") }
        else if (oct %in% c("kcde", "kcopula")) titlef <- "Distribution"
        else if (oct %in% "kde.loctest") { titlef <- "Signif difference" }
        else if (oct %in% "kms") { titlef <- "Cluster"; rev <- FALSE }
        #else if (oct %in% c("kdr", "kfs")) { titlef <- ""; titlec <- ""; rev <- FALSE }
        else if (oct %in% "kfs") { titlef <- "Signif curv" }
        else if (oct %in% "kdr") { titlef <- "Density ridge" }
        else if (oct %in% "ksupp") titlef <- "Support\nconvex hull" 
        titlec <- titlef

        if ((titlef=="" & titlec=="")) gu <- ggplot2::guides(fill=ggplot2::guide_legend(title=NULL, reverse=rev), colour=ggplot2::guide_legend(title=NULL, reverse=rev))
        else if (!is.null(titlef) & !is.null(titlec)) gu <- ggplot2::guides(fill=ggplot2::guide_legend(title=titlef, reverse=rev), colour=ggplot2::guide_legend(title=titlec, reverse=rev))
        else if (!is.null(titlef) & is.null(titlec)) gu <- ggplot2::guides(fill=ggplot2::guide_legend(title=titlef, reverse=rev))
        else if (is.null(titlef) & !is.null(titlec)) gu <- ggplot2::guides(colour=ggplot2::guide_legend(title=titlec, reverse=rev))
        else gu <- ggplot2::guides()

        if (oct %in% "kdde") gu <- c(ggplot2::guides(linetype="none"), gu)
    }

    return(gu)
}

## default colour scale
scale_colour_ks <- function(object)
{
    if (inherits(object, "ggplot")) object <- object$data
    oc <- object_class(object, type="ks")
    oct <- object_class(object, "tks")
       
    if (inherits(object, "sf_ks")) { object <- object$tidy_ks[1, ]; sfobject <- TRUE }
    else { sfobject <- FALSE }

    if (oc %in% c("kcde", "kcopula")) sc <- ggplot2::scale_colour_viridis_d()
    else if (oc %in% "kda") sc <- colorspace::scale_colour_discrete_qualitative(palette="Dark2", na.translate=FALSE) 
    else if (oc %in% "kda") sc <- colorspace::scale_colour_discrete_sequential()
    else if (oc %in% c("kdde", "kqdde")) sc <- scale_colour_discrete_diverging_breaks()
    else if (oc %in% c("kde", "kqde")) 
    {
        if (oct %in% "kcurv") sc <- colorspace::scale_colour_discrete_sequential(h1=30, c1=360, c2=30) 
        else sc <- colorspace::scale_colour_discrete_sequential(palette="Heat")
    }
    else if (oc %in% "kde.loctest") sc <- colorspace::scale_colour_discrete_qualitative(palette="Dark2", rev=TRUE)
    else if (oc %in% "kdr") sc <- ggplot2::scale_colour_manual(values=6)     
    else if (oc %in% "kfs") sc <- ggplot2::scale_colour_manual(values=7) 
    else if (oc %in% "kms") sc <- colorspace::scale_colour_discrete_qualitative(palette="Set2")  
    else sc <- NULL

    return (sc)
}
scale_color_ks <- scale_colour_ks

scale_colour_discrete_diverging_breaks <- function(..., breaks, extend1)
{
    if (missing(breaks))
        sc <- colorspace::scale_colour_discrete_diverging(...)
    else 
    {
        ## unblanced number of negative + positive contour breaks
        if (inherits(breaks, c("matrix","data.frame","tbl_df"))) breaks <- breaks$breaks
        breaks <- sort(breaks)       
        breaks.rle <- rle(sign(breaks))
        breaks.len <- max(breaks.rle$lengths)
        col.ind.neg <- floor(seq(1, breaks.len, len=length(which(breaks<=0))))
        col.ind.pos <- ceiling(seq(1, breaks.len, len=length(which(breaks>0)))) + breaks.len 
        col.ind <- c(col.ind.neg, col.ind.pos)
        nmax <- max(col.ind)
    
        ## for sf_ks objects, need to add 1 extra colour here 
        if (extend1)
            sc <- colorspace::scale_colour_discrete_diverging(..., nmax=nmax+1, order=c(col.ind, nmax+1))  
        else
            sc <- colorspace::scale_colour_discrete_diverging(..., nmax=nmax, order=col.ind)
    }

    return(sc)
}

## default fill colour scale
scale_fill_ks <- function(object, transp_neutral)
{
    if (inherits(object, "ggplot")) object <- object$data
    if (inherits(object, "sf_ks")) { object <- object$tidy_ks[1,]; sfobject <- TRUE;  } else {sfobject <- FALSE }
    oc <- object_class(object, type="ks")
    oct <- object_class(object, "tks")

    if (oc %in% c("kcde", "kcopula")) sc <- ggplot2::scale_fill_viridis_d()
    else if (oc %in% "kda") sc <- colorspace::scale_fill_discrete_qualitative(palette="Dark2", na.translate=FALSE) 
    else if (oc %in% c("kdde", "kqdde"))
    {
        if (missing(transp_neutral)) transp_neutral <- oc %in% "kdde"
        sc <- scale_fill_discrete_diverging_breaks(transp_neutral=transp_neutral)
    }
    else if (oc %in% c("kde", "kqde")) 
    {
        if (oct %in% "kcurv") sc <- colorspace::scale_fill_discrete_sequential(h1=30, c1=360, c2=30)
        else sc <- colorspace::scale_fill_discrete_sequential(palette="Heat") 
    }
    else if (oc %in% "kde.loctest")
    {
        if (sfobject) sc <- colorspace::scale_fill_discrete_qualitative(palette="Dark2", rev=TRUE, order=1:2, na.translate=FALSE)
        else sc <- colorspace::scale_fill_discrete_qualitative(palette="Dark2", rev=TRUE, order=1:2, na.translate=FALSE) 
    }
    else if (oc %in% "kdr") sc <- ggplot2::scale_fill_manual(values=6) 
    else if (oc %in% "kfs") sc <- ggplot2::scale_fill_manual(values=7)
    else if (oc %in% "kms") sc <- colorspace::scale_fill_discrete_qualitative(palette="Set2")  
    else sc <- NULL

    return (sc)
}

## ggplot scale version of ks::col.diverging for unblanced # negative + positive breaks
scale_fill_discrete_diverging_breaks <- function(..., breaks, transp_neutral, extend1)
{
    if (missing(breaks))
    { 
        sc <- colorspace::scale_fill_discrete_diverging(...)
        if (transp_neutral) sc <- scale_transp_neutral(sc, insert_neutral=TRUE)
    }
    else 
    {   
        ## unblanced number of negative + positive contour breaks
        if (inherits(breaks, c("matrix","data.frame","tbl_df"))) breaks <- breaks$breaks
        breaks <- sort(breaks)
        breaks.rle <- rle(sign(breaks))
        breaks.len <- max(breaks.rle$lengths)
        col.ind.neg <- floor(seq(1, breaks.len, len=length(which(breaks<=0))))
        col.ind.pos <- ceiling(seq(1, breaks.len, len=length(which(breaks>0)))) + breaks.len 
        col.ind <- c(col.ind.neg, col.ind.pos)
        nmax <- max(col.ind)
        
        ## for sf_ks objects, need to add 1 extra colour here 
        if (extend1)
            sc <- colorspace::scale_fill_discrete_diverging(..., nmax=nmax+1, order=c(col.ind, nmax+1))  
        else
            sc <- colorspace::scale_fill_discrete_diverging(..., nmax=nmax, order=col.ind)
        if (transp_neutral) sc <- scale_transp_neutral(sc, insert_neutral=FALSE, ind_neutral=length(col.ind.neg)+extend1)    
    }
    
    return(sc)
}

## make individual colours in colour scale to be transparent 
## insert_neutral=TRUE then add transparent after "ind"th colour
## insert_neutral=FALSE then replace "ind"th colour with transparent 
scale_transp_neutral <- function(x, ind_neutral=NULL, insert_neutral=FALSE)
{  
    Scale <- ggplot2::ggproto("Scale", x, 
        palette = function(self, n)
        { 
            ## palette some times appears to have 1 more colour for NA class
            n <- n-1
            pal <- x$palette(n) 
            if (is.null(ind_neutral)) ind_neutral <- ceiling(stats::median(seq_len(n)))
            if (ind_neutral<1) ind_neutral <- 1
            if (ind_neutral>n) ind_neutral <- n
            if (insert_neutral) pal <- pal[c(1:ind_neutral, ind_neutral:n)] 
            pal[ind_neutral] <- "transparent"
            return(pal)
        }
    )

    return(Scale)
}

scale_fill_remove <- function() scale_remove("fill")
scale_fill_continuous_remove <- function() scale_remove("fill", "continuous")
scale_fill_discrete_remove <- function() scale_remove("fill", "discrete")

scale_color_remove <- scale_colour_remove <- function() scale_remove("colour")
scale_color_continuous_remove <- scale_colour_continuous_remove <- function() scale_remove("colour", "continuous")
scale_color_discrete_remove <- scale_colour_discrete_remove <- function() scale_remove("colour", "discrete")

## 
scale_asymmetric <- function(object, breaks, transp_neutral)
{
    if (missing(object))
    {
        if (missing(transp_neutral)) transp_neutral <- FALSE
        sfobject <- TRUE 
    }
    else
    {
        if (inherits(object, "ggplot")) object <- object$data
        if (inherits(object, "sf_ks")) { object <- object$tidy_ks[1,]; sfobject <- TRUE;  } else { sfobject <- FALSE }
        oct <- object_class(object, type="tks")

        if (oct %in% "kdde")
        {
            if (missing(transp_neutral)) transp_neutral <- object_class(object,type="ks") %in% "kdde"
        }
    }
    if (missing(breaks))
    {
        scc <- scale_colour_discrete_diverging_breaks()
        scf <- scale_fill_discrete_diverging_breaks(transp_neutral=transp_neutral)
    }
    else 
    {
        scc <- scale_colour_discrete_diverging_breaks(breaks=breaks, extend1=sfobject)
        scf <- scale_fill_discrete_diverging_breaks(breaks=breaks, transp_neutral=transp_neutral, extend1=TRUE)
    }
    sc <- list(scale_colour_discrete_remove(), scc, scale_fill_discrete_remove(), scf)

    return(sc)
}

## remove scales from ggplot object to avoid multiple scales warning message 
ggplot_add.scale_rem <- function(object, plot, object_name) 
{
    if (length(plot$scales$scales)>0)
    {
        aes.ind <- which(sapply(plot$scales$scales, function(.) .$aesthetics %in% object[1] & inherits(., paste0("Scale", tools::toTitleCase(object[2])))))
        if (length(aes.ind)>0) plot$scales$scales[[aes.ind]] <- NULL
    }

    return(plot)
}

scale_remove <- function(scale1, scale2) 
{ 
    scale2 <- ifelse(missing(scale2), "",  sapply(scale2, function(.) match.arg(., c("continuous","discrete"))))
    rs <- structure(c(ggplot2::standardise_aes_names(scale1), scale2), class="scale_rem") 

    return(rs)
}
