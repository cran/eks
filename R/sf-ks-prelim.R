#####################################################################
## Auxiliary functions for st_* functions
#####################################################################

## plot function for sf_ks object
ggplot.sf_ks <- function(data=NULL, mapping=aes(), ..., which_geometry="sf")
{
    if (!is.null(data))
    {
        g <- match.arg(which_geometry, c("sf","grid"))
        data1 <- data[[g]] 
        transp_neutral <- object_class(data, type="ks") %in% "kdde"
    }
    else transp_neutral <- FALSE
    p <- ggplot2::ggplot(data=data1, mapping=mapping, ...) + labs_ks(data) + guides_ks(data)
    p <- p + scale_colour_ks(data) + scale_fill_ks(data, transp_neutral=transp_neutral)

    return(p)
}

plot.sf_ks <- function(x, ...) { plot_sf_ks(x=x, ...) }

plot_sf_ks <- function(x, which_geometry="sf", cont=c(25,50,75), abs_cont=breaks, breaks, which_deriv_ind=1, main="", pal, col, border, lty, transp_neutral, alpha, legend=TRUE, key.pos=NULL, legend.title, pos="bottomleft", digits, ...) 
{
    g <- match.arg(which_geometry, c("sf","grid"))
    y <- x[[g]]
    oct <- head(object_class(x, type="tks"), n=1)
    oc <- object_class(x, type="ks")
    if (missing(digits)) digits <- ifelse(oc %in% c("kcde", "kcopula", "kde.loctest", "kfs", "kdr"), 2, 4)

    if (missing(transp_neutral)) transp_neutral <- object_class(x, type="ks") %in% "kdde"
    if (missing(alpha)) { if (g=="sf") alpha <-1 else alpha <- 0.1 }
    
    if (missing(col))
    {
        if (oct %in% "kdr") col <- ggplot2::alpha(6, alpha=alpha)
        else if (oct %in% "kfs") col <- ggplot2::alpha(7, alpha=alpha)
        else if (oct %in% "ksupp") col <- NA
        else if (oct %in% "kquiver") col <- ggplot2::alpha(1, alpha=alpha)
    
        pal.missing <- FALSE
        if (missing(pal)) 
        {
            pal.missing <- TRUE
            if (oct %in% "kdde") 
            { 
                pal1 <- function(.) colorspace::diverging_hcl(n=., palette="Blue-Red", alpha=alpha)
                if (missing(breaks)) breaks1 <- sort(c(-(1:length(cont)), 1:length(cont)))
                else 
                {
                    if (inherits(breaks, c("matrix", "data.frame", "tbl_df"))) breaks1 <- breaks$breaks else breaks1 <- breaks
                    breaks1 <- sort(breaks1)
                }
                pal <- function(.) { ks::col.diverging(f=pal1, levels=breaks1, display="filled.contour", transp.neutral=transp_neutral) }
                col <- pal()
            }
            else if (oct %in% "kcurv") pal <- function(.) { colorspace::sequential_hcl(n=., h1=30, c1=360, c2=60, alpha=alpha, rev=TRUE) }
            else if (oct %in% "kde.loctest") pal <- function(.) { colorspace::qualitative_hcl(n=., palette="Dark2", alpha=alpha, rev=TRUE) }
            else if (oct %in% "kcde") pal <- function(.) { colorspace::sequential_hcl(n=., rev=FALSE, palette="Viridis", alpha=alpha) }
            else if (oct %in% c("kroc","kms")) pal <- function(.) { colorspace::qualitative_hcl(n=., palette="Set2", alpha=alpha) }
            else if (any(oct %in% "kda")) pal <- function(.) { colorspace::qualitative_hcl(n=., palette="Dark2", alpha=alpha) }
            else pal <- function(.) { colorspace::sequential_hcl(n=., rev=TRUE, palette="Heat", alpha=alpha) }
        }
        if ((oct %in% "kda") & (g=="sf"))
        {
            ng <- length(table(dplyr::group_indices(y)))
            cols <- pal(ng)
            if (missing(breaks)) nc <- length(cont) else nc <- length(breaks)
            col <- rep(cols, each=nc)
        }
    }
    else 
    {
        col <- ggplot2::alpha(col, alpha=alpha)
        pal.missing <- TRUE
    }

    if (g=="sf")
	{       
        if (oct %in% c("kdr", "kms"))
        {
            yd <- dplyr::select(y, dplyr::all_of("label"))
        }
        else if (oct %in% "ksupp")
        {
            yd <- dplyr::select(y, dplyr::all_of("contlabel")) 
        }
        else if (oct %in% "kdde")
        {
            y <- st_get_contour(x, cont=cont, breaks=abs_cont, digits=digits) 
            y <- y[y$deriv_ind == which_deriv_ind,]
            y$estimate <- droplevels(y$estimate)
            yd <- dplyr::select(y, dplyr::all_of("estimate")) 

            if (length(y$estimate)!=length(col)) { breaks1 <- tail(unfactor(y$estimate),n=-1); col <- pal() }
        }
        else
        {
            y <- st_get_contour(x, cont=cont, breaks=abs_cont, digits=digits)
            yd <- dplyr::select(y, dplyr::all_of("estimate"))
        }
        if (missing(lty)) 
        { 
            if (oct %in% "kdde") 
            { lty <- as.numeric(y$contline); lty[is.na(lty)] <- 0 }
            else lty <- rep(1, nrow(yd))
            if (any(names(yd) %in% "estimate")) lty[y$estimate==0] <- 0 
        }
        
        if (missing(col)) plot(yd, main=main, pal=pal, lty=lty, key.pos=key.pos, ...)
        else 
        { 
            if (oct %in% "kda") 
            {
                if (missing(border)) plot(yd, main=main, col=NA, border=col, lty=lty, key.pos=key.pos, ...)
                else plot(yd, main=main, col=NA, border=border, lty=lty, key.pos=key.pos, ...)
            }
            else 
            {
                if (missing(border)) plot(yd, main=main, col=col, lty=lty, key.pos=key.pos, ...)
                else plot(yd, main=main, col=col, border=border, lty=lty, key.pos=key.pos, ...)
            }
        }
	}
	else if (g=="grid") 
	{
        if (oct %in% "kdde") y <- y[y$deriv_ind == which_deriv_ind,]
        yd <- dplyr::select(y, dplyr::all_of("estimate"))
        if (missing(lty)) lty <- 1
        if (oct %in% "kda")
        {
            if (missing(col)) plot(yd, main=main, pal=pal, lty=lty, border=border, key.pos=key.pos, ...)
            else plot(yd, main=main, col=col, lty=lty, border=border, key.pos=key.pos, ...)
        }
        else 
        {
            if (missing(col)) plot(yd, main=main, pal=pal, lty=lty,, key.pos=key.pos, ...)
            else plot(yd, main=main, col=col, lty=lty, key.pos=key.pos, ...)
        }
    }

    ## add legend
    ## mapsf legends don't allow for line types in legend boxes so 
    ## st_ksupp plots aren't well supported
    if (legend)
    {
        if (!requireNamespace("mapsf", quietly=TRUE)) stop("Install the mapsf package as it is required.", call.=FALSE)       
        forms <- list(...)
        forms <- forms[names(forms) %in% c("pos", "val", "pal", "title", "title_cex", "val_cex", "col_na", "no_data", "no_data_txt", "frame", "border", "bg", "fg", "cex")]
        forms <- forms[!(names(forms) %in% "border")] 
        gu <- guides_ks(dplyr::add_row(dplyr::ungroup(x$tidy_ks), ks=list(2L)))
       
        if (is.list(gu)) gu <- list(guides=unlist(lapply(gu, function(.) unlist(.[["guides"]]))))
        if (missing(legend.title)) gu.title <- gu$guides$fill$params$title
        else gu.title <- legend.title
        if (g=="sf")
        {
            if (oct %in% c("kdde", "kdr")) 
            {
                contlabel <- sort(unique(y$contregion), decreasing=TRUE) 
                nc <- length(contlabel)
            }
            else 
            { 
                contlabel <- sort(unique(y$contlabel), decreasing=TRUE)
                if (missing(breaks)) contlabel <- paste0(contlabel,"%")
                nc <- length(contlabel)
            }
        }
        
        if (g=="grid")
        {
            if (oct %in% "kda")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), .data$label))
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=gv, pal=pal(length(gv)), pos=pos, title=gu.title), forms))
            }
            else
            {
                y2 <- st_get_contour(x, cont=cont, breaks=abs_cont, digits=digits)
                yd2 <- dplyr::select(y2, dplyr::all_of("estimate"))
                gv <- yd2$estimate
                do.call(mapsf::mf_legend, args=c(list(type="choro", val=gv, pal=pal(length(gv)), pos=pos, title=gu.title), forms))
            }
        }
        else if (!is.null(forms$border)) 
        {
            if (oct %in% "kda")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), dplyr::all_of(dplyr::group_vars(y))))
                ng <- length(gv)
                do.call(mapsf::mf_legend, args=c(list(type="symb", val=gv, pal=unique(forms$border), pos=pos, title=gu.title, cex=rep(3,ng), pch=rep("-", ng)), forms[!(names(forms) %in% c("cex","pch"))]))
            }
            else 
            {
                if (!all(is.na(forms$border)) & pal.missing) 
                {
                    forms <- forms[!(names(forms) %in% c("cex","pch"))]
                    do.call(mapsf::mf_legend, args=c(list(type="symb", val=contlabel, pal=rev(forms$border), pos=pos, title=gu.title, cex=3, pch="-"), forms[!(names(forms) %in% c("cex","pch"))]))
                }
                else do.call(mapsf::mf_legend, args=c(list(type="typo", val=contlabel, pal=rev(pal(nc)), pos=pos, title=gu.title), forms)) 
            }
        }
        else if (!missing(col))
        {
            if (oct %in% "kdr") do.call(mapsf::mf_legend, args=c(list(type="symb", val=levels(y$label), pal=col,  pos=pos, title=gu.title, cex=3, pch="-"), forms[!(names(forms) %in% c("cex","pch"))]))
            else if (oct %in% "kfs") do.call(mapsf::mf_legend, args=c(list(type="typo", val=y$label[1], pos=pos, pal=col, title=gu.title), forms))  
            ## mapsf::mf_legend doesn't accept lty argument
            else if (oct %in% "ksupp") do.call(mapsf::mf_legend, args=c(list(type="typo", val=contlabel, pos=pos, pal=col, title=gu.title), forms))
            else if (oct %in% "kda")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), dplyr::all_of(dplyr::group_vars(y))))
                ng <- length(gv)
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=gv, pal=pal(ng), pos=pos, title=gu.title), forms))
            }
            else 
            {
                if (!missing(breaks) | oct %in% "kdde") { col <- rev(col) }
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=contlabel, pal=col, pos=pos, title=gu.title), forms))
            }
        }
        else if (!missing(pal)) 
        {        
            if (oct %in% "kms") 
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), .data$label))
                nc <- length(gv)
                formsg <- list(...)
                if (is.null(formsg$pch)) pch <- 1 else pch <- unique(formsg$pch) 
                forms <- forms[!(names(forms) %in% c("cex","pch"))]
                do.call(mapsf::mf_legend, args=c(list(type="symb", val=gv, pal=pal(nc), pos=pos, title=gu.title, cex=rep(1,nc), pch=pch), forms[!(names(forms) %in% c("cex","pch"))]))
            }
            else if (oct %in% "kde.loctest")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), .data$label))
                ind <- order(gv, decreasing=TRUE)
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=gv[ind], pal=pal(length(gv))[ind], pos=pos, title=gu.title), forms))
            }
            else 
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=contlabel, pal=rev(pal(nc)), pos=pos, title=gu.title), forms))
        }
    }
}

## probability contour regions
contourLevels.sf_ks <- function(x, ..., cont=c(25,50,75)) { ks::contourLevels(x=x$tidy_ks, ..., cont=cont) }

## convert contour polygons in ks object to multipolygon geometry
## x = ks object
## output sf-geometry is named "geometry"   
.st_contourline <- function(x, cont, abs.cont, edge_zero=TRUE, digits)
{   
    ## flag to use values of abs.cont in contlabel rather than names(abs.cont)
    if (missing(abs.cont)) abs.cont.flag <- FALSE 
    else if (is.null(names(abs.cont))) abs.cont.flag <- TRUE
    else abs.cont.flag <- !all(grepl("[0-9]%", names(abs.cont)))
    
    if (missing(abs.cont) & missing(cont)) abs.cont <- x$cont
    if (!missing(cont)) 
    {
        if (class(x) %in% c("kfs", "kde.loctest")) abs.cont <- abs.cont
        else abs.cont <- contourLevels(x, cont=cont)
    }
    abs.cont <- sort(abs.cont)
   
    ## convert contourLines output to linestring
    cline2ls <- function(y) 
    { 
        sf::st_linestring(as.matrix(data.frame(x=y["x"], y=y["y"]))) 
    }
    
    ## add zero estimate values at edges of estimation grid 
    if (edge_zero)
    {
        x.zero <- x
        delta <- sapply(x$eval.points, function(.) head(diff(.), n=1))

        a <- 1
        x.zero$eval.points[[1]] <- c(min(x.zero$eval.points[[1]]) - delta[1]*(a:1), x.zero$eval.points[[1]], max(x.zero$eval.points[[1]]) + delta[1]*(1:a))
        x.zero$eval.points[[2]] <- c(min(x.zero$eval.points[[2]]) - delta[2]*(a:1), x.zero$eval.points[[2]], max(x.zero$eval.points[[2]]) + delta[2]*(1:a))
        x.zero$estimate <- matrix(0, ncol=ncol(x.zero$estimate)+2*a, nrow=nrow(x.zero$estimate)+2*a)
        x.zero$estimate[(1:nrow(x$estimate))+a, (1:ncol(x$estimate))+a] <- x$estimate
        x <- x.zero
    }

    j <- 0   
    for (i in 1:length(abs.cont))
    {
        if (is.list(x$estimate)) estimate <- x$estimate[[1]] else estimate <- x$estimate
    	cline <- contourLines(x=x$eval.points[[1]], y=x$eval.points[[2]], z=estimate, levels=abs.cont[i])
        cline.len <- lengths(lapply(cline, getElement, "x"))
        cline <- cline[cline.len>=4]
     	cont.polygon <- lapply(cline, cline2ls)
    	cont.polygon <- try(sf::st_polygon(cont.polygon), silent=TRUE)
        if (inherits(cont.polygon, "try-error")) cont.polygon <- sf::st_polygon()
        cont.polygon <- sf::st_sf(geometry=sf::st_geometry(cont.polygon))

        if (abs.cont.flag) 
            cont.polygon <- cbind(contlabel=round_signif(abs.cont[i], digits=digits), estimate=as.numeric(abs.cont[i]), cont.polygon) 
    	else 
            cont.polygon <- cbind(contlabel=as.numeric(gsub("%", "", names(abs.cont[i]))), estimate=as.numeric(abs.cont[i]), cont.polygon) 	
        j <- j + !sf::st_is_empty(cont.polygon)
    	if (j<=1) { cont.polygon.all <- cont.polygon ; j <- j+1 }
        else if (j>1) cont.polygon.all <- rbind(cont.polygon.all, cont.polygon)
    }

    cont.polygon.all <- cont.polygon.all[!sf::st_is_empty(cont.polygon.all),]
    cont.polygon.all <- dplyr::arrange(cont.polygon.all, .data$estimate, .data$contlabel)
    
    if (nrow(cont.polygon.all)>0)
    {
        cont.polygon.all <- sf::st_make_valid(cont.polygon.all, geos_method="valid_linework")
        cont.polygon.all <- dplyr::distinct(cont.polygon.all)
       
        if (abs.cont.flag) 
        {
            cont.polygon.all$contlabel <- factor(cont.polygon.all$contlabel)
        }
        else
        {
            cont.polygon.all$contlabel <- factor(cont.polygon.all$contlabel)
            levels(cont.polygon.all$contlabel) <- as.character(100-as.numeric(levels(cont.polygon.all$contlabel)))
        }
    }
   
	return(cont.polygon.all)
}

## x = tibble with ks column (output from first part of st_ks)
st_contourline <- function(x, cont, abs.cont, edge_zero=TRUE, digits, ...)
{
    ## convert contour polygons to multipolygon geometry
    if (missing(abs.cont)) cont.geom <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_contourline(untidy_ks(.x), cont=cont, edge_zero=edge_zero, digits=digits))) 
    else if (!missing(abs.cont)) cont.geom <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_contourline(untidy_ks(.x), abs.cont=abs.cont, edge_zero=edge_zero, digits=digits))) 
    
    return(cont.geom)
}

## create set of disjoint contour polygons
## x = multipolygon geometry of contour regions
st_contour_disjoint <- function(x)
{
    gv <- dplyr::group_vars(x)
    xc <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_contour_disjoint(.x))) 
    xc <- sf::st_sf(xc)
    xc <- sf::st_cast(xc, to = "MULTIPOLYGON")
    xc <- dplyr::ungroup(xc)
    if (length(gv)>0) xc <- dplyr::group_by(xc, dplyr::across({{gv}}))
    xc$.contlabel <- NULL

    return(xc)
}

.st_contour_disjoint <- function(x)
{
    if (any(as.numeric(as.character(x$contlabel))<0)) oc <- "kdde" else oc <- "kde"
    if (any(oc %in% "kdde"))
    {
        cp <- x
        cp <- dplyr::mutate(cp, estimate=as.numeric(as.character(.data$estimate)), contlabel=as.numeric(as.character(.data$contlabel)))
        if (!any(names(cp) %in% "deriv_ind")) cp <- dplyr::arrange(cp, .data$estimate) 
        else cp <- dplyr::arrange(cp, .data$deriv_ind, .data$estimate)

        cp1 <- dplyr::filter(cp, .data$contlabel<=0)
        cp1 <- sf::st_geometry(cp1)
        cp1 <- st_difference_sequence(cp1, headtail=FALSE)
        cp1 <- sf::st_sfc(cp1)
        cp1 <- c(head(sf::st_geometry(dplyr::filter(cp, .data$contlabel<=0)),n=1),cp1)
        cp1 <- sf::st_set_crs(cp1, sf::st_crs(x))
        
        cp2 <- dplyr::filter(cp, .data$contlabel>0)
        cp2 <- sf::st_geometry(cp2)
        cp2 <- st_difference_sequence(cp2, headtail=TRUE)
        cp2 <- sf::st_sfc(cp2)
        cp2 <- c(cp2, tail(sf::st_geometry(dplyr::filter(cp, .data$contlabel>0)),n=1))
        cp2 <- sf::st_set_crs(cp2, sf::st_crs(x))
        x <- sf::st_set_geometry(x, value=c(cp1,cp2))
    }
    else
    { 
        cp <- sf::st_geometry(x)
        cp <- st_difference_sequence(cp, headtail=TRUE)
        cp <- sf::st_sfc(cp)
        cp <- c(cp, tail(sf::st_geometry(x),n=1))
        cp <- sf::st_set_crs(cp, sf::st_crs(x))
        x <- sf::st_set_geometry(x, value=cp)

    }
    x <- sf::st_make_valid(x, geos_method="valid_linework") 
    
    return(x)
}

## compute st_difference on pairs of adjacent rows in sf object
## x = sf object
st_difference_sequence <- function(x, headtail=TRUE) 
{ 
    if (length(x)>1)
    {
        for (i in 1:(length(x)-1)) 
        {
            xh <- head(x, n=-1)
            xt <- tail(x, n=-1)
            if (headtail) 
            {
                xpt <- sf::st_difference(xh[i], xt[i])
                if (length(xpt)==0) xpt <- sf::st_sfc(sf::st_polygon(), crs=sf::st_crs(x))
                if (i==1) xp <- xpt else xp <- c(xp, xpt)
            }
            else
            {
                xpt <- sf::st_difference(xt[i], xh[i]) 
                if (length(xpt)==0) xpt <- sf::st_sfc(sf::st_polygon(), crs=sf::st_crs(x))
                if (i==1) xp <- xpt else xp <- c(xp, xpt)
            }
        }
    }
    else xp <- x[0,]
    
    return(xp)
}

## x = tibble with ks column (output from first part of st_kdde)
st_contourline_kdde <- function(x, cont, abs.cont, which_deriv_ind=1, digits)
{
    x.temp <- x
    if (missing(abs.cont)) 
    { 
        abs.cont <- contourLevels(x, cont=cont, which.deriv.ind=which_deriv_ind)
        abs.cont.flag <- FALSE
        abs.cont <- sort(c(abs.cont[1,], abs.cont[2,]))
    }
    else 
    {
        if (!is.vector(abs.cont) & !is.null(abs.cont$breaks)) abs.cont <- abs.cont$breaks
        abs.cont <- sort(abs.cont)
        abs.cont.flag <- TRUE
    }
    abs.cont <- abs.cont[!duplicated(abs.cont)]
    x.temp$estimate <- x$estimate[[which_deriv_ind]]

    cont.geom <- .st_contourline(x=x.temp, abs.cont=abs.cont, digits=digits) 

    if (!abs.cont.flag) 
    {
        cont.geom <- dplyr::mutate(cont.geom, contlabel=ifelse(.data$estimate<0, paste0("-", as.character(.data$contlabel)), as.character(.data$contlabel)))
        cont.geom$contlabel <- factor(cont.geom$contlabel, levels=c(-sort(cont), sort(cont, decreasing=TRUE)))
    }
    else cont.geom$contlabel <- factor(cont.geom$contlabel)
    
    return(cont.geom)
}

## x = tibble with ks column (output from first part of st_kda)
st_contourline_kda <- function(x, cont, digits)
{
    xj <- untidy_ks(x)
    ## replace KDE with weighted KDE
    j <- x$rn
    pp <- xj$prior.prob[j]
    abs.cont <- contourLevels(xj, cont=cont)[[j]]*pp
    xj$estimate <- xj$estimate[[j]]*pp
    xj$H <- xj$H[[j]]
    xj$cont <- xj$cont[[j]]
    class(xj) <- "kde"
    
    cont.geom <- .st_contourline(x=xj, abs.cont=abs.cont, digits=digits) 
    cont.geom <- dplyr::mutate(cont.geom, contlabel=factor(ifelse(.data$estimate<0, paste0("-", 100-as.numeric(as.character(.data$contlabel))), as.character(.data$contlabel))))
    cont.geom <- dplyr::arrange(cont.geom, sign(.data$estimate), abs(.data$estimate))
    levels(cont.geom) <- as.character(cont.geom$contlabel)
       
    return(cont.geom)
}

## convert eval.points in ks object to rectangle polygon/point geometry
## x = tibble with ks column (output from first part of st_ks)
st_evalpoints <- function(x)
{
    ## convert ks grid to rectangle polygon geometry
    if (inherits(untidy_ks(x[1,]), "kde.loctest"))
    {
        fhat.geom.neg <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_evalpoints(untidy_ks(.x)$fhat.diff.neg)))
        fhat.geom.pos <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_evalpoints(untidy_ks(.x)$fhat.diff.pos)))
        fhat.geom.neg <- dplyr::mutate(fhat.geom.neg, estimate=as.integer(.data$estimate>=0.5))
        fhat.geom.pos <- dplyr::mutate(fhat.geom.pos, estimate=as.integer(.data$estimate>=0.5))
        fhat.geom <- fhat.geom.pos
        fhat.geom$estimate <- fhat.geom.pos$estimate - fhat.geom.neg$estimate
    }
    else 
        fhat.geom <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_evalpoints(untidy_ks(.x))))
    fhat.geom <- sf::st_as_sf(fhat.geom, crs=sf::st_crs(x))
       
    return(fhat.geom)
}

## convert eval.points in ks object to rectangle polygon geometry
## x = ks object
## ouput sf-geometry column is named "geometry"
.st_evalpoints <- function(x)
{
    ## convert KDE to rectangle polygon geometry
    nx <- lengths(x$eval.points)[1]
    ny <- lengths(x$eval.points)[2]
    ep1 <- x$eval.points[[1]]
    ep2 <- x$eval.points[[2]]
    fhat.point <- sf::st_as_sf(x=expand.grid(x=ep1[c(1,nx)], y=ep2[c(1,ny)]), coords=1:2)
    fhat.geom <- sf::st_make_grid(fhat.point, n=c(nx-1,ny-1))
    fhat.geom.cent <- sf::st_centroid(sf::st_geometry(fhat.geom))
    fhat.geom.cent <- sf::st_coordinates(fhat.geom.cent)
    ## "geometry" is automatically given sf-geometry column when concatenating into data frame
    fhat.geom <- data.frame(estimate=predict(x, x=fhat.geom.cent), fhat.geom)
  
    return(fhat.geom)
}

## remove segments in linestring where distance between two consecutive points > len 
## (metres) useful to neaten appearance of polygons which are split to separate edges of 
## world map after projection  
## x = linestring or multipoint geometry
st_remove_long_segment <- function(x, len=500e3)
{
    x.geom <- dplyr::group_modify(x, .f=~dplyr::as_tibble(sf::st_coordinates(sf::st_geometry(.))))
    gvar <- setdiff(names(x.geom), c("X","Y")) 
    x.geom <- dplyr::group_by(x.geom, dplyr::across({{gvar}}))
    x.geom <- dplyr::mutate(x.geom, L1=cumsum(c(0, sqrt(diff(.data$X)^2 + diff(.data$Y)^2) > len)))
    x.geom <- dplyr::group_by(x.geom, dplyr::across({{gvar}}))
    valid.LS <- names(table(x.geom$L1))[table(x.geom$L1)>2]
    
    if (length(valid.LS)>0) 
    {
        x.geom <- dplyr::filter(x.geom, .data$L1 %in% valid.LS)
        x.geom <- dplyr::summarise(sf::st_as_sf(x.geom, coords=c("X","Y"), crs=sf::st_crs(x)), do_union=FALSE, .groups="drop")
        x.geom <- sf::st_cast(x.geom, to="LINESTRING")
        x.geom <- dplyr::filter(x.geom, as.numeric(sf::st_length(.data$geometry))>0)
        x.geom <- dplyr::group_by(x.geom, dplyr::across({{gvar}}))
        x.geom <- dplyr::summarise(sf::st_sf(x.geom), do_union=FALSE, .groups="drop")
        sfname <- attr(x.geom, "sf_column")
    
        x2 <- dplyr::right_join(sf::st_drop_geometry(x), x.geom, by=dplyr::group_vars(x))
        x2 <- dplyr::relocate(x2, !!sfname, .after=dplyr::last_col())
        x <- sf::st_sf(x2)
    }
    
    return(x)
}

## concatenate sf objects
c_sf <- function(..., labels)
{
    args <- list(...)
    nargs <- length(args)
    if (nargs>=1) 
    {
        sfname <- attr(args[[1]], "sf_column")
        tt <- dplyr::mutate(args[[1]], group=1L)
        if (nargs>=2)
        {
            for (i in 2:nargs) 
                tt <- dplyr::bind_rows(tt, dplyr::mutate(args[[i]], group=i))
        }
        if (!missing(labels)) tt$group <- factor(tt$group, labels=labels)
        tt <- dplyr::relocate(tt, !!sfname, .after=dplyr::last_col())
        tt <- sf::st_sf(tt)
    }
    else tt <- NULL

    return(tt)
}

## concatenate sf_ks objects
#c.sf_ks <- function(..., labels) { c_sf_ks(..., labels=labels) }

c.sf_ks <- function(..., labels)
{
    args <- list(...)
    nargs <- length(args)
    if (nargs>=1) 
    {
        csf <- list()
        csf$tidy_ks <- do.call(c.tidy_ks, lapply(args, getElement, "tidy_ks"))
        csf$grid <- do.call(c_sf, lapply(args, getElement, "grid"))
        csf$sf <- do.call(c_sf, lapply(args, getElement, "sf"))

        if (!missing(labels)) 
        {
            csf$tidy_ks$group <- factor(csf$tidy_ks$group, labels=labels)
            csf$grid$group <- factor(csf$grid$group, labels=labels)
            csf$sf$group <- factor(csf$sf$group, labels=labels)
        }

        csf$tidy_ks <- dplyr::group_by(csf$tidy_ks, .data$group)
        csf$grid <- dplyr::group_by(csf$grid, .data$group)
        csf$sf <- dplyr::group_by(csf$sf, .data$group)
        csf$tidy_ks <- as_tidy_ks(csf$tidy_ks)
    }
    else csf <- NULL
    class(csf) <- "sf_ks"

    return(csf)
}

rbind.sf_ks <- c.sf_ks

## add coordinates as attributes to point geometry
st_add_coordinates <- function(x, as_sf=TRUE, rename=TRUE)
{
    sfname <- attr(x, "sf_column")
    xc <- sf::st_cast(x, to="POINT", warn=FALSE)
    xcoord <- sf::st_coordinates(sf::st_geometry(xc))
    xc$X <- xcoord[,1]
    xc$Y <- xcoord[,2]
    
    if (!as_sf)  xc <- sf::st_drop_geometry(xc)
    else xc <- dplyr::relocate(xc, dplyr::all_of(c("X", "Y")), .before=!!sfname)
    if (inherits(x, "tbl")) xc <- dplyr::as_tibble(xc)
    if (as_sf)  xc <- sf::st_sf(xc)
    if (rename) xc <- dplyr::rename(xc, lon=dplyr::all_of("X"), lat=dplyr::all_of("Y"))

    return(xc)
}

## returns segments of bounding box of x
## h1 = (xmax, ymin) -> (xmin, ymin)
## v1 = (xmin, ymin) -> (xmin, ymax)
## h2 = (xmin, ymax) -> (xmax, ymax)
## v2 = (xmax, ymax) -> (xmax, ymin)
st_bbox_segments <- function(x)
{
    xbbox <- sf::st_bbox(x)
    xcrs <- sf::st_crs(x)

    ## vertical and horizontal lines of bounding box
    h1 <- sf::st_sfc(sf::st_linestring(rbind(xbbox[c(3,2)], xbbox[c(1,2)])), crs=xcrs)
    v1 <- sf::st_sfc(sf::st_linestring(rbind(xbbox[c(1,2)], xbbox[c(1,4)])), crs=xcrs)
    h2 <- sf::st_sfc(sf::st_linestring(rbind(xbbox[c(1,4)], xbbox[c(3,4)])), crs=xcrs)
    v2 <- sf::st_sfc(sf::st_linestring(rbind(xbbox[c(3,4)], xbbox[c(3,2)])), crs=xcrs)

    return(c(h1, v1, h2, v2))
}

## returns segments of rectangle
## first value = len of "horizontal" side (smallest Hausdorff dist to h1 of st_bbox_segments) 
## second value = len of "vertical" side (smallest Hausdorff dist to v1 of st_bbox_segments) 
st_rectangle_segments <- function(x)
{
    sfname <- attr(x, "sf_column")
    xbbox.seg <- st_bbox_segments(x)
    xseg <- suppressWarnings(sf::st_collection_extract(lwgeom::st_split(sf::st_cast(x[1,], to="LINESTRING"), sf::st_cast(x[1,], to="POINT", warn=FALSE)), "LINESTRING"))
    xseg.dist1 <- sf::st_distance(xseg, xbbox.seg, which="Hausdorff")
    xseg.dist2 <- sf::st_distance(xseg, xbbox.seg)
    units(xseg.dist1) <- units(xseg.dist2) <- NULL
    xseg.dist <- xseg.dist1 + xseg.dist2
    xseg.ind <- apply(xseg.dist, 2, which.min)
    if (any(class(xseg) %in% "sfc"))
        xseg <- sf::st_sf(.side=xseg.ind, geometry=xseg)
    else 
        xseg <- dplyr::mutate(xseg, .side=!!xseg.ind, .before=!!sfname)
    xseg <- dplyr::arrange(xseg, .data$.side)
    xseg <- sf::st_sf(xseg)

    return (xseg)
}

## add "%" suffix - deprecated but kept for backwards compatibility 
label_percent <- function(y)
{
    factor(y, labels=paste0(levels(y),"%"))
}

## read GPKG and renames "geom" as "geometry"
st_read_gpkg <- function(..., as_tibble=TRUE)
{
    sf::st_set_geometry(sf::st_read(..., as_tibble=as_tibble), "geometry")
}

## convenience function to write to progress geopackage for diagnostics
write_temp <- function(x, layer, dsn) 
{ 
    if (!inherits(x, "sf")) x <- sf::st_sf(x)
    if (missing(layer))
    {
        if (any(sf::st_is(x, "POINT"))) layer <- "temp3"
        else if (any(sf::st_is(x, "POLYGON"))) layer <- "temp4"
        else layer <- "temp"
    }
    epsg <- sf::st_crs(x)$epsg
    if (missing(dsn)) 
    { 
        if (is.na(epsg)) dsn <- "~/Downloads/flowmap_progress.gpkg"
        else dsn <- paste0("~/Downloads/flowmap_progress_", epsg,".gpkg")
    }
    x <- dplyr::select(x, !dplyr::where(is.list))

    sf::write_sf(x, dsn=dsn, layer=layer)
}

